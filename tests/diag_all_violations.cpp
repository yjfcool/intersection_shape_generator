// 全量审计：逐份数据集运行一次完整生成，再按“检查类别”逐条输出审计发现。
//
// 输出两路：
//   1) stdout —— 与历史格式兼容（`VIOLATION:` / `PRECURSOR:` / 每份汇总行 /
//      `TOTAL VIOLATIONS`），既有 awk 逐对差分脚本可继续使用；末尾新增“按检查
//      类别”和“按文件”的统计表。
//   2) CSV（`--csv <路径>`）—— 统一发现记录，列与 tools/test_report.py 约定一致：
//      source,category,check,dataset,subject,severity,metric,value,threshold,detail,location
//      由 tools/test_report.py 汇总成“按检查类别分 sheet”的 Excel 报表。
//
// 检查来源分两类：
//   A) 与生成侧同源的约束评估器（ConstraintEvaluator + ConstraintProfile）——
//      围栏、障碍、Boundary、路沿净距、自交、段间 G1、普通形态、掉头形态、
//      人行横道。输出侧一旦报违约，说明生成期用的同一把尺子在最终结果上不成立。
//   B) 审计专有口径 —— 同簇非端点相交（1.5m 业务判交，严格，不随生成侧豁免放宽）、
//      控制多边形互穿前兆、端点 G1、固定形态几何偏差、生成器自报状态、单路口耗时。
//
// 用法：
//   diag_all_violations [数据目录] [--csv 输出.csv] [--only 子串] [--quiet]
//                       [--road-edge-clearance 米]
#include "constraints/cluster_order.h"
#include "constraints/constraint_evaluator.h"
#include "constraints/fence_check.h"
#include "curve/curve_utils.h"
#include "domain/scene_context.h"
#include "generation/connectivity_generation_context.h"
#include "generation/uturn_shape.h"
#include "intersection_shape_generator.h"
#include "io/iodata_json.h"
#include "optimizer/sdf_field.h"
#include "preprocessing/uturn_family_builder.h"
#include "toolkits/toolkits.h"
#include "utils.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <dirent.h>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

using namespace isg;

namespace {

// ── 统一发现记录 ────────────────────────────────────────────────────────────
struct Finding {
    std::string category;  ///< 检查类别 id（Excel 按此分 sheet）
    std::string check;     ///< 具体规则 id（如 shape.uturn.leads）
    std::string dataset;   ///< 数据文件名（按文件统计的键）
    std::string subject;   ///< 审计对象：连接 id 或 "a|b" 连接对
    std::string severity;  ///< error / warn / info
    std::string metric;    ///< 数值指标名（可空）
    std::string detail;    ///< 说明（含生成侧给出的 reason 与坐标）
    double value = 0.0;
    double threshold = 0.0;
    bool numeric = false;
};

std::vector<Finding> g_findings;

void emit(const std::string& category, const std::string& check,
          const std::string& dataset, const std::string& subject,
          const char* severity, const std::string& detail,
          const char* metric = "", double value = 0.0, double threshold = 0.0,
          bool numeric = false) {
    Finding f;
    f.category = category;
    f.check = check;
    f.dataset = dataset;
    f.subject = subject;
    f.severity = severity;
    f.metric = metric;
    f.detail = detail;
    f.value = value;
    f.threshold = threshold;
    f.numeric = numeric;
    g_findings.push_back(f);
}

void emitNum(const std::string& category, const std::string& check,
             const std::string& dataset, const std::string& subject,
             const char* severity, const std::string& detail,
             const char* metric, double value, double threshold) {
    emit(category, check, dataset, subject, severity, detail, metric, value,
         threshold, true);
}

bool startsWith(const std::string& s, const char* prefix) {
    const size_t n = std::strlen(prefix);
    return s.size() >= n && s.compare(0, n, prefix) == 0;
}

/// 规则 id → 检查类别 id。类别是修复动作的组织单位，故按“同一处代码/同一条产品
/// 规则”归并，而不是按 severity 归并。
std::string categoryOf(const std::string& check) {
    if (startsWith(check, "cluster.")) return "cluster_cross";
    if (startsWith(check, "precursor.")) return "cluster_precursor";
    if (startsWith(check, "fence.")) return "fence";
    if (startsWith(check, "shape.crosswalk")) return "crosswalk";
    if (startsWith(check, "shape.uturn")) return "uturn_shape";
    if (startsWith(check, "shape.ordinary")) return "ordinary_shape";
    if (startsWith(check, "shape.curvature")) return "curvature";
    if (startsWith(check, "shape.g1") || startsWith(check, "g1.")) return "g1";
    if (startsWith(check, "shape.single_segment")) return "ordinary_shape";
    if (startsWith(check, "fixed_shape.")) return "fixed_shape";
    if (check == "physical.obstacle") return "obstacle";
    if (check == "physical.boundary") return "boundary";
    if (check == "physical.road_edge_clearance") return "road_edge";
    if (check == "geometry.self_intersection") return "self_intersect";
    if (startsWith(check, "geometry.")) return "generation";
    if (startsWith(check, "generator.")) return "generator_report";
    if (startsWith(check, "topology.")) return "topology";
    if (startsWith(check, "perf.")) return "perf";
    return "other";
}

std::string csvEscape(const std::string& s) {
    if (s.find_first_of(",\"\n\r") == std::string::npos) return s;
    std::string out = "\"";
    for (size_t i = 0; i < s.size(); ++i) {
        if (s[i] == '"') out += "\"\"";
        else out += s[i];
    }
    out += "\"";
    return out;
}

bool writeFindingsCsv(const std::string& path) {
    FILE* f = std::fopen(path.c_str(), "w");
    if (!f) return false;
    std::fprintf(f, "source,category,check,dataset,subject,severity,metric,"
                    "value,threshold,detail,location\n");
    for (size_t i = 0; i < g_findings.size(); ++i) {
        const Finding& fd = g_findings[i];
        char value[64] = "";
        char threshold[64] = "";
        if (fd.numeric) {
            std::snprintf(value, sizeof(value), "%.4f", fd.value);
            std::snprintf(threshold, sizeof(threshold), "%.4f", fd.threshold);
        }
        std::fprintf(f, "audit,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n",
                     csvEscape(fd.category).c_str(), csvEscape(fd.check).c_str(),
                     csvEscape(fd.dataset).c_str(), csvEscape(fd.subject).c_str(),
                     csvEscape(fd.severity).c_str(), csvEscape(fd.metric).c_str(),
                     value, threshold, csvEscape(fd.detail).c_str(),
                     "tests/diag_all_violations.cpp");
    }
    std::fclose(f);
    return true;
}

// ── 几何/口径辅助 ──────────────────────────────────────────────────────────
/// 掉头端点走廊宽度：进入端点与退出端点在掉头轴法向上的间距。
double uturnEndpointCorridor(const BezierCurve& u) {
    if (u.empty()) return 1e18;
    const Vec2d t0 = u.startTan(), t1 = u.endTan();
    if (t0.norm() < 1e-8 || t1.norm() < 1e-8) return 1e18;
    Vec2d axis = t0.normalized() - t1.normalized();
    if (axis.norm() < 1e-8) return 1e18;
    axis = axis.normalized();
    const Vec2d lateral(-axis.y(), axis.x());
    return std::abs((u.endPt() - u.startPt()).dot(lateral));
}

/// 审计用约束档位：与生成侧同源，但不打开任何豁免。
ConstraintProfile auditProfile(int mode, double road_edge_clearance_override) {
    ConstraintProfile profile;
    profile.require_single_segment = false;  // 掉头合法地取三段式
    profile.enforce_fence = true;
    profile.check_self_intersection = true;
    profile.check_obstacle = true;
    profile.check_boundary = true;
    profile.check_g1 = true;        // 段间 G1（默认 15°）
    profile.check_curvature = false;  // 曲率上限随横向间隙/弦长变化，统一上限会误报；
                                      // 掉头曲率由 shape.uturn.* 覆盖
    profile.check_ordinary_shape = true;
    profile.check_uturn_shape = true;
    profile.check_crosswalk = true;
    profile.check_cluster = false;  // 同簇交由 1.5m 严格口径的成对审计负责，避免与
                                    // cluster_endpoint_tolerance=1e-4 的口径重复计数
    profile.require_fixed_shape_preservation = false;  // 依赖候选侧标记，输出侧改用几何偏差
    profile.allow_structural_cross = false;
    profile.allow_obstacle_cross_exemption = false;
    // 绕障清距已经烘进 SDF 的 buffer（与生成侧同参构建），此处再叠加会双计。
    profile.obstacle_clearance = 0.0;
    profile.road_edge_clearance = road_edge_clearance_override >= 0.0
        ? road_edge_clearance_override
        : roadEdgeAvoidanceClearanceForMode(mode);
    profile.samples = 64;
    return profile;
}

/// 与 IntersectionShapeGenerator::generate 同源的 SDF 取景框。
BoundingBox2d outerBoundingBox(const IntersectionInput& input) {
    BoundingBox2d roi;
    if (!input.area.geometry.outer.empty()) {
        roi = input.area.geometry.bbox();
    } else {
        for (size_t i = 0; i < input.lanes.size(); ++i)
            for (size_t j = 0; j < input.lanes[i].geometry.points.size(); ++j)
                roi.expand(input.lanes[i].geometry.points[j]);
        roi.min_pt -= Vec2d(20, 20);
        roi.max_pt += Vec2d(20, 20);
    }
    return roi;
}

/// 端点 G1：曲线端点切向与车道端点切向的最小点积。几何掉头按业务口径豁免——
/// 家族错开站位会有意打破端点 G1（与回归套件 endpointG1Min 同口径）。
double endpointG1Min(const BezierCurve& curve, const IntersectionInput& input,
                     const ConnectivityCurve& cc) {
    if (curve.empty()) return 1.0;
    const std::pair<Vec2d, Vec2d> entry = input.entryPtDir(cc.entry_lane_id);
    const std::pair<Vec2d, Vec2d> exit = input.exitPtDir(cc.exit_lane_id);
    const Vec2d t0 = entry.second.norm() > 1e-8 ? entry.second.normalized() : Vec2d(1, 0);
    const Vec2d t1 = exit.second.norm() > 1e-8 ? exit.second.normalized() : Vec2d(1, 0);
    const Vec2d st = curve.startTan().norm() > 1e-8 ? curve.startTan().normalized() : Vec2d(1, 0);
    const Vec2d et = curve.endTan().norm() > 1e-8 ? curve.endTan().normalized() : Vec2d(1, 0);
    return std::min(st.dot(t0), et.dot(t1));
}

/// 采样曲线到参考折线的最大偏差，用于固定形态保持审计。
double maxDeviationToPolyline(const BezierCurve& curve, const LineString2d& line) {
    const std::vector<Vec2d> ref = toVec2dArray(line.points);
    if (ref.size() < 2 || curve.empty()) return 0.0;
    const std::vector<Vec2d> pts = curve.sampleByArcLength(96);
    double worst = 0.0;
    for (size_t i = 0; i < pts.size(); ++i) {
        double best = 1e18;
        for (size_t j = 0; j + 1 < ref.size(); ++j)
            best = std::min(best, pointToSegment(pts[i], ref[j], ref[j + 1]).first);
        worst = std::max(worst, best);
    }
    return worst;
}

/// 几何掉头判定：与 shape_constraint 内的 isGeometricUTurn 同口径。
bool isGeometricUTurn(const std::pair<Vec2d, Vec2d>& entry,
                      const std::pair<Vec2d, Vec2d>& exit) {
    return entry.second.norm() > 1e-8 && exit.second.norm() > 1e-8 &&
           entry.second.normalized().dot(exit.second.normalized()) < -0.5;
}

std::string formatPoint(const Vec2d& pt) {
    char buf[64];
    std::snprintf(buf, sizeof(buf), "(%.3f, %.3f)", pt.x(), pt.y());
    return buf;
}

const char* severityOf(const ConstraintResult& result) {
    if (result.state == ConstraintState::Exempt) return "info";
    return result.severity == ConstraintSeverity::Hard ? "error" : "warn";
}

// ── 运行选项与每份数据集的汇总 ─────────────────────────────────────────────
struct Options {
    std::string dir;
    std::string csv;
    std::string only;
    bool quiet = false;
    double road_edge_clearance = -1.0;  ///< <0 表示按 mode 规则取值
};

struct DatasetSummary {
    std::string name;
    bool generated = false;
    double gen_ms = 0.0;
    int violations = 0;
    int precursors = 0;
    int constrained_pairs = 0;
    int exempt_pairs = 0;
    int errors = 0;
    int warns = 0;
};

/// 同簇非端点相交与控制多边形前兆。审计口径始终严格：结构性汇入豁免只作用于生成侧
/// 三段式候选搜索的兜底档（见 searchSegmentedUTurnTwoPass），输出审计不放宽，否则
/// 无法用同一把尺子对比历次改动。命中“恰好一条掉头 + 该掉头端点走廊退化”的违约会
/// 额外记录 merge_funnel 半径，便于把结构性汇入与真违约分开看。
void auditClusterPairs(const std::string& name, const IntersectionInput& input,
                       const IntersectionOutput& output, bool quiet,
                       DatasetSummary& summary) {
    std::unordered_map<ConnId, const ConnectivityCurve*> cmap;
    for (size_t i = 0; i < output.connectivity_curves.size(); ++i)
        cmap[output.connectivity_curves[i].id] = &output.connectivity_curves[i];

    ClusterOrderSolver solver;
    solver.build(input.connectivities, input.lanes, input.lane_groups, input.crosswalks);

    const std::vector<CurvePair>& pairs = solver.pairs();
    for (size_t i = 0; i < pairs.size(); ++i) {
        const CurvePair& p = pairs[i];
        if (p.exempt == CrossExemption::StructuralCross) {
            ++summary.exempt_pairs;
            continue;
        }
        ++summary.constrained_pairs;
        std::unordered_map<ConnId, const ConnectivityCurve*>::const_iterator ia =
            cmap.find(p.id_a);
        std::unordered_map<ConnId, const ConnectivityCurve*>::const_iterator ib =
            cmap.find(p.id_b);
        if (ia == cmap.end() || ib == cmap.end()) continue;
        if (!ia->second->curve || !ib->second->curve) continue;
        const BezierCurve& ca = *ia->second->curve;
        const BezierCurve& cb = *ib->second->curve;
        const std::string subject = p.id_a + "|" + p.id_b;
        if (curvesIntersectBusiness(ca, cb, 1.5)) {
            ++summary.violations;
            const bool a_ut = curveLooksUTurnForClusterExemption(ca);
            const bool b_ut = curveLooksUTurnForClusterExemption(cb);
            const bool degenerate = (a_ut != b_ut) &&
                uturnEndpointCorridor(a_ut ? ca : cb) < 0.60;
            const double funnel =
                degenerate ? sharedEndpointMergeFunnelRadius(ca, cb) : 0.0;
            char detail[256];
            std::snprintf(detail, sizeof(detail),
                          "同簇非端点相交 exempt=%d shared_ep=%d merge_funnel=%.3f",
                          (int)p.exempt, (int)p.shared_endpoint, funnel);
            emitNum("cluster_cross", "cluster.intersection", name, subject, "error",
                    detail, "merge_funnel", funnel, 0.0);
            if (!quiet)
                std::printf("  VIOLATION: %s <-> %s (exempt=%d, shared_ep=%d,"
                            " merge_funnel=%.3f)\n",
                            p.id_a.c_str(), p.id_b.c_str(), (int)p.exempt,
                            (int)p.shared_endpoint, funnel);
        } else if (p.shared_endpoint &&
                   sharedEndpointControlPolylinesCross(ca, cb, 0.30)) {
            // 控制多边形折线互穿是中段相交的几何前兆：采样判定尚未报违约，
            // 但控制多边形已不嵌套，属于需要关注的隐患。
            ++summary.precursors;
            emit("cluster_precursor", "precursor.control_polyline", name, subject,
                 "warn", "共享端点控制多边形折线互穿（中段相交前兆）");
            if (!quiet)
                std::printf("  PRECURSOR: %s <-> %s control polylines cross\n",
                            p.id_a.c_str(), p.id_b.c_str());
        }
    }
}

/// 输入固有形态是否穿障：穿障的固有形态生成侧允许重新生成，不能按“未保持固定
/// 形态”记违约（与 ConnectivityGenerationSession 的 fixedGeometryHitsObstacle 同义）。
bool polylineHitsObstacles(const LineString2d& line, const SDFField& sdf) {
    if (!sdf.valid()) return false;
    const std::vector<Vec2d> pts = toVec2dArray(line.points);
    for (size_t i = 0; i + 1 < pts.size(); ++i)
        for (int k = 0; k <= 8; ++k) {
            const Vec2d p = pts[i] + (pts[i + 1] - pts[i]) * (k / 8.0);
            if (sdf.queryWithGrad(p).first < 0.0) return true;
        }
    return false;
}

/// 逐条输出曲线跑一遍与生成侧同源的约束评估器，再补三项审计专有检查：
/// 端点 G1、固定形态几何偏差、生成器自报状态。
void auditCurves(const std::string& name, const IntersectionInput& raw_input,
                 const IntersectionOutput& output, const Options& opts) {
    // 与生成侧同源的规范化：先补全组引用，再统一同组端点切向。缺这一步，端点切向
    // 与生成期不一致，形态与 G1 类判定会整体误报。
    IntersectionInput input = InputNormalizer(raw_input);
    const ConnectivityDirectionConfig dir_cfg;  // 与生成器默认档一致
    ConnectivityDirectionNormalizer(input, dir_cfg);

    SDFField sdf;
    if (!input.obstacles.empty())
        sdf.build(outerBoundingBox(input), input.obstacles, 0.2,
                  obstacleAvoidanceClearanceForMode(input.mode));

    ClusterOrderSolver solver;
    solver.build(input.connectivities, input.lanes, input.lane_groups, input.crosswalks);
    const ClusterTopology topology = solver.topology();

    SceneContext scene(input);
    scene.bindSdf(&sdf);
    scene.bindClusterTopology(&topology);

    std::unordered_map<ConnId, const ConnectivityCurve*> cmap;
    GenerationState state;
    for (size_t i = 0; i < output.connectivity_curves.size(); ++i) {
        const ConnectivityCurve& cc = output.connectivity_curves[i];
        cmap[cc.id] = &cc;
        if (cc.curve) state.accepted_curves[cc.id] = *cc.curve;
    }

    const UTurnFamilyBuilder family_builder;
    const ConstraintEvaluator evaluator;
    const ConstraintProfile profile =
        auditProfile(input.mode, opts.road_edge_clearance);

    for (size_t i = 0; i < input.connectivities.size(); ++i) {
        const Connectivity& conn = input.connectivities[i];
        std::unordered_map<ConnId, const ConnectivityCurve*>::const_iterator it =
            cmap.find(conn.id);
        if (it == cmap.end() || !it->second->curve || it->second->curve->empty()) {
            emit("generation", "geometry.empty", name, conn.id, "error",
                 "连接没有生成曲线");
            if (!opts.quiet)
                std::printf("  ERROR   [generation] %s: 连接没有生成曲线\n",
                            conn.id.c_str());
            continue;
        }
        const ConnectivityCurve& cc = *it->second;
        const BezierCurve& curve = *cc.curve;

        const UTurnFamilyInfo family = family_builder.build(
            conn, input, dir_cfg.uturn_alignment_scope, &solver, false);
        CurveGenerationContext context = CurveGenerationContextBuilder().build(
            scene, conn, dir_cfg.uturn_alignment_scope, &solver, &family);
        context.profile = profile;

        const ConstraintReport report = evaluator.evaluate(curve, context, state);
        for (size_t r = 0; r < report.results.size(); ++r) {
            const ConstraintResult& result = report.results[r];
            if (result.state == ConstraintState::Satisfied ||
                result.state == ConstraintState::NotApplicable)
                continue;
            std::string detail = result.reason;
            if (!result.locations.empty())
                detail += " @ " + formatPoint(result.locations.front());
            const char* severity = severityOf(result);
            emitNum(categoryOf(result.id), result.id, name, cc.id, severity, detail,
                    "violation", result.violation, 0.0);
            if (!opts.quiet)
                std::printf("  %-7s [%s] %s: %s\n",
                            std::strcmp(severity, "error") == 0 ? "ERROR" :
                            (std::strcmp(severity, "warn") == 0 ? "WARN" : "INFO"),
                            categoryOf(result.id).c_str(), cc.id.c_str(),
                            detail.c_str());
        }

        // 端点 G1（业务口径）：几何掉头豁免，家族错开站位会有意打破端点 G1。
        if (!isGeometricUTurn(context.entry, context.exit)) {
            const double g1 = endpointG1Min(curve, input, cc);
            if (g1 < 0.999) {
                const bool hard = g1 < 0.99;
                emitNum("g1", "g1.endpoint", name, cc.id, hard ? "error" : "warn",
                        "端点切向与车道端点切向不一致", "cos", g1, 0.99);
                if (!opts.quiet)
                    std::printf("  %-7s [g1] %s: 端点 G1 cos=%.5f\n",
                                hard ? "ERROR" : "WARN", cc.id.c_str(), g1);
            }
        }

        // 固定形态保持：几何掉头按最新掉头形态约束统一重生成，不参与此项；固有形态
        // 穿障时生成侧允许重生成，只记 info。
        if (conn.geometry.points.size() >= 2 &&
            !isGeometricUTurn(context.entry, context.exit)) {
            const double deviation = maxDeviationToPolyline(curve, conn.geometry);
            const double tolerance = 0.10;
            if (deviation > tolerance) {
                const bool obstacle_regen = polylineHitsObstacles(conn.geometry, sdf);
                const char* severity =
                    obstacle_regen ? "info" : (conn.fixed_shape ? "error" : "warn");
                const std::string detail = obstacle_regen
                    ? "固有形态穿障，生成侧允许重新生成"
                    : (conn.fixed_shape ? "声明 fixed_shape 的连接偏离输入固有形态"
                                        : "带固有形态的连接偏离输入几何（未声明 fixed_shape）");
                emitNum("fixed_shape", "fixed_shape.geometry_deviation", name, cc.id,
                        severity, detail, "max_deviation", deviation, tolerance);
                if (!opts.quiet)
                    std::printf("  %-7s [fixed_shape] %s: 偏离固有形态 %.3fm (%s)\n",
                                std::strcmp(severity, "error") == 0 ? "ERROR" :
                                (std::strcmp(severity, "warn") == 0 ? "WARN" : "INFO"),
                                cc.id.c_str(), deviation, detail.c_str());
            }
        }

        // 生成器自报状态：与独立复算互为交叉验证，二者不一致本身就是线索。
        // 分级依据「是否存在具体违规原因」：Degraded 且带 reason 说明确有规则被破坏；
        // Degraded 但 reason 为空、仅记录了豁免交叉，属设计允许的结果，记为 info，
        // 否则统计表会被这类记录淹没而掩盖真正的违规。
        if (cc.status != CurveStatus::OK) {
            const char* status_name = cc.status == CurveStatus::Infeasible ? "Infeasible"
                : (cc.status == CurveStatus::Degraded ? "Degraded" : "WarnA2");
            const bool exempt_only =
                cc.violation.reason.empty() && !cc.violation.exempt_crosses.empty();
            const char* severity =
                cc.status == CurveStatus::Infeasible ? "error"
                : (!cc.violation.reason.empty() ? "error" : (exempt_only ? "info" : "warn"));
            emitNum("generator_report", "generator.status", name, cc.id, severity,
                    std::string("生成器自报 status=") + status_name +
                        (cc.violation.reason.empty()
                             ? (exempt_only ? "（仅记录豁免交叉，无具体违规原因）" : "（无具体违规原因）")
                             : "，reason=" + cc.violation.reason),
                    "exempt_crosses", (double)cc.violation.exempt_crosses.size(), 0.0);
        }
        if (cc.violation.max_obstacle_penetration > 1e-3)
            emitNum("generator_report", "generator.obstacle_penetration", name, cc.id,
                    "error", "生成器自报障碍穿透", "penetration",
                    cc.violation.max_obstacle_penetration, 0.0);
        if (cc.violation.max_fence_overflow > 1e-3) {
            // 与 fence.containment 同口径：面本身不含两连接点之间的直线弦时，外溢是
            // 输入面画小了而非曲线画歪了，降为 info 并把弦外溢一并报出以便定标。
            const double chord = cc.curve
                ? fenceChordOverflow(input.area.geometry, cc.curve->startPt(),
                                     cc.curve->endPt(), 160)
                : 0.0;
            const bool forced =
                fenceOverflowForcedByFace(cc.violation.max_fence_overflow, chord);
            char detail[192];
            std::snprintf(detail, sizeof(detail),
                          "生成器自报围栏外溢（弦外溢=%.4fm，%s）", chord,
                          forced ? "面不含该弦，几何强制" : "超出弦强制下限");
            emitNum("generator_report", "generator.fence_overflow", name, cc.id,
                    forced ? "info" : "error", detail, "overflow",
                    cc.violation.max_fence_overflow, chord);
        }
        if (cc.violation.type != ViolationInfo::InfeasibilityType::None)
            emit("generator_report", "generator.infeasibility", name, cc.id, "error",
                 std::string("生成器自报不可行类型=") +
                     std::to_string((int)cc.violation.type) +
                     (cc.violation.reason.empty() ? "" : "，reason=" + cc.violation.reason));
    }
}

/// 审计一份数据集：生成 + 拓扑报告 + 耗时 + 同簇成对审计 + 逐曲线约束审计。
DatasetSummary auditDataset(const std::string& path, const std::string& name,
                            const Options& opts) {
    DatasetSummary summary;
    summary.name = name;
    const size_t begin = g_findings.size();

    IntersectionInput input = IntersectionIO::loadFromFile(path);
    if (input.connectivities.empty()) {
        emit("generation", "geometry.no_connectivity", name, "-", "info",
             "数据集无连通关系，跳过");
        return summary;
    }

    IntersectionShapeGenerator gen;
    IntersectionOutput output;
    const std::chrono::steady_clock::time_point t0 = std::chrono::steady_clock::now();
    const bool ok = gen.generate(input, output);
    summary.gen_ms = std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - t0).count();
    summary.generated = ok;

    const ValidationReport& validation = gen.lastReport();
    for (size_t i = 0; i < validation.errors.size(); ++i)
        emit("topology", "topology.error", name, "-", "error", validation.errors[i]);
    for (size_t i = 0; i < validation.warnings.size(); ++i)
        emit("topology", "topology.warning", name, "-", "warn", validation.warnings[i]);

    if (!ok) {
        emit("generation", "geometry.generation_failed", name, "-", "error",
             "generate() 返回 false");
        if (!opts.quiet) std::fprintf(stderr, "  GENERATION FAILED\n");
    } else {
        // 单路口生成耗时：15s 是产品上限，始终记录一行以便逐版本对比。
        const double limit = 15000.0;
        emitNum("perf", "perf.generation_ms", name, "-",
                summary.gen_ms > limit ? "error" : "info", "单路口生成耗时（含审计外的生成全流程）",
                "ms", summary.gen_ms, limit);
        auditClusterPairs(name, input, output, opts.quiet, summary);
        auditCurves(name, input, output, opts);
    }

    for (size_t i = begin; i < g_findings.size(); ++i) {
        if (g_findings[i].severity == "error") ++summary.errors;
        else if (g_findings[i].severity == "warn") ++summary.warns;
    }
    if (!opts.quiet)
        std::printf("  constrained_pairs=%d, exempt_pairs=%d, violations=%d,"
                    " precursors=%d, errors=%d, warns=%d, gen=%.0fms\n",
                    summary.constrained_pairs, summary.exempt_pairs,
                    summary.violations, summary.precursors, summary.errors,
                    summary.warns, summary.gen_ms);
    return summary;
}

/// 检查类别的中文名，仅用于 stdout 表格可读性；Excel 侧的类别名由 tools/test_report.py 统一维护。
const char* categoryLabel(const std::string& id) {
    if (id == "cluster_cross") return "同簇相交";
    if (id == "cluster_precursor") return "同簇相交前兆";
    if (id == "fence") return "围栏越界";
    if (id == "crosswalk") return "人行横道";
    if (id == "uturn_shape") return "调头形态";
    if (id == "ordinary_shape") return "常规形态";
    if (id == "curvature") return "曲率";
    if (id == "g1") return "G1连续";
    if (id == "fixed_shape") return "固有形态保持";
    if (id == "obstacle") return "障碍避让";
    if (id == "boundary") return "边界避让";
    if (id == "road_edge") return "路缘避让";
    if (id == "self_intersect") return "自交";
    if (id == "generation") return "生成失败";
    if (id == "generator_report") return "生成器自报";
    if (id == "topology") return "拓扑校验";
    if (id == "perf") return "性能";
    return "其他";
}

/// 统计一格：error / warn / info 三档计数。
struct Tally {
    int error = 0, warn = 0, info = 0;
    void add(const std::string& severity) {
        if (severity == "error") ++error;
        else if (severity == "warn") ++warn;
        else ++info;
    }
    int total() const { return error + warn + info; }
};

/// 打印“按检查类别”“按文件”两张错误统计表，便于按类别/按文件组织修复。
void printStatistics() {
    std::map<std::string, Tally> by_category;
    std::map<std::string, Tally> by_dataset;
    for (size_t i = 0; i < g_findings.size(); ++i) {
        by_category[g_findings[i].category].add(g_findings[i].severity);
        by_dataset[g_findings[i].dataset].add(g_findings[i].severity);
    }

    std::printf("\n==== 按检查类别统计 ====\n");
    std::printf("%-20s %-16s %8s %8s %8s %8s\n", "category", "类别", "error",
                "warn", "info", "total");
    for (std::map<std::string, Tally>::const_iterator it = by_category.begin();
         it != by_category.end(); ++it)
        std::printf("%-20s %-16s %8d %8d %8d %8d\n", it->first.c_str(),
                    categoryLabel(it->first), it->second.error, it->second.warn,
                    it->second.info, it->second.total());

    std::printf("\n==== 按文件统计 ====\n");
    std::printf("%-40s %8s %8s %8s %8s\n", "dataset", "error", "warn", "info", "total");
    for (std::map<std::string, Tally>::const_iterator it = by_dataset.begin();
         it != by_dataset.end(); ++it)
        std::printf("%-40s %8d %8d %8d %8d\n", it->first.c_str(), it->second.error,
                    it->second.warn, it->second.info, it->second.total());
}

/// 解析命令行；无法识别的参数按“数据目录”对待（保持历史用法 diag_all_violations datas）。
bool parseOptions(int argc, char** argv, Options& opts) {
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--csv" && i + 1 < argc) opts.csv = argv[++i];
        else if (a == "--only" && i + 1 < argc) opts.only = argv[++i];
        else if (a == "--quiet") opts.quiet = true;
        else if (a == "--road-edge-clearance" && i + 1 < argc)
            opts.road_edge_clearance = std::atof(argv[++i]);
        else if (startsWith(a, "--")) {
            std::fprintf(stderr, "未知参数: %s\n", a.c_str());
            return false;
        } else opts.dir = a;
    }
    if (opts.dir.empty()) opts.dir = std::string(PROJECT_ROOT_DIR) + "/datas";
    return true;
}

}  // namespace

int main(int argc, char** argv) {
    Options opts;
    if (!parseOptions(argc, argv, opts)) return 2;

    DIR* d = opendir(opts.dir.c_str());
    if (!d) {
        std::fprintf(stderr, "Cannot open %s\n", opts.dir.c_str());
        return 1;
    }
    std::vector<std::string> files;
    struct dirent* ent = nullptr;
    while ((ent = readdir(d)) != nullptr) {
        const std::string name = ent->d_name;
        if (name.size() < 6) continue;
        if (name.substr(name.size() - 5) != ".json") continue;
        if (!opts.only.empty() && name.find(opts.only) == std::string::npos) continue;
        files.push_back(name);
    }
    closedir(d);
    std::sort(files.begin(), files.end());  // 固定顺序，便于逐版本 diff

    int total_violations = 0, total_precursors = 0;
    for (size_t i = 0; i < files.size(); ++i) {
        if (!opts.quiet) std::printf("\n==== %s ====\n", files[i].c_str());
        const DatasetSummary s =
            auditDataset(opts.dir + "/" + files[i], files[i], opts);
        total_violations += s.violations;
        total_precursors += s.precursors;
    }
    std::printf("\n==== TOTAL VIOLATIONS: %d ====\n", total_violations);
    std::printf("==== TOTAL PRECURSORS: %d ====\n", total_precursors);
    printStatistics();

    if (!opts.csv.empty()) {
        if (!writeFindingsCsv(opts.csv)) {
            std::fprintf(stderr, "写入 CSV 失败: %s\n", opts.csv.c_str());
            return 1;
        }
        std::printf("\nfindings CSV -> %s (%zu 条)\n", opts.csv.c_str(),
                    g_findings.size());
    }
    return total_violations > 0 ? 1 : 0;
}
