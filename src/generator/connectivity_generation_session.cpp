#include "generator/connectivity_generation_session.h"
#include "curve/hermite_init.h"
#include "curve/curve_utils.h"
#include "constraints/boundary_safety.h"
#include "constraints/fence_check.h"
#include "constraints/road_edge_clearance.h"
#include "utils.h"
#include "optimizer/sdf_field.h"
#include "constraints/infeasibility_detector.h"
#include "constraints/uturn_envelope_constraint.h"
#include "ordering/generation_planner.h"
#include "initialization/fixed_shape_initializer.h"
#include "initialization/ordinary_curve_initializer.h"
#include "initialization/uturn_curve_initializer.h"
#include "initialization/avoidance_candidate_generator.h"
#include "initialization/curve_initializer_registry.h"
#include "generation/candidate_generation.h"
#include "generation/curve_sampling.h"
#include "generation/uturn_multi_constraint_search.h"
#include "generation/uturn_shape.h"
#include "generation/segmented_uturn_candidate_search.h"
#include "generation/curve_processing.h"
#include "generation/repair_and_audit.h"
#include "toolkits/toolkits.h"
#include "preprocessing/uturn_family_builder.h"
#include <chrono>
#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <cstdlib>
#include <cstdio>
#include <functional>
#include <queue>
#include <string>

namespace isg {

// U型调头硬形态约束：进入/退出两侧都没有人行横道约束时，首尾直行段才使用2m保底距离进入中间单段掉头弧。
static constexpr double kUTurnNoCrosswalkMinLead = 2.0;
// 同配置家族U-turn的首/尾平齐点沿轴向法线相向微移，避免直行段完全重叠/贴合。
// 该量必须与需求中的家族分档一致：按短径到长径逆序 N，沿两个
// 平齐点连线相向移动 N*0.01m。不能使用分米级偏移再由单条曲线
// 独立截断，否则会把短径中弧挤成尖弧并破坏家族嵌套关系。
static constexpr double kUTurnSharedEndpointLeadStagger = 0.01;
static constexpr double kUTurnSharedEndpointPairLeadStagger = 0.01;
// 同簇曲线只允许在真实连接点相遇；连接点外的贴合/重叠/相交都按违规处理。
static constexpr double kClusterEndpointTol = 0.30;
// 多段严格避让回退的同簇端点容差。该口径只豁免真实连接点的浮点误差，
// 防止两段候选在 0.15m 之外发生短距离真实相交仍被 1.5m 近端窗口掩盖。
static constexpr double kStrictSameClusterEndpointTol = 0.15;
// "中段互穿"判据的端点容差：与最终同簇审计一致，用于把连接点附近不可
// 避免的汇聚/贴合与真正远离端点的互穿区分开。仅用于增量守卫的复核。
static constexpr double kMidSpanCrossEndpointTol = 1.5;
// 共享连接点的同簇配对在配对修复里的相交口径。
//
// 两条曲线从同一连接点以同一切向扇出时，端点邻域的横向间距按
// `lat(s) ≈ Δκ0·s²/2` 增长：即使两条曲线的起点曲率差已足以让它们在
// 10m 外分开 0.3m 以上，在离端点 1m 处的间距也只有毫米级，远低于采样
// 间距，采样折线必然在那里互相穿插。若配对修复继续用 0.15m 的端点
// 容差，这些**纯数值**的近端点穿插会被记成真实穿越，把审计口径下干净
// 的候选整批剔除（110003285 的 81 只有两个审计干净的单段形态，二者都
// 因此被 0.15m 判为"仍与 80 相交"而无法录取，`80|81` 因此长期无解）。
// 因此共享连接点的配对统一采用与最终审计一致的 1.5m 端点容差；
// 连接点邻域的"不得贴合"由 `curvesHaveForbiddenSameClusterIntersection`
// 的贴合项和文档 §6.5 的连接点容差规则负责，不靠相交计数表达。
static constexpr double kSharedEndpointPairCrossTol = kMidSpanCrossEndpointTol;
// 共享连接点上"近直行 × 近直行"配对是否计入配对修复的同簇交叉账本。
//
// 配对修复本身只处理"同侧转向 × 同侧转向"和"直行 × 转向"两类对，近直行
// 兄弟之间从来不是修复目标；但它们同样受同簇非端点不相交约束（最终审计
// 只豁免 StructuralCross）。若账本也把它们排除在外，为躲开同族转向而把
// 一条近直行压成极端非对称形态（h0 由 23 压到 3.7）时，它与其它近直行
// 兄弟的转入率次序反转对本阶段完全不可见：100000012-nu 的 16 为躲
// 43285424 把 κ0 抬到 0.0590，越过 14 的 0.0554，凭空造出 14|16。
// 因此把这类对纳入计数（只计数，不新增修复目标），口径与审计一致。
static bool sharedEndpointStraightPairParticipates(double turn_a, double turn_b,
                                                  CrossExemption exempt) {
    return turn_a < 0.25 && turn_b < 0.25 &&
           exempt != CrossExemption::StructuralCross;
}
// 同入口/同出口三段式U-turn首尾直行段从同一点扇出时，平齐点之间至少
// 需要达到该间距，才能在最终同簇标注里只检查后续中弧/尾段。
static constexpr double kUTurnSharedLeadFanMinSeparation = 0.005;


////////////////////////////////////////////////////////////
// 调试辅助函数: 当环境变量ISG_DEBUG_UTURN设置时返回true。
////////////////////////////////////////////////////////////
// 返回是否启用 U-turn 初始生成和修复阶段的诊断输出。
// 结果在进程生命周期内缓存；仅检查环境变量是否存在，不解析其内容。
static bool isgDebugUTurn() {
    static bool v = (std::getenv("ISG_DEBUG_UTURN") != nullptr);
    return v;
}

// 返回是否启用 Boundary 候选搜索的详细诊断输出。
// 该开关只影响 stderr 日志，不改变候选生成、约束判定或最终曲线。
static bool isgDebugBoundaryRepair() {
    static bool v = (std::getenv("ISG_DEBUG_BOUNDARY_REPAIR") != nullptr);
    return v;
}

// 返回是否启用同簇成对修复的诊断输出。
// 环境变量存在即开启，并通过静态缓存避免重复访问进程环境。
static bool isgDebugPairRepair() {
    static bool v = (std::getenv("ISG_DEBUG_PAIR_REPAIR") != nullptr);
    return v;
}

// 返回是否启用阶段和单连接生成耗时统计。
// 开启后仅增加 profile 日志，不改变算法选择和结果。
static bool isgProfile() {
    static bool v = (std::getenv("ISG_PROFILE") != nullptr);
    return v;
}

// 共享族回溯会对同一候选执行大量成对相交检查。使用与
// curvesIntersectBusiness 相同的采样密度预采样，避免每次检查都重新构造
// 两条曲线的采样折线；最终输出审计仍调用原始高精度判定。
static SampledCurve sampleStrictPairCurve(const BezierCurve& curve) {
    SampledCurve sampled;
    sampled.start = curve.startPt();
    sampled.end = curve.endPt();
    const int count = std::max(
        48, static_cast<int>(std::ceil(curve.arcLength() / 0.20)) + 1);
    sampled.pts = curve.sampleByArcLength(std::min(count, 240));
    buildSampledCurveIndex(sampled);
    return sampled;
}

struct StrictFamilyProfileScope {
    const char* name;
    std::chrono::steady_clock::time_point started =
        std::chrono::steady_clock::now();

    explicit StrictFamilyProfileScope(const char* family_name)
        : name(family_name) {}

    ~StrictFamilyProfileScope() {
        if (!isgProfile())
            return;
        fprintf(stderr, "[ISG_PROFILE] strict family %s %.3f ms\n", name,
                std::chrono::duration<double, std::milli>(
                    std::chrono::steady_clock::now() - started).count());
    }
};


////////////////////////////////////////////////////////////
// 编译辅助函数
////////////////////////////////////////////////////////////
// 判别输入是否只包含连接、车道等几何数据，而不包含会触发物理约束的
// 障碍物、边界、人行横道或停车线。纯几何输入使用更轻量的拓扑收口路径。
static bool isPureGeometricInput(const IntersectionInput& input) {
    return input.obstacles.empty() && input.boundaries.empty() &&
           input.crosswalks.empty() && input.stop_lines.empty();
}

// ── 端点粘连 boundary 过滤（已移除） ───────────────────────────────────────
// 历史上这里有一个 filterEndpointAdherentBoundaries：只要 boundary 的某个端点
// 贴在曲线端点附近、且折线"首尾弦方向"与该端点侧车道切向对齐，就整条剔除。
// 该判据有两个缺陷，已由 boundary_safety.h 的逐段豁免（见
// endpointGrazeExemptSegmentMasks）取代：
//   1. 粒度过粗——整条折线被剔除，发夹形鼻端折回的另一条腿也不再判违；
//   2. 方向量取错——中央分隔带/渠化岛鼻端是两条近平行腿加一个急弯，首尾弦
//      几乎垂直于两条腿（实测 |dot T1| = 0.015 / 0.006），远低于 0.70 阈值，
//      于是真正需要豁免的场景反而命中不了。
// 新判据改用"从贴合端点沿折线走出的连续共线段区间"，并要求曲线相对该区间的
// 最大横向偏移不超过 graze_tol，只豁免厘米级的共享端点数值擦碰。

// 使用 BoundarySafety 的道路中心和采样规则检查曲线是否穿越边界或越出道路边缘。
// 空边界集合直接视为安全；返回值只表示物理边界风险，不包含障碍物或同簇风险。
// 端点贴合 RoadEdge 的共线擦碰由 boundary_safety.h 的豁免判据剔除：那类接触
// 只是曲线端点与折线端点重合造成的厘米级数值退化，不是穿越路缘。
static bool curveIntersectsBoundaries(
        const BezierCurve& curve, const IntersectionInput& input) {
    if (input.boundaries.empty())
        return false;
    auto safety = curveBoundarySafetyForInput(curve, input, 128, 0.15);
    return safety.intersects || safety.outside_road_edge;
}

// 对曲线进行采样后的线段级边界相交检查。
// 只允许交点落在曲线真实首点或尾点的浮点误差内；Boundary 端点和端点后的
// 共线贴合/重叠不属于连接点豁免。返回 true 表示存在一个非端点真实接触。
static bool curveRawIntersectsBoundariesImpl(
        const BezierCurve& curve, const std::vector<Boundary>& boundaries,
        bool road_edges_only, double curve_endpoint_tol, double boundary_endpoint_tol,
        const std::vector<std::vector<bool>>* exempt_masks = nullptr) {
    (void)boundary_endpoint_tol;
    (void)exempt_masks;
    if (curve.empty())
        return false;
    const double endpoint_tol = std::min(
        std::max(0.0, curve_endpoint_tol), kConnectionPointTolerance);
    BoundingBox2d curve_box = curve.bbox();
    double bbox_pad = endpoint_tol + 0.02;
    curve_box.min_pt -= Vec2d(bbox_pad, bbox_pad);
    curve_box.max_pt += Vec2d(bbox_pad, bbox_pad);
    auto pts = curve.sampleByArcLength(std::max(
            64, std::min(240, (int)std::ceil(curve.arcLength() / 0.18) + 1)));
    if (pts.size() < 2)
        return false;
    // 采样折线 × 边界折线是 O(Ns×Nb) 的双层循环：240 个采样点配一条 50 段的
    // 边界就有上万次 segmentHasForbiddenBoundaryContact，而绝大多数配对的包围盒
    // 根本不相交。这里先算出曲线分段包围盒，再对每条边界预算分段包围盒，
    // 用包围盒拒绝替代逐对求交。
    // segmentHasForbiddenBoundaryContact 的第一道判据就是 segmentsIntersect，
    // 包围盒不相交时它必然返回 false，因此剪枝不改变任何结论；1e-9 的余量
    // 对应 segmentsIntersect 共线重叠判据里的参数容差。
    const double kBoxSlack = 1e-9;
    const std::size_t nseg = pts.size() - 1;
    std::vector<BoundingBox2d> curve_seg_box(nseg);
    for (std::size_t i = 0; i < nseg; ++i) {
        curve_seg_box[i].expand(pts[i]);
        curve_seg_box[i].expand(pts[i + 1]);
    }
    const auto boxesApart = [kBoxSlack](const BoundingBox2d& x,
                                        const BoundingBox2d& y) {
        return x.max_pt[0] < y.min_pt[0] - kBoxSlack ||
               y.max_pt[0] < x.min_pt[0] - kBoxSlack ||
               x.max_pt[1] < y.min_pt[1] - kBoxSlack ||
               y.max_pt[1] < x.min_pt[1] - kBoxSlack;
    };
    std::vector<BoundingBox2d> bnd_seg_box;
    for (size_t bi = 0; bi < boundaries.size(); ++bi) {
        const auto& bnd = boundaries[bi];
        if ((road_edges_only && bnd.type != Boundary::Type::RoadEdge) ||
            bnd.geometry.points.size() < 2)
            continue;
        const BoundingBox2d bnd_box = bnd.geometry.bbox();
        if (!curve_box.intersects(bnd_box))
            continue;
        const auto& bpts = bnd.geometry.points;
        const std::size_t nbseg = bpts.size() - 1;
        bnd_seg_box.assign(nbseg, BoundingBox2d());
        for (std::size_t j = 0; j < nbseg; ++j) {
            bnd_seg_box[j].expand(bpts[j]);
            bnd_seg_box[j].expand(bpts[j + 1]);
        }
        for (std::size_t i = 0; i < nseg; ++i) {
            const BoundingBox2d& abox = curve_seg_box[i];
            if (boxesApart(abox, bnd_box))
                continue;
            for (std::size_t j = 0; j < nbseg; ++j) {
                if (boxesApart(abox, bnd_seg_box[j]))
                    continue;
                if (!segmentHasForbiddenBoundaryContact(
                        pts[i], pts[i + 1], bpts[j], bpts[j + 1],
                        pts.front(), pts.back(), endpoint_tol))
                    continue;
                // 共享连接点的"离开楔形"不判违（见 boundary_safety.h）。
                // 判据是全局的（要求所有接触都属于楔形），因此在第一处接触上
                // 求值一次即可定论；严格路径因此不付任何额外代价。
                if (curveBoundaryContactsAreDepartureWedges(
                        curve, boundaries, road_edges_only, curve_endpoint_tol))
                    return false;
                return true;
            }
        }
    }
    return false;
}

// 检查曲线是否与任意类型的边界发生非端点真实接触。
// 这是生成期与最终审计共用的严格包装器；参数中的历史宽松容差会在底层收敛到
// kConnectionPointTolerance，避免调用者重新打开 Boundary 端点或端点贴合豁免。
static bool curveRawIntersectsAnyBoundary(
        const BezierCurve& curve, const std::vector<Boundary>& boundaries,
        double curve_endpoint_tol, double boundary_endpoint_tol) {
    return curveRawIntersectsBoundariesImpl(
            curve, boundaries, false, curve_endpoint_tol, boundary_endpoint_tol);
}

// 三段式 U-turn 的首尾直行段可能位于短 RoadEdge 的中心反侧，但仍在
// 合法进入/退出车道内。这里仅判真实非端点相交；RoadEdge 的 1m 净距
// 由 roadEdgeClearanceViolation 独立检查，不能用短边延长线的 outside 判定替代真实穿越。
// 对三段式 U-turn 使用真实线段穿越判定，对其它曲线复用常规 BoundarySafety
// 规则；这样既允许合法首尾直行段，也不会放宽中间掉头弧的边界约束。
static bool curveIntersectsBoundariesUTurn(
        const BezierCurve& curve, const IntersectionInput& input) {
    if (curve.numSegments() != 3)
        return curveIntersectsBoundaries(curve, input);
    return curveRawIntersectsAnyBoundary(curve, input.boundaries, 0.15, 0.10);
}

////////////////////////////////////////////////////////////
// ConnectivityGenerationSession 成员函数
////////////////////////////////////////////////////////////
// 初始化单路口生成会话，保存 LBFGS 求解器配置和方向归一化配置，供后续
// 单曲线生成、物理修复、拓扑修复及最终审计共享。
ConnectivityGenerationSession::ConnectivityGenerationSession(
    const LBFGSConfig& config, const ConnectivityDirectionConfig& direction_config)
    : solver_(config), direction_cfg_(direction_config) {}

// ─────────────────────────────────────────────────────────────────────────────
//  buildSiblings
//
//  对每条已生成曲线:
//  - 判定簇归属(同进入组 OR 同退出组)
//  - 设置 exempt_a1: 仅当 enterGroup 与 exitGroup 均不匹配时为true
//    (= 不同道路arm的交叉交通 = 结构性交叉)
//  - 设置 expected_side: 来自 ClusterOrderSolver::expectedSideOf(cid, id)
//    (排序为降序, 故符号正确)
//  - 设置 ref_perp: 来自 ClusterOrderSolver::refPerpOf()
// ─────────────────────────────────────────────────────────────────────────────
std::vector<SiblingCurve> ConnectivityGenerationSession::buildSiblings(
    const ConnId& id, const std::unordered_map<ConnId,BezierCurve>& done,
    const ClusterOrderSolver& cs, const std::vector<Connectivity>& conns,
    bool constrained_only, const std::unordered_set<ConnId>* fixed_shape_ids) const {
    std::vector<SiblingCurve> sibs;
    for (auto& kv : done) {
        auto& cid = kv.first;
        auto& curve = kv.second;
        if (cid == id)
            continue;

        SiblingCurve s;
        s.id = cid;
        s.curve = curve;
        s.fixed_shape = fixed_shape_ids && fixed_shape_ids->count(cid) > 0;

        auto ex = cs.exemptionOf(id, cid);
        bool in_cluster = cs.pairExists(id, cid);
        if (!in_cluster) {
            // 无簇关系 → 不同arm → 结构性交叉交通
            s.exempt_a1 = true;
        } else if (ex == CrossExemption::StructuralCross) {
            // 拓扑倒置 → 必须交叉 → 豁免
            s.exempt_a1 = true;
        } else {
            // 同簇约束配对 → 施加惩罚
            s.exempt_a1 = false;
        }
        if (constrained_only && s.exempt_a1)
            continue;

        s.exempt_a2_radius = (ex == CrossExemption::ObstacleCross) ? 1.5 : 0.0;

        // expected_side 以 cid 视角相对 id:
        // cs.expectedSideOf(cid, id):
        //   +1 → 兄弟(cid)在当前(id)左侧
        //   -1 → 兄弟(cid)在当前(id)右侧
        s.expected_side = cs.expectedSideOf(cid, id);

        // evalCluster的固定横向参考轴
        s.ref_perp = cs.refPerpOf(id, cid);

        // 标识共享首/尾连接点；只用于连接点判定和侧向排序，不允许端点外贴合/重叠。
        s.shared_endpoint = cs.isSharedEndpoint(id, cid);

        sibs.push_back(std::move(s));
    }
    return sibs;
}

/// 沿曲线自适应采样查询最小SDF值, 性能优化: 上限从320降至100, 不调用queryWithGrad(仅用query)
// 沿曲线按弧长自适应采样并返回最小 signed distance field 值。
// 采样数量限制在固定范围内，以在障碍物安全性和单连接生成耗时之间取得平衡；
// 无效 SDF 返回正的大值，表示该检查无法发现障碍物穿透。
static double minSDFAlongCurveAdaptive(const BezierCurve& curve, const SDFField& sdf) {
    if (!sdf.valid())
        return 1e18;
    // 性能优化: 上限从320降至100
    int n = std::max(20, std::min(100, (int)std::ceil(curve.arcLength() / 0.30) + 1));
    double m = 1e18;
    for (const auto& pt : curve.sampleByArcLength(n))
        m = std::min(m, sdf.queryWithGrad(pt).first);
    return m;
}
// 前置声明：使用障碍物硬几何检查曲线是否命中障碍物，并可返回命中点。
static bool curveIntersectsObstacles(
    const BezierCurve& curve, const std::vector<Obstacle>& obstacles, Vec2d* out);
// 前置声明：计算曲线相对当前 mode 的 RoadEdge 净距违规量。
static double roadEdgeClearanceViolation(
    const BezierCurve& curve, const std::vector<Boundary>& boundaries, int mode);

// 汇总单条输出曲线的障碍物、围栏、边界、道路边缘、形态和自交状态，
// 并通过 FinalCurveAuditor 写回统一的状态与违规原因。该函数只更新 cc 审计字段。
void ConnectivityGenerationSession::validate(
    ConnectivityCurve& cc, const IntersectionInput& input, const SDFField& sdf) const {
    if (!cc.curve) return;
    auto& c = *cc.curve;
    FinalCurveAuditSnapshot audit;
    double ms = minSDFAlongCurveAdaptive(c, sdf);
    audit.max_obstacle_penetration = std::max(0.0, -ms);
    Vec2d obstacle_hit;
    audit.obstacle_intersection =
        curveIntersectsObstacles(c, input.obstacles, &obstacle_hit);
    if (!input.area.geometry.outer.empty()) {
        // 性能优化: 上限从240降至80,降低采样数
        int n = std::max(20, std::min(80, (int)std::ceil(c.arcLength() / 0.25) + 1));
        audit.update_fence_overflow = true;
        audit.max_fence_overflow = curveFenceOverflow(c, input.area.geometry, n);
    }
    audit.self_intersection = curveSelfIntersectsBusiness(c, 1.0);

    // 普通转向与U型调头都使用完整Boundary集合：curveBoundarySafety只裁掉首尾
    // 连续贴行段，而端点贴合RoadEdge的共线擦碰由 boundary_safety.h 的逐段豁免
    // 判据剔除（见 endpointGrazeExemptSegmentMasks），因此不再需要整条 boundary
    // 级别的端点粘连过滤——那个判据用折线首尾弦方向，对发夹形鼻端会失效。
    auto _entry = input.entryPtDir(cc.entry_lane_id);
    auto _exit  = input.exitPtDir(cc.exit_lane_id);
    Vec2d p0_v = _entry.first;
    Vec2d t0_v = _entry.second;
    Vec2d p1_v = _exit.first;
    Vec2d t1_v = _exit.second;

    // 检测是否为U-turn连接
    bool is_conn_uturn = (t0_v.norm() > 1e-8 && t1_v.norm() > 1e-8 &&
                          t0_v.normalized().dot(t1_v.normalized()) < -0.5);
    audit.uturn = is_conn_uturn;

    bool boundary_cross = false;
    if (is_conn_uturn) {
        boundary_cross = curveIntersectsBoundariesUTurn(c, input);
    } else {
        boundary_cross =
            curveIntersectsBoundaries(c, input) ||
            curveRawIntersectsAnyBoundary(c, input.boundaries, 0.15, 0.10);
    }
    audit.boundary_intersection = boundary_cross;

    bool road_edge_clearance_violation =
        roadEdgeClearanceViolation(c, input.boundaries, input.mode) > 0.05;
    audit.road_edge_clearance_violation = road_edge_clearance_violation;
    bool ordinary_single_axis_violation = false;
    bool ordinary_single_nonphysical =
        !is_conn_uturn && c.numSegments() == 1 &&
        ms >= -0.05 && cc.violation.max_fence_overflow <= 0.05 &&
        !boundary_cross && !road_edge_clearance_violation;
    if (ordinary_single_nonphysical &&
        !ordinarySingleCubicControlsValid(
            c, p0_v, t0_v, p1_v, t1_v, 1e-5, false)) {
        ordinary_single_axis_violation = true;
    }
    audit.ordinary_single_axis_violation = ordinary_single_axis_violation;
    FinalCurveAuditor().apply(audit, cc);
}

// 建立连接 ID 到结果数组下标的索引，供批量修复阶段稳定定位曲线。
// 重复 ID 时后出现的结果覆盖先出现的结果。
static std::unordered_map<ConnId, size_t> resultIndexById(
    const std::vector<ConnectivityCurve>& results) {
    std::unordered_map<ConnId, size_t> idx;
    for (size_t i = 0; i < results.size(); ++i)
        idx[results[i].id] = i;
    return idx;
}

// 从已有结果提取非空曲线快照，形成连接 ID 到 Bezier 曲线的快速查询表。
// 快照按值复制，后续候选替换不会反向修改 results。
static std::unordered_map<ConnId, BezierCurve> curveMapFromResults(
    const std::vector<ConnectivityCurve>& results) {
    std::unordered_map<ConnId, BezierCurve> curves;
    for (auto& cc : results)
        if (cc.curve)
            curves[cc.id] = *cc.curve;
    return curves;
}

// 判断连接是否携带至少两个点的固有几何形态，可用于构造固定形态曲线。
static bool hasFixedGeometry(const Connectivity& conn) {
    return conn.geometry.points.size() >= 2;
}

/// 将固有形态转为BezierCurve
// 将连接输入中的固有几何转换为输出曲线，并预先记录固定曲线的障碍物风险。
// 固有几何本身保留在 geometry 字段；若转换失败则返回无曲线结果。
static ConnectivityCurve makeFixedGeometryCurve(
    const Connectivity& conn, const IntersectionInput& input, const SDFField& sdf) {
    ConnectivityCurve cc;
    cc.id = conn.id;
    cc.entry_lane_id = conn.entry_lane_id;
    cc.exit_lane_id = conn.exit_lane_id;
    cc.turn_type = conn.turn_type;
    cc.geometry = conn.geometry;
    cc.fixed_shape = conn.fixed_shape;
    BezierCurve curve = FixedShapeInitializer().build(conn);
    if (!curve.empty())
        cc.curve = std::make_shared<BezierCurve>(curve);

    if (cc.curve) {
        double ms = minSDFAlongCurveAdaptive(*cc.curve, sdf);
        cc.violation.max_obstacle_penetration = std::max(0.0, -ms);
        Vec2d obstacle_hit;
        if (curveIntersectsObstacles(*cc.curve, input.obstacles, &obstacle_hit)) {
            cc.violation.max_obstacle_penetration =
                std::max(cc.violation.max_obstacle_penetration, 0.01);
            cc.violation.reason = "fixed geometry intersects obstacle";
            cc.status = CurveStatus::Degraded;
        }
    }
    return cc;
}

// 构造固定几何的临时曲线并检查其是否穿越任一障碍物。
// 用于决定该 fixed shape 是否必须降级为可重新生成的普通连接。
static bool fixedGeometryHitsObstacle(
    const Connectivity& conn, const IntersectionInput& input) {
    if (!hasFixedGeometry(conn))
        return false;
    BezierCurve curve = FixedShapeInitializer().build(conn);
    if (curve.empty())
        return false;
    return curveIntersectsObstacles(curve, input.obstacles, nullptr);
}

// 前置声明：判断两条曲线是否存在必须修复的同簇相交或贴合。
// allow_merge_funnel 保留在签名中用于兼容已有调用点，但严格同簇判定不接受
// 任何 merge-funnel 豁免。
static bool curvesHaveForbiddenSameClusterIntersection(
    const BezierCurve& a, const BezierCurve& b, double endpoint_tol = 0.15,
    bool allow_merge_funnel = false);

// 判断点是否落在多边形外环任一线段的容差范围内。
// 该结果与 polygonContains 组合使用，用于把障碍物边界视为硬几何区域。
static bool pointOnPolygonBoundary(const Vec2d& pt, const Polygon2d& poly, double tol = 1e-6) {
    const auto& ring = poly.outer;
    if (ring.size() < 2)
        return false;
    for (int i = 0; i < (int)ring.size(); ++i) {
        const Vec2d& a = ring[i];
        const Vec2d& b = ring[(i + 1) % ring.size()];
        if (pointToSegment(pt, a, b).first <= tol)
            return true;
    }
    return false;
}

// 判断点位于多边形内部或其边界上，统一处理内部点和边界点两种障碍物命中。
static bool pointInsideOrOnPolygon(const Vec2d& pt, const Polygon2d& poly, double tol = 1e-6) {
    return polygonContains(poly, pt) || pointOnPolygonBoundary(pt, poly, tol);
}

// 检查线段是否与多边形外环相交，并可返回第一处相交点。
// 只检查外环边，不负责判断线段完全位于多边形内部的情形。
static bool polygonSegmentIntersects(
    const Vec2d& a, const Vec2d& b, const Polygon2d& poly, Vec2d* out = nullptr) {
    const auto& ring = poly.outer;
    if (ring.size() < 2)
        return false;
    for (int i = 0; i < (int)ring.size(); ++i) {
        Vec2d isect;
        if (segmentsIntersect(a, b, ring[i], ring[(i + 1) % ring.size()], &isect)) {
            if (out)
                *out = isect;
            return true;
        }
    }
    return false;
}

// 返回障碍物的硬判定几何：优先使用原始 geometry，缺失时回退到 buffered_geometry。
// 返回引用的生命周期由输入障碍物保证。
static const Polygon2d& obstacleHardGeometry(const Obstacle& obs) {
    return obs.geometry.outer.empty() ? obs.buffered_geometry : obs.geometry;
}

// 通过曲线采样点和采样线段与障碍物硬多边形的相交，判断曲线是否进入障碍物。
// 可选输出命中位置；端点不做特殊豁免，因为障碍物命中本身属于物理风险。
static bool curveIntersectsObstacles(
    const BezierCurve& curve, const std::vector<Obstacle>& obstacles, Vec2d* out = nullptr) {
    if (obstacles.empty())
        return false;

    auto sampled = sampleCurveForIntersections(curve, 64);
    if (sampled.pts.size() < 2)
        return false;

    for (const auto& obs : obstacles) {
        const Polygon2d& poly = obstacleHardGeometry(obs);
        if (poly.outer.size() < 3)
            continue;

        BoundingBox2d obs_box = poly.bbox();
        if (!sampled.bbox.intersects(obs_box))
            continue;

        for (int i = 1; i + 1 < (int)sampled.pts.size(); ++i) {
            if (pointInsideOrOnPolygon(sampled.pts[i], poly)) {
                if (out)
                    *out = sampled.pts[i];
                return true;
            }
        }

        for (int i = 0; i + 1 < (int)sampled.pts.size(); ++i) {
            Vec2d isect;
            if (polygonSegmentIntersects(sampled.pts[i], sampled.pts[i + 1], poly, &isect)) {
                if (out)
                    *out = isect;
                return true;
            }
        }
    }
    return false;
}

// 返回点到两条采样曲线四个端点的最小距离，用于过滤共享连接点附近的合法相遇。
static double distToSampledEndpoints(const Vec2d& pt, const SampledCurve& a, const SampledCurve& b) {
    double d = (pt - a.start).norm();
    d = std::min(d, (pt - a.end).norm());
    d = std::min(d, (pt - b.start).norm());
    d = std::min(d, (pt - b.end).norm());
    return d;
}

// 对两条采样曲线执行带业务端点豁免的线段相交检查，并可选检测近距离平行贴合。
// 返回 true 表示存在端点容差之外的交叉或重叠，并可写出冲突位置。
//
// 分块剪枝说明：内层用 b 的分块粗包围盒整段跳过与当前 a 段必然分离的区间。
// 判据与逐段的 max<min 分离判断完全一致（BoundingBox2d::intersects 用闭区间），
// 所以跳过的都是原本也会被 continue 的配对；i、j 仍按升序访问，写出的首个
// 冲突位置与逐段遍历一致。
static bool sampledCurvesIntersectBusiness(
    const SampledCurve& a, const SampledCurve& b, double endpoint_tol,
    Vec2d* out = nullptr, bool detect_near_overlap = false) {
    if (a.pts.size() < 2 || b.pts.size() < 2 || !a.bbox.intersects(b.bbox))
        return false;
    const int na = static_cast<int>(a.pts.size()) - 1;
    const int nb = static_cast<int>(b.pts.size()) - 1;
    const bool has_boxes = a.segment_boxes.size() == a.pts.size() - 1 &&
        b.segment_boxes.size() == b.pts.size() - 1;
    const bool has_blocks = has_boxes &&
        b.block_boxes.size() ==
            static_cast<std::size_t>((nb + kSampleBlockSize - 1) / kSampleBlockSize) &&
        a.block_boxes.size() ==
            static_cast<std::size_t>((na + kSampleBlockSize - 1) / kSampleBlockSize);
    const int nblock = has_blocks ? static_cast<int>(b.block_boxes.size()) : 1;
    for (int ai = 0; ai < na; ++ai) {
        if (has_blocks) {
            // a 侧整块与 b 的总包围盒分离时，跳过该块的全部 a 段。
            if ((ai % kSampleBlockSize) == 0 &&
                !a.block_boxes[ai / kSampleBlockSize].intersects(b.bbox)) {
                ai += kSampleBlockSize - 1;
                continue;
            }
            if (!a.segment_boxes[ai].intersects(b.bbox))
                continue;
        }
        for (int blk = 0; blk < nblock; ++blk) {
            int bi_begin = 0;
            int bi_end = nb;
            if (has_blocks) {
                if (!a.segment_boxes[ai].intersects(b.block_boxes[blk]))
                    continue;
                bi_begin = blk * kSampleBlockSize;
                bi_end = std::min(nb, bi_begin + kSampleBlockSize);
            }
            for (int bi = bi_begin; bi < bi_end; ++bi) {
                if (has_boxes &&
                    !a.segment_boxes[ai].intersects(b.segment_boxes[bi]))
                    continue;
                const Vec2d& a0 = a.pts[ai];
                const Vec2d& a1 = a.pts[ai + 1];
                const Vec2d& b0 = b.pts[bi];
                const Vec2d& b1 = b.pts[bi + 1];
                if (std::max(a0.x(), a1.x()) < std::min(b0.x(), b1.x()) ||
                    std::max(b0.x(), b1.x()) < std::min(a0.x(), a1.x()) ||
                    std::max(a0.y(), a1.y()) < std::min(b0.y(), b1.y()) ||
                    std::max(b0.y(), b1.y()) < std::min(a0.y(), a1.y()))
                    continue;
                const Vec2d amid = a.segment_midpoints.size() == a.pts.size() - 1
                    ? a.segment_midpoints[ai]
                    : 0.5 * (a0 + a1);
                const Vec2d bmid = b.segment_midpoints.size() == b.pts.size() - 1
                    ? b.segment_midpoints[bi]
                    : 0.5 * (b0 + b1);
                if ((amid - bmid).squaredNorm() > 900.0)
                    continue;
                Vec2d isect;
                if (!segmentsIntersect(a0, a1, b0, b1, &isect)) {
                    if (!detect_near_overlap)
                        continue;
                    Vec2d ad = a1 - a0;
                    Vec2d bd = b1 - b0;
                    if (ad.norm() < 1e-8 || bd.norm() < 1e-8)
                        continue;
                    if (std::abs(ad.normalized().dot(bd.normalized())) < 0.96)
                        continue;
                    double near_dist = std::min({
                        pointToSegment(a0, b0, b1).first,
                        pointToSegment(a1, b0, b1).first,
                        pointToSegment(b0, a0, a1).first,
                        pointToSegment(b1, a0, a1).first
                    });
                    if (near_dist > 0.18)
                        continue;
                    isect = 0.5 * (amid + bmid);
                }
                if (distToSampledEndpoints(isect, a, b) <= endpoint_tol)
                    continue;
                if (out)
                    *out = isect;
                return true;
            }
        }
    }
    return false;
}

// 统计一条采样曲线与兄弟采样曲线的约束交叉数，可跳过结构性交叉豁免项。
// detect_near_overlap 开启时把近距离同向贴合也计为冲突。
static int sampledSiblingCrossCount(
    const SampledCurve& curve, const std::vector<SampledSiblingCurve>& siblings,
    bool constrained_only, double endpoint_tol = kClusterEndpointTol,
    bool detect_near_overlap = false) {
    int count = 0;
    for (const auto& sib : siblings) {
        if (constrained_only && sib.exempt_a1)
            continue;
        if (sampledCurvesIntersectBusiness(
                curve, sib.sampled, endpoint_tol, nullptr, detect_near_overlap))
            ++count;
    }
    return count;
}

static int sampledSiblingCrossCount(
    const BezierCurve& curve, const std::vector<SampledSiblingCurve>& siblings,
    bool constrained_only, double endpoint_tol = kClusterEndpointTol,
    bool detect_near_overlap = false, bool allow_merge_funnel = false) {
    int count = 0;
    SampledCurve sample = sampleCurveForIntersections(curve);
    for (const auto& sib : siblings) {
        if (constrained_only && sib.exempt_a1)
            continue;
        if (!sampledCurvesIntersectBusiness(
                sample, sib.sampled, endpoint_tol, nullptr, detect_near_overlap))
            continue;
        if (!sib.curve.empty() &&
            !curvesHaveForbiddenSameClusterIntersection(
                curve, sib.curve, endpoint_tol, allow_merge_funnel))
            continue;
        ++count;
    }
    return count;
}

// 统计候选曲线与固定形态兄弟之间的真实同簇冲突数。
// 只计入非结构性交叉且 fixed_shape 为真的兄弟，用于保护固定几何优先级。
static int rawFixedSiblingCrossCount(
    const BezierCurve& curve, const std::vector<SiblingCurve>& siblings,
        double endpoint_tol = kClusterEndpointTol) {
    int count = 0;
    for (const auto& sib : siblings) {
        if (!sib.fixed_shape || sib.exempt_a1)
            continue;
        if (curvesHaveForbiddenSameClusterIntersection(
                curve, sib.curve, endpoint_tol))
            ++count;
    }
    return count;
}

// 量化共享端点曲线相对参考法线的侧向排序违反程度。
// 只检查共享端点后有限距离内的采样点，返回值越大表示扇出顺序越不符合预期。
static double sharedEndpointSideViolation(
    const SampledCurve& curve, const std::vector<SampledSiblingCurve>& siblings,
    double endpoint_tol = kClusterEndpointTol) {
    if (curve.pts.size() < 4)
        return 0.0;

    double violation = 0.0;
    constexpr double SIDE_MARGIN = 0.25;
    for (const auto& sib : siblings) {
        if (sib.exempt_a1 || !sib.shared_endpoint)
            continue;
        if (sib.expected_side == 0 || sib.ref_perp.norm() < 1e-9)
            continue;
        if (sib.sampled.pts.size() < 4)
            continue;

        bool share_start = (curve.start - sib.sampled.start).norm() <= endpoint_tol;
        bool share_end = (curve.end - sib.sampled.end).norm() <= endpoint_tol;
        if (!share_start && !share_end)
            continue;

        for (const auto& pt : curve.pts) {
            Vec2d anchor = share_start ? curve.start : curve.end;
            double d = (pt - anchor).norm();
            if (d <= endpoint_tol || d > 8.0)
                continue;

            double best_d2 = std::numeric_limits<double>::infinity();
            double best_lat = 0.0;
            for (const auto& sp : sib.sampled.pts) {
                Vec2d sib_anchor = share_start ? sib.sampled.start : sib.sampled.end;
                double sd = (sp - sib_anchor).norm();
                if (sd <= endpoint_tol || sd > 8.0)
                    continue;
                double d2 = (sp - pt).squaredNorm();
                if (d2 < best_d2) {
                    best_d2 = d2;
                    best_lat = sp.dot(sib.ref_perp);
                }
            }
            if (!std::isfinite(best_d2) || best_d2 > 36.0)
                continue;

            double diff = pt.dot(sib.ref_perp) - best_lat;
            if (sib.expected_side == +1)
                violation += std::max(0.0, diff + SIDE_MARGIN);
            else
                violation += std::max(0.0, -diff + SIDE_MARGIN);
        }
    }
    return violation;
}

// 从 RepairImpactClosure 关联对中找出与指定种子相连且实际发生禁止交叉的连接。
// 结果用于非纯几何场景的定向重生成，避免无关连接被纳入修复预算。
static std::unordered_set<ConnId> crossingIdsTouchingSeeds(
    const std::vector<ConnectivityCurve>& results,
    const ClusterOrderSolver& cs, const std::unordered_set<ConnId>& seeds,
    double endpoint_tol = kClusterEndpointTol) {
    std::unordered_set<ConnId> bad;
    if (seeds.empty())
        return bad;
    auto idx = resultIndexById(results);
    std::unordered_map<ConnId, SampledCurve> samples;
    samples.reserve(results.size());
    for (const auto& cc : results)
        if (cc.curve)
            samples.emplace(cc.id, sampleCurveForIntersections(*cc.curve));
    const RepairImpactClosure closure =
        RepairImpactClosureBuilder().build(cs.pairs(), seeds);
    for (const auto& p : closure.pairs) {
        auto ia = idx.find(p.id_a), ib = idx.find(p.id_b);
        if (ia == idx.end() || ib == idx.end())
            continue;
        auto sa = samples.find(p.id_a);
        auto sb = samples.find(p.id_b);
        if (sa == samples.end() || sb == samples.end())
            continue;
        double pair_tol = endpoint_tol;
        if (curvesHaveForbiddenSameClusterIntersection(
                *results[ia->second].curve, *results[ib->second].curve,
                pair_tol)) {
            if (seeds.count(p.id_a))
                bad.insert(p.id_a);
            if (seeds.count(p.id_b))
                bad.insert(p.id_b);
        }
    }
    return bad;
}

// 检查两条曲线是否在共享端点之外发生业务相交或重叠。
// 与普通业务判定相比，该包装器明确表达“连接点外不得相遇”的调用意图。
static bool curvesIntersectBeyondSharedEndpointOverlap(
    const BezierCurve& a, const BezierCurve& b, double endpoint_tol = 0.15) {
    return curvesIntersectBusiness(a, b, endpoint_tol);
}

// 在采样线段层面检测两条曲线是否存在非端点的同向近距离贴合。
// 该检查独立于真正的线段相交，用来禁止重叠、贴行和共享中段。
static bool sampledCurvesHaveForbiddenAdherence(
        const SampledCurve& a, const SampledCurve& b,
        double endpoint_tol = 0.15, double adherent_tol = 0.18) {
    if (a.pts.size() < 2 || b.pts.size() < 2)
        return false;
    BoundingBox2d abox = a.bbox;
    abox.min_pt -= Vec2d(adherent_tol, adherent_tol);
    abox.max_pt += Vec2d(adherent_tol, adherent_tol);
    if (!abox.intersects(b.bbox))
        return false;
    for (int ai = 0; ai + 1 < (int)a.pts.size(); ++ai) {
        Vec2d a0 = a.pts[ai];
        Vec2d a1 = a.pts[ai + 1];
        Vec2d ad = a1 - a0;
        double alen = ad.norm();
        if (alen < 1e-8)
            continue;
        Vec2d au = ad / alen;
        for (int bi = 0; bi + 1 < (int)b.pts.size(); ++bi) {
            Vec2d b0 = b.pts[bi];
            Vec2d b1 = b.pts[bi + 1];
            Vec2d bd = b1 - b0;
            double blen = bd.norm();
            if (blen < 1e-8)
                continue;
            Vec2d bu = bd / blen;
            if (std::abs(au.dot(bu)) < 0.96)
                continue;
            if (std::abs(cross2d(au, b0 - a0)) > adherent_tol ||
                std::abs(cross2d(au, b1 - a0)) > adherent_tol)
                continue;

            double b0s = (b0 - a0).dot(au);
            double b1s = (b1 - a0).dot(au);
            double lo = std::max(0.0, std::min(b0s, b1s));
            double hi = std::min(alen, std::max(b0s, b1s));
            if (lo > hi + adherent_tol)
                continue;
            Vec2d witness = a0 + 0.5 * (lo + hi) * au;
            if (distToSampledEndpoints(witness, a, b) <= endpoint_tol)
                continue;
            return true;
        }
    }
    return false;
}

// 对 Bezier 曲线进行统一采样后检查非端点贴合，隐藏采样参数并复用采样实现。
static bool curvesHaveForbiddenAdherence(
        const BezierCurve& a, const BezierCurve& b,
        double endpoint_tol = 0.15) {
    return sampledCurvesHaveForbiddenAdherence(
        sampleCurveForIntersections(a, 64),
        sampleCurveForIntersections(b, 64),
        endpoint_tol);
}

// 左/右转和掉头的几何穿越可豁免；贴合/重叠仍禁止。
// 判断同簇曲线之间的几何交叉是否属于允许的 U-turn 与左右转结构性交叉。
// 任何非端点贴合优先判为禁止；两条曲线形态不满足豁免组合时也返回 false。
// 结构性汇入/分流的端点容差：掉头与"共享同一端点的非掉头连接"之间，收敛段内的
// 贴行与穿越由车道拓扑决定，任何造型都消不掉，必须按收敛段长度而不是固定容差豁免。
//
// 三重限定，缺一不可：
//
// ①**恰好一条是掉头**。同家族的掉头两两之间同样共享端点，但它们的左右次序由
//   家族横向阶梯与嵌套规则单独承载（小径必须整体落在大径控制面内），一旦按收敛段
//   豁免，110000703-u 的 41|44 这类"外层腿穿进内层两腿之间"的真实违约就会重新放行。
//
// ②**该掉头的端点走廊退化**（进出车道端点横向间距 < kDegenerateUTurnCorridor）。
//   这是"冲突不可消除"的充要几何条件：走廊只有几十厘米时，掉头被人行横道逼出的
//   首尾直段几乎与出口车道共线，汇入同一出口车道的兄弟必然贴着它跑并至少穿一次；
//   走廊有米级宽度的掉头则完全有重新造型的余地，必须继续按严格规则修复。
//   放开这一条的代价实测很大：全量 21 份里有 38 对同簇交叉会被一并放行，其中
//   110003449 的 12|52、13|53 正是上一轮靠二段式横向偏置刚修好的对。
//
// ③**冲突落在收敛段内**（见 sharedEndpointMergeFunnelRadius）。
//
// 110000703-u 的 41 是唯一命中三条的成员：端点走廊 0.3099m，人行横道要求首尾直段
// 9.32/9.24m，而汇入同一条出口车道的近直行 13 在这 9.3m 内与 41 的入口直段横向只差
// 0.005~0.31m，且 13 的终点就落在 41 的走廊内部——13 必须从走廊外进入走廊内，
// 穿越在拓扑上强制发生。此前正是这一对把 41 逼回单段拱形，破坏了掉头的
// shape.crosswalk.segments 约束。
static const double kDegenerateUTurnCorridor = 0.60;

static double uturnEndpointCorridor(const BezierCurve& u) {
    if (u.empty())
        return std::numeric_limits<double>::infinity();
    const Vec2d t0 = u.startTan();
    const Vec2d t1 = u.endTan();
    if (t0.norm() < 1e-8 || t1.norm() < 1e-8)
        return std::numeric_limits<double>::infinity();
    Vec2d axis = (t0.normalized() - t1.normalized());
    if (axis.norm() < 1e-8)
        return std::numeric_limits<double>::infinity();
    axis = axis.normalized();
    const Vec2d lateral(-axis.y(), axis.x());
    return std::abs((u.endPt() - u.startPt()).dot(lateral));
}

static double sameClusterMergeFunnelTol(const BezierCurve& a,
                                        const BezierCurve& b) {
    // 诊断开关：ISG_MERGE_FUNNEL=0 关闭本豁免，用于源码级 A/B 消融。
    static const bool enabled = [] {
        const char* v = std::getenv("ISG_MERGE_FUNNEL");
        return v == nullptr || v[0] != '0';
    }();
    if (!enabled)
        return 0.0;
    const bool a_uturn = curveLooksUTurnForClusterExemption(a);
    const bool b_uturn = curveLooksUTurnForClusterExemption(b);
    if (a_uturn == b_uturn)
        return 0.0;
    if (uturnEndpointCorridor(a_uturn ? a : b) >= kDegenerateUTurnCorridor)
        return 0.0;
    return sharedEndpointMergeFunnelRadius(a, b);
}

// 共享同一出/入口端点的“明显转向 + 近直行”在汇入端点后的极短距离内
// 允许存在一阶收敛贴合：两条曲线随后必须分离，且不能发生真实业务相交。
// 这不是一般性的同簇放宽，而是避免把道路拓扑必然产生的端点汇入误判为
// 造型冲突，进而把规范转向把手压短成尖钩。
static bool sharedTurnStraightConvergenceOnly(
    const BezierCurve& a, const BezierCurve& b, double endpoint_tol) {
    if (a.empty() || b.empty() || curveLooksUTurnForClusterExemption(a) ||
        curveLooksUTurnForClusterExemption(b))
        return false;
    const bool shared_start = (a.startPt() - b.startPt()).norm() <= endpoint_tol;
    const bool shared_end = (a.endPt() - b.endPt()).norm() <= endpoint_tol;
    if (shared_start == shared_end)
        return false;
    const Vec2d ta = shared_start ? a.startTan() : a.endTan();
    const Vec2d tb = shared_start ? b.startTan() : b.endTan();
    if (ta.norm() < 1e-8 || tb.norm() < 1e-8 ||
        ta.normalized().dot(tb.normalized()) < 0.98)
        return false;
    auto turn_strength = [](const BezierCurve& curve) {
        const Vec2d chord = curve.endPt() - curve.startPt();
        const Vec2d tangent = curve.startTan();
        if (chord.norm() < 1e-8 || tangent.norm() < 1e-8)
            return 0.0;
        return std::abs(cross2d(tangent.normalized(), chord.normalized()));
    };
    const double ta_strength = turn_strength(a);
    const double tb_strength = turn_strength(b);
    if (!((ta_strength > 0.35 && tb_strength < 0.25) ||
          (tb_strength > 0.35 && ta_strength < 0.25)))
        return false;
    // 真实相交优先于收敛豁免；这里只处理端点邻域的近距离贴合。
    return !curvesIntersectBusiness(a, b, endpoint_tol);
}

// 汇总同簇相交、接触和重叠规则，返回是否存在必须修复的非豁免冲突。
// 先排除允许的结构性交叉，再检查业务相交，因此是同簇约束的统一入口。
//
// `allow_merge_funnel` 默认关闭：全部既有调用点保持严格判定，逐字节不改变行为。
// 只有三段式掉头候选搜索的"兜底档"（见 searchSegmentedUTurnTwoPass）会打开它，
// 从而把结构性汇入豁免限制在"人行横道净空与同簇非交无法同时满足"的死局内。
// 一律打开会外溢：实测 100000643 的 125-113、115-103 会被一并放行。
static bool curvesHaveForbiddenSameClusterIntersection(
    const BezierCurve& a, const BezierCurve& b, double endpoint_tol,
    bool allow_merge_funnel) {
    // StructuralCross 只在调用方根据 CurvePair 拓扑关系跳过；几何层不再根据
    // 曲线形态推断任何结构性豁免。除真实首/尾连接点误差外，相交与重叠均违规。
    (void)allow_merge_funnel;
    return curvesIntersectBusiness(a, b, endpoint_tol);
}

// 扫描全部非结构性交叉约束，收集当前结果中发生禁止相交且未受固定形态保护的
// 连接 ID。该集合是批量重生成阶段的冲突入口。
static std::unordered_set<ConnId> allConstrainedCrossingIds(
    const std::vector<ConnectivityCurve>& results,
    const ClusterOrderSolver& cs,
    const std::unordered_set<ConnId>& preserved_fixed_ids,
    double endpoint_tol = kClusterEndpointTol) {
    std::unordered_set<ConnId> bad;
    auto idx = resultIndexById(results);
    std::unordered_map<ConnId, SampledCurve> samples;
    samples.reserve(results.size());
    for (const auto& cc : results)
        if (cc.curve)
            samples.emplace(cc.id, sampleCurveForIntersections(*cc.curve));

    for (auto& p : cs.pairs()) {
        if (p.exempt == CrossExemption::StructuralCross)
            continue;
        auto ia = idx.find(p.id_a), ib = idx.find(p.id_b);
        if (ia == idx.end() || ib == idx.end())
            continue;
        auto sa = samples.find(p.id_a);
        auto sb = samples.find(p.id_b);
        if (sa == samples.end() || sb == samples.end())
            continue;
        double pair_tol = endpoint_tol;
        if (!curvesHaveForbiddenSameClusterIntersection(
                *results[ia->second].curve, *results[ib->second].curve, pair_tol))
            continue;

        if (!preserved_fixed_ids.count(p.id_a))
            bad.insert(p.id_a);
        if (!preserved_fixed_ids.count(p.id_b))
            bad.insert(p.id_b);
    }
    return bad;
}

// 对最终结果执行一次完整同簇审计，并把冲突位置、状态和原因写回输出结果。
// 此函数用于最终标注，不负责重新生成曲线。
static void annotateClusterCrossings(
    std::vector<ConnectivityCurve>& results,
    const ClusterOrderSolver& cs,
    double endpoint_tol = kClusterEndpointTol) {
    // 同簇标注可能在纯拓扑收口和最终收口阶段重复执行；该字段表示本次
    // 最终标注的发现，不能把前一阶段已经替换掉的旧曲线交叉继续带入状态。
    for (auto& cc : results)
        cc.violation.exempt_crosses.clear();
    std::unordered_map<ConnId, SampledCurve> samples;
    samples.reserve(results.size());
    // 标注只用于输出告警，24个采样点可兼顾速度与稳定性。
    const int N_ANNOTATE = 24;
    for (const auto& cc : results)
        if (cc.curve)
            samples.emplace(cc.id, sampleCurveForIntersections(*cc.curve, N_ANNOTATE));
    const FinalPairViolationDetector detector =
        [&](const CurvePair& pair, const BezierCurve& curve_a,
            const BezierCurve& curve_b, Vec2d& location) {
            const auto sa = samples.find(pair.id_a);
            const auto sb = samples.find(pair.id_b);
            if (sa == samples.end() || sb == samples.end())
                return false;
            const double pair_tol = endpoint_tol;
            BezierCurve check_a = curve_a;
            BezierCurve check_b = curve_b;
            const bool both_segmented_uturn =
                curveLooksUTurnForClusterExemption(check_a) &&
                curveLooksUTurnForClusterExemption(check_b) &&
                check_a.numSegments() == 3 && check_b.numSegments() == 3 &&
                segmentLooksStraight(check_a.segs.front()) &&
                segmentLooksStraight(check_b.segs.front()) &&
                segmentLooksStraight(check_a.segs.back()) &&
                segmentLooksStraight(check_b.segs.back());
            bool trimmed_shared_lead = false;
            if (both_segmented_uturn &&
                (check_a.startPt() - check_b.startPt()).norm() <= endpoint_tol &&
                (check_a.segs.front().ctrl[3] -
                 check_b.segs.front().ctrl[3]).norm() >
                    kUTurnSharedLeadFanMinSeparation) {
                check_a.segs.erase(check_a.segs.begin());
                check_b.segs.erase(check_b.segs.begin());
                trimmed_shared_lead = true;
            }
            if (both_segmented_uturn &&
                (curve_a.endPt() - curve_b.endPt()).norm() <= endpoint_tol &&
                !check_a.empty() && !check_b.empty() &&
                (curve_a.segs.back().ctrl[0] -
                 curve_b.segs.back().ctrl[0]).norm() >
                    kUTurnSharedLeadFanMinSeparation) {
                check_a.segs.pop_back();
                check_b.segs.pop_back();
                trimmed_shared_lead = true;
            }
            const bool interior = trimmed_shared_lead
                ? curvesIntersectBusiness(check_a, check_b, 1.5)
                : curvesHaveForbiddenSameClusterIntersection(
                      curve_a, curve_b, pair_tol);
            if (!interior)
                return false;
            if (isgDebugPairRepair()) {
                fprintf(stderr,
                        "[ISG_PROFILE] annotate cluster pair %s|%s trimmed=%d segmented=%d start=%.3f q0=%.3f end=%.3f q1=%.3f\n",
                        pair.id_a.c_str(), pair.id_b.c_str(),
                        trimmed_shared_lead ? 1 : 0,
                        both_segmented_uturn ? 1 : 0,
                        (curve_a.startPt() - curve_b.startPt()).norm(),
                        (curve_a.segs.front().ctrl[3] -
                         curve_b.segs.front().ctrl[3]).norm(),
                        (curve_a.endPt() - curve_b.endPt()).norm(),
                        (curve_a.segs.back().ctrl[0] -
                         curve_b.segs.back().ctrl[0]).norm());
            }
            sampledCurvesIntersectBusiness(
                sa->second, sb->second, pair_tol, &location);
            return true;
        };
    const FinalPairAuditor auditor;
    auditor.apply(auditor.audit(cs.pairs(), results, detector), results);
}

// 构造自然单段 G1 三次曲线，快速预判其是否会进入障碍物或 SDF 禁区。
static bool naturalCubicHitsObstacle(
    const Vec2d& p0, const Vec2d& t0, const Vec2d& p1, const Vec2d& t1,
    const SDFField& sdf, const std::vector<Obstacle>& obstacles) {
    Vec2d T0 = t0.norm() > 1e-8 ? t0.normalized() : Vec2d(1, 0);
    Vec2d T1 = t1.norm() > 1e-8 ? t1.normalized() : Vec2d(1, 0);
    BezierCurve trial;
    trial.segs.push_back(makeCubicG1(p0, T0, p1, T1, 0.4));
    if (curveIntersectsObstacles(trial, obstacles))
        return true;
    return sdf.valid() && minSDFAlongCurveAdaptive(trial, sdf) < 0.0;
}

// 采样曲线内部点并测量其到所有 RoadEdge 的净距，逐点按"该点可达净距下限"折算亏欠。
// 仅跳过曲线真实首/尾连接点的浮点误差，连接点之后的贴合、重叠和近边缘段都必须参与
// 净距评估；可达下限由连接点自身净距与端点切向共同决定，详见
// `constraints/road_edge_clearance.h`。
static RoadEdgeClearanceMeasure measureRoadEdgeClearanceAlongCurve(
    const BezierCurve& curve, const std::vector<Boundary>& boundaries,
    double clearance, double endpoint_skip = kConnectionPointTolerance) {
    RoadEdgeClearanceMeasure measure;
    if (curve.empty() || boundaries.empty())
        return measure;
    const double endpoint_tol = std::min(
        std::max(0.0, endpoint_skip), kConnectionPointTolerance);
    int n = std::max(32, std::min(160, (int)std::ceil(curve.arcLength() / 0.20) + 1));
    const std::vector<Vec2d> samples = curve.sampleByArcLength(n);
    std::vector<BoundingBox2d> road_edge_boxes(boundaries.size());
    for (std::size_t bi = 0; bi < boundaries.size(); ++bi) {
        if (boundaries[bi].type == Boundary::Type::RoadEdge)
            road_edge_boxes[bi] = boundaries[bi].geometry.bbox();
    }
    // 逐点求到 RoadEdge 的最小距离；bound 为已知上界，用于 bbox 剪枝。
    auto nearestEdge = [&](const Vec2d& pt, double bound) -> double {
        double best = bound;
        for (std::size_t boundary_index = 0;
             boundary_index < boundaries.size(); ++boundary_index) {
            const auto& bnd = boundaries[boundary_index];
            if (bnd.type != Boundary::Type::RoadEdge)
                continue;
            const BoundingBox2d& boundary_box = road_edge_boxes[boundary_index];
            const double dx = pt.x() < boundary_box.min_pt.x()
                ? boundary_box.min_pt.x() - pt.x()
                : (pt.x() > boundary_box.max_pt.x()
                    ? pt.x() - boundary_box.max_pt.x() : 0.0);
            const double dy = pt.y() < boundary_box.min_pt.y()
                ? boundary_box.min_pt.y() - pt.y()
                : (pt.y() > boundary_box.max_pt.y()
                    ? pt.y() - boundary_box.max_pt.y() : 0.0);
            if (std::isfinite(best) && dx * dx + dy * dy >= best * best)
                continue;
            for (int bi = 0; bi + 1 < (int)bnd.geometry.points.size(); ++bi) {
                best = std::min(
                    best, pointToSegment(pt, bnd.geometry.points[bi],
                                         bnd.geometry.points[bi + 1]).first);
            }
        }
        return best;
    };
    const double inf = std::numeric_limits<double>::infinity();
    measure.start_distance = nearestEdge(curve.startPt(), inf);
    measure.end_distance = nearestEdge(curve.endPt(), inf);
    const Vec2d start_tan = curve.startTan().norm() > 1e-12
        ? curve.startTan().normalized() : Vec2d(1, 0);
    const Vec2d end_tan = curve.endTan().norm() > 1e-12
        ? curve.endTan().normalized() : Vec2d(1, 0);
    for (auto& pt : samples) {
        const double to_start = (pt - curve.startPt()).norm();
        const double to_end = (pt - curve.endPt()).norm();
        if (to_start <= endpoint_tol || to_end <= endpoint_tol)
            continue;
        const double distance = nearestEdge(pt, inf);
        if (!std::isfinite(distance))
            continue;
        measure.valid = true;
        if (distance < measure.minimum) {
            measure.minimum = distance;
            measure.location = pt;
        }
        if (distance >= clearance)
            continue;
        const bool near_start = to_start <= to_end;
        const double span = near_start ? to_start : to_end;
        double ray_distance = -1.0;
        if (roadEdgeClearanceNeedsFloor(span)) {
            const Vec2d origin = near_start ? curve.startPt() : curve.endPt();
            const Vec2d direction = near_start ? start_tan : -end_tan;
            ray_distance = nearestEdge(origin + direction * span, inf);
        }
        const double floor_value = roadEdgeClearanceFloor(
            clearance, span, ray_distance);
        if (floor_value - distance > measure.deficit) {
            measure.deficit = floor_value - distance;
            measure.deficit_location = pt;
        }
    }
    return measure;
}

// 将曲线到 RoadEdge 的净距转换为当前 mode 下可归责于曲线的净距违规量。
// 无需净距约束或没有有效边界时返回零。端点邻域的亏欠若完全由连接点位置强制，
// 则不计入违规——否则生成侧会为一笔还不掉的账把每个候选都判成物理风险。
static double roadEdgeClearanceViolation(
    const BezierCurve& curve, const std::vector<Boundary>& boundaries, int mode) {
    double clearance = roadEdgeAvoidanceClearanceForMode(mode);
    if (clearance <= 0.0)
        return 0.0;
    const RoadEdgeClearanceMeasure measure =
        measureRoadEdgeClearanceAlongCurve(curve, boundaries, clearance);
    return roadEdgeClearanceDeficit(measure, clearance);
}

// 计算候选曲线的边界避让代价，合并道路边缘净距和真实边界穿越惩罚。
// 该值供 U-turn 搜索器排序使用，不直接修改曲线。
static double boundaryAvoidancePenalty(
    const BezierCurve& curve, const IntersectionInput& input,
    bool is_uturn = false) {
    double penalty = roadEdgeClearanceViolation(curve, input.boundaries, input.mode);
    // 所有mode都禁止与Boundary贴合、重叠或相交；mode=2在此基础上
    // 额外要求RoadEdge非端点1m净距。
    if ((is_uturn ? curveIntersectsBoundariesUTurn(curve, input)
                  : curveIntersectsBoundaries(curve, input)) ||
        curveRawIntersectsAnyBoundary(curve, input.boundaries, 0.15, 0.10))
        penalty = std::max(penalty, 0.10);
    return penalty;
}

// 仅检查曲线与 RoadEdge 的非端点真实相交，忽略其它 Boundary 类型。
static bool curveRawIntersectsRoadEdge(
    const BezierCurve& curve, const std::vector<Boundary>& boundaries,
    double curve_endpoint_tol = 0.15, double boundary_endpoint_tol = 0.10) {
    return curveRawIntersectsBoundariesImpl(
        curve, boundaries, true, curve_endpoint_tol, boundary_endpoint_tol);
}

// 检查一条直线段是否真实穿越 RoadEdge；只忽略曲线侧真实首/尾连接点。
static bool segmentRawIntersectsRoadEdge(
    const Vec2d& a, const Vec2d& b, const std::vector<Boundary>& boundaries,
    double endpoint_tol = 0.15, double boundary_endpoint_tol = 0.10) {
    (void)boundary_endpoint_tol;
    const double strict_endpoint_tol = std::min(
        std::max(0.0, endpoint_tol), kConnectionPointTolerance);
    for (const auto& bnd : boundaries) {
        if (bnd.type != Boundary::Type::RoadEdge ||
            bnd.geometry.points.size() < 2)
            continue;
        std::vector<Vec2d> pts = toVec2dArray(bnd.geometry.points);
        for (int i = 0; i + 1 < (int)pts.size(); ++i) {
            if (segmentHasForbiddenBoundaryContact(
                    a, b, pts[i], pts[i + 1], a, b, strict_endpoint_tol))
                return true;
        }
    }
    return false;
}

// 判断线段及其指定宽度走廊是否触碰 RoadEdge。
// 用于在自然弦线本身未相交但候选搜索可能撞边时提前识别修复需求。
static bool segmentCorridorTouchesRoadEdge(
    const Vec2d& a, const Vec2d& b, const std::vector<Boundary>& boundaries,
    double corridor_tol = 1.20) {
    BoundingBox2d seg_box;
    seg_box.expand(a);
    seg_box.expand(b);
    seg_box.min_pt -= Vec2d(corridor_tol, corridor_tol);
    seg_box.max_pt += Vec2d(corridor_tol, corridor_tol);
    for (const auto& bnd : boundaries) {
        if (bnd.type != Boundary::Type::RoadEdge ||
            bnd.geometry.points.size() < 2 ||
            !seg_box.intersects(bnd.geometry.bbox()))
            continue;
        std::vector<Vec2d> pts = toVec2dArray(bnd.geometry.points);
        for (int i = 0; i + 1 < (int)pts.size(); ++i) {
            if (segmentsIntersect(a, b, pts[i], pts[i + 1]))
                return true;
            double d = std::min({
                pointToSegment(pts[i], a, b).first,
                pointToSegment(pts[i + 1], a, b).first,
                pointToSegment(a, pts[i], pts[i + 1]).first,
                pointToSegment(b, pts[i], pts[i + 1]).first
            });
            if (d <= corridor_tol)
                return true;
        }
    }
    return false;
}

// 检查曲线是否完全位于粗围栏内。生成期必须与 ConstraintEvaluator 的
// curveInsideFence 使用同一严格口径，不能让低密度采样或 0.10m 边界宽容
// 把候选放过去，最终审计再判定为围栏越界。
static bool curveLeavesFence(const BezierCurve& curve, const Polygon2d& fence) {
    if (fence.outer.empty())
        return false;
    return !curveInsideFence(curve, fence, 64);
}

// 检测曲线离散曲率符号是否发生有效翻转，以识别不期望的 S 形路径。
// 小于 eps 的曲率视为数值噪声，不参与符号判断。
static bool curveHasCurvatureSignFlip(const BezierCurve& curve, double eps = 0.10) {
    int sign = 0;
    for (const auto& seg : curve.segs) {
        for (int i = 1; i < 20; ++i) {
            double t = (double)i / 20.0;
            Vec2d d1 = seg.evalDeriv1(t);
            Vec2d d2 = seg.evalDeriv2(t);
            double den = std::pow(d1.squaredNorm(), 1.5);
            if (den < 1e-12)
                continue;
            double k = cross2d(d1, d2) / den;
            if (std::abs(k) <= eps)
                continue;
            int s = k > 0.0 ? 1 : -1;
            if (sign != 0 && s != sign)
                return true;
            sign = s;
        }
    }
    return false;
}

// 非 U-turn 转向允许的单向绕转跨度上限（弧度）：200°。
// 跨度定义与取舍理由见 curve_utils.h 的 curveTurningSpan：有符号累加的 max-min，
// 因此合法的 S 形换道弧只贡献较大的单侧弯角，不会被当成折返。
// 实测正常转向的跨度 ≤ ~150°（110003285 的换道弧 14 为 60.2°、左转扇出族为 0°），
// 而折返型病态曲线同向绕满一整圈到 448°，阈值取在两者中间。
//
// 它补上了 curveHasCurvatureSignFlip 的盲区：后者每段只取 t = i/20（i∈[1,20)）
// 的 20 个内部点，且用 eps=0.10 过滤，因此绕圈过程中曲率始终保持同号时完全
// 看不到。共享端点转向对修复曾因此产出过总转角 448°、首尾切向差只有 88°
// 的两段曲线（局部曲率半径 0.5 m，实际不可行驶）。
static const double kMaxNonUTurnTurningSpan = 200.0 * M_PI / 180.0;

// 校验非 U-turn 单段转向的弧长、曲率和曲率符号约束。
// 短急弯采用受限例外；其它转向必须保持可见单拱形态且不能形成尖钩或 S 弯。
// 直行/转向分界必须与最终形态审计一致：`evaluateOrdinaryShape` 与全部
// `straight_like` 判定都用 0.25，此前本函数内部却用 0.35 放行，导致
// turn_strength 落在 (0.25, 0.35] 的换道弧在运行时门禁被无条件通过、
// 却在最终审计里按转向判定——`evaluateOrdinaryShape` 已经按 0.25 分派到
// 本函数的那一支因此完全失效。统一到 0.25。
static bool isNonUTurnTurnShapeAcceptable(
    const BezierCurve& curve, double chord_len, double turn_strength) {
    if (chord_len < 1e-6 || turn_strength < 0.25)
        return true;
    // 绕转判据必须在短急弯例外之前生效：同向绕过大半圈的曲线不论长短都不可行驶。
    if (!isgLegacyExcess() &&
        curveTurningSpan(curve) > kMaxNonUTurnTurningSpan)
        return false;
    double arc_chord = curve.arcLength() / chord_len;
    // 极短的单侧右转可以用局部曲率峰值较高的单段三次曲线合法表达。
    // 该例外仅覆盖短曲线，较长转向仍必须形成可见的平滑弧，不能坍缩成尖钩。
    if (curve.numSegments() == 1 && chord_len <= 10.0 &&
        curve.arcLength() <= 10.5 && arc_chord >= 1.005 &&
        arc_chord <= 1.08 && curve.maxCurvature(40) <= 6.0 &&
        !curveHasCurvatureSignFlip(curve))
        return true;
    // 普通左/右转不能为避障或排序被压成近直线。中长距离转向要求
    // 保留更明显的弧长余量，短小转向则放宽以避免过约束。
    // 下限按端点转角折算（见 ordinaryTurnArcChordFloor）：常数是按 90° 写的，
    // 浅转时超过同转角圆弧的理论比例，会把合规单拱判成"被压平"。
    double min_arc = ordinaryTurnArcChordFloor(
        curveEndpointTurnAngle(curve), chord_len);
    // A short turn next to a Boundary chain can need a slightly longer smooth
    // arch to remain on the legal side of the edge. Long turns keep the
    // original 1.35 limit and all turns still pass curvature and sign checks.
    const double max_arc = chord_len <= 15.0 ? 1.40 : 1.35;
    if (arc_chord < min_arc || arc_chord > max_arc)
        return false;
    if (curve.maxCurvature(20) > 2.5)
        return false;
    return !curveHasCurvatureSignFlip(curve);
}

// 校验非 U-turn 多段转向是否仍是首尾直行、单中弧的连续拱形。
// 除首尾段直线性和 G1 方向外，还要求中弧具有有效曲率且全曲线无符号翻转。
static bool isNonUTurnArchShapeAcceptable(
    const BezierCurve& curve, double chord_len, double turn_strength,
    const Vec2d& entry_tan, const Vec2d& exit_tan) {
    if (curve.empty() || chord_len < 1e-6 || turn_strength <= 0.35 ||
        curveHasCurvatureSignFlip(curve))
        return false;
    if (curve.numSegments() == 1)
        return isNonUTurnTurnShapeAcceptable(curve, chord_len, turn_strength);
    if (curve.numSegments() != 3)
        return false;

    const BezierSegment& first = curve.segs.front();
    const BezierSegment& arc = curve.segs[1];
    const BezierSegment& last = curve.segs.back();
    Vec2d first_dir = first.ctrl[3] - first.ctrl[0];
    Vec2d last_dir = last.ctrl[3] - last.ctrl[0];
    if (first_dir.norm() < 1e-6 || last_dir.norm() < 1e-6 ||
        !segmentLooksStraight(first) || !segmentLooksStraight(last) ||
        first_dir.normalized().dot(entry_tan.normalized()) < 0.98 ||
        last_dir.normalized().dot(exit_tan.normalized()) < 0.98 ||
        arc.maxCurvature(30) <= 0.03)
        return false;
    return true;
}

// 计算曲线采样点相对端点弦线的最大横向偏移与弦长之比。
// 结果用于区分直行样式和明显偏转样式。
static double maxLateralChordDeviationRatio(
    const BezierCurve& curve, double chord_len) {
    if (chord_len < 1e-6 || curve.empty())
        return std::numeric_limits<double>::infinity();
    Vec2d start = curve.startPt();
    Vec2d chord = curve.endPt() - start;
    if (chord.norm() < 1e-6)
        return std::numeric_limits<double>::infinity();
    Vec2d dir = chord.normalized();
    double max_dev = 0.0;
    for (const auto& pt : curve.sampleByArcLength(96))
        max_dev = std::max(max_dev, std::abs(cross2d(dir, pt - start)));
    return max_dev / chord_len;
}

// 判断曲线是否满足宽松直行形态，包括弧弦比、横向偏移和最大曲率限制。
static bool isStraightLikeShapeAcceptable(
    const BezierCurve& curve, double chord_len) {
    if (chord_len < 1e-6 || curve.empty())
        return false;
    double arc_chord = curve.arcLength() / chord_len;
    return arc_chord <= 1.08 &&
           maxLateralChordDeviationRatio(curve, chord_len) <= 0.08 &&
           curve.maxCurvature(40) <= 3.0;
}

// 严格同簇闭包在单段域无解时允许的最小多段普通表达。它只用于长距离
// 近直行连接的有限路点绕行；端点切向、内部 G1、弧长、曲率和单向形态仍
// 必须满足，正常普通连接继续使用单段 cubic。
//
// 直行分支的曲率上限取 R>=4m（口径与 shape_constraint.cpp 的
// kTwoSegmentStraightMaxCurvature 一致）：原值 3.0 允许 0.33m 半径，对
// 需求 1.4① 等于没有约束，100000385-u 的直行 41 就是这样折出 R=1.37m 接点的。
static constexpr double kTwoSegmentStraightMaxCurvature = 0.25;

// 物理修复搜索里"折形"的打分权重（乘 curveChordBudgetOvershoot 的返回比例）。
//
// Boundary / Obstacle / Fence 修复的候选网格里含 alpha 0.55~1.30，这些候选的
// 控制多边形沿弦倒序、曲线被折成 Z 形/L 形。修复搜索原先只按"物理是否安全 +
// 交叉数 + 弧长 + 曲率"打分，完全看不见折形，于是只要折形候选物理安全就可能胜出
// （100000385-u 的直行 55 折了 55.5% 弦长、左转 21 的两段各折 24.1%/16.2%）。
//
// 折形不能一律淘汰：左转 21 的折形两段候选是唯一能绕开 Boundary 的表达，淘汰它
// 只是把形态问题换成更严重的物理违规。所以折形按幅度计惩罚——权重取 50，
// 使它压过同一打分里的弧长（0.03/m，长曲线约 1.5）与曲率（约 0.1~1）项，
// 但远低于交叉数项（1000/条）：能不折就不折，物理上必须折才折。
static constexpr double kChordBudgetOvershootPenalty = 50.0;

static bool isTwoSegmentOrdinaryWaypointShapeAcceptable(
    const BezierCurve& curve, const Vec2d& p0, const Vec2d& t0,
    const Vec2d& p1, const Vec2d& t1, double chord_len,
    double turn_strength) {
    if (curve.empty() || curve.numSegments() != 2 || chord_len < 1e-6 ||
        t0.norm() < 1e-8 || t1.norm() < 1e-8)
        return false;
    const Vec2d start_tan = curve.startTan();
    const Vec2d end_tan = curve.endTan();
    if (start_tan.norm() < 1e-8 || end_tan.norm() < 1e-8 ||
        start_tan.normalized().dot(t0.normalized()) < 0.95 ||
        end_tan.normalized().dot(t1.normalized()) < 0.95)
        return false;
    const BezierSegment& first = curve.segs.front();
    const BezierSegment& second = curve.segs.back();
    const Vec2d left_tan = first.evalDeriv1(1.0);
    const Vec2d right_tan = second.evalDeriv1(0.0);
    if (left_tan.norm() < 1e-8 || right_tan.norm() < 1e-8 ||
        left_tan.normalized().dot(right_tan.normalized()) < 0.98 ||
        (first.ctrl[3] - first.ctrl[0]).norm() < 0.30 ||
        (second.ctrl[3] - second.ctrl[0]).norm() < 0.30)
        return false;
    // 每段控制多边形必须沿自身弦单调，理由同单段的弦方向联合预算
    // （见 cubicControlPolygonMonotone）：两段各自把手吃满自身弦时同样
    // 会折出 Z 形控制多边形。
    if (!cubicControlPolygonMonotone(first) ||
        !cubicControlPolygonMonotone(second))
        return false;
    const double ratio = curve.arcLength() / chord_len;
    double max_lateral = 0.0;
    const Vec2d chord_dir = (p1 - p0).normalized();
    for (const Vec2d& point : curve.sampleByArcLength(96))
        max_lateral = std::max(
            max_lateral, std::abs(cross2d(chord_dir, point - p0)));
    if (turn_strength <= 0.25)
        return ratio <= 1.12 && max_lateral / chord_len <= 0.12 &&
               curve.maxCurvature(40) <= kTwoSegmentStraightMaxCurvature;
    return ratio >= ordinaryTurnArcChordFloor(
                        curveEndpointTurnAngle(curve), chord_len) &&
           ratio <= 1.45 && max_lateral / chord_len <= 0.45 &&
           curve.maxCurvature(40) <= 2.5 &&
           curveTurningSpan(curve) <= 200.0 * M_PI / 180.0 &&
           !curveHasCurvatureSignFlip(curve);
}

// 判断曲线是否满足基础直行的严格形态限制。
// 在宽松直行判定之上进一步限制弧弦比和最大曲率，供自然直行恢复使用。
static bool isStrictStraightBaseShapeAcceptable(
    const BezierCurve& curve, double chord_len) {
    if (!isStraightLikeShapeAcceptable(curve, chord_len))
        return false;
    return curve.arcLength() / chord_len < 1.02 &&
           curve.maxCurvature(40) < 0.10;
}

// 判断非 U-turn 曲线是否存在需要形态恢复的风险，而不是简单判断是否不合格。
// 返回 true 表示弧长不足、曲率过高或曲率符号翻转等高风险情形。
static bool isNonUTurnTurnShapeRisk(
    const BezierCurve& curve, double chord_len, double turn_strength) {
    if (chord_len < 1e-6 || turn_strength <= 0.35)
        return false;
    if (isNonUTurnTurnShapeAcceptable(curve, chord_len, turn_strength))
        return false;
    double arc_chord = curve.arcLength() / chord_len;
    double min_arc = ordinaryTurnArcChordRestoreFloor(
        curveEndpointTurnAngle(curve), chord_len);
    return arc_chord < min_arc || curve.maxCurvature(20) > 2.0 ||
           curveHasCurvatureSignFlip(curve);
}

// 将二维向量绕原点旋转指定弧度并返回旋转后的向量。
static Vec2d rotateVec(const Vec2d& v, double radians) {
    double c = std::cos(radians);
    double s = std::sin(radians);
    return Vec2d(v.x() * c - v.y() * s, v.x() * s + v.y() * c);
}

struct CurveRisk {
    bool obstacle = false;
    bool boundary = false;
    bool fence = false;
    int sibling_crosses = 0;

    // 返回是否存在至少一种物理风险：障碍物、边界或围栏违规。
    // 同簇交叉不属于 physical，单独由 sibling_crosses 表示。
    bool physical() const {
        return obstacle || boundary || fence;
    }
};

// 将 Bezier 曲线写入 ConnectivityCurve，并按弧长重新生成输出折线几何。
// 提供固定几何时保留该几何，否则清空旧 geometry 后从曲线采样填充。
static void setConnectivityCurveGeometry(
    ConnectivityCurve& cc, const BezierCurve& curve, const LineString2d* fixed_geometry = nullptr) {
    cc.curve = std::make_shared<BezierCurve>(curve);
    cc.geometry.points.clear();
    if (fixed_geometry && fixed_geometry->points.size() >= 2) {
        cc.geometry = *fixed_geometry;
        return;
    }
    if (!curve.empty()) {
        int n = std::max(2, std::min(240, (int)std::ceil(curve.arcLength() / 0.3) + 1));
        cc.geometry.points = toVec3dArray(curve.sampleByArcLength(n));
    }
}

// 统一评估候选曲线的障碍物、边界、围栏和同簇交叉风险。
// U-turn 使用专门边界规则；include_fence 控制围栏检查，结果仅描述风险不修改曲线。
// allow_merge_funnel 只由三段式掉头兜底档传入 true，其余调用点保持严格判定。
static CurveRisk assessCurveRisk(
    const BezierCurve& curve, const IntersectionInput& input, const SDFField& sdf,
    const std::vector<SampledSiblingCurve>& sampled_siblings, bool include_fence = true,
    bool is_uturn = false, bool allow_merge_funnel = false) {
    CurveRisk risk;
    double ms = minSDFAlongCurveAdaptive(curve, sdf);
    risk.obstacle = curveIntersectsObstacles(curve, input.obstacles) || (ms < 0.0);

    // 对于U-turn，使用特殊的边界检查逻辑，优先考虑弧形边界穿越
    if (is_uturn) {
        // 使用特殊的U-turn边界检查，优先穿越弧形边界
        risk.boundary = curveIntersectsBoundariesUTurn(curve, input);
    } else {
        risk.boundary = curveIntersectsBoundaries(curve, input) ||
            curveRawIntersectsAnyBoundary(curve, input.boundaries, 0.15, 0.10);
    }

    if (roadEdgeClearanceViolation(curve, input.boundaries, input.mode) > 0.0)
        risk.boundary = true;
    risk.fence = include_fence && (!input.area.is_rough && curveLeavesFence(curve, input.area.geometry));
    risk.sibling_crosses =
        sampledSiblingCrossCount(curve, sampled_siblings, true, kClusterEndpointTol,
                                 false, allow_merge_funnel);
    return risk;
}

// 当当前曲线存在物理风险时，在有限把手长度集合中寻找安全的单段 G1 三次曲线。
// 候选必须消除障碍物、边界和围栏风险，并在 U-turn 拓扑场景下保持同簇无交叉；
// 成功时原地替换 curve，失败时保持输入不变并返回 false。
//
// 本函数是**捷径**而非修复搜索：它的语义是"一条朴素的单段三次曲线已经把物理问题
// 解决了，那就直接采用、跳过后续昂贵搜索"。因此折形候选必须在此被排除——沿弦倒序
// 的 Z/L 形控制多边形本身就不是"朴素的单段三次曲线"，而是审计会直接报
// `ordinary single cubic control points leave endpoint direction axes` 的形态违约。
// 100000385-u 的直行 55 正是这样产生的：grid 里 alpha=0.80 使两个把手各占弦长 0.8
// （弦 53.958、把手 43.167、弦向用量 83.923，超支 29.964 即弦长的 55.5%），它以
// "同簇交叉数最少"胜出（录取是 sibling_crosses 的字典序，打分惩罚只在同交叉数内
// 比较，因此加权惩罚对它无效）。
//
// 排除折形并不等于"在物理路径上一律禁止折形"：本函数返回 false 只是放弃捷径，
// 后续 tryObstacleBypassCandidate / 优化器 / tryBoundarySafeCandidate 依旧可以
// 输出折形曲线——它们的表达力更强（两侧把手可不等长、可用两段），实测 55 由此
// 得到单调且 Boundary 安全的非对称把手形态（h0 10.25 / h1 43.17，弦向富余
// 2.035m，maxk 0.0971）。真正"折形是唯一物理安全形态"的场景仍然保留：左转 21 在
// tryBoundarySafeCandidate 的 4649 个候选中只有一个 full_safe 候选，且它是折形的
// （两段分别超支弦长的 24.1% / 16.2%），那条路径按 curveChordBudgetOvershoot 打分
// 惩罚而不淘汰，因此不受本处影响。
static bool tryPhysicalSafeSingleCubic(
    const Vec2d& p0, const Vec2d& t0, const Vec2d& p1, const Vec2d& t1,
    const IntersectionInput& input, const SDFField& sdf,
    const std::vector<SampledSiblingCurve>& sampled_siblings, BezierCurve& curve) {
    CurveRisk current = assessCurveRisk(curve, input, sdf, sampled_siblings);
    if (!current.physical())
        return false;

    Vec2d T0 = t0.norm() > 1e-8 ? t0.normalized() : Vec2d(1, 0);
    Vec2d T1 = t1.norm() > 1e-8 ? t1.normalized() : Vec2d(1, 0);
    bool require_topology_clear = (T0.dot(T1) < -0.5);
    BezierCurve best;
    bool have_best = false;
    int best_cross = std::numeric_limits<int>::max();
    double best_shape = std::numeric_limits<double>::max();

    // Boundary 穿越修复有时需要比自然区间 0.30 至 0.50 更短的保守把手；
    // 较小的单段三次曲线更贴近车道间弦线，可避免切过短 RoadEdge 端帽。
    for (double alpha : {0.24, 0.18, 0.12, 0.30, 0.34, 0.38, 0.42,
                         0.46, 0.50, 0.55, 0.60, 0.70, 0.80}) {
        BezierCurve candidate;
        candidate.segs.push_back(makeCubicG1(p0, T0, p1, T1, alpha));
        // 弦方向联合预算：先做这一步（O(1)，无采样），既排除折形捷径，
        // 也省下近直行几何上必然折形的大 alpha 候选的风险评估开销。
        if (!cubicControlPolygonMonotone(candidate.segs.front()))
            continue;
        if (curveSelfIntersectsBusiness(candidate, 1.0))
            continue;
        CurveRisk risk = assessCurveRisk(candidate, input, sdf, sampled_siblings);
        if (risk.physical() || (require_topology_clear && risk.sibling_crosses > 0))
            continue;
        double shape_score = current.boundary
            ? std::abs(alpha - 0.24)
            : std::abs(alpha - 0.40);
        if (!have_best ||
            risk.sibling_crosses < best_cross ||
            (risk.sibling_crosses == best_cross && shape_score < best_shape)) {
            best = candidate;
            have_best = true;
            best_cross = risk.sibling_crosses;
            best_shape = shape_score;
        }
    }

    if (!have_best)
        return false;
    curve = best;
    return true;
}

// 当单段曲线形态发散、曲率过高或弧长异常时，搜索满足控制点轴向和形态约束的
// 单段候选。候选还需保持物理风险不变差以及同簇交叉数不增加。
static bool tryShapeSafeSingleCubic(
    const Vec2d& p0, const Vec2d& t0, const Vec2d& p1, const Vec2d& t1,
    const IntersectionInput& input, const SDFField& sdf,
    const std::vector<SampledSiblingCurve>& sampled_siblings, bool include_fence,
    BezierCurve& curve) {
    double chord_len = (p1 - p0).norm();
    if (chord_len < 1e-6)
        return false;
    bool debug = isgDebugUTurn();
    Vec2d T0 = t0.norm() > 1e-8 ? t0.normalized() : (p1 - p0).normalized();
    Vec2d chord_dir = (p1 - p0).normalized();
    double turn_strength = std::abs(cross2d(T0, chord_dir));
    double arc_chord = curve.arcLength() / chord_len;
    bool shape_bad = curve.maxCurvature(20) > 2.0 ||
        curve.arcLength() > std::max(chord_len * 1.8, chord_len + 12.0) ||
        (turn_strength > 0.35 && (arc_chord < 1.05 || curveHasCurvatureSignFlip(curve)));
    if (!shape_bad)
        return false;

    // 发散曲线允许用更保守的单段候选直接替换。
    bool severely_divergent = isCurveSeverelyDivergent(curve, p0, p1, chord_len);

    Vec2d T1 = t1.norm() > 1e-8 ? t1.normalized() : (p1 - p0).normalized();
    BezierCurve best = curve;
    CurveRisk best_risk = assessCurveRisk(best, input, sdf, sampled_siblings, include_fence);
    double best_score = 1000.0 * best_risk.sibling_crosses + best.maxCurvature(20) + 0.02 * best.arcLength();
    bool improved = false;

    const OrdinaryCurveInitializer ordinary_initializer;
    auto consider_candidate = [&](const BezierCurve& candidate) {
        if (candidate.empty())
            return;
        if (curveSelfIntersectsBusiness(candidate, 1.0))
            return;
        if (!ordinarySingleCubicControlsValid(
                candidate, p0, T0, p1, T1, 1e-5, true))
            return;
        if (turn_strength > 0.35 &&
            !isNonUTurnTurnShapeAcceptable(candidate, chord_len, turn_strength))
            return;
        CurveRisk risk = assessCurveRisk(candidate, input, sdf, sampled_siblings, include_fence);
        if (!severely_divergent && risk.physical())
            return;
        double penalty = severely_divergent ? 0.0 : 1000.0;
        double score = penalty * risk.sibling_crosses + candidate.maxCurvature(20) + 0.02 * candidate.arcLength();
        // 发散曲线没有保留价值，首个可用单段候选即可作为改进。
        if (severely_divergent && !improved) {
            best = candidate;
            best_score = score;
            improved = true;
            return;
        }
        if (score + 1e-6 < best_score) {
            best = candidate;
            best_score = score;
            improved = true;
            if (risk.sibling_crosses == 0 && candidate.maxCurvature(20) < 1.0)
                return;
        }
    };

    // 先尝试方向交点 2/3 升阶候选，允许首尾把手不对称；这是短急弯
    // 在方向交点很近时仍保持单拱形态的关键候选。
    consider_candidate(ordinary_initializer.buildPreferredSingleCubic(
        p0, T0, p1, T1));
    // 短急弯常出现一侧方向交点很近、另一侧仍有充足把手空间。
    // 对称 alpha 无法表达这种几何，补充有界的非对称首尾组合。
    const OrdinarySingleCubicHandleBounds handle_bounds =
        ordinarySingleCubicHandleBounds(p0, T0, p1, T1, true);
    const double preferred_start = handle_bounds.has_direction_intersection
        ? (2.0 / 3.0) * handle_bounds.direction_intersection_start
        : chord_len / 3.0;
    const double preferred_end = handle_bounds.has_direction_intersection
        ? (2.0 / 3.0) * handle_bounds.direction_intersection_end
        : chord_len / 3.0;
    const std::vector<double> start_handles = {
        std::min(handle_bounds.start_max, preferred_start),
        handle_bounds.start_max};
    const std::vector<double> end_handles = {
        std::min(handle_bounds.end_max, preferred_end),
        0.80 * handle_bounds.end_max, handle_bounds.end_max};
    for (double start_handle : start_handles) {
        for (double end_handle : end_handles) {
            BezierSegment segment;
            segment.ctrl[0] = p0;
            segment.ctrl[1] = p0 + T0 * std::max(
                handle_bounds.start_min, start_handle);
            segment.ctrl[2] = p1 - T1 * std::max(
                handle_bounds.end_min, end_handle);
            segment.ctrl[3] = p1;
            BezierCurve asymmetric;
            asymmetric.segs.push_back(segment);
            consider_candidate(asymmetric);
        }
    }
    // 再从常用把手长度到短把手递减搜索，作为有界回退。
    for (double alpha : {0.50, 0.42, 0.34, 0.26, 0.22, 0.18, 0.14, 0.12}) {
        consider_candidate(ordinary_initializer.buildSingleCubic(
            p0, T0, p1, T1, alpha));
    }

    if (!improved)
        return false;
    curve = best;
    return true;
}

// 普通左右转在无固定形态且单段候选无物理风险时，统一恢复为
// 单段自然弧；同簇风险只参与排序，不再作为保留多段S弯的豁免理由。
// 尝试把无固定形态的多段普通转向恢复为自然单段 G1 三次曲线。
// 仅在非 U-turn、无物理风险且转向明显时执行，并选择不增加同簇风险的最佳把手长度。
static bool tryNaturalSingleCubicIfSafe(
    const Vec2d& p0, const Vec2d& t0, const Vec2d& p1, const Vec2d& t1,
    const IntersectionInput& input, const SDFField& sdf,
    const std::vector<SampledSiblingCurve>& sampled_siblings, bool include_fence,
    BezierCurve& curve) {
    Vec2d chord = p1 - p0;
    double chord_len = chord.norm();
    if (chord_len < 1e-6)
        return false;
    Vec2d T0 = t0.norm() > 1e-8 ? t0.normalized() : chord.normalized();
    Vec2d T1 = t1.norm() > 1e-8 ? t1.normalized() : chord.normalized();
    if (T0.dot(T1) < -0.5)
        return false;

    CurveRisk current_risk = assessCurveRisk(
        curve, input, sdf, sampled_siblings, include_fence);
    if (current_risk.obstacle || current_risk.fence)
        return false;

    double turn_strength = std::abs(cross2d(T0, chord.normalized()));
    if (turn_strength <= 0.35)
        return false;

    BezierCurve best = curve;
    bool improved = false;
    double best_score = std::numeric_limits<double>::infinity();
    const OrdinaryCurveInitializer initializer;
    std::vector<BezierCurve> candidates;
    // 方向射线交点存在时，preferred candidate 将控制点放在交点前 1/3
    // （即从端点量取交点距离的 2/3），避免优化器把右转压成近直线。
    candidates.push_back(initializer.buildPreferredSingleCubic(p0, T0, p1, T1));
    for (double alpha : {0.42, 0.46, 0.50, 0.38, 0.34}) {
        BezierCurve candidate;
        candidate.segs.push_back(makeCubicG1(p0, T0, p1, T1, alpha));
        candidates.push_back(candidate);
    }
    for (std::size_t candidate_index = 0; candidate_index < candidates.size();
         ++candidate_index) {
        const BezierCurve& candidate = candidates[candidate_index];
        if (candidate.empty() || curveSelfIntersectsBusiness(candidate, 1.0))
            continue;
        if (!ordinarySingleCubicControlsValid(
                candidate, p0, T0, p1, T1, 1e-5, true))
            continue;
        if (!isNonUTurnTurnShapeAcceptable(candidate, chord_len, turn_strength))
            continue;
        double pure_min_arc = ordinaryTurnArcChordRestoreFloor(
            curveEndpointTurnAngle(candidate), chord_len);
        if (candidate.arcLength() / chord_len < pure_min_arc)
            continue;
        CurveRisk risk = assessCurveRisk(
            candidate, input, sdf, sampled_siblings, include_fence);
        if (risk.physical())
            continue;
        if (risk.sibling_crosses > current_risk.sibling_crosses)
            continue;
        double arc_chord = candidate.arcLength() / chord_len;
        double score = 1000.0 * risk.sibling_crosses +
                       10.0 * std::abs(arc_chord - 1.10) +
                       candidate.maxCurvature(40) + 0.02 * candidate.arcLength();
        // 稳定的方向交点候选是普通转向的规范形态：控制点从端点量取
        // 交点距离的 2/3（交点前 1/3）。只要它通过全部物理、形态和同簇
        // 门禁，就不能被弧长评分偏好的对称短把手覆盖。
        if (candidate_index == 0) {
            best = candidate;
            curve = candidate;
            return true;
        }
        if (!improved || score + 1e-6 < best_score) {
            best = candidate;
            best_score = score;
            improved = true;
        }
    }
    if (!improved)
        return false;
    curve = best;
    return true;
}

// 尝试把直行或近直行曲线恢复为严格直行样式的自然单段三次曲线。
// 候选必须满足端点切向、严格直行形态、物理安全和同簇交叉不恶化条件。
static bool tryNaturalStraightSingleCubicIfSafe(
    const Vec2d& p0, const Vec2d& t0, const Vec2d& p1, const Vec2d& t1,
    const IntersectionInput& input, const SDFField& sdf,
    const std::vector<SampledSiblingCurve>& sampled_siblings, bool include_fence,
    BezierCurve& curve) {
    Vec2d chord = p1 - p0;
    double chord_len = chord.norm();
    if (chord_len < 1e-6)
        return false;
    Vec2d T0 = t0.norm() > 1e-8 ? t0.normalized() : chord.normalized();
    Vec2d T1 = t1.norm() > 1e-8 ? t1.normalized() : chord.normalized();
    // 声明直行在斜交路口中入口切向可能与端点弦有中等夹角；允许进入
    // 候选搜索，但最终候选仍必须满足严格直行形态且不能增加同簇交叉。
    if (std::abs(cross2d(T0, chord.normalized())) >= 0.40)
        return false;

    CurveRisk current_risk = assessCurveRisk(
        curve, input, sdf, sampled_siblings, include_fence);
    if (current_risk.physical())
        return false;
    if (curve.numSegments() <= 1 &&
        isStraightLikeShapeAcceptable(curve, chord_len))
        return false;

    BezierCurve best;
    bool have_best = false;
    double best_score = std::numeric_limits<double>::infinity();
    bool debug_straight_repair =
        std::getenv("ISG_DEBUG_STRAIGHT_REPAIR") != nullptr;
    int dbg_total = 0;
    int dbg_shape = 0;
    int dbg_phys = 0;
    int dbg_cross = 0;
    int dbg_obstacle = 0;
    int dbg_boundary = 0;
    int dbg_fence = 0;
    auto consider_single = [&](const BezierCurve& candidate) {
        ++dbg_total;
        if (candidate.empty() ||
            !isStrictStraightBaseShapeAcceptable(candidate, chord_len)) {
            ++dbg_shape;
            return;
        }
        CurveRisk candidate_risk = assessCurveRisk(
            candidate, input, sdf, sampled_siblings, include_fence);
        if (candidate_risk.physical()) {
            ++dbg_phys;
            dbg_obstacle += candidate_risk.obstacle ? 1 : 0;
            dbg_boundary += candidate_risk.boundary ? 1 : 0;
            dbg_fence += candidate_risk.fence ? 1 : 0;
            return;
        }
        if (candidate_risk.sibling_crosses > current_risk.sibling_crosses) {
            ++dbg_cross;
            return;
        }
        double score = 1000.0 * candidate_risk.sibling_crosses +
                       candidate.maxCurvature(40) + 0.02 * candidate.arcLength();
        if (!have_best || score < best_score) {
            best = candidate;
            best_score = score;
            have_best = true;
        }
    };
    if (curve.numSegments() > 1) {
        BezierSegment seg;
        seg.ctrl[0] = p0;
        seg.ctrl[1] = curve.segs.front().ctrl[1];
        seg.ctrl[2] = curve.segs.back().ctrl[2];
        seg.ctrl[3] = p1;
        BezierCurve candidate;
        candidate.segs.push_back(seg);
        Vec2d st = candidate.startTan();
        Vec2d et = candidate.endTan();
        if (st.norm() > 1e-8 && et.norm() > 1e-8 &&
            st.normalized().dot(T0) > 0.99 &&
            et.normalized().dot(T1) > 0.99)
            consider_single(candidate);
    }
    for (double alpha : {0.08, 0.12, 0.18, 0.24, 0.30, 0.36, 0.40, 0.46, 0.52, 0.60}) {
        BezierCurve candidate;
        candidate.segs.push_back(makeCubicG1(p0, T0, p1, T1, alpha));
        consider_single(candidate);
    }
    if (!have_best) {
        if (debug_straight_repair) {
            fprintf(stderr,
                    "[STRAIGHT-NATURAL] reject p0=(%.2f,%.2f) p1=(%.2f,%.2f) current_cross=%d total=%d shape=%d phys=%d cross=%d phys_obs=%d phys_bnd=%d phys_fence=%d\n",
                    p0.x(), p0.y(), p1.x(), p1.y(),
                    current_risk.sibling_crosses,
                    dbg_total, dbg_shape, dbg_phys, dbg_cross,
                    dbg_obstacle, dbg_boundary, dbg_fence);
        }
        return false;
    }
    curve = best;
    if (debug_straight_repair) {
        fprintf(stderr,
                "[STRAIGHT-NATURAL] accept p0=(%.2f,%.2f) p1=(%.2f,%.2f) current_cross=%d score=%.3f\n",
                p0.x(), p0.y(), p1.x(), p1.y(),
                current_risk.sibling_crosses, best_score);
    }
    return true;
}

// 针对已确认的 Boundary 违规曲线执行有界候选搜索。
// 搜索包含端点开口、保守折点和边界侧方候选；每个候选均经过形态、物理风险、
// 同簇交叉和真实边界相交门禁，成功后写回 curve。
static bool tryBoundarySafeCandidate(
    const Vec2d& p0, const Vec2d& t0, const Vec2d& p1, const Vec2d& t1,
    const IntersectionInput& input, const SDFField& sdf,
    const std::vector<SampledSiblingCurve>& sampled_siblings, bool include_fence,
    BezierCurve& curve, bool force_search = false) {
    CurveRisk current = assessCurveRisk(curve, input, sdf, sampled_siblings, include_fence);
    if (!current.physical() && !force_search)
        return false;

    Vec2d chord = p1 - p0;
    double chord_len = chord.norm();
    if (chord_len < 1e-6)
        return false;
    Vec2d chord_dir = chord.normalized();
    Vec2d T0 = t0.norm() > 1e-8 ? t0.normalized() : chord_dir;
    Vec2d T1 = t1.norm() > 1e-8 ? t1.normalized() : chord_dir;
    Vec2d center = boundarySafetyCenter(input);
    Vec2d mid0 = 0.5 * (p0 + p1);
    Vec2d inward = center - mid0;
    if (inward.norm() < 1e-8)
        return false;
    inward.normalize();

    double turn_strength = std::abs(cross2d(T0, chord_dir));
    bool straight_like = turn_strength < 0.25;
    BezierCurve best = curve;
    int best_cross = std::numeric_limits<int>::max();
    double best_score = std::numeric_limits<double>::infinity();
    bool have_best = false;
    int current_cross =
        sampledSiblingCrossCount(curve, sampled_siblings, true, kClusterEndpointTol);
    int debug_considered = 0;
    int debug_endpoint_hook = 0;
    int debug_shape_reject = 0;
    int debug_cross_reject = 0;
    int debug_min_cross_reject = std::numeric_limits<int>::max();
    int debug_endpoint_physical_safe = 0;
    int debug_endpoint_min_cross = std::numeric_limits<int>::max();
    double debug_endpoint_min_boundary_score = std::numeric_limits<double>::infinity();
    int debug_full_safe = 0;
    std::vector<Boundary> violated_boundaries;
    BoundingBox2d repair_box = curve.bbox();
    repair_box.min_pt -= Vec2d(2.0, 2.0);
    repair_box.max_pt += Vec2d(2.0, 2.0);
    for (const auto& bnd : input.boundaries) {
        if (!bnd.geometry.bbox().intersects(repair_box))
            continue;
        std::vector<Boundary> one_boundary{bnd};
        bool raw_through = curveRawIntersectsAnyBoundary(
            curve, one_boundary, 0.15, 0.10);
        BoundarySafetyResult bs = curveBoundarySafety(
            curve, one_boundary, center, 128);
        if (raw_through || bs.intersects || bs.outside_road_edge ||
            (force_search &&
             (curveRawIntersectsRoadEdge(curve, one_boundary) ||
              segmentRawIntersectsRoadEdge(p0, p1, one_boundary) ||
              segmentCorridorTouchesRoadEdge(p0, p1, one_boundary))))
            violated_boundaries.push_back(bnd);
    }
    // 仅当当前曲线确有 Boundary 违规时才启动 Boundary 修复搜索。
    // Obstacle、Fence 和纯形态风险各有专用修复路径；送入宽范围 Boundary 网格
    // 只会增加数千个高成本候选，无法改善约束结果。
    if (violated_boundaries.empty())
        return false;
    auto boundaryViolationScore = [&](const BezierCurve& c,
                                      const std::vector<Boundary>& bnds) {
        double score = 0.0;
        for (const auto& b : bnds) {
            const std::vector<Boundary> one_boundary{b};
            BoundarySafetyResult bs = curveBoundarySafety(
                c, one_boundary, center, 128, 0.15);
            if (curveRawIntersectsAnyBoundary(
                    c, one_boundary, 0.15, 0.10))
                score += 10000.0;
            if (bs.intersects)
                score += 10000.0;
            if (bs.outside_road_edge)
                score += 1.0 + bs.outside_penalty;
        }
        return score;
    };
    auto debugViolatedIds = [&]() {
        std::string ids;
        for (const auto& b : violated_boundaries) {
            if (!ids.empty())
                ids += ",";
            ids += b.id;
        }
        return ids;
    };
    auto consider = [&](const BezierCurve& candidate, double extra_bias,
                        bool endpoint_opening_hook = false) {
        ++debug_considered;
        if (endpoint_opening_hook)
            ++debug_endpoint_hook;
        if (candidate.empty() || curveSelfIntersectsBusiness(candidate, 1.0))
            return;
        double arc_chord = candidate.arcLength() / chord_len;
        double straight_arc_limit = endpoint_opening_hook ? 1.12 : 1.06;
        if (straight_like && arc_chord > straight_arc_limit) {
            ++debug_shape_reject;
            return;
        }
        if (!straight_like && (arc_chord < 1.02 || arc_chord > 1.45)) {
            ++debug_shape_reject;
            return;
        }
        double maxk = candidate.maxCurvature(40);
        double maxk_limit = straight_like ? 1.0 : 2.5;
        if (maxk > maxk_limit) {
            ++debug_shape_reject;
            return;
        }
        CurveRisk risk = assessCurveRisk(candidate, input, sdf, sampled_siblings, include_fence);
        int allowed_cross = (endpoint_opening_hook && !straight_like)
            ? std::max(current_cross, 30)
            : current_cross;
        if (risk.physical()) {
            if (endpoint_opening_hook && risk.boundary && !violated_boundaries.empty()) {
                debug_endpoint_min_boundary_score = std::min(
                    debug_endpoint_min_boundary_score,
                    boundaryViolationScore(candidate, violated_boundaries));
            }
            return;
        }
        if (endpoint_opening_hook) {
            ++debug_endpoint_physical_safe;
            debug_endpoint_min_cross =
                std::min(debug_endpoint_min_cross, risk.sibling_crosses);
        }
        if (risk.sibling_crosses > allowed_cross) {
            ++debug_cross_reject;
            debug_min_cross_reject = std::min(debug_min_cross_reject, risk.sibling_crosses);
            return;
        }
        if (!violated_boundaries.empty() &&
            boundaryViolationScore(candidate, violated_boundaries) > 1e-6)
            return;
        ++debug_full_safe;
        double target_arc = straight_like ? 1.03 : 1.12;
        double shape_k = endpoint_opening_hook && straight_like
            ? std::min(maxk, 0.12)
            : maxk;
        double score = 1000.0 * risk.sibling_crosses
                     + 8.0 * std::abs(arc_chord - target_arc)
                     + shape_k + 0.03 * candidate.arcLength()
                     + kChordBudgetOvershootPenalty *
                           curveChordBudgetOvershoot(candidate)
                     + 0.05 * extra_bias;
        if (!have_best ||
            risk.sibling_crosses < best_cross ||
            (risk.sibling_crosses == best_cross && score < best_score)) {
            best = candidate;
            best_cross = risk.sibling_crosses;
            best_score = score;
            have_best = true;
        }
    };

    auto add_endpoint_opening_candidates = [&]() {
        constexpr double endpoint_match_tol = 0.35;
        std::vector<Vec2d> entry_open_knots;
        std::vector<Vec2d> exit_open_knots;
        for (const auto& bnd : violated_boundaries) {
            std::vector<Vec2d> pts = toVec2dArray(bnd.geometry.points);
            if (pts.size() < 2)
                continue;
            const Vec2d& first = pts.front();
            const Vec2d& last = pts.back();
            auto emit_exit_opening = [&](const Vec2d& other_endpoint) {
                Vec2d opening = other_endpoint - p1;
                if (opening.norm() < 1e-6)
                    return;
                for (double scale : {1.15, 1.35, 1.60, 2.00}) {
                    Vec2d knot = p1 + opening * scale;
                    exit_open_knots.push_back(knot);
                    Vec2d to_knot = knot - p0;
                    Vec2d to_exit = p1 - knot;
                    std::vector<Vec2d> mid_tans;
                    if (to_exit.norm() > 1e-8)
                        mid_tans.push_back(to_exit.normalized());
                    if (to_knot.norm() > 1e-8 && to_exit.norm() > 1e-8) {
                        Vec2d bis = to_knot.normalized() + to_exit.normalized();
                        if (bis.norm() > 1e-8)
                            mid_tans.push_back(bis.normalized());
                    }
                    mid_tans.push_back(T1);
                    for (const auto& mt : mid_tans) {
                        for (double alpha : {0.005, 0.01, 0.02, 0.04, 0.08, 0.12, 0.20}) {
                            BezierCurve candidate;
                            candidate.segs.push_back(makeCubicG1(p0, T0, knot, mt, alpha));
                            candidate.segs.push_back(makeCubicG1(knot, mt, p1, T1, alpha));
                            consider(candidate, 2.0 + scale, true);
                        }
                    }
                }
            };
            auto emit_entry_opening = [&](const Vec2d* other_endpoint) {
                if (other_endpoint) {
                    Vec2d opening = *other_endpoint - p0;
                    if (opening.norm() > 1e-6) {
                        for (double scale : {1.00, 1.15, 1.35, 1.60}) {
                            Vec2d knot = p0 + opening * scale;
                            entry_open_knots.push_back(knot);
                            Vec2d to_knot = knot - p0;
                            Vec2d to_exit = p1 - knot;
                            std::vector<Vec2d> mid_tans;
                            if (to_exit.norm() > 1e-8)
                                mid_tans.push_back(to_exit.normalized());
                            if (to_knot.norm() > 1e-8 && to_exit.norm() > 1e-8) {
                                Vec2d bis = to_knot.normalized() + to_exit.normalized();
                                if (bis.norm() > 1e-8)
                                    mid_tans.push_back(bis.normalized());
                            }
                            mid_tans.push_back(chord_dir);
                            for (const auto& mt : mid_tans) {
                                for (double alpha : {0.005, 0.01, 0.02, 0.04, 0.08, 0.12}) {
                                    BezierCurve candidate;
                                    candidate.segs.push_back(makeCubicG1(p0, T0, knot, mt, alpha));
                                    candidate.segs.push_back(makeCubicG1(knot, mt, p1, T1, alpha));
                                    consider(candidate, 1.5 + scale, true);
                                }
                            }
                        }
                    }
                }
                Vec2d normal(-T0.y(), T0.x());
                for (double lead : {1.5, 2.5, 3.5, 5.0}) {
                    for (double side : {-0.8, 0.0, 0.8}) {
                        Vec2d knot = p0 + T0 * lead + normal * side;
                        entry_open_knots.push_back(knot);
                        Vec2d to_knot = knot - p0;
                        Vec2d to_exit = p1 - knot;
                        std::vector<Vec2d> mid_tans;
                        if (to_exit.norm() > 1e-8)
                            mid_tans.push_back(to_exit.normalized());
                        if (to_knot.norm() > 1e-8 && to_exit.norm() > 1e-8) {
                            Vec2d bis = to_knot.normalized() + to_exit.normalized();
                            if (bis.norm() > 1e-8)
                                mid_tans.push_back(bis.normalized());
                        }
                        mid_tans.push_back(chord_dir);
                        for (const auto& mt : mid_tans) {
                            for (double alpha : {0.005, 0.01, 0.02, 0.04, 0.08, 0.12}) {
                                BezierCurve candidate;
                                candidate.segs.push_back(makeCubicG1(p0, T0, knot, mt, alpha));
                                candidate.segs.push_back(makeCubicG1(knot, mt, p1, T1, alpha));
                                consider(candidate, 1.0 + lead + std::abs(side), true);
                            }
                        }
                    }
                }
            };

            if ((p1 - first).norm() <= endpoint_match_tol)
                emit_exit_opening(last);
            if ((p1 - last).norm() <= endpoint_match_tol)
                emit_exit_opening(first);
            if ((p0 - first).norm() <= endpoint_match_tol)
                emit_entry_opening(&last);
            if ((p0 - last).norm() <= endpoint_match_tol)
                emit_entry_opening(&first);
        }
        for (const auto& entry_knot : entry_open_knots) {
            for (const auto& exit_knot : exit_open_knots) {
                if ((exit_knot - entry_knot).norm() < 1e-6)
                    continue;
                Vec2d mid_tan = (exit_knot - entry_knot).normalized();
                std::vector<Vec2d> pts = {p0, entry_knot, exit_knot, p1};
                std::vector<Vec2d> tans = {T0, mid_tan, mid_tan, T1};
                for (double alpha : {0.005, 0.01, 0.02, 0.04, 0.08}) {
                    consider(makeCurveFromKnots(pts, tans, alpha),
                             4.0 + (entry_knot - p0).norm() + (exit_knot - p1).norm(),
                             true);
                }
            }
        }
        if (straight_like && !exit_open_knots.empty()) {
            for (const auto& exit_knot : exit_open_knots) {
                for (double lead : {3.0, 5.0, 8.0, 11.0}) {
                    for (double offset : {0.0, 1.5, 3.0, 5.0, 7.0}) {
                        Vec2d entry_knot = p0 + T0 * lead + inward * offset;
                        if ((entry_knot - p0).norm() < 1e-6 ||
                            (exit_knot - entry_knot).norm() < 1e-6)
                            continue;
                        Vec2d mid_tan = (exit_knot - entry_knot).normalized();
                        std::vector<Vec2d> pts = {p0, entry_knot, exit_knot, p1};
                        std::vector<Vec2d> tans = {T0, mid_tan, mid_tan, T1};
                        for (double alpha : {0.005, 0.01, 0.02, 0.04, 0.08}) {
                            consider(makeCurveFromKnots(pts, tans, alpha),
                                     3.0 + lead + offset + (exit_knot - p1).norm(),
                                     true);
                        }
                    }
                }
            }
        }
    };

    auto add_conservative_two_segment_candidates = [&]() {
        Vec2d chord_perp(-chord_dir.y(), chord_dir.x());
        std::vector<Vec2d> dirs = {inward, chord_perp, -chord_perp, T0, -T1};
        for (double frac : {0.15, 0.25, 0.35, 0.50, 0.65, 0.75, 0.85}) {
            Vec2d base = p0 + chord * frac;
            for (const Vec2d& raw_dir : dirs) {
                if (raw_dir.norm() < 1e-8)
                    continue;
                Vec2d dir = raw_dir.normalized();
                for (double offset : {-4.0, -2.0, 0.0, 2.0, 4.0}) {
                    Vec2d knot = base + dir * offset;
                    if ((knot - p0).norm() < 0.30 ||
                        (p1 - knot).norm() < 0.30)
                        continue;
                    Vec2d mt = p1 - knot;
                    if (mt.norm() < 1e-8)
                        mt = chord_dir;
                    mt.normalize();
                    for (double alpha : {0.12, 0.18, 0.24, 0.30}) {
                        BezierCurve candidate;
                        candidate.segs.push_back(makeCubicG1(p0, T0, knot, mt, alpha));
                        candidate.segs.push_back(makeCubicG1(knot, mt, p1, T1, alpha));
                        consider(candidate, 0.6 + std::abs(offset) * 0.1, true);
                    }
                }
            }
        }
    };

    // 根据实际命中的 Boundary 折线构造局部法向绕行候选。全局 x/y 网格对斜向
    // 道路边缘的有效采样非常稀疏，尤其在边缘由多段折线组成时会反复落回原边缘；
    // 这里沿命中段的两侧法向取有限偏移，并只保留命中段相邻的局部顶点链，既能
    // 绕过鼻端/折点，又不会把单连接搜索扩成全局穷举。
    auto add_local_boundary_detour_candidates = [&]() {
        const bool dense_boundary_scene = input.boundaries.size() >= 20;
        const std::vector<double> offsets = dense_boundary_scene
            ? std::vector<double>{0.35, 0.60, 0.90, 1.30, 1.80}
            : std::vector<double>{0.30, 0.55, 0.85, 1.20, 1.70};
        auto vertex_normal = [](const std::vector<Vec2d>& points, int index,
                                int side) {
            Vec2d tangent(0, 0);
            if (index > 0)
                tangent += points[index] - points[index - 1];
            if (index + 1 < (int)points.size())
                tangent += points[index + 1] - points[index];
            if (tangent.norm() < 1e-8)
                tangent = Vec2d(1, 0);
            tangent.normalize();
            Vec2d normal(-tangent.y(), tangent.x());
            return side * normal;
        };
        for (const auto& boundary : violated_boundaries) {
            std::vector<Vec2d> points = toVec2dArray(boundary.geometry.points);
            if (points.size() < 2)
                continue;
            const int last = (int)points.size() - 1;
            for (int hit_segment = 0; hit_segment < last; ++hit_segment) {
                const Vec2d edge = points[hit_segment + 1] - points[hit_segment];
                if (edge.norm() < 1e-8)
                    continue;
                const int first_vertex = std::max(0, hit_segment - 1);
                const int last_vertex = std::min(last, hit_segment + 2);
                for (int side : {-1, 1}) {
                    for (double offset : offsets) {
                        std::vector<Vec2d> detour;
                        detour.reserve(last_vertex - first_vertex + 1);
                        for (int index = first_vertex; index <= last_vertex; ++index)
                            detour.push_back(
                                points[index] + offset * vertex_normal(points, index, side));
                        if (detour.empty())
                            continue;
                        std::vector<Vec2d> knots;
                        knots.reserve(detour.size() + 2);
                        knots.push_back(p0);
                        knots.insert(knots.end(), detour.begin(), detour.end());
                        knots.push_back(p1);
                        std::vector<Vec2d> tangents(knots.size(), chord_dir);
                        tangents.front() = T0;
                        tangents.back() = T1;
                        for (size_t index = 1; index + 1 < knots.size(); ++index) {
                            Vec2d tangent = knots[index + 1] - knots[index - 1];
                            tangents[index] = tangent.norm() > 1e-8
                                ? tangent.normalized() : chord_dir;
                        }
                        for (double alpha : {0.05, 0.10, 0.16, 0.22}) {
                            consider(makeCurveFromKnots(knots, tangents, alpha),
                                     offset + (double)(last_vertex - first_vertex), true);
                        }
                    }
                }
            }
        }
    };
    add_conservative_two_segment_candidates();
    if (have_best && best_cross == 0) {
        curve = best;
        return true;
    }
    add_endpoint_opening_candidates();
    add_local_boundary_detour_candidates();
    if (have_best && best_cross == 0) {
        curve = best;
        return true;
    }

    auto add_boundary_side_candidates = [&]() {
        for (const auto& bnd : violated_boundaries) {
            if (bnd.geometry.points.size() < 2)
                continue;
            BoundingBox2d box = bnd.geometry.bbox();
            if (box.empty())
                continue;

            std::vector<double> side_y = {
                box.max_pt.y() + 0.5, box.max_pt.y() + 1.0,
                box.max_pt.y() + 2.0, box.max_pt.y() + 3.5,
                box.max_pt.y() + 5.0, box.max_pt.y() + 7.0,
                box.min_pt.y() - 0.5, box.min_pt.y() - 1.0,
                box.min_pt.y() - 2.0, box.min_pt.y() - 3.5,
                box.min_pt.y() - 5.0, box.min_pt.y() - 7.0,
            };
            std::vector<double> xs = {
                box.min_pt.x() - 0.8,
                box.min_pt.x() + 0.8,
                box.min_pt.x() + 2.4,
                box.min_pt.x() + 3.4,
                0.5 * (box.min_pt.x() + box.max_pt.x()),
                box.max_pt.x() - 0.8,
                box.max_pt.x() + 0.8,
            };
            Vec2d mid_tan = chord_dir;
            std::vector<double> alphas = straight_like
                ? std::vector<double>{0.18, 0.24, 0.30, 0.36}
                : std::vector<double>{0.24, 0.30, 0.36, 0.42, 0.50};

            for (double y : side_y) {
                for (double x : xs) {
                    Vec2d knot(x, y);
                    for (double alpha : alphas) {
                        std::vector<Vec2d> pts = {p0, knot, p1};
                        std::vector<Vec2d> tans = {T0, mid_tan, T1};
                        consider(makeCurveFromKnots(pts, tans, alpha),
                                 std::abs(y - mid0.y()) + std::abs(x - mid0.x()));
                    }
                }

                std::vector<Vec2d> pair_knots = {
                    Vec2d(box.min_pt.x() - 0.8, y),
                    Vec2d(box.min_pt.x() + 0.8, y),
                    Vec2d(box.min_pt.x() + 2.4, y),
                    Vec2d(box.min_pt.x() + 3.4, y),
                    Vec2d(0.5 * (box.min_pt.x() + box.max_pt.x()), y),
                    Vec2d(box.max_pt.x() - 0.8, y),
                    Vec2d(box.max_pt.x() + 0.8, y),
                };
                std::sort(pair_knots.begin(), pair_knots.end(),
                          [&](const Vec2d& a, const Vec2d& b) {
                              return (a - p0).dot(chord_dir) < (b - p0).dot(chord_dir);
                          });
                for (int i = 0; i + 1 < (int)pair_knots.size(); ++i) {
                    Vec2d k0 = pair_knots[i];
                    Vec2d k1 = pair_knots[i + 1];
                    if ((k1 - k0).norm() < 1e-6)
                        continue;
                    Vec2d kt = (k1 - k0).normalized();
                    for (double alpha : alphas) {
                        std::vector<Vec2d> pts = {p0, k0, k1, p1};
                        std::vector<Vec2d> tans = {T0, kt, kt, T1};
                        consider(makeCurveFromKnots(pts, tans, alpha),
                                 2.0 + std::abs(y - mid0.y()) +
                                 (k0 - p0).norm() * 0.05 + (k1 - p1).norm() * 0.05);
                    }
                }
            }
        }
    };

    // 密集边界场景使用上面的局部折线候选即可；旧的全局网格在单路口内会
    // 枚举数千个候选，却无法针对斜向链提供新的几何信息。
    if (input.boundaries.size() < 20)
        add_boundary_side_candidates();

    for (double alpha : {0.30, 0.34, 0.38, 0.42, 0.46, 0.52, 0.60}) {
        BezierCurve single;
        single.segs.push_back(makeCubicG1(p0, T0, p1, T1, alpha));
        consider(single, 0.0);
    }

    Vec2d mid_tan_base = T0 + T1;
    if (mid_tan_base.norm() < 1e-8)
        mid_tan_base = chord_dir;

    for (double frac : {0.25, 0.35, 0.65}) {
        for (double offset : {-1.0, 0.0, 1.0, 2.0}) {
            Vec2d knot = p0 + chord * frac + inward * offset;
            Vec2d to_exit = p1 - knot;
            std::vector<Vec2d> mid_tans = {chord_dir, mid_tan_base.normalized()};
            if (to_exit.norm() > 1e-8)
                mid_tans.push_back(to_exit.normalized());
            for (const auto& mt : mid_tans) {
                std::vector<double> knot_alphas = straight_like
                    ? std::vector<double>{0.24, 0.30, 0.36}
                    : std::vector<double>{0.24, 0.30, 0.36, 0.45, 0.55, 0.65};
                for (double alpha : knot_alphas) {
                    std::vector<Vec2d> pts = {p0, knot, p1};
                    std::vector<Vec2d> tans = {T0, mt, T1};
                    consider(makeCurveFromKnots(pts, tans, alpha), offset);
                }
            }
        }
    }

    if (straight_like) {
        for (double tail : {0.20, 0.40, 0.60}) {
            Vec2d q_tail = p1 - T1 * tail;
            for (double frac : {0.70, 0.80, 0.88}) {
                for (double offset : {-6.0, -4.0, -2.0, 0.0, 2.0}) {
                    Vec2d knot = p0 + chord * frac + inward * offset;
                    for (double alpha0 : {0.16, 0.20, 0.24}) {
                        for (double alpha1 : {0.18, 0.30, 0.42}) {
                            BezierCurve candidate;
                            candidate.segs.push_back(makeCubicG1(p0, T0, knot, chord_dir, alpha0));
                            candidate.segs.push_back(makeCubicG1(knot, chord_dir, q_tail, T1, alpha0));
                            candidate.segs.push_back(makeCubicG1(q_tail, T1, p1, T1, alpha1));
                            consider(candidate, std::abs(offset) + tail);
                        }
                    }
                }
            }
        }
    }

    if (!straight_like) {
        Vec2d chord_perp(-chord_dir.y(), chord_dir.x());
        Vec2d entry_perp(-T0.y(), T0.x());
        Vec2d turn_tan = T0 + T1;
        if (turn_tan.norm() < 1e-8)
            turn_tan = chord_dir;
        turn_tan.normalize();
        for (double lead_frac : {0.22, 0.27, 0.32}) {
            double lead = std::max(3.0, chord_len * lead_frac);
            for (double side : {-1.0, 1.0}) {
                for (double lateral_frac : {0.04, 0.06, 0.08}) {
                    Vec2d knot = p0 + T0 * lead + entry_perp * (side * chord_len * lateral_frac);
                    Vec2d to_exit = p1 - knot;
                    std::vector<Vec2d> mid_tans = {chord_dir, turn_tan};
                    if (to_exit.norm() > 1e-8)
                        mid_tans.push_back(to_exit.normalized());
                    for (const auto& mt : mid_tans) {
                        for (double alpha0 : {0.30, 0.36, 0.42}) {
                            for (double alpha1 : {0.30, 0.36, 0.42}) {
                                BezierCurve candidate;
                                candidate.segs.push_back(makeCubicG1(p0, T0, knot, mt, alpha0));
                                candidate.segs.push_back(makeCubicG1(knot, mt, p1, T1, alpha1));
                                consider(candidate, lead_frac * 10.0 + lateral_frac * 10.0);
                            }
                        }
                    }
                }
            }
        }
        for (double frac : {0.30, 0.35, 0.42, 0.55, 0.65}) {
            for (double side : {-1.0, 1.0}) {
                for (double offset : {3.0, 5.0, 7.0, 9.0}) {
                    Vec2d knot = p0 + chord * frac + chord_perp * (side * offset);
                    std::vector<Vec2d> mid_tans = {turn_tan, chord_dir};
                    Vec2d to_exit = p1 - knot;
                    if (to_exit.norm() > 1e-8)
                        mid_tans.push_back(to_exit.normalized());
                    for (const auto& mt : mid_tans) {
                        for (double alpha : {0.36, 0.45, 0.55, 0.65}) {
                            std::vector<Vec2d> pts = {p0, knot, p1};
                            std::vector<Vec2d> tans = {T0, mt, T1};
                            consider(makeCurveFromKnots(pts, tans, alpha),
                                     offset + std::abs(frac - 0.45));
                        }
                    }
                }
            }
        }
    }

    std::vector<double> offsets = straight_like
        ? std::vector<double>{1.0, 2.0}
        : std::vector<double>{1.8, 4.0, 7.5};
    for (double offset : offsets) {
        Vec2d mid = mid0 + inward * offset;
        for (double alpha : {0.24, 0.34, 0.45}) {
            std::vector<Vec2d> pts = {p0, mid, p1};
            std::vector<Vec2d> tans = {T0, mid_tan_base.normalized(), T1};
            consider(makeCurveFromKnots(pts, tans, alpha), offset);

            double lead = std::max(2.0, std::min(8.0, chord_len * 0.20));
            Vec2d q0 = p0 + T0 * lead + inward * offset;
            Vec2d q1 = p1 - T1 * lead + inward * offset;
            std::vector<Vec2d> pts2 = {p0, q0, q1, p1};
            std::vector<Vec2d> tans2 = {T0, chord_dir, chord_dir, T1};
            consider(makeCurveFromKnots(pts2, tans2, alpha), offset);

            Vec2d exit_knot = p1 - T1 * lead + inward * offset;
            std::vector<Vec2d> pts3 = {p0, exit_knot, p1};
            std::vector<Vec2d> tans3 = {T0, chord_dir, T1};
            consider(makeCurveFromKnots(pts3, tans3, alpha), offset);

            Vec2d mid_exit = 0.55 * p0 + 0.45 * p1 + inward * (0.7 * offset);
            std::vector<Vec2d> pts4 = {p0, mid_exit, exit_knot, p1};
            std::vector<Vec2d> tans4 = {T0, chord_dir, chord_dir, T1};
            consider(makeCurveFromKnots(pts4, tans4, alpha), offset);
        }
    }
    if (!have_best) {
        if (isgDebugBoundaryRepair()) {
            fprintf(stderr,
                    "[BOUNDARY-REPAIR] failed p0=(%.2f,%.2f) p1=(%.2f,%.2f) considered=%d endpoint=%d shape_reject=%d cross_reject=%d full_safe=%d current_cross=%d violated=%zu ids=%s\n",
                    p0.x(), p0.y(), p1.x(), p1.y(),
                    debug_considered, debug_endpoint_hook, debug_shape_reject,
                    debug_cross_reject, debug_full_safe,
                    current_cross, violated_boundaries.size(), debugViolatedIds().c_str());
            if (debug_cross_reject > 0) {
                fprintf(stderr,
                        "[BOUNDARY-REPAIR] min rejected sibling crosses=%d\n",
                        debug_min_cross_reject);
            }
            if (debug_endpoint_hook > 0) {
                fprintf(stderr,
                        "[BOUNDARY-REPAIR] endpoint physical safe=%d min endpoint sibling crosses=%d\n",
                        debug_endpoint_physical_safe,
                        debug_endpoint_min_cross == std::numeric_limits<int>::max()
                            ? -1 : debug_endpoint_min_cross);
                fprintf(stderr,
                        "[BOUNDARY-REPAIR] min endpoint boundary score=%.6f\n",
                        std::isfinite(debug_endpoint_min_boundary_score)
                            ? debug_endpoint_min_boundary_score : -1.0);
            }
        }
        return false;
    }
    if (isgDebugBoundaryRepair()) {
        fprintf(stderr,
                "[BOUNDARY-REPAIR] accepted full p0=(%.2f,%.2f) p1=(%.2f,%.2f) considered=%d endpoint=%d shape_reject=%d cross_reject=%d full_safe=%d current_cross=%d best_cross=%d best_arc=%.3f best_maxk=%.3f violated=%zu ids=%s\n",
                p0.x(), p0.y(), p1.x(), p1.y(),
                debug_considered, debug_endpoint_hook, debug_shape_reject,
                debug_cross_reject, debug_full_safe,
                current_cross, best_cross, best.arcLength() / chord_len,
                best.maxCurvature(40), violated_boundaries.size(), debugViolatedIds().c_str());
    }
    curve = best;
    return true;
}

// 根据连接入口切向与端点弦线计算带符号的转向强度。
// 正负表示左右转方向，绝对值表示偏转程度；连接或弦线无效时返回零。
static double signedTurnStrengthOfConnId(const SceneView& scene, const ConnId& id) {
    const Connectivity* conn = scene.connectivity(id);
    if (!conn)
        return 0.0;
    const IntersectionInput& input = scene.input();
    auto entry = scene.entryFrame(conn->entry_lane_id);
    auto exit_ = scene.exitFrame(conn->exit_lane_id);
    Vec2d chord = exit_.first - entry.first;
    if (chord.norm() < 1e-8)
        return 0.0;
    Vec2d T0 = entry.second.norm() > 1e-8
        ? entry.second.normalized() : chord.normalized();
    return cross2d(T0, chord.normalized());
}

// ── 共享端点扇出族的内/外侧判定 ──────────────────────────────────
//
// 历史实现用 pair.ref_perp 上的 |横向偏移| 比较内外侧。但 pairRefPerp 对
// 共享端点对取的是"共享端点 → 两个自由端点中点"方向的左法向，两个自由
// 端点在该轴上的投影恒为等大反号（中点在轴上投影为零），于是
// |lateral_a| > |lateral_b| 完全由浮点噪声决定：同一个扇出族内会得出互相
// 矛盾的内外侧结论（如 4|6 判 4 为外、5|6 判 5 为外），把手顺序修复因此
// 可能把中间成员的共享侧把手抬到最外侧成员之上，反而制造非端点相交。
//
// 改用共享端点自身的参考系度量：自由端点到"共享端点 + 共享切向"射线的
// 垂距越大者为外侧。该量对整族一致、与配对枚举顺序无关，且与物理含义
// 吻合——垂距更大意味着转弯半径更大，共享侧把手必须更长才能保持嵌套。
static double sharedEndpointFanOutwardness(
    const Vec2d& shared_pt, const Vec2d& shared_dir, const Vec2d& free_pt) {
    if (shared_dir.norm() < 1e-8)
        return 0.0;
    return std::abs(cross2d(shared_dir.normalized(), free_pt - shared_pt));
}

// 共享侧把手长度顺序只在"自由端方向一致"的成员之间等价于扇出嵌套顺序。
// 自由端指向不同 arm 的两条曲线靠方向而不是把手长度分离，对它们施加把手
// 顺序会把无关成员（例如同一进入车道上的近直行）拉进排序链，破坏真正扇出
// 族的单调性。要求两个自由端切向近似同向后再排序。
static bool sharedEndpointFanDirectionsAligned(
    const Vec2d& free_dir_a, const Vec2d& free_dir_b, double min_dot = 0.90) {
    if (free_dir_a.norm() < 1e-8 || free_dir_b.norm() < 1e-8)
        return false;
    return free_dir_a.normalized().dot(free_dir_b.normalized()) >= min_dot;
}

// 为当前已有结果建立同簇约束邻接表，过滤掉没有曲线的连接。
// 该表供交叉计数和局部修复快速查找可能影响当前连接的兄弟曲线。
static std::unordered_map<ConnId, std::vector<ConnId>> constrainedNeighborMap(
    const std::vector<ConnectivityCurve>& results, const ClusterOrderSolver& cs) {
    std::unordered_set<ConnId> available_ids;
    for (const auto& result : results)
        if (result.curve)
            available_ids.insert(result.id);
    return RepairImpactClosureBuilder().neighbors(cs.pairs(), available_ids);
}

// 统计指定连接候选曲线与其所有约束邻居的禁止相交数量。
// 结构性交叉已由邻接构建规则排除，端点容差用于允许真实连接点相遇。
static int constrainedCrossCountForId(
    const ConnId& id, const BezierCurve& curve,
    const std::vector<ConnectivityCurve>& results,
    const std::unordered_map<ConnId, size_t>& result_idx,
    const std::unordered_map<ConnId, std::vector<ConnId>>& neighbors,
    double endpoint_tol = kClusterEndpointTol) {
    int count = 0;
    auto nit = neighbors.find(id);
    if (nit == neighbors.end())
        return 0;
    for (const auto& other : nit->second) {
        auto it = result_idx.find(other);
        if (it == result_idx.end())
            continue;
        const auto& other_cc = results[it->second];
        if (!other_cc.curve)
            continue;
        if (curvesHaveForbiddenSameClusterIntersection(
                curve, *other_cc.curve, endpoint_tol))
            ++count;
    }
    return count;
}

// 统计候选曲线与约束邻居中固定形态曲线的禁止相交数量。
// 该计数用于候选排序和固定形态保护，不改变邻居曲线。
static int constrainedFixedCrossCountForId(
    const ConnId& id, const BezierCurve& curve,
    const std::vector<ConnectivityCurve>& results,
    const std::unordered_map<ConnId, size_t>& result_idx,
    const std::unordered_map<ConnId, std::vector<ConnId>>& neighbors,
    double endpoint_tol = kClusterEndpointTol) {
    int count = 0;
    auto nit = neighbors.find(id);
    if (nit == neighbors.end())
        return 0;
    for (const auto& other : nit->second) {
        auto it = result_idx.find(other);
        if (it == result_idx.end())
            continue;
        const auto& other_cc = results[it->second];
        if (!other_cc.fixed_shape || !other_cc.curve)
            continue;
        if (curvesHaveForbiddenSameClusterIntersection(
                curve, *other_cc.curve, endpoint_tol))
            ++count;
    }
    return count;
}

// 只统计候选曲线与共享端点约束邻居的禁止相交数量。
// 结构性交叉不计入，适合判断共享入口/出口的局部扇出冲突。
static int constrainedSharedEndpointCrossCountForId(
    const ConnId& id, const BezierCurve& curve,
    const std::vector<ConnectivityCurve>& results,
    const std::unordered_map<ConnId, size_t>& result_idx,
    const std::unordered_map<ConnId, std::vector<ConnId>>& neighbors,
    const ClusterOrderSolver& cs, double endpoint_tol = kClusterEndpointTol) {
    int count = 0;
    auto nit = neighbors.find(id);
    if (nit == neighbors.end())
        return 0;
    for (const auto& other : nit->second) {
        if (!cs.isSharedEndpoint(id, other))
            continue;
        if (cs.exemptionOf(id, other) == CrossExemption::StructuralCross)
            continue;
        auto it = result_idx.find(other);
        if (it == result_idx.end())
            continue;
        const auto& other_cc = results[it->second];
        if (!other_cc.curve)
            continue;
        if (curvesHaveForbiddenSameClusterIntersection(
                curve, *other_cc.curve, endpoint_tol))
            ++count;
    }
    return count;
}

// 判断候选曲线是否会给当前连接带来"新的"禁止相交：只统计当前曲线与该
// 邻居尚未冲突、而候选曲线开始冲突的邻居。把手顺序修复不能用总交叉数
// 比较来把关（目标对本身可能已计入旧数），但必须禁止修复一个对的同时
// 在同族其它成员上引入新的违约。
static bool introducesNewConstrainedCrossForId(
    const ConnId& id, const BezierCurve& current, const BezierCurve& candidate,
    const std::vector<ConnectivityCurve>& results,
    const std::unordered_map<ConnId, size_t>& result_idx,
    const std::unordered_map<ConnId, std::vector<ConnId>>& neighbors,
    const std::unordered_set<ConnId>& ignored_ids,
    double endpoint_tol = kClusterEndpointTol) {
    auto nit = neighbors.find(id);
    if (nit == neighbors.end())
        return false;
    for (const auto& other : nit->second) {
        if (ignored_ids.count(other))
            continue;
        auto it = result_idx.find(other);
        if (it == result_idx.end())
            continue;
        const auto& other_cc = results[it->second];
        if (!other_cc.curve)
            continue;
        if (curvesHaveForbiddenSameClusterIntersection(
                current, *other_cc.curve, endpoint_tol)) {
            // 旧曲线已经与该邻居冲突时不能直接放过候选：端点带内的贴合
            // (endpoint_tol 只有 0.30m，同一连接点出发的同族成员几乎必然
            // 命中) 会把整个邻居屏蔽掉，于是候选可以在远离端点的中段公然
            // 穿过它。110003285 的 81 正是靠这个缺口越过 76。因此额外用
            // 审计级端点容差复核一次"中段互穿"是否被新引入。
            if (curvesIntersectBusiness(
                    candidate, *other_cc.curve, kMidSpanCrossEndpointTol) &&
                !curvesIntersectBusiness(
                    current, *other_cc.curve, kMidSpanCrossEndpointTol))
                return true;
            continue;
        }
        if (curvesHaveForbiddenSameClusterIntersection(
                candidate, *other_cc.curve, endpoint_tol))
            return true;
    }
    return false;
}

// 检查候选曲线是否与同簇约束邻居中任一单段曲线出现"控制多边形折线互穿"。
// 该判定只对共享端点的单段表达生效，是同簇非端点相交的几何前兆，用于在
// 把手顺序修复阶段提前剔除会破坏扇出嵌套的候选。
// 覆盖范围是完整控制折线 P0→P1→P2→P3，而不只是两条 P1→P2 中间连线：
// 一条曲线的中间连线穿过另一条曲线的首/尾把手连线同样代表互穿。
static bool breaksSharedEndpointControlPolygonNesting(
    const ConnId& id, const BezierCurve& candidate,
    const std::vector<ConnectivityCurve>& results,
    const std::unordered_map<ConnId, size_t>& result_idx,
    const std::unordered_map<ConnId, std::vector<ConnId>>& neighbors,
    const ClusterOrderSolver& cs,
    const std::unordered_set<ConnId>& ignored_ids) {
    if (candidate.numSegments() != 1)
        return false;
    auto nit = neighbors.find(id);
    if (nit == neighbors.end())
        return false;
    for (const auto& other : nit->second) {
        if (ignored_ids.count(other))
            continue;
        if (!cs.isSharedEndpoint(id, other))
            continue;
        if (cs.exemptionOf(id, other) == CrossExemption::StructuralCross)
            continue;
        auto it = result_idx.find(other);
        if (it == result_idx.end())
            continue;
        const auto& other_cc = results[it->second];
        if (!other_cc.curve)
            continue;
        if (sharedEndpointControlPolylinesCross(
                candidate, *other_cc.curve, kClusterEndpointTol))
            return true;
    }
    return false;
}

// 从当前连接的冲突邻居中提取第一个有效的主导横向参考法线。
// 无冲突或没有有效参考轴时返回零向量，供定向偏移修复决定是否放弃搜索。
static Vec2d dominantConstrainedRefPerpForId(
    const ConnId& id, const BezierCurve& curve,
    const std::vector<ConnectivityCurve>& results,
    const std::unordered_map<ConnId, size_t>& result_idx,
    const std::unordered_map<ConnId, std::vector<ConnId>>& neighbors,
    const ClusterOrderSolver& cs, double endpoint_tol = kClusterEndpointTol) {
    auto nit = neighbors.find(id);
    if (nit == neighbors.end())
        return Vec2d(0, 0);
    for (const auto& other : nit->second) {
        auto it = result_idx.find(other);
        if (it == result_idx.end())
            continue;
        const auto& other_cc = results[it->second];
        if (!other_cc.curve)
            continue;
        if (!curvesHaveForbiddenSameClusterIntersection(
                curve, *other_cc.curve, endpoint_tol))
            continue;
        Vec2d ref = cs.refPerpOf(id, other);
        if (ref.norm() > 1e-8)
            return ref.normalized();
    }
    return Vec2d(0, 0);
}

// 收集纯几何场景中需要重修的连接 ID：仅包括非结构性交叉且未被固定形态保护的
// 一方或双方，确保后续拓扑修复预算只作用于真实冲突集合。
static std::unordered_set<ConnId> collectPureGeometryRepairIds(
    const std::vector<ConnectivityCurve>& results, const ClusterOrderSolver& cs,
    const std::unordered_set<ConnId>& preserved_fixed_ids, double endpoint_tol = kClusterEndpointTol) {
    std::unordered_set<ConnId> repair_ids;
    auto idx = resultIndexById(results);
    for (const auto& p : cs.pairs()) {
        if (p.exempt == CrossExemption::StructuralCross)
            continue;
        auto ia = idx.find(p.id_a), ib = idx.find(p.id_b);
        if (ia == idx.end() || ib == idx.end())
            continue;
        const auto& ca = results[ia->second];
        const auto& cb = results[ib->second];
        if (!ca.curve || !cb.curve)
            continue;
        double pair_tol = endpoint_tol;
        if (!curvesHaveForbiddenSameClusterIntersection(
                *ca.curve, *cb.curve, pair_tol))
            continue;
        bool fixed_a = preserved_fixed_ids.count(p.id_a) > 0;
        bool fixed_b = preserved_fixed_ids.count(p.id_b) > 0;
        if (!fixed_a && !fixed_b) {
            repair_ids.insert(p.id_a);
            repair_ids.insert(p.id_b);
        } else if (!fixed_a) {
            repair_ids.insert(p.id_a);
        } else if (!fixed_b) {
            repair_ids.insert(p.id_b);
        }
    }
    return repair_ids;
}

// 找出彼此冲突的固定形态连接，并按端点弦长保留较长者、移除较短者的保护资格。
// 该策略为自然曲线恢复释放必要的修复空间。
static std::unordered_set<ConnId> collectConflictingPreservedFixedIds(
    const std::vector<ConnectivityCurve>& results, const ClusterOrderSolver& cs,
    const std::unordered_set<ConnId>& preserved_fixed_ids,
    double endpoint_tol = kClusterEndpointTol) {
    std::unordered_set<ConnId> drop_ids;
    auto idx = resultIndexById(results);
    for (const auto& p : cs.pairs()) {
        if (p.exempt == CrossExemption::StructuralCross)
            continue;
        if (!preserved_fixed_ids.count(p.id_a) ||
            !preserved_fixed_ids.count(p.id_b))
            continue;
        auto ia = idx.find(p.id_a), ib = idx.find(p.id_b);
        if (ia == idx.end() || ib == idx.end())
            continue;
        const auto& ca = results[ia->second];
        const auto& cb = results[ib->second];
        if (!ca.curve || !cb.curve)
            continue;
        double pair_tol = endpoint_tol;
        if (!curvesHaveForbiddenSameClusterIntersection(
                *ca.curve, *cb.curve, pair_tol))
            continue;
        double span_a = (ca.curve->endPt() - ca.curve->startPt()).norm();
        double span_b = (cb.curve->endPt() - cb.curve->startPt()).norm();
        drop_ids.insert(span_a <= span_b ? p.id_b : p.id_a);
    }
    return drop_ids;
}

// 为普通非 U-turn 连接构造不依赖环境约束的自然单段候选。
// 仅返回切向有效、无自交且弧长形态合理的候选，供固定形态阻塞预判使用。
//
// 候选还必须**物理可达**（不穿障、不穿 Boundary/RoadEdge、不出精围栏）。
// 本函数唯一的用途是回答"保留这条固有形态，是否会挡住另一条连接的自然单段形态"，
// 而牺牲固有形态的理由（架构设计文档 6.8.4：非 U 型固有形态原则上保留）只有
// "被牺牲后另一条真的能拿到自然单段形态"。若该自然候选本身就因物理原因不可用，
// 那条连接无论如何都会走避让路径、永远不会采用这个形态，此时撤销固有形态保护
// 是纯粹的损失：100000385-u 的右转 929910（fixed_shape=1，18 点输入折线）就是因为
// 同入的直行 59 的"自然直行候选"与它相交而被撤保护并重新生成，偏离输入形态 1.829m
// 且被压成 maxk 0.385（R 2.6m）的两段畸形；而 59 的自然单段候选实测 11 个全部因
// 物理风险被拒（tryNaturalStraightSingleCubicIfSafe 报 phys=11/11），59 最终仍是
// 两段曲线，`59|929910` 的同簇相交也依然存在——牺牲没有换到任何东西。
static bool naturalSingleSegmentCandidate(
    const Connectivity& conn, const IntersectionInput& input, const SDFField& sdf,
    BezierCurve& out_curve) {
    if (hasFixedGeometry(conn))
        return false;
    auto entry = input.entryPtDir(conn.entry_lane_id);
    auto exit_ = input.exitPtDir(conn.exit_lane_id);
    Vec2d p0 = entry.first;
    Vec2d p1 = exit_.first;
    Vec2d chord = p1 - p0;
    double chord_len = chord.norm();
    if (chord_len < 1e-6 || entry.second.norm() < 1e-8 || exit_.second.norm() < 1e-8)
        return false;
    Vec2d T0 = entry.second.normalized();
    Vec2d T1 = exit_.second.normalized();
    if (T0.dot(T1) < -0.5)
        return false;

    double turn_strength = std::abs(cross2d(T0, chord.normalized()));
    std::vector<double> alphas = turn_strength < 0.25
        ? std::vector<double>{0.30, 0.34, 0.38, 0.42}
        : std::vector<double>{0.34, 0.38, 0.42, 0.46};

    BezierCurve best;
    double best_score = std::numeric_limits<double>::infinity();
    for (double alpha : alphas) {
        BezierCurve candidate;
        candidate.segs.push_back(makeCubicG1(p0, T0, p1, T1, alpha));
        if (candidate.empty() || curveSelfIntersectsBusiness(candidate, 1.0))
            continue;
        double arc_chord = candidate.arcLength() / chord_len;
        if (turn_strength < 0.25) {
            if (arc_chord > 1.04)
                continue;
        } else if (arc_chord < 1.05 || arc_chord > 1.35) {
            continue;
        }
        // 物理可达性：不穿障、不穿 Boundary/RoadEdge、不出精围栏。
        // 放在形态过滤之后，避免对必然被形态淘汰的候选做 SDF 采样。
        // 同簇交叉不在此判定（siblings 传空）：本函数问的正是"若不被这条固有形态
        // 挡住，它能否取到自然形态"，把同簇交叉算进来会自我循环。
        if (assessCurveRisk(candidate, input, sdf, {}).physical())
            continue;
        double score = candidate.maxCurvature(20) + 0.02 * candidate.arcLength();
        if (score < best_score) {
            best = candidate;
            best_score = score;
        }
    }
    if (best.empty())
        return false;
    out_curve = best;
    return true;
}

// 检查固定形态是否会阻塞其它连接的自然单段形态，并收集需要取消保护的固定 ID。
// 只处理非结构性交叉，且不会直接修改结果数组。
//
// 自然候选必须物理可达（见 naturalSingleSegmentCandidate）：撤销固有形态保护的
// 唯一收益是让被挡的连接取到自然单段形态，若那条连接因物理原因根本取不到，
// 撤保护就只剩形态损失。
static std::unordered_set<ConnId> collectFixedIdsBlockingNaturalShapes(
    const IntersectionInput& input, const SDFField& sdf,
    const std::vector<ConnectivityCurve>& fixed_results,
    const ClusterOrderSolver& cs, const std::unordered_set<ConnId>& preserved_fixed_ids,
    double endpoint_tol = kClusterEndpointTol) {
    std::unordered_set<ConnId> drop_ids;
    auto idx = resultIndexById(fixed_results);
    for (const auto& conn : input.connectivities) {
        BezierCurve natural;
        if (!naturalSingleSegmentCandidate(conn, input, sdf, natural))
            continue;
        for (const auto& p : cs.pairs()) {
            if (p.exempt == CrossExemption::StructuralCross)
                continue;
            ConnId fixed_id;
            if (p.id_a == conn.id && preserved_fixed_ids.count(p.id_b)) {
                fixed_id = p.id_b;
            } else if (p.id_b == conn.id && preserved_fixed_ids.count(p.id_a)) {
                fixed_id = p.id_a;
            } else {
                continue;
            }
            auto fit = idx.find(fixed_id);
            if (fit == idx.end() || !fixed_results[fit->second].curve)
                continue;
            if (curvesHaveForbiddenSameClusterIntersection(
                    natural, *fixed_results[fit->second].curve, endpoint_tol))
                drop_ids.insert(fixed_id);
        }
    }
    return drop_ids;
}

// 根据输入连接和入口/出口方向判断其是否为几何意义上的 U-turn。
// 该判断独立于 turn_type 字段，用于选择三段式调头生成路径。
static bool isGeometricUTurnConn(const Connectivity& conn, const IntersectionInput& input) {
    return UTurnFamilyBuilder().isGeometricUTurn(conn, input);
}

// 返回 U-turn 家族排序所用的几何半径键，半径越大表示应越偏向外层深拱。
static double uturnRadiusKey(const Connectivity& conn, const IntersectionInput& input) {
    return UTurnFamilyBuilder().radiusKey(conn, input);
}

// 根据 U-turn 家族排名计算共享入口或出口平齐点的横向错开量。
// 家族不足两条时返回零；较大的半径排名获得更大的错开，以避免中弧贴合。
// out_family_step 回传家族统一分档步长（相邻半径名次的横向站位之差）。
// out_entry_stagger / out_exit_stagger 回传入口侧与出口侧各自裁定的站位：
// 家族是一条交替由共享入口、共享出口连接的链，两侧的单调次序不是同一个次序，
// 必须分别下发给 buildSegmented（见 UTurnFamilyLadder）。返回值是两侧较大者，
// 供只关心量级的旧口径（例如是否启用分档）使用。
static double uturnSharedEndpointStagger(
    const Connectivity& conn, const IntersectionInput& input,
    const ClusterOrderSolver& cs, UTurnAlignmentScope scope,
    double* out_family_step = nullptr,
    double* out_entry_stagger = nullptr,
    double* out_exit_stagger = nullptr) {
    if (out_family_step)
        *out_family_step = 0.0;
    if (out_entry_stagger)
        *out_entry_stagger = 0.0;
    if (out_exit_stagger)
        *out_exit_stagger = 0.0;
    // 分档必须与平齐站位使用相同的同入/同出传递连通家族。若只在当前
    // entry 或 exit 的一跳家族中排序，连接两个家族的 U-turn 会在中弧
    // 起点/终点拿到同一分档，重新产生非连接点交叉。
    const bool require_shared_endpoint_pair =
        scope == UTurnAlignmentScope::LaneEndpoint;
    const UTurnFamilyRank rank = UTurnFamilyBuilder().radiusRank(
        conn, input, scope, &cs, require_shared_endpoint_pair);
    if (rank.family_size < 2) {
        // 单成员 U-turn 也可能与普通直行/转向共享真实连接点。若不做
        // 厘米级扇出，其强制 2m 直段会与普通曲线的首/尾把手重叠。
        // 这里只在确有共享端点的普通兄弟时启用，不改变独立 U-turn。
        for (const auto& other : input.connectivities) {
            if (other.id == conn.id || isGeometricUTurnConn(other, input))
                continue;
            if (other.entry_lane_id == conn.entry_lane_id ||
                other.exit_lane_id == conn.exit_lane_id) {
                if (out_entry_stagger)
                    *out_entry_stagger = kUTurnSharedEndpointLeadStagger;
                if (out_exit_stagger)
                    *out_exit_stagger = kUTurnSharedEndpointLeadStagger;
                return kUTurnSharedEndpointLeadStagger;
            }
        }
        return 0.0;
    }
    // family按短径到长径排列；错开量按逆序给N，使内层曲线更紧、
    // 外层曲线更开，同时保留同入口直行段的扇出间距。
    // 纯几何输入严格使用需求规定的厘米级分档。存在 Boundary、Obstacle
    // 或 Crosswalk 时，候选还必须为物理避让保留扇出空间；此时保留旧的
    // 分米级基线，避免障碍绕行后的多条外层 U-turn 重新穿回内层。
    const bool physical_scene = !isPureGeometricInput(input);
    const double nominal_step = physical_scene && rank.family_size >= 3
        ? 0.25 : kUTurnSharedEndpointLeadStagger;
    // 分米级基线只是"名义"步长：家族里若有走廊被压到厘米级的退化成员，
    // 名义步长会让它的分档被自身走廊裁掉，而与它共享端点的外层邻居拿到完整
    // 分档，阶梯反号后外层入口直段横扫内层走廊（100000412 的 34|36、16|18、
    // 8|10）。因此最终站位改由家族统一裁定：先按名义步长排阶梯，再按各成员
    // 走廊与共享端点拓扑压掉冲突的档位，见 familyLateralLadder。
    const UTurnFamilyLadder ladder = UTurnFamilyBuilder().familyLateralLadder(
        conn, input, scope, nominal_step, &cs, require_shared_endpoint_pair);
    if (out_family_step)
        *out_family_step = ladder.step;
    if (out_entry_stagger)
        *out_entry_stagger = ladder.entry_stagger;
    if (out_exit_stagger)
        *out_exit_stagger = ladder.exit_stagger;
    return ladder.stagger;
}

// 纯几何形态恢复：当直行/普通转向被固有形态或拓扑修复拉成异常多段、
// 过短弧或S形时，搜索更自然的单段/少段候选，并确保同簇交叉不变差。
// 在纯几何拓扑修复后，为异常多段或 S 形普通转向搜索自然形态恢复候选。
// 录取候选必须保持当前同簇交叉和固定形态交叉不增加，并满足转向弧形约束。
static bool tryPureGeometryShapeRestore(
    const Connectivity& conn, const IntersectionInput& input,
    const std::vector<ConnectivityCurve>& results,
    const std::unordered_map<ConnId, size_t>& result_idx,
    const std::unordered_map<ConnId, std::vector<ConnId>>& neighbors,
    BezierCurve& out_curve) {
    auto it = result_idx.find(conn.id);
    if (it == result_idx.end() || !results[it->second].curve)
        return false;
    const BezierCurve& current = *results[it->second].curve;

    auto entry = input.entryPtDir(conn.entry_lane_id);
    auto exit_ = input.exitPtDir(conn.exit_lane_id);
    Vec2d p0 = entry.first;
    Vec2d t0 = entry.second;
    Vec2d p1 = exit_.first;
    Vec2d t1 = exit_.second;
    Vec2d chord = p1 - p0;
    double chord_len = chord.norm();
    if (chord_len < 1e-6)
        return false;

    Vec2d T0 = t0.norm() > 1e-8 ? t0.normalized() : chord.normalized();
    Vec2d T1 = t1.norm() > 1e-8 ? t1.normalized() : chord.normalized();
    double turn_strength = std::abs(cross2d(T0, chord.normalized()));
    bool turn_shape_risk = isNonUTurnTurnShapeRisk(current, chord_len, turn_strength);
    if (turn_strength <= 0.35 || !turn_shape_risk)
        return false;

    int current_cross = constrainedCrossCountForId(
        conn.id, current, results, result_idx, neighbors, kClusterEndpointTol);
    int current_fixed_cross = constrainedFixedCrossCountForId(
        conn.id, current, results, result_idx, neighbors, kClusterEndpointTol);

    std::vector<BezierCurve> candidates;
    candidates.reserve(96);
    auto addCandidate = [&](BezierCurve&& candidate) {
        if (!candidate.empty())
            candidates.push_back(std::move(candidate));
    };

    // 此阶段位于拓扑修复之后，必须保持交互性能，因此搜索范围有意收紧。
    // 下方录取门禁仍保证恢复后的转向形态不会恶化同簇或 fixed shape 交叉。
    Vec2d chord_dir = chord.normalized();
    Vec2d chord_perp{-chord_dir.y(), chord_dir.x()};
    double natural_side = cross2d(T0, T1) >= 0.0 ? 1.0 : -1.0;
    for (double alpha : {0.42, 0.38, 0.46, 0.34, 0.50, 0.30}) {
        BezierCurve single;
        single.segs.push_back(makeCubicG1(p0, T0, p1, T1, alpha));
        addCandidate(std::move(single));
    }

    Vec2d base_mid_tan = T0 + T1;
    if (base_mid_tan.norm() < 1e-8)
        base_mid_tan = chord_dir;
    base_mid_tan.normalize();
    for (double frac : {0.42, 0.50, 0.58}) {
        for (double side : {natural_side, -natural_side}) {
            for (double offset : {1.5, 2.8, 4.2, 6.0, 8.0}) {
                Vec2d mid = p0 + chord * frac + chord_perp * (side * offset);
                for (double deg : {0.0, side * 30.0, side * 55.0}) {
                    Vec2d mid_tan = rotateVec(base_mid_tan, deg * M_PI / 180.0);
                    if (mid_tan.norm() < 1e-8)
                        continue;
                    mid_tan.normalize();
                    for (double alpha0 : {0.24, 0.30, 0.42, 0.50}) {
                        for (double alpha1 : {0.20, 0.24, 0.36, 0.46}) {
                            BezierCurve shaped;
                            shaped.segs.push_back(makeCubicG1(p0, T0, mid, mid_tan, alpha0));
                            shaped.segs.push_back(makeCubicG1(mid, mid_tan, p1, T1, alpha1));
                            addCandidate(std::move(shaped));
                        }
                    }
                }
            }
        }
    }

    std::vector<CurveCandidate> audited_candidates;
    audited_candidates.reserve(candidates.size());
    for (const auto& candidate : candidates) {
        if (candidate.empty() || curveSelfIntersectsBusiness(candidate, 1.0))
            continue;
        if (!isNonUTurnTurnShapeAcceptable(candidate, chord_len, turn_strength))
            continue;
        double arc_chord = candidate.arcLength() / chord_len;
        double pure_min_arc = ordinaryTurnArcChordRestoreFloor(
            curveEndpointTurnAngle(candidate), chord_len);
        if (arc_chord < pure_min_arc)
            continue;
        int cross_count = constrainedCrossCountForId(
            conn.id, candidate, results, result_idx, neighbors, kClusterEndpointTol);
        int fixed_cross_count = constrainedFixedCrossCountForId(
            conn.id, candidate, results, result_idx, neighbors, kClusterEndpointTol);
        if (fixed_cross_count > current_fixed_cross || cross_count > current_cross)
            continue;
        double score = 0.5 * std::abs(arc_chord - 1.12) +
            candidate.maxCurvature(20) + 0.02 * candidate.arcLength();
        CurveCandidate accepted;
        accepted.curve = candidate;
        accepted.origin = CandidateOrigin::Repaired;
        accepted.hard_violation_count = fixed_cross_count;
        accepted.new_cluster_crosses = cross_count;
        accepted.quality_score = score;
        accepted.acceptReport(ConstraintReport());
        audited_candidates.push_back(std::move(accepted));
    }

    const std::size_t best = CandidateSelector().selectBest(audited_candidates);
    if (best == CandidateSelector::npos)
        return false;
    out_curve = audited_candidates[best].curve;
    return true;
}

// 判断端点弦线是否进入障碍物硬多边形，包括中点在内部和与外环相交两种情况。
static bool chordIntersectsObstacle(
    const Vec2d& p0, const Vec2d& p1, const Obstacle& obs) {
    const Polygon2d& poly = obstacleHardGeometry(obs);
    if (poly.outer.size() < 3)
        return false;
    if (pointInsideOrOnPolygon(0.5 * (p0 + p1), poly))
        return true;
    return polygonSegmentIntersects(p0, p1, poly);
}

// 围绕阻塞端点弦线的障碍物生成有限路点绕行候选。
// 候选按障碍物两侧、净距和同簇交叉排序，并必须通过物理风险及 U-turn 拓扑门禁。
static bool tryObstacleBypassCandidate(
    const Vec2d& p0, const Vec2d& t0, const Vec2d& p1, const Vec2d& t1,
    const IntersectionInput& input, const SDFField& sdf,
    const std::vector<SampledSiblingCurve>& sampled_siblings, BezierCurve& curve) {
    if (input.obstacles.empty())
        return false;

    Vec2d along = p1 - p0;
    double len = along.norm();
    if (len < 1e-6)
        return false;
    along /= len;
    Vec2d perp{-along[1], along[0]};
    Vec2d T0 = t0.norm() > 1e-8 ? t0.normalized() : Vec2d(1, 0);
    Vec2d T1 = t1.norm() > 1e-8 ? t1.normalized() : Vec2d(1, 0);
    bool require_topology_clear = (T0.dot(T1) < -0.5);

    BezierCurve best;
    bool have_best = false;
    int best_cross = std::numeric_limits<int>::max();
    double best_score = std::numeric_limits<double>::max();
    Vec2d mid = 0.5 * (p0 + p1);
    Vec2d center = input.area.geometry.outer.empty() ? Vec2d(0, 0) : Vec2d(0, 0);
    if (!input.area.geometry.outer.empty()) {
        center.setZero();
        int cnt = 0;
        for (const auto& p : input.area.geometry.outer) {
            center += xyOf(p);
            ++cnt;
        }
        if (cnt > 0)
            center /= cnt;
    }
    double center_lat = (center - mid).dot(perp);

    for (const auto& obs : input.obstacles) {
        const Polygon2d& poly = obstacleHardGeometry(obs);
        if (poly.outer.size() < 3)
            continue;
        if (!chordIntersectsObstacle(p0, p1, obs) &&
            !curveIntersectsObstacles(curve, std::vector<Obstacle>{obs}))
            continue;

        double lon_min = 1e18, lon_max = -1e18;
        double lat_min = 1e18, lat_max = -1e18;
        for (const auto& p : poly.outer) {
            double lon = (p - p0).dot(along);
            double lat = (p - p0).dot(perp);
            lon_min = std::min(lon_min, lon);
            lon_max = std::max(lon_max, lon);
            lat_min = std::min(lat_min, lat);
            lat_max = std::max(lat_max, lat);
        }
        double lon_mid = std::max(len * 0.15, std::min(len * 0.85, 0.5 * (lon_min + lon_max)));
        double obs_lat_mid = 0.5 * (lat_min + lat_max);
        std::vector<int> sides;
        if (std::abs(center_lat) > 0.25 && std::abs(obs_lat_mid) > 0.25 &&
            ((center_lat > 0) != (obs_lat_mid > 0))) {
            sides = {(center_lat > 0) ? +1 : -1, (center_lat > 0) ? -1 : +1};
        } else {
            sides = {+1, -1};
        }

        for (int side : sides) {
            double edge_lat = side > 0 ? lat_max : lat_min;
            for (double clearance : {0.65, 1.0, 1.4, 1.9, 2.5, 3.2}) {
                double apex_lat = edge_lat + side * clearance;
                Vec2d apex = p0 + lon_mid * along + apex_lat * perp;
                double spread = std::max(2.0, std::min(len * 0.28, (lon_max - lon_min) * 0.8 + 1.5));
                double lon_a = std::max(len * 0.08, lon_mid - spread);
                double lon_b = std::min(len * 0.92, lon_mid + spread);
                Vec2d q0 = p0 + lon_a * along + apex_lat * 0.55 * perp;
                Vec2d q1 = p0 + lon_b * along + apex_lat * 0.55 * perp;
                std::vector<Vec2d> pts;
                if ((q0 - p0).norm() > 1.0 && (q1 - p1).norm() > 1.0 && lon_b > lon_a + 0.5)
                    pts = {p0, q0, apex, q1, p1};
                else
                    pts = {p0, apex, p1};

                BezierCurve candidate = AvoidanceCandidateGenerator().buildWaypointCurve(
                    pts, t0, t1);
                if (candidate.empty() || curveSelfIntersectsBusiness(candidate, 1.0))
                    continue;
                CurveRisk risk = assessCurveRisk(candidate, input, sdf, sampled_siblings);
                if (risk.physical() || (require_topology_clear && risk.sibling_crosses > 0))
                    continue;
                double shape_score = clearance + 0.02 * candidate.arcLength() +
                    (side * center_lat >= 0.0 ? 0.0 : 0.4);
                if (!have_best ||
                    risk.sibling_crosses < best_cross ||
                    (risk.sibling_crosses == best_cross && shape_score < best_score)) {
                    best = candidate;
                    have_best = true;
                    best_cross = risk.sibling_crosses;
                    best_score = shape_score;
                }
            }
        }
    }

    if (!have_best)
        return false;
    curve = best;
    return true;
}

// 使用有界 U-turn 搜索器寻找相对 reference 的安全调头候选。
// auditor 统一检查物理风险和兄弟交叉，min_lead 参数约束首尾直行段最小长度。
static bool tryBoundedUTurnCandidate(
    const Vec2d& p0, const Vec2d& t0, const Vec2d& p1, const Vec2d& t1,
    const IntersectionInput& input, const SDFField& sdf,
    const std::vector<SampledSiblingCurve>& sampled_siblings,
    const BezierCurve& reference, BezierCurve& out_curve,
    double min_lead0 = 0.0, double min_lead1 = 0.0) {
    const BoundedUTurnCandidateAuditor auditor =
        [&](const BezierCurve& candidate) {
            const CurveRisk risk = assessCurveRisk(
                candidate, input, sdf, sampled_siblings, true, true);
            BoundedUTurnCandidateAudit audit;
            audit.physical_violation = risk.physical();
            audit.sibling_crosses = risk.sibling_crosses;
            return audit;
        };
    return BoundedUTurnCandidateSearch().search(
        p0, t0, p1, t1, reference, auditor, out_curve,
        min_lead0, min_lead1);
}

// 为 U-turn 多约束搜索器组装采样、同簇交叉、障碍物、边界和共享端点侧向惩罚回调。
// 返回的后端捕获当前输入和 SDF，只读服务于候选评分。
static UTurnSearchBackend makeUTurnSearchBackend(
    const IntersectionInput& input, const SDFField& sdf,
    const std::vector<SampledSiblingCurve>& sampled_siblings) {
    UTurnSearchBackend backend;
    backend.sample = [](const BezierCurve& curve) {
        return sampleCurveForIntersections(curve);
    };
    backend.sibling_cross_count = [&](const BezierCurve& curve) {
        return sampledSiblingCrossCount(
            curve, sampled_siblings, true, kClusterEndpointTol,
            input.obstacles.empty());
    };
    backend.obstacle_penalty = [&](const BezierCurve& curve) {
        if (input.obstacles.empty())
            return 0.0;
        double penalty = sdf.valid()
            ? std::max(0.0, -minSDFAlongCurveAdaptive(curve, sdf)) : 0.0;
        if (curveIntersectsObstacles(curve, input.obstacles))
            penalty = std::max(penalty, 0.10);
        return penalty;
    };
    backend.boundary_penalty = [&](const BezierCurve& curve) {
        return boundaryAvoidancePenalty(curve, input, true);
    };
    backend.endpoint_side_violation = [&](const SampledCurve& curve) {
        return sharedEndpointSideViolation(
            curve, sampled_siblings, kClusterEndpointTol);
    };
    return backend;
}

// 创建三段式 U-turn 候选审计器，统一执行物理风险和同簇交叉检查。
// include_fence 由调用阶段决定是否把围栏纳入硬门禁。
// allow_merge_funnel 只在两阶段搜索的兜底档为 true，见 searchSegmentedUTurnTwoPass。
static SegmentedUTurnAuditor makeSegmentedUTurnAuditor(
    const IntersectionInput& input, const SDFField& sdf,
    const std::vector<SampledSiblingCurve>& sampled_siblings,
    bool allow_merge_funnel = false) {
    return [&input, &sdf, &sampled_siblings, allow_merge_funnel](
               const BezierCurve& curve, bool include_fence) {
        const CurveRisk risk = assessCurveRisk(
            curve, input, sdf, sampled_siblings, include_fence, true,
            allow_merge_funnel);
        SegmentedUTurnAudit audit;
        audit.physical_violation = risk.physical();
        audit.sibling_crosses = risk.sibling_crosses;
        return audit;
    };
}

// 走廊反转判据：三段式掉头的首尾直段沿 U 轴法线各自漂移
// `lead × |T0·lateral|`，两段漂移之和一旦达到端点走廊宽度
// `|(p1-p0)·lateral|`，平齐站位处的有向走廊就会翻符号——首尾直段必须换侧，
// 三段式在几何上无法避免与"汇入同一条出口车道的兄弟"相交。
//
// 这是区分"真死局"与"只是没搜到更好候选"的关键量：
//   110000703-u 的 41：走廊 0.3099m，人行横道要求首尾直段 9.318+9.237m，
//                      |T0·lateral|=0.01726 ⇒ 漂移 0.3203 ≥ 0.3099，反转；
//   100000643 的 125：走廊 0.1643m，无横道只有 2m 底限 ⇒ 漂移 0.1136 < 0.1643，
//                      不反转，必须继续按严格规则修复（放开会让 125-113 相交）。
static bool segmentedUTurnCorridorInverts(
    const Vec2d& p0, const Vec2d& t0, const Vec2d& p1, const Vec2d& t1,
    double min_lead0, double min_lead1) {
    if (t0.norm() < 1e-8 || t1.norm() < 1e-8)
        return false;
    const Vec2d u0 = t0.normalized();
    Vec2d axis = u0 - t1.normalized();
    if (axis.norm() < 1e-8)
        return false;
    axis = axis.normalized();
    const Vec2d lateral(-axis.y(), axis.x());
    const double corridor = std::abs((p1 - p0).dot(lateral));
    if (corridor >= kDegenerateUTurnCorridor)
        return false;
    const double drift =
        (std::max(0.0, min_lead0) + std::max(0.0, min_lead1)) *
        std::abs(u0.dot(lateral));
    return drift >= corridor;
}

// 三段式掉头候选搜索的两阶段包装。
//
// 第一阶段用严格同簇判定，与本次改动前逐字节一致：任何当前能搜到合规候选的
// 掉头，形态与耗时都不受影响。
//
// 第二阶段（结构性汇入豁免兜底）需要同时满足三个前提，缺一不可：
//   ① 严格档在 832 个几何候选 + 偏置候选里一无所获；
//   ② **当前形态确实违反了人行横道分段要求**——即入参曲线不是三段式、或中间
//      弧没能整体避开净空人行横道集合；
//   ③ **端点走廊在平齐站位处反转**（见 segmentedUTurnCorridorInverts），即
//      "人行横道净空要求的首尾直段"与"同簇非端点不相交"在该车道几何下确实
//      无法同时满足。
//
// ②③ 是范围控制的关键。只按 ① 兜底会外溢：100000643 里 125、115 的严格档同样
// 返回 false（只是没能改进现状，并非死局），一旦放行豁免，125-113、115-103 就由
// 通过转为相交——这正是"每修一次掉头就触发一处旧违约"的机制。附带收益是耗时：
// 只有真正的死局才付两遍搜索的代价。
static bool searchSegmentedUTurnTwoPass(
    const Vec2d& entry, const Vec2d& entry_tangent, const Vec2d& exit,
    const Vec2d& exit_tangent, const IntersectionInput& input,
    const SDFField& sdf,
    const std::vector<SampledSiblingCurve>& sampled_siblings,
    bool include_fence, BezierCurve& curve, double min_lead0, double min_lead1,
    const std::vector<Crosswalk>* crosswalks_for_clearance,
    double aligned_point_stagger, double base_lead0_extra_after_align,
    double base_lead1_extra_after_align, double aligned_family_station,
    double aligned_entry_stagger, double aligned_exit_stagger,
    double family_stagger_step) {
    const bool crosswalk_shape_violated =
        crosswalks_for_clearance != nullptr &&
        !crosswalks_for_clearance->empty() &&
        (curve.empty() || curve.numSegments() != 3 ||
         !segmentedUTurnMiddleArcClearsCrosswalks(
             curve, *crosswalks_for_clearance));
    const bool deadlocked = crosswalk_shape_violated &&
        segmentedUTurnCorridorInverts(entry, entry_tangent, exit, exit_tangent,
                                      min_lead0, min_lead1);
    // 严格需求下不再使用第二遍 merge-funnel 宽松搜索；保留一次有限候选搜索，
    // 避免在无法满足形态/物理约束时重复消耗整张候选网格。
    const int passes = 1;
    for (int pass = 0; pass < passes; ++pass) {
        BezierCurve working = curve;
        if (SegmentedUTurnCandidateSearch().search(
                entry, entry_tangent, exit, exit_tangent, input,
                sampled_siblings,
                makeSegmentedUTurnAuditor(input, sdf, sampled_siblings,
                                          pass == 1),
                include_fence, working, min_lead0, min_lead1,
                crosswalks_for_clearance, aligned_point_stagger,
                base_lead0_extra_after_align, base_lead1_extra_after_align,
                aligned_family_station, aligned_entry_stagger,
                aligned_exit_stagger, family_stagger_step)) {
            // 诊断开关：ISG_DEBUG_MERGE_FALLBACK=1 打印兜底档真正生效的位置，
            // 用于确认结构性汇入豁免没有在别的路口悄悄改变形态。
            static const bool debug_fallback = [] {
                const char* v = std::getenv("ISG_DEBUG_MERGE_FALLBACK");
                return v != nullptr && v[0] == '1';
            }();
            if (debug_fallback && pass == 1) {
                fprintf(stderr,
                        "[merge-fallback] p0=(%.3f,%.3f) p1=(%.3f,%.3f)"
                        " segs=%d\n",
                        entry.x(), entry.y(), exit.x(), exit.y(),
                        (int)working.numSegments());
            }
            curve = working;
            return true;
        }
    }
    return false;
}

class SingleCurveGenerationPipeline {
public:
    typedef std::function<void(ConnectivityCurve&, const IntersectionInput&, const SDFField&)> Validator;

    SingleCurveGenerationPipeline(
        LBFGSSolver& solver, const ClusterOrderSolver& cluster_solver,
        const ConnectivityDirectionConfig& direction_config, const Validator& validator)
        : solver_(solver), cluster_solver_(cluster_solver),
          direction_cfg_(direction_config), validator_(validator) {}

    // 按连接类型选择初始曲线、候选修复、优化器和最终审计，生成一条完整输出曲线。
    // siblings 提供已生成兄弟约束；out_physical_risk 报告最终是否仍有物理风险。
    ConnectivityCurve run(
        const CurveGenerationContext& context, const std::vector<SiblingCurve>& siblings, bool* out_physical_risk);

private:
    void validate(ConnectivityCurve& curve, const IntersectionInput& input, const SDFField& sdf) const {
        validator_(curve, input, sdf);
    }

    LBFGSSolver& solver_;
    const ClusterOrderSolver& cluster_solver_;
    const ConnectivityDirectionConfig& direction_cfg_;
    Validator validator_;
};

// 执行单连接生成流水线：解析端点与 U-turn 家族约束，构造初解，按物理/拓扑/形态
// 风险选择修复或优化路径，最后写入几何并完成统一审计。该函数不修改其它连接结果。
ConnectivityCurve SingleCurveGenerationPipeline::run(
    const CurveGenerationContext& context, const std::vector<SiblingCurve>& siblings, bool* out_physical_risk) {
    if (out_physical_risk)
        *out_physical_risk = false;
    if (!context.connectivity || !context.scene || !context.scene->sdf ||
        !context.scene->coarse_sdf)
        return ConnectivityCurve();
    const Connectivity& conn = *context.connectivity;
    const IntersectionInput& input = context.scene->view.input();
    const SDFField& sdf = *context.scene->sdf;
    const SDFField& sdf_coarse = *context.scene->coarse_sdf;
    ConnectivityCurve cc;
    cc.id = conn.id;
    cc.entry_lane_id = conn.entry_lane_id;
    cc.exit_lane_id = conn.exit_lane_id;
    cc.turn_type = context.turn;
    cc.fixed_shape = conn.fixed_shape;

    const std::pair<Vec2d, Vec2d>& _entry = context.entry;
    Vec2d p0 = _entry.first;
    Vec2d t0 = _entry.second;
    const std::pair<Vec2d, Vec2d>& _exit = context.exit;
    Vec2d p1 = _exit.first;
    Vec2d t1 = _exit.second;
    const bool is_uturn_geom = context.uturn_family.geometric_uturn;

    double lw = 3.5;
    if (auto* l = input.findLane(conn.entry_lane_id))
        lw = l->width;
    bool enforce_fence = !(input.obstacles.empty() && input.boundaries.empty());
    // Boundary 密集场景统一限制重复候选搜索。mode=2 的 RoadEdge 相交和
    // 1m 净距仍由每个候选的 assessCurveRisk/Boundary 门禁完整校验。
    bool complex_boundary_fast_path =
        input.obstacles.empty() && input.boundaries.size() >= 20;

    // ── 构造初始曲线：同时考虑同簇拓扑与障碍绕行方向 ─────────────────────
    // 对障碍绕行只纳入非豁免兄弟曲线，避免初始绕行方向与已生成同簇曲线冲突。
    BezierCurve initial;
    // U型调头专用：首尾直行段最小长度(米)。
    // 有近旁人行横道的一侧取跨越远边距离；缺失侧仅随U型轴线平齐补长。
    // 只有两侧都没有人行横道时才强制使用2m保底。
    double uturn_min_lead0 = 0.0;
    double uturn_min_lead1 = 0.0;
    double uturn_shared_endpoint_stagger = 0.0;
    double uturn_shared_entry_stagger = 0.0;
    double uturn_shared_exit_stagger = 0.0;
    double uturn_lead0_extra_after_align = 0.0;
    double uturn_lead1_extra_after_align = 0.0;
    double uturn_alignment_station =
        std::numeric_limits<double>::quiet_NaN();
    double uturn_family_stagger_step = 0.0;
    const UTurnFamilyInfo& uturn_family_info = context.uturn_family;
    std::vector<Crosswalk> uturn_clearance_crosswalks;
    if (is_uturn_geom) {
        uturn_min_lead0 = uturn_family_info.lead0;
        uturn_min_lead1 = uturn_family_info.lead1;
        // 只有多成员平齐家族才拥有不可越过的共同站位。单条 U-turn
        // 仍要求 q0/q1 自身平齐和最小 lead，但必须允许候选对称加长，
        // 才能在近邻直行簇外找到不相交的中弧。
        if (uturn_family_info.rank.family_size > 1)
            uturn_alignment_station = uturn_family_info.aligned_station;
        uturn_clearance_crosswalks =
            uturn_family_info.clearance_crosswalks;
        uturn_shared_endpoint_stagger = uturnSharedEndpointStagger(
            conn, input, cluster_solver_,
            direction_cfg_.uturn_alignment_scope,
            &uturn_family_stagger_step,
            &uturn_shared_entry_stagger,
            &uturn_shared_exit_stagger);
        // 同配置家族的多条U-turn按两侧最大顺序号计算统一错开距离，
        // 首/尾两个已轴向平齐的点沿轴向法线同步相向微移；
        // 移动后仍保持 q0/q1 的轴向平齐关系。
        // 入口侧与出口侧的分档量由 familyLateralLadder 分别裁定：两侧共享的
        // 是不同的端点链，压小一侧不应连带压小另一侧。
    }

    // 兄弟曲线采样缓存供U型调头多约束求解器和后续交叉检测复用。
    std::vector<SampledSiblingCurve> sampled_siblings;
    auto ensure_sampled_siblings = [&]() -> const std::vector<SampledSiblingCurve>& {
        if (sampled_siblings.empty() && !siblings.empty())
            sampled_siblings = sampleSiblingsForIntersections(
                siblings, 32, !complex_boundary_fast_path);
        return sampled_siblings;
    };

    // ── 修复 74|76: 拱弧深度偏好 (depth_preference) ──────────────────────
    // +1 = 当前U型调头应为"外层深拱"
    // -1 = 当前U型调头应为"内层浅拱"
    //  0 = 无偏好 (普通搜索)
    // 声明在 is_uturn_geom 分支外, 以便后续在优化器跳过逻辑中使用。
    //
    // 计算方法: 遍历所有共进入车道U型调头配对，累计更深/更浅计数。
    // net = deeper_count - shallower_count。
    // 若净值大于0则偏好更深，小于0则偏好更浅，等于0则使用默认搜索。
    // 这能处理"中间"U型调头: 它同时比一侧深、比另一侧浅时净偏好为0。
    int depth_preference = 0;      // lead0 (共进入) 深度偏好
    int lead1_depth_preference = 0; // lead1 (共退出) 深度偏好, 修复 74|73
    // lateral_preference: 横向偏置偏好
    // +1.0 = 当前U型调头应横向偏移到 +lat_dir 方向
    // -1.0 = 当前U型调头应横向偏移到 -lat_dir 方向
    // 0.0 = 无偏好
    // 这使内外层U型调头的弧线横向分离，降低中段相交风险。
    double lateral_preference = 0.0;

    if (is_uturn_geom) {
        // ── 同端点U型调头家族深度偏好 ───────────────────────────
        // 大径U型调头应作为外层深拱包裹小径U型调头。此前这里从
        // expected_side/ref_perp/axis 符号反推"深/浅", 在不同arm方向下会翻转,
        // 例如 100000643 中 75(turn_gap=11.3m) 被判成浅拱,
        // 76(turn_gap=8.0m) 被判成深拱, 导致二者中段相交。
        //
        // 现在直接使用几何半径键(turn_gap)排序:
        // - 共进入实际车道: 大径给 lead0 更深, 小径更浅;
        // - 共退出实际车道: 大径给 lead1 更深, 小径更浅。
        // 中间半径会同时有更大/更小兄弟, 净偏好自然为0。
        {
            int deeper_count = 0, shallower_count = 0;
            int lead1_deeper = 0, lead1_shallower = 0;
            double current_radius = uturnRadiusKey(conn, input);

            for (const auto& pair : cluster_solver_.pairs()) {
                if (pair.exempt != CrossExemption::None) continue;
                ConnId other_id;
                if (pair.id_a == conn.id) other_id = pair.id_b;
                else if (pair.id_b == conn.id) other_id = pair.id_a;
                else continue;

                // 查找兄弟连通关系。
                const Connectivity* other_conn = nullptr;
                for (const auto& c : input.connectivities) {
                    if (c.id == other_id) { other_conn = &c; break; }
                }
                if (!other_conn) continue;

                // 兄弟必须也是U型调头(几何判断)。
                if (!isGeometricUTurnConn(*other_conn, input)) continue;
                double other_radius = uturnRadiusKey(*other_conn, input);
                if (std::abs(current_radius - other_radius) < 0.5) continue;
                int dp = (current_radius > other_radius) ? +1 : -1;

                // 区分共进入 vs 共退出
                bool shared_entry = (other_conn->entry_lane_id == conn.entry_lane_id &&
                                     other_conn->exit_lane_id != conn.exit_lane_id);
                bool shared_exit = (other_conn->exit_lane_id == conn.exit_lane_id &&
                                    other_conn->entry_lane_id != conn.entry_lane_id);

                if (shared_entry) {
                    if (dp > 0) ++deeper_count;
                    else if (dp < 0) ++shallower_count;
                } else if (shared_exit) {
                    // 共退出: lead1 深度偏好
                    if (dp > 0) ++lead1_deeper;
                    else if (dp < 0) ++lead1_shallower;
                }

                // ── 横向偏好计算 ───────────────────────────────────────
                // expected_side 表示兄弟在当前曲线的哪一侧；求解器的
                // lateral_bias 沿调头局部 lat_dir，因此要把簇横向参考轴
                // ref_perp 投影到 lat_dir 上，得到可用于候选搜索的有符号偏置。
                {
                    int es = cluster_solver_.expectedSideOf(other_id, conn.id);
                    Vec2d rp = cluster_solver_.refPerpOf(conn.id, other_id);
                    if (es == 0 || rp.norm() < 1e-9)
                        continue;
                    Vec2d T0_n = t0.norm() > 1e-8 ? t0.normalized() : Vec2d(1, 0);
                    Vec2d T1_n = t1.norm() > 1e-8 ? t1.normalized() : Vec2d(1, 0);
                    Vec2d ot_axis = (T0_n - T1_n).normalized();
                    if (ot_axis.dot(T0_n) < 0) ot_axis = -ot_axis;
                    Vec2d lat_dir{-ot_axis[1], ot_axis[0]};
                    double lat_ref_proj = lat_dir.dot(rp.normalized());
                    if (std::abs(lat_ref_proj) > 0.1) {
                        // 转换规则: 期望侧 × (lat_dir 与 ref_perp 的同/反向)。
                        double lp = es * (lat_ref_proj > 0 ? 1.0 : -1.0);
                        // 累加: 取最强的偏好
                        if (std::abs(lp) > std::abs(lateral_preference))
                            lateral_preference = lp;
                    }
                }
            }
            int net = deeper_count - shallower_count;
            depth_preference = (net > 0) ? +1 : (net < 0 ? -1 : 0);
            int lead1_net = lead1_deeper - lead1_shallower;
            lead1_depth_preference = (lead1_net > 0) ? +1 : (lead1_net < 0 ? -1 : 0);
            if (isgDebugUTurn()) {
                fprintf(stderr,
                        "[UTURN-PREF] %s depth=%d lead1_depth=%d stagger=%.2f\n",
                        conn.id.c_str(), depth_preference,
                        lead1_depth_preference,
                        uturn_shared_endpoint_stagger);
            }
        }

        // ── U型调头多约束初解 ─────────────────────────────────────────
        // 直接把采样兄弟曲线传入多约束搜索器，在尺度、横向偏置、
        // 首尾直行延长和把手偏置网格中搜索联合代价最小候选。这样初始曲线
        // 一开始就能同时满足G1、同簇不交叉、障碍/围栏和调头包络约束。
        const auto& sampled = ensure_sampled_siblings();

        // 包络检查使用默认对齐调头弧作为参考。
        BezierCurve ref_arch = UTurnCurveInitializer().buildAligned(p0, t0, p1, t1,
                                                  Vec2d{-t0.normalized().y(), t0.normalized().x()},
                                                  0.0, 1.0, uturn_min_lead0, uturn_min_lead1);
        double uturn_chord_len = (p1 - p0).norm();
        auto uturn_arch_shape_ok = [&](const BezierCurve& curve) {
            if (uturn_chord_len < 10.0)
                return true;
            double arc_chord = curve.arcLength() / uturn_chord_len;
            return arc_chord > 1.40 && curve.maxCurvature(20) < 0.50;
        };

        UTurnSolverResult uturn_sol = UTurnMultiConstraintSearch().search(
            p0, t0, p1, t1, input, sampled, ref_arch,
            makeUTurnSearchBackend(input, sdf, sampled),
            uturn_min_lead0, uturn_min_lead1,
            depth_preference, lead1_depth_preference, lateral_preference,
            &uturn_clearance_crosswalks);

        if (std::isfinite(uturn_sol.cost)) {
            initial = uturn_sol.curve;
        } else {
            // 极端几何下无候选通过硬过滤时，先回退默认弧，后续优化仍可修复。
            initial = ref_arch;
        }
        bool pure_geometric_input = isPureGeometricInput(input);
        bool consider_segmented_arch =
            !pure_geometric_input || !std::isfinite(uturn_sol.cost) ||
            initial.numSegments() < 3;
        if (consider_segmented_arch) {
            BezierCurve segmented_arch =
                UTurnCurveInitializer().buildSegmented(
                    p0, t0, p1, t1, uturn_min_lead0, uturn_min_lead1,
                    2.0 / 3.0, uturn_shared_endpoint_stagger,
                    uturn_lead0_extra_after_align,
                    uturn_lead1_extra_after_align,
                    uturn_shared_entry_stagger,
                    uturn_shared_exit_stagger);
            CurveRisk segmented_risk = assessCurveRisk(
                segmented_arch, input, sdf, sampled, enforce_fence, true);
            if (uturn_arch_shape_ok(segmented_arch) &&
                !segmented_risk.physical() &&
                !curveSelfIntersectsBusiness(segmented_arch, 1.0)) {
                int initial_cross = sampledSiblingCrossCount(
                    initial, sampled, true, kClusterEndpointTol);
                double initial_shape = initial.maxCurvature(20) + 0.02 * initial.arcLength();
                double segmented_shape = segmented_arch.maxCurvature(20) + 0.02 * segmented_arch.arcLength();
                bool solver_failed = !std::isfinite(uturn_sol.cost);
                if ((solver_failed && segmented_risk.sibling_crosses <= initial_cross) ||
                    segmented_risk.sibling_crosses < initial_cross ||
                    (segmented_risk.sibling_crosses == initial_cross &&
                     segmented_arch.numSegments() > initial.numSegments() &&
                     segmented_arch.maxCurvature(20) < 3.0) ||
                    (segmented_risk.sibling_crosses == initial_cross &&
                     segmented_shape + 1e-6 < initial_shape)) {
                    initial = segmented_arch;
                }
            }
        }
    } else {
        std::vector<std::vector<Vec2d>> sib_polys_for_init;
        if (!siblings.empty() &&
            !input.obstacles.empty() &&
            naturalCubicHitsObstacle(p0, t0, p1, t1, sdf, input.obstacles)) {
            sib_polys_for_init.reserve(siblings.size());
            for (auto& sib : siblings) {
                if (!sib.exempt_a1)
                    sib_polys_for_init.push_back(sib.curve.sampleByArcLength(12));
            }
        }
        if (input.obstacles.empty()) {
            CurveInitializationOptions options;
            options.allow_fixed_shape = false;
            options.ordinary_alphas.assign(1, 0.4);
            const std::vector<CurveCandidate> candidates =
                CurveInitializerRegistry().build(context, options);
            if (!candidates.empty())
                initial = candidates.front().curve;
        } else {
            initial = buildInitialCurve(
                p0, t0, p1, t1, sdf, input.area.geometry,
                sib_polys_for_init);
        }
    }

    if (!is_uturn_geom) {
        int current_cross = sampledSiblingCrossCount(
            initial, ensure_sampled_siblings(), true, kClusterEndpointTol);
        double turn_strength = 0.0;
        bool has_shared_endpoint_sibling = false;
        for (const auto& sib : siblings) {
            if (!sib.exempt_a1 && sib.shared_endpoint) {
                has_shared_endpoint_sibling = true;
                break;
            }
        }
        Vec2d chord = p1 - p0;
        if (chord.norm() > 1e-8 && t0.norm() > 1e-8)
            turn_strength = std::abs(cross2d(t0.normalized(), chord.normalized()));
        if (turn_strength > 0.35 && !has_shared_endpoint_sibling) {
            for (double alpha : {0.50, 0.46, 0.42}) {
                BezierCurve candidate;
                candidate.segs.push_back(makeCubicG1(
                    p0, t0.norm()>1e-8 ? t0.normalized() : Vec2d(1,0),
                    p1, t1.norm()>1e-8 ? t1.normalized() : Vec2d(1,0),
                    alpha));
                if (curveSelfIntersectsBusiness(candidate, 1.0))
                    continue;
                Vec2d chord = p1 - p0;
                double chord_len = chord.norm();
                if (chord_len > 1e-8 &&
                    !isNonUTurnTurnShapeAcceptable(candidate, chord_len, turn_strength))
                    continue;
                if (assessCurveRisk(candidate, input, sdf, ensure_sampled_siblings(), enforce_fence).physical())
                    continue;
                int cross_count = sampledSiblingCrossCount(
                    candidate, ensure_sampled_siblings(), true, kClusterEndpointTol);
                if (cross_count < current_cross ||
                    (current_cross == 0 && cross_count == 0)) {
                    initial = candidate;
                    current_cross = cross_count;
                    break;
                }
            }
        }
    }

    // ── 自适应把手长度：初解若与兄弟曲线交叉，优先尝试更短G1把手 ─────────
    // buildInitialCurve 在无障碍时可能直接返回 alpha=0.4 的自然弧，未必考虑
    // 同簇拓扑；较短把手常能从源头避开同出口/同入口兄弟曲线的局部交叉。
    if (!is_uturn_geom) {
        bool has_shared_endpoint_sibling = false;
        for (const auto& sib : siblings) {
            if (!sib.exempt_a1 && sib.shared_endpoint) {
                has_shared_endpoint_sibling = true;
                break;
            }
        }
        // 性能优化: alpha扫描从6个降至3个。
        std::vector<double> try_alphas = {0.38, 0.30, 0.24};
        int best_cross = sampledSiblingCrossCount(
            initial, ensure_sampled_siblings(), true, kClusterEndpointTol);
        BezierCurve best_curve = initial;
        double best_alpha = 0.4;
        for (double alpha : try_alphas) {
            if (best_cross == 0) break;
            // 用更短把手重建单段曲线，同时保持端点G1。
            BezierCurve c_short = OrdinaryCurveInitializer().buildSingleCubic(
                p0, t0, p1, t1, alpha);
            if (c_short.empty())
                continue;
            BezierSegment shorter = c_short.segs.front();
            // 确认缩短后的弧线仍不穿越障碍。
            bool arc_ok = true;
            if (sdf.valid()) {
                for (int k=1; k<16; ++k) {
                    auto kv = sdf.queryWithGrad(shorter.evaluate((double)k/16));
                    if (kv.first < 0) { arc_ok = false; break; }
                }
            }
            if (arc_ok) {
                Vec2d chord = p1 - p0;
                double chord_len = chord.norm();
                double short_turn_strength = 0.0;
                if (chord_len > 1e-8 && t0.norm() > 1e-8)
                    short_turn_strength = std::abs(cross2d(t0.normalized(), chord.normalized()));
                if (!isNonUTurnTurnShapeAcceptable(c_short, chord_len, short_turn_strength))
                    continue;
                if (assessCurveRisk(c_short, input, sdf, ensure_sampled_siblings(), enforce_fence).physical())
                    continue;
                int cross_count = sampledSiblingCrossCount(
                    c_short, ensure_sampled_siblings(), true, kClusterEndpointTol);
                if (cross_count < best_cross ||
                    (cross_count == best_cross && alpha > best_alpha)) {
                    best_cross = cross_count;
                    best_curve = c_short;
                    best_alpha = alpha;
                }
            }
        }
        initial = best_curve;
        if (best_cross > 0 && has_shared_endpoint_sibling) {
            for (double alpha : {0.18, 0.14, 0.10, 0.06}) {
                BezierCurve c_short;
                c_short.segs.push_back(makeCubicG1(
                    p0, t0.norm()>1e-8 ? t0.normalized() : Vec2d(1,0),
                    p1, t1.norm()>1e-8 ? t1.normalized() : Vec2d(1,0),
                    alpha));
                Vec2d chord = p1 - p0;
                double chord_len = chord.norm();
                double short_turn_strength = 0.0;
                if (chord_len > 1e-8 && t0.norm() > 1e-8)
                    short_turn_strength = std::abs(cross2d(t0.normalized(), chord.normalized()));
                if (!isNonUTurnTurnShapeAcceptable(c_short, chord_len, short_turn_strength))
                    continue;
                int cross_count = sampledSiblingCrossCount(
                    c_short, ensure_sampled_siblings(), true, kClusterEndpointTol);
                if (cross_count < best_cross) {
                    initial = c_short;
                    best_cross = cross_count;
                    if (best_cross == 0)
                        break;
                }
            }
        }
    } else if (is_uturn_geom &&
               sampledSiblingCrossCount(
                   initial, ensure_sampled_siblings(), true, kClusterEndpointTol) > 0) {
        // ── U型调头残余形态收口 ─────────────────────────────────────
        // 复杂 Boundary 链场景中初始多约束求解已使用当前完整兄弟集合；
        // 普通/纯几何场景仍保留二次搜索以维持既有U型三段形态回归。
        const auto& sampled = ensure_sampled_siblings();
        double uturn_chord_len = (p1 - p0).norm();
        auto uturn_arch_shape_ok = [&](const BezierCurve& curve) {
            if (uturn_chord_len < 10.0)
                return true;
            double arc_chord = curve.arcLength() / uturn_chord_len;
            return arc_chord > 1.40 && curve.maxCurvature(20) < 0.50;
        };
        if (!complex_boundary_fast_path) {
            BezierCurve ref_arch = UTurnCurveInitializer().buildAligned(
                p0, t0, p1, t1,
                Vec2d{-t0.normalized().y(), t0.normalized().x()},
                0.0, 1.0, uturn_min_lead0, uturn_min_lead1);
            UTurnSolverResult uturn_sol = UTurnMultiConstraintSearch().search(
                p0, t0, p1, t1, input, sampled, ref_arch,
                makeUTurnSearchBackend(input, sdf, sampled),
                uturn_min_lead0, uturn_min_lead1,
                0, 0, 0.0,  // 修复阶段不额外施加深度/横向偏好。
                &uturn_clearance_crosswalks);
            if (std::isfinite(uturn_sol.cost) &&
                uturn_sol.sibling_crosses <
                    sampledSiblingCrossCount(initial, sampled, true, kClusterEndpointTol)) {
                initial = uturn_sol.curve;
            }
        }
        if (!isPureGeometricInput(input)) {
            BezierCurve segmented_arch =
                UTurnCurveInitializer().buildSegmented(
                    p0, t0, p1, t1, uturn_min_lead0, uturn_min_lead1,
                    2.0 / 3.0, uturn_shared_endpoint_stagger,
                    uturn_lead0_extra_after_align,
                    uturn_lead1_extra_after_align,
                    uturn_shared_entry_stagger,
                    uturn_shared_exit_stagger);
            CurveRisk segmented_risk = assessCurveRisk(
                segmented_arch, input, sdf, sampled, enforce_fence, true);
            int initial_cross = sampledSiblingCrossCount(
                initial, sampled, true, kClusterEndpointTol);
            if (uturn_arch_shape_ok(segmented_arch) &&
                !segmented_risk.physical() &&
                !curveSelfIntersectsBusiness(segmented_arch, 1.0) &&
                (segmented_risk.sibling_crosses < initial_cross ||
                 (segmented_risk.sibling_crosses == initial_cross &&
                  segmented_arch.numSegments() > initial.numSegments() &&
                  segmented_arch.maxCurvature(20) < 3.0) ||
                 (segmented_risk.sibling_crosses == initial_cross &&
                  segmented_arch.maxCurvature(20) + 0.02 * segmented_arch.arcLength() <
                      initial.maxCurvature(20) + 0.02 * initial.arcLength()))) {
                initial = segmented_arch;
            }
        }
    }

    const auto& sampled_for_gate = ensure_sampled_siblings();
    CurveRisk risk = assessCurveRisk(
        initial, input, sdf, sampled_for_gate, enforce_fence, is_uturn_geom);
    // 单条生成阶段不再使用 waypoint 调头拓扑候选；它会把无避让的大型调头
    // 压成两段折弧。后续显式三段候选在同一风险门禁下完成形态收口。
    bool fixed_shape_cross = false;
    if (risk.sibling_crosses > 0) {
        for (const auto& sib : siblings) {
            if (!sib.fixed_shape || sib.exempt_a1)
                continue;
            if (curvesHaveForbiddenSameClusterIntersection(
                    initial, sib.curve, kClusterEndpointTol)) {
                fixed_shape_cross = true;
                break;
            }
        }
    }
    // 只要存在物理风险、固有形态冲突或同簇交叉，就进入优化/修复流程；
    // 该门禁不区分转向类型，避免U型调头在初次生成时跳过同簇约束。
    bool needs_optimization = risk.physical() || fixed_shape_cross ||
        risk.sibling_crosses > 0;
    if (out_physical_risk)
        *out_physical_risk = false;
    PreCheckResult pre;
    if (needs_optimization && sdf_coarse.valid() && !input.area.is_rough) {
        pre = preCheck(sdf_coarse, input.area.geometry, p0, p1, lw, input.boundaries);
        if (pre.type == ViolationInfo::InfeasibilityType::TopologicalBlock)
            return makeFallbackCurve(pre, conn, p0, p1);
    }

    if (!needs_optimization) {
        if (is_uturn_geom) {
            BezierCurve segmented = initial;
            searchSegmentedUTurnTwoPass(
                p0, t0, p1, t1, input, sdf, sampled_for_gate, enforce_fence,
                segmented, uturn_min_lead0, uturn_min_lead1,
                &uturn_clearance_crosswalks, uturn_shared_endpoint_stagger,
                uturn_lead0_extra_after_align, uturn_lead1_extra_after_align,
                uturn_alignment_station,
                uturn_shared_entry_stagger, uturn_shared_exit_stagger,
                uturn_family_stagger_step);
            initial = segmented;
        }
        if (!is_uturn_geom)
            tryShapeSafeSingleCubic(p0, t0, p1, t1, input, sdf, sampled_for_gate, enforce_fence, initial);
        if (!is_uturn_geom) {
            tryNaturalSingleCubicIfSafe(
                p0, t0, p1, t1, input, sdf, sampled_for_gate, enforce_fence, initial);
            tryNaturalStraightSingleCubicIfSafe(
                p0, t0, p1, t1, input, sdf, sampled_for_gate, enforce_fence, initial);
        } else if (initial.numSegments() != 3 || initial.maxCurvature(40) >= 3.0) {
            BezierCurve segmented_arch =
                UTurnCurveInitializer().buildSegmented(
                    p0, t0, p1, t1, uturn_min_lead0, uturn_min_lead1,
                    2.0 / 3.0, uturn_shared_endpoint_stagger,
                    uturn_lead0_extra_after_align,
                    uturn_lead1_extra_after_align,
                    uturn_shared_entry_stagger,
                    uturn_shared_exit_stagger);
            CurveRisk current_uturn_risk =
                assessCurveRisk(
                    initial, input, sdf, sampled_for_gate, enforce_fence, true);
            CurveRisk segmented_risk =
                assessCurveRisk(
                    segmented_arch, input, sdf, sampled_for_gate,
                    enforce_fence, true);
            if (!segmented_arch.empty() &&
                !curveSelfIntersectsBusiness(segmented_arch, 1.0) &&
                !segmented_risk.physical() &&
                segmented_risk.sibling_crosses <= current_uturn_risk.sibling_crosses &&
                segmented_arch.maxCurvature(40) < 3.0) {
                initial = segmented_arch;
            }
        }
        setConnectivityCurveGeometry(cc, initial);
        validate(cc, input, sdf);
        return cc;
    }

    if (is_uturn_geom) {
        // 性能优化: 移除冗余的二次多约束求解
        // (原代码再次调用多约束搜索器，但首次调用已穷举搜索，
        //  二次调用极少改进且耗时显著)
        int current_cross = sampledSiblingCrossCount(
            initial, sampled_for_gate, true, kClusterEndpointTol);
        if (current_cross == 0 && !risk.physical() && !fixed_shape_cross) {
            // 初解已无同簇交叉且无物理风险时，仍优先把大型调头表达成
            // 首直行 + 中间弧 + 尾直行，避免纯 de Casteljau 拆分保留畸形单弧。
            BezierCurve segmented = initial;
            searchSegmentedUTurnTwoPass(
                p0, t0, p1, t1, input, sdf, sampled_for_gate, enforce_fence,
                segmented, uturn_min_lead0, uturn_min_lead1,
                &uturn_clearance_crosswalks, uturn_shared_endpoint_stagger,
                uturn_lead0_extra_after_align, uturn_lead1_extra_after_align,
                uturn_alignment_station,
                uturn_shared_entry_stagger, uturn_shared_exit_stagger,
                uturn_family_stagger_step);
            initial = segmented;
            setConnectivityCurveGeometry(cc, initial);
            validate(cc, input, sdf);
            return cc;
        }
    }

    BezierCurve safe_single = initial;
    if (!is_uturn_geom &&
        tryPhysicalSafeSingleCubic(p0, t0, p1, t1, input, sdf, sampled_for_gate, safe_single)) {
        setConnectivityCurveGeometry(cc, safe_single);
        validate(cc, input, sdf);
        return cc;
    }

    BezierCurve bypass_candidate = initial;
    if (!is_uturn_geom && tryObstacleBypassCandidate(
            p0, t0, p1, t1, input, sdf, sampled_for_gate, bypass_candidate)) {
        setConnectivityCurveGeometry(cc, bypass_candidate);
        validate(cc, input, sdf);
        return cc;
    }

    // ── 优化 ────────────────────────────────────────────────────────────────
    // 普通单段曲线在无物理避让时只允许两个内部控制点沿端点方向轴移动；
    // U型调头和已进入物理避让路径的曲线保留例外自由度。
    CurveOptimizationOptions optimization_options;
    optimization_options.enforce_fence = enforce_fence;
    optimization_options.constrain_single_cubic_axes =
        !is_uturn_geom && !risk.physical();
    BezierCurve opt = CurveOptimizer(solver_).optimize(
        initial, input, sdf, siblings, t0, t1, optimization_options);

    OptimizationResultOptions result_options;
    result_options.geometric_uturn = is_uturn_geom;
    BezierCurve final_c = OptimizationResultProcessor(solver_).process(
        initial, opt, sdf, input.area.geometry,
        p0, t0, p1, t1, result_options);

    const UTurnEnvelopeConstraint envelope_constraint;
    bool shape_risk = is_uturn_geom &&
        (envelope_constraint.exceeds(final_c, initial, p0, p1) ||
         envelope_constraint.collapses(final_c, initial, p0, p1));
    // 调头包络坍缩必须回退，否则优化器可能把调头压扁并重新引入形态违规。
    if (is_uturn_geom && shape_risk) {
        BezierCurve fallback = initial;
        if (tryBoundedUTurnCandidate(p0, t0, p1, t1, input, sdf, sampled_for_gate,
                                     initial, fallback, uturn_min_lead0, uturn_min_lead1)) {
            final_c = fallback;
        } else {
            final_c = initial;
        }
        shape_risk = false;
    }
    if (is_uturn_geom) {
        searchSegmentedUTurnTwoPass(
            p0, t0, p1, t1, input, sdf, sampled_for_gate, enforce_fence,
            final_c, uturn_min_lead0, uturn_min_lead1,
            &uturn_clearance_crosswalks, uturn_shared_endpoint_stagger,
            uturn_lead0_extra_after_align, uturn_lead1_extra_after_align,
            uturn_alignment_station,
            uturn_shared_entry_stagger, uturn_shared_exit_stagger,
            uturn_family_stagger_step);
    }
    const PhysicalRiskAuditor physical_auditor =
        [&](const BezierCurve& candidate) {
            const CurveRisk risk = assessCurveRisk(
                candidate, input, sdf, sampled_for_gate,
                enforce_fence, is_uturn_geom);
            PhysicalRiskSnapshot snapshot;
            snapshot.obstacle = risk.obstacle;
            snapshot.boundary = risk.boundary;
            snapshot.fence = risk.fence;
            snapshot.sibling_crosses = risk.sibling_crosses;
            return snapshot;
        };
    const PhysicalRepairAttempt boundary_repair =
        [&](const BezierCurve& current, BezierCurve& repaired) {
            repaired = current;
            if (!tryBoundarySafeCandidate(
                    p0, t0, p1, t1, input, sdf, sampled_for_gate,
                    enforce_fence, repaired))
                return false;
            return rawFixedSiblingCrossCount(
                       repaired, siblings, kClusterEndpointTol) <=
                   rawFixedSiblingCrossCount(
                       current, siblings, kClusterEndpointTol);
        };
    const PhysicalRepairAttempt obstacle_repair =
        [&](const BezierCurve& current, BezierCurve& repaired) {
            repaired = current;
            return tryObstacleBypassCandidate(
                p0, t0, p1, t1, input, sdf, sampled_for_gate, repaired);
        };
    const PhysicalRepairResult physical_repair =
        PhysicalRepairCoordinator().repair(
            initial, final_c, is_uturn_geom, physical_auditor,
            boundary_repair, obstacle_repair);
    final_c = physical_repair.curve;
    bool boundary_repaired = physical_repair.boundary_repaired;

    // ── U型调头最终形态/Crosswalk保护 ───────────────────────────────
    // 不再只检查单段曲线顶点附近3m；最终输出必须显式为首直行+单中弧+
    // 尾直行，并逐点确认整个中弧都位于相关 crosswalk 集合外。
    if (is_uturn_geom &&
        (final_c.numSegments() != 3 ||
         !segmentedUTurnMiddleArcClearsCrosswalks(
             final_c, uturn_clearance_crosswalks))) {
        BezierCurve segmented = final_c;
        if (searchSegmentedUTurnTwoPass(
                p0, t0, p1, t1, input, sdf, sampled_for_gate,
                enforce_fence, segmented,
                uturn_min_lead0, uturn_min_lead1,
                &uturn_clearance_crosswalks, uturn_shared_endpoint_stagger,
                uturn_lead0_extra_after_align, uturn_lead1_extra_after_align,
                uturn_alignment_station,
                uturn_shared_entry_stagger, uturn_shared_exit_stagger,
                uturn_family_stagger_step)) {
            final_c = segmented;
        }
    }

    if (!is_uturn_geom)
        tryShapeSafeSingleCubic(p0, t0, p1, t1, input, sdf, sampled_for_gate, enforce_fence, final_c);
    if (!is_uturn_geom && !boundary_repaired) {
        Vec2d chord = p1 - p0;
        double chord_len2 = chord.norm();
        Vec2d T0 = t0.norm() > 1e-8 ? t0.normalized() : (chord_len2 > 1e-8 ? chord.normalized() : Vec2d(1, 0));
        Vec2d T1 = t1.norm() > 1e-8 ? t1.normalized() : (chord_len2 > 1e-8 ? chord.normalized() : Vec2d(1, 0));
        double turn_strength = chord_len2 > 1e-8 ? std::abs(cross2d(T0, chord.normalized())) : 0.0;
        double arc_chord = chord_len2 > 1e-8 ? final_c.arcLength() / chord_len2 : 1.0;
        bool straight_like = turn_strength < 0.25;
        auto constrainedSiblingCrosses = [&](const BezierCurve& candidate, bool fixed_only) {
            int count = 0;
            for (const auto& sib : siblings) {
                if (sib.exempt_a1 || (fixed_only && !sib.fixed_shape))
                    continue;
                if (curvesHaveForbiddenSameClusterIntersection(
                        candidate, sib.curve, kClusterEndpointTol))
                    ++count;
            }
            return count;
        };

        int current_cross = constrainedSiblingCrosses(final_c, false);
        int current_fixed_cross = constrainedSiblingCrosses(final_c, true);
        bool flattened_turn = turn_strength > 0.35 && arc_chord < 1.05;
        bool turn_shape_risk = isNonUTurnTurnShapeRisk(final_c, chord_len2, turn_strength);
        bool fixed_shape_cross_final = current_fixed_cross > 0;
        if (turn_shape_risk || fixed_shape_cross_final) {
            BezierCurve best = final_c;
            int best_cross = current_cross;
            int best_fixed_cross = current_fixed_cross;
            double best_score = std::numeric_limits<double>::infinity();
            bool improved = false;

            std::vector<BezierCurve> candidates;
            candidates.reserve(96);
            auto addCandidate = [&](const BezierCurve& candidate) {
                if (!candidate.empty())
                    candidates.push_back(candidate);
            };
            Vec2d chord_dir = chord.normalized();
            Vec2d chord_perp{-chord_dir.y(), chord_dir.x()};
            double turn_sign = cross2d(T0, chord_dir);
            double natural_side = turn_sign >= 0.0 ? -1.0 : 1.0;
            std::vector<double> sides = straight_like
                ? std::vector<double>{natural_side, -natural_side}
                : std::vector<double>{natural_side};
            std::vector<double> alphas = straight_like
                ? std::vector<double>{0.30, 0.36, 0.42}
                : std::vector<double>{0.42, 0.46, 0.50, 0.55};
            for (double alpha : alphas) {
                BezierCurve single;
                single.segs.push_back(makeCubicG1(p0, T0, p1, T1, alpha));
                addCandidate(single);

                for (double side : sides) {
                    std::vector<double> offsets = straight_like
                        ? std::vector<double>{0.75, 1.25, 1.8, 2.5}
                        : std::vector<double>{2.0, 3.5, 5.0, 6.5};
                    for (double offset : offsets) {
                        Vec2d off = chord_perp * (side * offset);
                        for (double frac : {0.42, 0.50, 0.58}) {
                            Vec2d mid = p0 + chord * frac + off;
                            Vec2d mid_tan = T0 + T1;
                            if (mid_tan.norm() < 1e-8)
                                mid_tan = chord_dir;
                            mid_tan.normalize();
                            BezierCurve shaped;
                            shaped.segs.push_back(makeCubicG1(p0, T0, mid, mid_tan, alpha));
                            shaped.segs.push_back(makeCubicG1(mid, mid_tan, p1, T1, alpha));
                            addCandidate(shaped);
                        }

                        double lead = std::max(2.0, std::min(7.0, chord_len2 * 0.20));
                        Vec2d q0 = p0 + T0 * lead + off;
                        Vec2d q1 = p1 - T1 * lead + off;
                        std::vector<Vec2d> pts2 = {p0, q0, q1, p1};
                        std::vector<Vec2d> tans2 = {T0, chord_dir, chord_dir, T1};
                        addCandidate(makeCurveFromKnots(pts2, tans2, alpha));
                    }
                }
            }

            for (const auto& candidate : candidates) {
                if (candidate.empty() || curveSelfIntersectsBusiness(candidate, 1.0))
                    continue;
                double cand_arc_chord = candidate.arcLength() / chord_len2;
                double cand_maxk40 = candidate.maxCurvature(40);
                if (straight_like && (cand_arc_chord > 1.08 || cand_maxk40 > 1.0))
                    continue;
                if (!straight_like &&
                    !isNonUTurnTurnShapeAcceptable(candidate, chord_len2, turn_strength))
                    continue;
                if (!flattened_turn && !straight_like && cand_arc_chord < 1.03)
                    continue;
                if (!flattened_turn && cand_arc_chord > (straight_like ? 1.20 : 1.35))
                    continue;
                CurveRisk candidate_risk = assessCurveRisk(
                    candidate, input, sdf, sampled_for_gate, enforce_fence);
                if (candidate_risk.physical())
                    continue;
                int fixed_cross_count = constrainedSiblingCrosses(candidate, true);
                int cross_count = constrainedSiblingCrosses(candidate, false);
                if (fixed_cross_count > current_fixed_cross || cross_count > current_cross)
                    continue;
                double target_arc = (flattened_turn || !straight_like) ? 1.10 : 1.02;
                double score = 0.5 * std::abs(cand_arc_chord - target_arc) +
                    cand_maxk40 + 0.02 * candidate.arcLength();
                if (!improved ||
                    fixed_cross_count < best_fixed_cross ||
                    (fixed_cross_count == best_fixed_cross && cross_count < best_cross) ||
                    (fixed_cross_count == best_fixed_cross && cross_count == best_cross &&
                     score < best_score)) {
                    best = candidate;
                    best_cross = cross_count;
                    best_fixed_cross = fixed_cross_count;
                    best_score = score;
                    improved = true;
                }
            }
            if (improved)
                final_c = best;
        }
    }
    if (!is_uturn_geom) {
        Vec2d chord = p1 - p0;
        double chord_len2 = chord.norm();
        if (chord_len2 > 1e-8) {
            Vec2d T0 = t0.norm() > 1e-8 ? t0.normalized() : chord.normalized();
            double turn_strength = std::abs(cross2d(T0, chord.normalized()));
            bool straight_like = turn_strength < 0.25;
            double arc_chord = final_c.arcLength() / chord_len2;
            bool shape_bad = straight_like
                ? (arc_chord > 1.08 || final_c.maxCurvature(40) > 1.0)
                : isNonUTurnTurnShapeRisk(final_c, chord_len2, turn_strength);
            if (shape_bad) {
                BezierCurve repaired = initial;
                if (tryBoundarySafeCandidate(
                        p0, t0, p1, t1, input, sdf, sampled_for_gate,
                        enforce_fence, repaired, true)) {
                    CurveRisk repaired_risk = assessCurveRisk(
                        repaired, input, sdf, sampled_for_gate, enforce_fence);
                    double repaired_arc_chord = repaired.arcLength() / chord_len2;
                    bool repaired_shape_ok = straight_like
                        ? (repaired_arc_chord <= 1.08 &&
                           repaired.maxCurvature(40) <= 1.0)
                        : isNonUTurnTurnShapeAcceptable(
                              repaired, chord_len2, turn_strength);
                    if (!repaired_risk.physical() && repaired_shape_ok)
                        final_c = repaired;
                }
            }
        }
    }
    CurveRisk exit_risk = assessCurveRisk(
        final_c, input, sdf, sampled_for_gate, enforce_fence, is_uturn_geom);
    if (!is_uturn_geom && exit_risk.boundary) {
        BezierCurve repaired = final_c;
        if (tryBoundarySafeCandidate(
                p0, t0, p1, t1, input, sdf, sampled_for_gate, enforce_fence,
                repaired)) {
            int current_fixed_raw_cross =
                rawFixedSiblingCrossCount(final_c, siblings, kClusterEndpointTol);
            int repaired_fixed_raw_cross =
                rawFixedSiblingCrossCount(repaired, siblings, kClusterEndpointTol);
            if (repaired_fixed_raw_cross <= current_fixed_raw_cross) {
                final_c = repaired;
                exit_risk = assessCurveRisk(
                    final_c, input, sdf, sampled_for_gate, enforce_fence, false);
            }
        }
    }
    if (!is_uturn_geom && !exit_risk.physical()) {
        tryNaturalSingleCubicIfSafe(
            p0, t0, p1, t1, input, sdf, sampled_for_gate, enforce_fence, final_c);
        tryNaturalStraightSingleCubicIfSafe(
            p0, t0, p1, t1, input, sdf, sampled_for_gate, enforce_fence, final_c);
        exit_risk = assessCurveRisk(
            final_c, input, sdf, sampled_for_gate, enforce_fence, false);
    }
    if (is_uturn_geom) {
        if (searchSegmentedUTurnTwoPass(
                p0, t0, p1, t1, input, sdf, sampled_for_gate, enforce_fence,
                final_c, uturn_min_lead0, uturn_min_lead1,
                &uturn_clearance_crosswalks, uturn_shared_endpoint_stagger,
                uturn_lead0_extra_after_align, uturn_lead1_extra_after_align,
                uturn_alignment_station,
                uturn_shared_entry_stagger, uturn_shared_exit_stagger,
                uturn_family_stagger_step))
            exit_risk = assessCurveRisk(
                final_c, input, sdf, sampled_for_gate, enforce_fence, true);
    }
    if (out_physical_risk) {
        *out_physical_risk = exit_risk.physical() || shape_risk;
    }
    setConnectivityCurveGeometry(cc, final_c);
    validate(cc, input, sdf);
    // 三段式 U-turn 的窄通道预警不能覆盖已通过实际 Boundary、横道、
    // 曲率和自交审计的结果：预检查沿端点弦估计宽度，而 U-turn 会先沿
    // 首尾直段通过该区域，直接把它标成 WarnA2 会把合规掉头误报为未修复。
    if (pre.narrow_passage && cc.status == CurveStatus::OK &&
        !(is_uturn_geom && final_c.numSegments() == 3))
        cc.status = CurveStatus::WarnA2;
    return cc;
}

// 构造单曲线流水线并生成指定连接，注入会话级求解器、簇排序器和审计器。
// 这是 run 内部批量生成单条连接的统一入口。
ConnectivityCurve ConnectivityGenerationSession::generateOne(
    const CurveGenerationContext& context, const std::vector<SiblingCurve>& siblings, bool* out_physical_risk) {
    const SingleCurveGenerationPipeline::Validator validator =
        [&](ConnectivityCurve& curve, const IntersectionInput& input,
            const SDFField& sdf) { validate(curve, input, sdf); };
    return SingleCurveGenerationPipeline(
        solver_, cluster_solver_, direction_cfg_, validator).run(
            context, siblings, out_physical_risk);
}

// 运行完整路口曲线生成流程：规范化输入，构建簇拓扑和 U-turn 家族快照，按计划
// 生成连接，处理固定形态、物理风险、同簇交叉、形态恢复及最终标注，并返回结果。
// out_ms 可选接收总耗时；函数不修改调用者传入的 raw_input 和 SDF 对象内容。
std::vector<ConnectivityCurve> ConnectivityGenerationSession::run(
    const IntersectionInput& raw_input, SDFField& sdf, double* out_ms) {
    auto t0 = std::chrono::steady_clock::now();

    // 初始化曲线首尾方向
    IntersectionInput input = raw_input;
    ConnectivityDirectionNormalizer(input, direction_cfg_);

    bool pure_geometry = input.mode != 2 && input.obstacles.empty() && input.boundaries.empty();
    // mode=2 也可复用密集 Boundary 快速路径；该路径只缩减候选枚举，不放宽 RoadEdge 的相交与净距硬约束。
    bool complex_boundary_fast_path = input.obstacles.empty() && input.boundaries.size() >= 20;

    SDFField sdf_coarse;
    auto roi = input.area.geometry.empty() ? BoundingBox2d{} : input.area.geometry.bbox();
    if (roi.width() < 1) {
        for (auto& l : input.lanes) {
            for (auto& p : l.geometry.points)
                roi.expand(p);
        }
        roi.min_pt -= Vec2d(20, 20);
        roi.max_pt += Vec2d(20, 20);
    }
    if (!input.obstacles.empty())
        sdf_coarse.build(roi, input.obstacles, 0.5, obstacleAvoidanceClearanceForMode(input.mode));

    // 构建同簇顺序求解器，供兄弟曲线约束、结构性交叉豁免和修复阶段复用。
    cluster_solver_.build(input.connectivities, input.lanes, input.lane_groups, input.crosswalks);
    const ClusterTopology cluster_topology = cluster_solver_.topology();
    SceneContext scene(input);
    scene.bindSdf(&sdf, &sdf_coarse);
    scene.bindClusterTopology(&cluster_topology);

    // 家族几何只依赖本次规范化输入、拓扑和对齐配置。
    // 一次构建后供初始生成、重生成及修复阶段共享，避免各阶段独立求站位造成重复计算或快照漂移。
    std::unordered_map<ConnId, UTurnFamilyInfo> uturn_family_snapshots;
    std::unordered_map<ConnId, UTurnFamilyInfo> paired_uturn_family_snapshots;
    const UTurnFamilyBuilder family_builder;
    for (const auto& conn : input.connectivities) {
        uturn_family_snapshots[conn.id] = family_builder.build(
            conn, input, direction_cfg_.uturn_alignment_scope,
            &cluster_solver_, false);
        paired_uturn_family_snapshots[conn.id] = family_builder.build(
            conn, input, direction_cfg_.uturn_alignment_scope,
            &cluster_solver_, direction_cfg_.uturn_alignment_scope == UTurnAlignmentScope::LaneEndpoint);
    }
    const auto familySnapshot = [&](const ConnId& id) -> const UTurnFamilyInfo& {
        return uturn_family_snapshots.at(id);
    };
    const auto pairedFamilySnapshot = [&](const ConnId& id) -> const UTurnFamilyInfo& {
        return paired_uturn_family_snapshots.at(id);
    };

    // 按直行、左右转、U型调头及内外侧关系构建生成顺序。
    const GenerationPlan plan = GenerationPlanner().build(scene, cluster_solver_);
    const RepairBudget repair_budget;

    std::vector<ConnectivityCurve> results;
    results.reserve(input.connectivities.size());
    std::unordered_map<ConnId, BezierCurve> done;
    std::unordered_set<ConnId> physical_risk_ids;
    std::unordered_set<ConnId> preserved_fixed_ids;
    std::unordered_set<ConnId> roadedge_boundary_repair_ids;
    for (const auto& conn : input.connectivities) {
        // 最新掉头形态约束优先于输入 fixed_shape: 几何掉头必须统一
        // 重新生成为“首端直行段 + 单段掉头弧 + 尾端直行段”。
        if (!hasFixedGeometry(conn) || isGeometricUTurnConn(conn, input) ||
            fixedGeometryHitsObstacle(conn, input))
            continue;
        auto cc = makeFixedGeometryCurve(conn, input, sdf);
        if (cc.curve) {
            done[conn.id] = *cc.curve;
            preserved_fixed_ids.insert(conn.id);
        }
        results.push_back(std::move(cc));
    }

    auto conflicting_fixed_ids = collectConflictingPreservedFixedIds(
        results, cluster_solver_, preserved_fixed_ids, kClusterEndpointTol);
    auto shape_blocking_fixed_ids = collectFixedIdsBlockingNaturalShapes(
        input, sdf, results, cluster_solver_, preserved_fixed_ids, kClusterEndpointTol);
    conflicting_fixed_ids.insert(
        shape_blocking_fixed_ids.begin(), shape_blocking_fixed_ids.end());
    if (!conflicting_fixed_ids.empty()) {
        for (const auto& cid : conflicting_fixed_ids) {
            preserved_fixed_ids.erase(cid);
            done.erase(cid);
        }
        results.erase(
            std::remove_if(results.begin(), results.end(),
                           [&](const ConnectivityCurve& cc) {
                               return conflicting_fixed_ids.count(cc.id) > 0;
                           }),
            results.end());
    }

    for (const auto& group : plan.batches) {
        for (auto& cid : group.conn_ids) {
            if (preserved_fixed_ids.count(cid))
                continue;
            const Connectivity* conn = scene.view.connectivity(cid);
            if (!conn) continue;

            auto sibs = buildSiblings(
                cid, done, cluster_solver_, input.connectivities, pure_geometry, &preserved_fixed_ids);

            bool phys_risk = false;
            auto conn_t0 = std::chrono::steady_clock::now();
            const CurveGenerationContext generation_context =
                CurveGenerationContextBuilder().build(
                    scene, *conn, direction_cfg_.uturn_alignment_scope,
                    &cluster_solver_, &familySnapshot(conn->id));
            auto cc = generateOne(generation_context, sibs, &phys_risk);
            if (isgProfile()) {
                double ms = std::chrono::duration<double, std::milli>(
                    std::chrono::steady_clock::now() - conn_t0).count();
                fprintf(stderr, "[ISG_PROFILE] conn %s generated in %.3f ms\n",
                        cid.c_str(), ms);
            }
            if (hasFixedGeometry(*conn) && !preserved_fixed_ids.count(cid))
                cc.fixed_shape = false;
            if (cc.curve) {
                done[cid] = *cc.curve;
                if (phys_risk)
                    physical_risk_ids.insert(cid);
            }

            results.push_back(std::move(cc));
        }
        // 每个优先级组完成后，把贴近障碍的交叉标为软约束，避免后续过度修复。
        cluster_solver_.checkAndMarkA2(done, sdf, 1.5);
    }
    if (isgProfile()) {
        fprintf(stderr, "[ISG_PROFILE] initial generation stage %.3f ms\n",
                std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count());
    }

    // 同端点同簇掉头交叉必须在最终标注前成对修复。即使没有 Crosswalk、Obstacle 或 Boundary，纯几何输入也必须执行该收口；
    // 否则只能依靠单条候选，无法同时保证同簇非交叉和出口侧顺序。
    auto repairMode2UTurnPairs = [&]() {
        bool any_changed = false;
        for (int pass = 0; pass < repair_budget.mode2_uturn_pair_passes; ++pass) {
            auto result_idx = resultIndexById(results);
            auto all_done = curveMapFromResults(results);
            bool changed = false;
            int repaired_pairs = 0;
            auto other_constrained_cross_count = [&](
                    const ConnId& cid, const ConnId& other_id, const BezierCurve& candidate) {
                int count = 0;
                for (const auto& kv : all_done) {
                    if (kv.first == cid || kv.first == other_id)
                        continue;
                    if (cluster_solver_.exemptionOf(cid, kv.first) ==
                        CrossExemption::StructuralCross)
                        continue;
                    if (curvesHaveForbiddenSameClusterIntersection(
                            candidate, kv.second, kClusterEndpointTol))
                        ++count;
                }
                return count;
            };
            auto introduces_new_constrained_cross = [&](
                    const ConnId& cid, const ConnId& paired_id,
                    const BezierCurve& current_curve,
                    const BezierCurve& candidate) {
                for (const auto& kv : all_done) {
                    if (kv.first == cid || kv.first == paired_id)
                        continue;
                    if (cluster_solver_.exemptionOf(cid, kv.first) ==
                        CrossExemption::StructuralCross)
                        continue;
                    bool currently_crosses =
                        curvesHaveForbiddenSameClusterIntersection(current_curve, kv.second, kClusterEndpointTol);
                    bool candidate_crosses =
                        curvesHaveForbiddenSameClusterIntersection(candidate, kv.second, kClusterEndpointTol);
                    if (!currently_crosses && candidate_crosses)
                        return true;
                }
                return false;
            };
            auto preservesSharedExitUTurnOrder = [&](
                    const ConnId& a_id, const BezierCurve& a,
                    const ConnId& b_id, const BezierCurve& b) {
                if (a.numSegments() != 3 || b.numSegments() != 3 ||
                    (a.endPt() - b.endPt()).norm() > 0.20)
                    return true;
                int expected_side = cluster_solver_.expectedSideOf(a_id, b_id);
                Vec2d ref_perp = cluster_solver_.refPerpOf(a_id, b_id);
                if (expected_side == 0 || ref_perp.norm() < 1e-8)
                    return true;
                ref_perp.normalize();
                // 对同出口 U-turn，middle arc 的 ctrl[2] 是靠近退出直行段的
                // 控制点。它必须保留簇排序指定的一侧，避免消交候选把内外层
                // 翻转；小于 1cm 的数值噪声不视作有效排序。
                double side = (a.segs[1].ctrl[2] - b.segs[1].ctrl[2]).dot(ref_perp);
                return expected_side * side >= 0.01;
            };
            auto segmented_pair_candidates = [&](const Connectivity& conn) {
                std::vector<BezierCurve> candidates;
                auto entry = scene.view.entryFrame(conn.entry_lane_id);
                auto exit_ = scene.view.exitFrame(conn.exit_lane_id);
                Vec2d p0 = entry.first;
                Vec2d p1 = exit_.first;
                double chord_len = (p1 - p0).norm();
                double segmented_max_curvature = segmentedUTurnMaxCurvatureLimit(p0, entry.second, p1, exit_.second);
                Vec2d T0 = entry.second.norm() > 1e-8
                    ? entry.second.normalized() : Vec2d(1, 0);
                Vec2d T1 = exit_.second.norm() > 1e-8
                    ? exit_.second.normalized() : -T0;
                Vec2d exit_back = -T1;
                Vec2d axis = T0 + exit_back;
                if (axis.norm() < 1e-8)
                    axis = T0;
                axis.normalize();
                if (axis.dot(T0) < 0.0)
                    axis = -axis;
                const UTurnFamilyInfo& family_info =
                    pairedFamilySnapshot(conn.id);
                const double alignment_station = family_info.aligned_station;
                double shared_entry_stagger = 0.0;
                double shared_exit_stagger = 0.0;
                double shared_endpoint_stagger = uturnSharedEndpointStagger(
                    conn, input, cluster_solver_,
                    direction_cfg_.uturn_alignment_scope, nullptr,
                    &shared_entry_stagger, &shared_exit_stagger);
                // 家族平齐站位是唯一站位，pair repair 不得单独轴向前推。
                const std::vector<double> extra_options = {0.0};
                const std::vector<double> alpha_options = input.mode == 2
                    ? std::vector<double>{2.0 / 3.0, 0.85, 1.0, 1.25, 1.50,
                                           1.75, 2.0, 0.50, 0.38, 0.28, 0.16}
                    : std::vector<double>{2.0 / 3.0, 0.50, 0.38, 0.28};
                // 禁止额外分米级错开覆盖 N*0.01m 家族分档。
                const std::vector<double> stagger_options = {0.0};
                const size_t candidate_cap = input.mode == 2 ? 96 : 32;
                for (double extra0 : extra_options) {
                    for (double extra1 : extra_options) {
                        for (double arc_alpha : alpha_options) {
                            for (double extra_stagger : stagger_options) {
                            // 即使只共享一侧端点，q0/q1 也必须使用相同总错开量。
                            // 继续传旧的 entry/exit 值会覆盖
                            // buildSegmented 的通用错开量，使
                            // extra_stagger 实际失效。
                            double total_stagger = std::max(
                                0.0, shared_endpoint_stagger + extra_stagger);
                            BezierCurve c = UTurnCurveInitializer().buildSegmented(
                                p0, entry.second, p1, exit_.second,
                                family_info.lead0, family_info.lead1, arc_alpha,
                                total_stagger, extra0, extra1,
                                total_stagger, total_stagger);
                            if (c.empty() || c.numSegments() != 3 ||
                                curveSelfIntersectsBusiness(c, 1.0) ||
                                c.maxCurvature(40) >= segmented_max_curvature ||
                                !segmentLooksStraight(c.segs.front()) ||
                                !segmentLooksStraight(c.segs.back()) ||
                                (chord_len >= 1.0 && !segmentedUTurnMiddleArcLooksRound(c, axis)))
                                continue;
                            Vec2d q0 = c.segs.front().ctrl[3];
                            Vec2d q1 = c.segs.back().ctrl[0];
                            if (std::abs((q0 - q1).dot(axis)) > 0.05)
                                continue;
                            if (std::isfinite(alignment_station) &&
                                (std::abs(q0.dot(axis) - alignment_station) > 0.05 ||
                                 std::abs(q1.dot(axis) - alignment_station) > 0.05))
                                continue;
                            if (c.arcLength() / std::max(1e-6, (p1 - p0).norm()) < 1.35)
                                continue;
                            if (assessCurveRisk(c, input, sdf, {}, true).physical())
                                continue;
                            candidates.push_back(std::move(c));
                            if (candidates.size() >= candidate_cap)
                                return candidates;
                            }
                        }
                    }
                }
                return candidates;
            };

            for (size_t ia_pos = 0; ia_pos < results.size(); ++ia_pos) {
                for (size_t ib_pos = ia_pos + 1; ib_pos < results.size(); ++ib_pos) {
                    auto ia = result_idx.find(results[ia_pos].id);
                    auto ib = result_idx.find(results[ib_pos].id);
                    if (ia == result_idx.end() || ib == result_idx.end())
                        continue;
                    if (!results[ia->second].curve || !results[ib->second].curve)
                        continue;
                    const Connectivity* ca_conn =
                        scene.view.connectivity(results[ia->second].id);
                    const Connectivity* cb_conn =
                        scene.view.connectivity(results[ib->second].id);
                    if (!ca_conn || !cb_conn ||
                        !isGeometricUTurnConn(*ca_conn, input) ||
                        !isGeometricUTurnConn(*cb_conn, input))
                        continue;
                    bool shared_endpoint =
                        (results[ia->second].curve->startPt() -
                         results[ib->second].curve->startPt()).norm() <= 0.20 ||
                        (results[ia->second].curve->endPt() -
                         results[ib->second].curve->endPt()).norm() <= 0.20;
                    if (!shared_endpoint)
                        continue;
                    if (cluster_solver_.exemptionOf(
                            results[ia->second].id,
                            results[ib->second].id) == CrossExemption::StructuralCross)
                        continue;
                    if (!curvesIntersectBusiness(
                            *results[ia->second].curve,
                            *results[ib->second].curve, kClusterEndpointTol))
                        continue;

                    std::vector<BezierCurve> cand_a =
                        segmented_pair_candidates(*ca_conn);
                    std::vector<BezierCurve> cand_b =
                        segmented_pair_candidates(*cb_conn);
                    if (isgDebugPairRepair()) {
                        fprintf(stderr,
                                "[PAIR-UTURN] %s-%s cand_a=%zu cand_b=%zu\n",
                                results[ia->second].id.c_str(),
                                results[ib->second].id.c_str(),
                                cand_a.size(), cand_b.size());
                    }
                    BezierCurve best_a;
                    BezierCurve best_b;
                    bool repair_only_a = false;
                    bool repair_only_b = false;
                    bool have_pair_repair = false;
                    double best_pair_score = std::numeric_limits<double>::infinity();
                    int current_other_cross =
                        other_constrained_cross_count(
                            results[ia->second].id, results[ib->second].id, *results[ia->second].curve) +
                        other_constrained_cross_count(
                            results[ib->second].id, results[ia->second].id, *results[ib->second].curve);
                    int dbg_pair_hit = 0;
                    int dbg_other_excess = 0;
                    int dbg_viable_pair = 0;
                    auto consider_single = [&](
                            bool change_a, const BezierCurve& candidate) {
                        const ConnId& change_id = change_a
                            ? results[ia->second].id : results[ib->second].id;
                        const ConnId& fixed_id = change_a
                            ? results[ib->second].id : results[ia->second].id;
                        const BezierCurve& fixed_curve = change_a
                            ? *results[ib->second].curve : *results[ia->second].curve;
                        if (!preservesSharedExitUTurnOrder(change_id, candidate, fixed_id, fixed_curve))
                            return;
                        if (curvesHaveForbiddenSameClusterIntersection(candidate, fixed_curve, kClusterEndpointTol)) {
                            ++dbg_pair_hit;
                            return;
                        }
                        const BezierCurve& current_curve = change_a
                            ? *results[ia->second].curve : *results[ib->second].curve;
                        if (introduces_new_constrained_cross(change_id, fixed_id, current_curve, candidate)) {
                            ++dbg_other_excess;
                            return;
                        }
                        ++dbg_viable_pair;
                        int candidate_other_cross =
                            other_constrained_cross_count(change_id, fixed_id, candidate) +
                            other_constrained_cross_count(fixed_id, change_id, fixed_curve);
                        if (candidate_other_cross > current_other_cross) {
                            ++dbg_other_excess;
                            return;
                        }
                        double score =
                            1000.0 * candidate_other_cross +
                            candidate.maxCurvature(40) +
                            0.02 * candidate.arcLength();
                        if (!have_pair_repair || score < best_pair_score) {
                            if (change_a) {
                                best_a = candidate;
                                best_b = fixed_curve;
                                repair_only_a = true;
                                repair_only_b = false;
                            } else {
                                best_a = fixed_curve;
                                best_b = candidate;
                                repair_only_a = false;
                                repair_only_b = true;
                            }
                            best_pair_score = score;
                            have_pair_repair = true;
                        }
                    };
                    if (cand_b.empty()) {
                        for (const auto& a_candidate : cand_a)
                            consider_single(true, a_candidate);
                    }
                    if (cand_a.empty()) {
                        for (const auto& b_candidate : cand_b)
                            consider_single(false, b_candidate);
                    }
                    for (const auto& a_candidate : cand_a) {
                        for (const auto& b_candidate : cand_b) {
                            if (!preservesSharedExitUTurnOrder(
                                    results[ia->second].id, a_candidate,
                                    results[ib->second].id, b_candidate)) {
                                continue;
                            }
                            // 共享连接点仅允许端点本身重合；不能再用 1.5m
                            // 的宽豁免把已分离的候选错判为相交或贴合。
                            if (curvesHaveForbiddenSameClusterIntersection(
                                    a_candidate, b_candidate, kClusterEndpointTol)) {
                                ++dbg_pair_hit;
                                continue;
                            }
                            ++dbg_viable_pair;
                            int candidate_other_cross =
                                other_constrained_cross_count(
                                    results[ia->second].id, results[ib->second].id, a_candidate) +
                                other_constrained_cross_count(
                                    results[ib->second].id, results[ia->second].id, b_candidate);
                            if (introduces_new_constrained_cross(
                                    results[ia->second].id, results[ib->second].id,
                                    *results[ia->second].curve, a_candidate) ||
                                introduces_new_constrained_cross(
                                    results[ib->second].id, results[ia->second].id,
                                    *results[ib->second].curve, b_candidate)) {
                                ++dbg_other_excess;
                                continue;
                            }
                            if (candidate_other_cross > current_other_cross) {
                                ++dbg_other_excess;
                                continue;
                            }
                            double score =
                                1000.0 * candidate_other_cross +
                                a_candidate.maxCurvature(40) +
                                b_candidate.maxCurvature(40) +
                                0.02 * (a_candidate.arcLength() + b_candidate.arcLength());
                            if (!have_pair_repair || score < best_pair_score) {
                                best_a = a_candidate;
                                best_b = b_candidate;
                                best_pair_score = score;
                                have_pair_repair = true;
                            }
                        }
                    }
                    if (!have_pair_repair) {
                        if (isgDebugPairRepair()) {
                            fprintf(stderr,
                                    "[PAIR-UTURN] %s-%s no pair repair current_other=%d\n",
                                    results[ia->second].id.c_str(),
                                    results[ib->second].id.c_str(),
                                    current_other_cross);
                            fprintf(stderr,
                                    "[PAIR-UTURN] reject stats pair_hit=%d viable_pair=%d other_excess=%d\n",
                                    dbg_pair_hit, dbg_viable_pair, dbg_other_excess);
                        }
                        continue;
                    }
                    std::vector<CurvePatchEntry> patch_entries;
                    if (!repair_only_b)
                        patch_entries.push_back(CurvePatchEntry{results[ia->second].id, best_a});
                    if (!repair_only_a)
                        patch_entries.push_back(CurvePatchEntry{results[ib->second].id, best_b});
                    const bool committed = AtomicCurvePatch().apply(
                        patch_entries, results,
                        [&](ConnectivityCurve& target, const BezierCurve& replacement) {
                            setConnectivityCurveGeometry(target, replacement);
                            validate(target, input, sdf);
                        });
                    if (!committed)
                        continue;
                    all_done[results[ia->second].id] = best_a;
                    all_done[results[ib->second].id] = best_b;
                    changed = true;
                    any_changed = true;
                    ++repaired_pairs;
                    if (repaired_pairs >= repair_budget.mode2_uturn_pairs_per_pass)
                        break;
                }
                if (repaired_pairs >= repair_budget.mode2_uturn_pairs_per_pass)
                    break;
            }
            if (!changed)
                break;
        }
        return any_changed;
    };

    // 声明直行与同簇U型掉头共起点/共切向时, 自然单段cubic可能沿掉头lead走得过久并穿过弧段边缘。
    // 这里只对已发生非豁免掉头交叉的直行做小预算分段偏移, 保持直行形态并降低同簇交叉后才接受。
    auto repairStraightUTurnCrossings = [&](int max_repairs) {
        if (input.crosswalks.empty())
            return false;
        auto result_idx = resultIndexById(results);
        auto neighbors = constrainedNeighborMap(results, cluster_solver_);
        int repaired = 0;
        for (const auto& group : plan.batches) {
            for (const auto& cid : group.conn_ids) {
                if (repaired >= max_repairs || preserved_fixed_ids.count(cid))
                    continue;
                auto ri = result_idx.find(cid);
                if (ri == result_idx.end() || !results[ri->second].curve)
                    continue;
                const Connectivity* conn = scene.view.connectivity(cid);
                if (!conn || conn->turn_type != ConnTurnType::Straight ||
                    isGeometricUTurnConn(*conn, input))
                    continue;

                auto entry = scene.view.entryFrame(conn->entry_lane_id);
                auto exit_ = scene.view.exitFrame(conn->exit_lane_id);
                Vec2d chord = exit_.first - entry.first;
                double chord_len = chord.norm();
                if (chord_len < 1e-8)
                    continue;
                Vec2d T0 = entry.second.norm() > 1e-8
                    ? entry.second.normalized() : chord.normalized();
                Vec2d T1 = exit_.second.norm() > 1e-8
                    ? exit_.second.normalized() : chord.normalized();
                Vec2d chord_dir = chord.normalized();
                double turn_strength = std::abs(cross2d(T0, chord_dir));
                if (turn_strength >= 0.40)
                    continue;

                const BezierCurve& current = *results[ri->second].curve;
                int current_cross = constrainedCrossCountForId(
                    cid, current, results, result_idx, neighbors, kClusterEndpointTol);
                bool current_shape_bad =
                    !isStraightLikeShapeAcceptable(current, chord_len) ||
                    current.maxCurvature(40) >= 0.10;
                auto uturn_cross_count = [&](const BezierCurve& curve) {
                    int count = 0;
                    auto nit = neighbors.find(cid);
                    if (nit == neighbors.end())
                        return count;
                    for (const auto& other_id : nit->second) {
                        auto oi = result_idx.find(other_id);
                        if (oi == result_idx.end() || !results[oi->second].curve)
                            continue;
                        const Connectivity* other_conn = scene.view.connectivity(other_id);
                        if (!other_conn || !isGeometricUTurnConn(*other_conn, input))
                            continue;
                        if (curvesHaveForbiddenSameClusterIntersection(
                                curve, *results[oi->second].curve, kClusterEndpointTol))
                            ++count;
                    }
                    return count;
                };
                bool has_shared_start_uturn_neighbor = false;
                auto nit_for_uturn_neighbor = neighbors.find(cid);
                if (nit_for_uturn_neighbor != neighbors.end()) {
                    for (const auto& other_id : nit_for_uturn_neighbor->second) {
                        auto oi = result_idx.find(other_id);
                        if (oi == result_idx.end() || !results[oi->second].curve)
                            continue;
                        const Connectivity* other_conn = scene.view.connectivity(other_id);
                        if (!other_conn || !isGeometricUTurnConn(*other_conn, input))
                            continue;
                        if ((results[oi->second].curve->startPt() -
                             entry.first).norm() <= 0.30) {
                            has_shared_start_uturn_neighbor = true;
                            break;
                        }
                    }
                }
                int current_uturn_cross = uturn_cross_count(current);
                if (current_shape_bad && current_uturn_cross <= 0 &&
                    !has_shared_start_uturn_neighbor)
                    continue;
                if (current_cross <= 0 && !current_shape_bad)
                    continue;
                if (current_uturn_cross <= 0 && !current_shape_bad)
                    continue;

                Vec2d ref_perp = dominantConstrainedRefPerpForId(
                    cid, current, results, result_idx, neighbors,
                    cluster_solver_, kClusterEndpointTol);
                if (ref_perp.norm() < 1e-8 && current_shape_bad) {
                    auto nit = neighbors.find(cid);
                    if (nit != neighbors.end()) {
                        for (const auto& other_id : nit->second) {
                            const Connectivity* other_conn = scene.view.connectivity(other_id);
                            if (!other_conn ||
                                !isGeometricUTurnConn(*other_conn, input))
                                continue;
                            Vec2d rp = cluster_solver_.refPerpOf(cid, other_id);
                            if (rp.norm() > 1e-8) {
                                ref_perp = rp.normalized();
                                break;
                            }
                        }
                    }
                }
                if (ref_perp.norm() < 1e-8)
                    continue;
                ref_perp.normalize();

                BezierCurve best;
                bool have_best = false;
                int best_cross = current_cross;
                int best_uturn_cross = current_uturn_cross;
                double best_score = std::numeric_limits<double>::infinity();
                auto consider = [&](const BezierCurve& candidate, double offset) {
                    if (candidate.empty() ||
                        curveSelfIntersectsBusiness(candidate, 1.0) ||
                        !isStraightLikeShapeAcceptable(candidate, chord_len) ||
                        candidate.maxCurvature(40) >= 0.10)
                        return;
                    int cross_count = constrainedCrossCountForId(
                        cid, candidate, results, result_idx, neighbors, kClusterEndpointTol);
                    if (cross_count > current_cross)
                        return;
                    int cand_uturn_cross = uturn_cross_count(candidate);
                    if (cand_uturn_cross > current_uturn_cross)
                        return;
                    if (!current_shape_bad &&
                        cross_count >= current_cross &&
                        cand_uturn_cross >= current_uturn_cross)
                        return;
                    if (current_shape_bad &&
                        cross_count > current_cross)
                        return;
                    CurveRisk risk = assessCurveRisk(candidate, input, sdf, {}, true);
                    if (risk.physical())
                        return;
                    double score = 1000.0 * cross_count +
                        500.0 * cand_uturn_cross + offset +
                        50.0 * candidate.maxCurvature(40) +
                        kChordBudgetOvershootPenalty *
                            curveChordBudgetOvershoot(candidate) +
                        0.02 * candidate.arcLength();
                    if (!have_best || cross_count < best_cross ||
                        (cross_count == best_cross &&
                         cand_uturn_cross < best_uturn_cross) ||
                        (cross_count == best_cross &&
                         cand_uturn_cross == best_uturn_cross &&
                         score < best_score)) {
                        best = candidate;
                        best_cross = cross_count;
                        best_uturn_cross = cand_uturn_cross;
                        best_score = score;
                        have_best = true;
                    }
                };

                for (double side : {-1.0, 1.0}) {
                    for (double offset : {0.45, 0.70, 1.00, 1.35, 1.80, 2.40,
                                          3.00, 4.25}) {
                        Vec2d off = ref_perp * (side * offset);
                        for (double lead : {18.0, 20.0, 24.0, 28.0}) {
                            if (lead > chord_len * 0.45)
                                continue;
                            for (double alpha : {0.34, 0.42, 0.55}) {
                                Vec2d q0 = entry.first + T0 * lead + off;
                                Vec2d q1 = exit_.first - T1 * lead + off;
                                consider(makeCurveFromKnots(
                                             {entry.first, q0, q1, exit_.first},
                                             {T0, chord_dir, chord_dir, T1},
                                             alpha),
                                         offset);

                                Vec2d r0 = entry.first + chord * 0.18 + off * 0.65;
                                Vec2d rm = entry.first + chord * 0.50 + off;
                                Vec2d r1 = entry.first + chord * 0.82 + off * 0.65;
                                consider(makeCurveFromKnots(
                                             {entry.first, r0, rm, r1, exit_.first},
                                             {T0, chord_dir, chord_dir, chord_dir, T1},
                                             alpha),
                                         offset);
                            }
                        }
                    }
                }

                if (!have_best)
                    continue;
                setConnectivityCurveGeometry(results[ri->second], best);
                validate(results[ri->second], input, sdf);
                ++repaired;
            }
        }
        return repaired > 0;
    };

    // 同入/同出同向普通转向在共用G1首尾段后仍可能在分岔/汇合区交叉。
    // 这里只对已确认非豁免交叉的非U型转向对做小预算两段候选搜索。
    auto repairSharedEndpointTurnPairs = [&](int max_repairs) {
        if (input.mode == 2)
            return false;
        auto result_idx = resultIndexById(results);
        auto neighbors = constrainedNeighborMap(results, cluster_solver_);
        int repaired = 0;
        auto try_repair_one = [&](const ConnId& cid, const ConnId& other_id) {
            auto ri = result_idx.find(cid);
            auto oi = result_idx.find(other_id);
            if (ri == result_idx.end() || oi == result_idx.end() ||
                !results[ri->second].curve || !results[oi->second].curve ||
                preserved_fixed_ids.count(cid))
                return false;
            const Connectivity* conn = scene.view.connectivity(cid);
            if (!conn || conn->turn_type == ConnTurnType::Straight ||
                isGeometricUTurnConn(*conn, input))
                return false;
            auto entry = scene.view.entryFrame(conn->entry_lane_id);
            auto exit_ = scene.view.exitFrame(conn->exit_lane_id);
            Vec2d chord = exit_.first - entry.first;
            double chord_len = chord.norm();
            if (chord_len < 1e-8)
                return false;
            Vec2d T0 = entry.second.norm() > 1e-8
                ? entry.second.normalized() : chord.normalized();
            Vec2d T1 = exit_.second.norm() > 1e-8
                ? exit_.second.normalized() : chord.normalized();
            Vec2d chord_dir = chord.normalized();
            double turn_strength = std::abs(cross2d(T0, chord_dir));
            if (turn_strength <= 0.35)
                return false;

            const BezierCurve& current = *results[ri->second].curve;
            bool current_pair_strict_cross = curvesIntersectBusiness(
                current, *results[oi->second].curve, 0.15);
            int current_shared_cross = constrainedSharedEndpointCrossCountForId(
                cid, current, results, result_idx, neighbors, cluster_solver_, kClusterEndpointTol);
            int current_cross = constrainedCrossCountForId(
                cid, current, results, result_idx, neighbors, kClusterEndpointTol);
            int current_fixed_cross = constrainedFixedCrossCountForId(
                cid, current, results, result_idx, neighbors, kClusterEndpointTol);
            if (current_cross <= 0 && current_shared_cross <= 0 &&
                !current_pair_strict_cross)
                return false;

            BezierCurve best;
            bool have_best = false;
            int best_cross = current_cross;
            int best_fixed_cross = current_fixed_cross;
            int best_segments = std::numeric_limits<int>::max();
            double best_score = std::numeric_limits<double>::infinity();
            auto consider = [&](const BezierCurve& candidate, double offset) {
                if (candidate.empty() ||
                    curveSelfIntersectsBusiness(candidate, 1.0) ||
                    !isNonUTurnTurnShapeAcceptable(
                        candidate, chord_len, turn_strength))
                    return;
                // 本阶段只处理普通（非直行、非掉头）转向，而 shape.ordinary.single_segment
                // 是无豁免的硬性形态要求：普通曲线必须是单段 cubic。同簇相交是可计数、
                // 可豁免、且会被后续阶段重新洗牌的软指标，形态是不可交易的硬约束，
                // 所以不允许"拆成两段换少一两个交叉"。这里直接把多段候选挡在门外，
                // 而不是只在交叉数并列时才让段数参与比较——并列判据永远轮不到，
                // 因为多段候选正是靠"交叉数更低"进门的。
                // ISG_LEGACY_SEGTIE=1 恢复旧行为（含下方两段候选网格）以便 A/B。
                if (!isgLegacySegTie() && candidate.numSegments() != 1)
                    return;
                if (assessCurveRisk(candidate, input, sdf, {}, true).physical())
                    return;
                if (curvesHaveForbiddenSameClusterIntersection(
                        candidate, *results[oi->second].curve, kClusterEndpointTol))
                    return;
                bool pair_strict_cross = curvesIntersectBusiness(
                    candidate, *results[oi->second].curve, 0.15);
                if (pair_strict_cross)
                    return;
                int cross_count = constrainedCrossCountForId(
                    cid, candidate, results, result_idx, neighbors, kClusterEndpointTol);
                int shared_cross = constrainedSharedEndpointCrossCountForId(
                    cid, candidate, results, result_idx, neighbors, cluster_solver_, kClusterEndpointTol);
                int fixed_cross = constrainedFixedCrossCountForId(
                    cid, candidate, results, result_idx, neighbors, kClusterEndpointTol);
                if (fixed_cross > current_fixed_cross ||
                    (cross_count >= current_cross &&
                     shared_cross >= current_shared_cross))
                    return;
                double score = 1000.0 * std::max(cross_count, shared_cross) +
                    0.25 * offset +
                    2.0 * std::abs(candidate.arcLength() / chord_len - 1.10) +
                    candidate.maxCurvature(40) +
                    kChordBudgetOvershootPenalty *
                        curveChordBudgetOvershoot(candidate) +
                    0.02 * candidate.arcLength();
                // 段数在交叉数、固定形态交叉数之后、综合分之前参与比较。新行为下
                // 候选全是单段，这个键恒相等；它只在 ISG_LEGACY_SEGTIE=1 的旧行为
                // 里生效——保留它是为了让 A/B 的两侧都可复现，而不是为了兜底：
                // 旧行为的教训正是"并列才比段数"不足以阻止多段候选进门。
                const int segments = candidate.numSegments();
                const int cross_metric = std::max(cross_count, shared_cross);
                const int seg_metric = isgLegacySegTie() ? 0 : segments;
                const int best_seg_metric =
                    isgLegacySegTie() ? 0 : best_segments;
                if (!have_best ||
                    cross_metric < best_cross ||
                    (cross_metric == best_cross &&
                     fixed_cross < best_fixed_cross) ||
                    (cross_metric == best_cross &&
                     fixed_cross == best_fixed_cross &&
                     seg_metric < best_seg_metric) ||
                    (cross_metric == best_cross &&
                     fixed_cross == best_fixed_cross &&
                     seg_metric == best_seg_metric &&
                     score < best_score)) {
                    best = candidate;
                    best_cross = cross_metric;
                    best_fixed_cross = fixed_cross;
                    best_segments = segments;
                    best_score = score;
                    have_best = true;
                }
            };

            for (double a0 : {0.001, 0.005, 0.01, 0.02, 0.04, 0.08, 0.12,
                              0.16, 0.22, 0.30, 0.42, 0.55}) {
                for (double a1 : {0.18, 0.25, 0.34, 0.42, 0.55, 0.70,
                                  0.90, 1.10, 1.30}) {
                    BezierSegment seg;
                    seg.ctrl[0] = entry.first;
                    seg.ctrl[1] = entry.first + T0 * (chord_len * a0);
                    seg.ctrl[2] = exit_.first - T1 * (chord_len * a1);
                    seg.ctrl[3] = exit_.first;
                    BezierCurve candidate;
                    candidate.segs.push_back(seg);
                    consider(candidate, std::abs(a1 - a0));
                }
            }

            // 两段候选网格：新行为下 consider() 会整体拒收，直接跳过 10368 个候选，
            // 既省掉本阶段绝大部分 assessCurveRisk / 交叉计数开销，也保证不会产生
            // 多段普通曲线。只在 ISG_LEGACY_SEGTIE=1 的旧行为下仍然枚举。
            if (isgLegacySegTie()) {
            Vec2d chord_perp{-chord_dir.y(), chord_dir.x()};
            for (double side : {-1.0, 1.0}) {
                for (double offset : {1.0, 2.0, 3.5, 5.0, 7.0, 10.0}) {
                    for (double frac : {0.15, 0.22, 0.30, 0.40, 0.55, 0.65}) {
                        Vec2d mid = entry.first + chord * frac +
                            chord_perp * (side * offset);
                        for (double deg : {-120.0, -90.0, -60.0, -30.0,
                                           0.0, 30.0, 60.0, 90.0, 120.0}) {
                            Vec2d mid_tan = rotateVec(
                                chord_dir, deg * M_PI / 180.0);
                            if (mid_tan.norm() < 1e-8)
                                continue;
                            mid_tan.normalize();
                            for (double a0 : {0.18, 0.24, 0.34, 0.44}) {
                                for (double a1 : {0.18, 0.24, 0.34, 0.44}) {
                                    BezierCurve candidate;
                                    candidate.segs.push_back(makeCubicG1(
                                        entry.first, T0, mid, mid_tan, a0));
                                    candidate.segs.push_back(makeCubicG1(
                                        mid, mid_tan, exit_.first, T1, a1));
                                    consider(candidate, offset);
                                }
                            }
                        }
                    }
                }
            }
            }
            if (!have_best)
                return false;
            setConnectivityCurveGeometry(results[ri->second], best);
            validate(results[ri->second], input, sdf);
            return true;
        };

        for (const auto& pair : cluster_solver_.pairs()) {
            if (repaired >= max_repairs ||
                pair.exempt == CrossExemption::StructuralCross)
                continue;
            auto ia = result_idx.find(pair.id_a);
            auto ib = result_idx.find(pair.id_b);
            if (ia == result_idx.end() || ib == result_idx.end() ||
                !results[ia->second].curve || !results[ib->second].curve)
                continue;
            bool pair_forbidden =
                curvesIntersectBusiness(
                    *results[ia->second].curve,
                    *results[ib->second].curve, 0.15);
            if (!pair.shared_endpoint || !pair_forbidden)
                continue;
            const Connectivity* ca =
                scene.view.connectivity(pair.id_a);
            const Connectivity* cb =
                scene.view.connectivity(pair.id_b);
            if (!ca || !cb || ca->turn_type == ConnTurnType::Straight ||
                cb->turn_type == ConnTurnType::Straight ||
                isGeometricUTurnConn(*ca, input) || isGeometricUTurnConn(*cb, input))
                continue;
            bool shared_start =
                (results[ia->second].curve->startPt() -
                 results[ib->second].curve->startPt()).norm() <= 0.30;
            bool shared_end =
                (results[ia->second].curve->endPt() -
                 results[ib->second].curve->endPt()).norm() <= 0.30;
            if (!shared_start && !shared_end)
                continue;
            bool changed = try_repair_one(pair.id_a, pair.id_b);
            if (!changed)
                changed = try_repair_one(pair.id_b, pair.id_a);
            if (changed) {
                ++repaired;
                result_idx = resultIndexById(results);
                neighbors = constrainedNeighborMap(results, cluster_solver_);
            }
        }
        return repaired > 0;
    };

    // Boundary-only 路口仍需要纯拓扑阶段的最终同簇收口；但其中的普通转向大
    // 网格与 Boundary 物理候选相互独立，后面按 complex_boundary_fast_path 跳过。
    bool topology_repair = input.obstacles.empty();
    if (topology_repair) {
        auto topology_t0 = std::chrono::steady_clock::now();
        auto result_idx = resultIndexById(results);
        auto neighbors = constrainedNeighborMap(results, cluster_solver_);
        // 固有掉头可能同时影响多条同簇直行/左转；一轮修复后其它曲线的
        // 位置变化会改变交叉关系，必须迭代复查到稳定。
        // 密集 Boundary 路口的掉头联合收口紧接着在下方执行；这里若再对每条
        // U-turn 重复三轮完整候选搜索，只会重算同一组物理约束，且不提供新的
        // 约束覆盖。普通拓扑场景仍保留原有多轮稳定化。
        const int pure_topology_passes = complex_boundary_fast_path
            ? 0 : repair_budget.pure_topology_passes;
        for (int pass = 0; pass < pure_topology_passes; ++pass) {
            auto pass_t0 = std::chrono::steady_clock::now();
            auto repair_ids = collectPureGeometryRepairIds(
                results, cluster_solver_, preserved_fixed_ids, kClusterEndpointTol);
            if (repair_ids.empty())
                break;
            bool changed = false;
            int repaired_this_pass = 0;
            for (const auto& group : plan.batches) {
                for (const auto& cid : group.conn_ids) {
                if (!repair_ids.count(cid))
                    continue;
                if (repaired_this_pass >= repair_budget.pure_topology_repairs_per_pass)
                    break;
                auto ri = result_idx.find(cid);
                if (ri == result_idx.end())
                    continue;
                if (constrainedSharedEndpointCrossCountForId(
                        cid, *results[ri->second].curve, results, result_idx,
                        neighbors, cluster_solver_, kClusterEndpointTol) <= 0)
                    continue;
                const Connectivity* conn = scene.view.connectivity(cid);
                if (!conn)
                    continue;
                if (isGeometricUTurnConn(*conn, input)) {
                    auto full_done = curveMapFromResults(results);
                    auto sibs = buildSiblings(
                        cid, full_done, cluster_solver_, input.connectivities,
                        pure_geometry, &preserved_fixed_ids);
                    auto sampled = sampleSiblingsForIntersections(sibs);
                    BezierCurve repaired = *results[ri->second].curve;
                    auto entry = scene.view.entryFrame(conn->entry_lane_id);
                    auto exit_ = scene.view.exitFrame(conn->exit_lane_id);
                    const UTurnFamilyInfo& family_info =
                        pairedFamilySnapshot(conn->id);
                    const double alignment_station = family_info.aligned_station;
                    double family_stagger_step = 0.0;
                    double shared_entry_stagger = 0.0;
                    double shared_exit_stagger = 0.0;
                    double shared_endpoint_stagger = uturnSharedEndpointStagger(
                        *conn, input, cluster_solver_,
                        direction_cfg_.uturn_alignment_scope,
                        &family_stagger_step,
                        &shared_entry_stagger, &shared_exit_stagger);
                    if (!searchSegmentedUTurnTwoPass(
                            entry.first, entry.second, exit_.first, exit_.second,
                            input, sdf, sampled, true, repaired,
                            family_info.lead0, family_info.lead1, nullptr,
                            shared_endpoint_stagger, 0.0, 0.0,
                            alignment_station, shared_entry_stagger,
                            shared_exit_stagger, family_stagger_step))
                        continue;
                    setConnectivityCurveGeometry(results[ri->second], repaired);
                    validate(results[ri->second], input, sdf);
                    changed = true;
                    ++repaired_this_pass;
                    continue;
                }
                // 非U型曲线的残余共享端点冲突由后续小预算定向修复处理；
                // 这里不进入高成本通用拓扑搜索，避免纯几何场景被少数
                // 冲突放大到分钟级。左/右转与掉头的允许相交已在
                // isAllowedSameClusterCrossing 中统一豁免。
                continue;
                }
                if (repaired_this_pass >= repair_budget.pure_topology_repairs_per_pass)
                    break;
            }
            if (isgProfile()) {
                fprintf(stderr, "[ISG_PROFILE] pure uturn repair pass %d stage %.3f ms\n",
                        pass,
                        std::chrono::duration<double, std::milli>(
                            std::chrono::steady_clock::now() - pass_t0).count());
            }
            if (!changed)
                break;
        }

        auto shape_restore_t0 = std::chrono::steady_clock::now();
        result_idx = resultIndexById(results);
        auto restore_neighbors = constrainedNeighborMap(results, cluster_solver_);
        int restored_count = 0;
        for (const auto& group : plan.batches) {
            for (auto& cid : group.conn_ids) {
                if (restored_count >= repair_budget.shape_restores ||
                    preserved_fixed_ids.count(cid))
                    continue;
                const Connectivity* conn = scene.view.connectivity(cid);
                if (!conn || isGeometricUTurnConn(*conn, input))
                    continue;
                auto ri = result_idx.find(cid);
                if (ri == result_idx.end() || !results[ri->second].curve)
                    continue;
                if (results[ri->second].curve->numSegments() <= 1)
                    continue;
                auto entry = scene.view.entryFrame(conn->entry_lane_id);
                auto exit_ = scene.view.exitFrame(conn->exit_lane_id);
                Vec2d chord = exit_.first - entry.first;
                double chord_len = chord.norm();
                if (chord_len < 1e-8)
                    continue;
                Vec2d T0 = entry.second.norm() > 1e-8 ? entry.second.normalized() : chord.normalized();
                double turn_strength = std::abs(cross2d(T0, chord.normalized()));
                bool straight_like = turn_strength < 0.25;
                const BezierCurve& current = *results[ri->second].curve;
                double arc_chord_current = current.arcLength() / chord_len;
                double pure_turn_min_arc = ordinaryTurnArcChordRestoreFloor(
                    curveEndpointTurnAngle(current), chord_len);
                bool shape_bad = straight_like
                    ? !isStraightLikeShapeAcceptable(current, chord_len)
                    : (isNonUTurnTurnShapeRisk(current, chord_len, turn_strength) ||
                       arc_chord_current < pure_turn_min_arc ||
                       curveHasCurvatureSignFlip(current));
                if (!shape_bad)
                    continue;
                BezierCurve restored;
                if (!tryPureGeometryShapeRestore(
                        *conn, input, results, result_idx,
                        restore_neighbors, restored))
                    continue;
                bool restored_shape_ok = straight_like
                    ? isStraightLikeShapeAcceptable(restored, chord_len)
                    : (isNonUTurnTurnShapeAcceptable(restored, chord_len, turn_strength) &&
                       restored.arcLength() / chord_len >= pure_turn_min_arc &&
                       !curveHasCurvatureSignFlip(restored));
                if (!restored_shape_ok)
                    continue;
                setConnectivityCurveGeometry(results[ri->second], restored);
                validate(results[ri->second], input, sdf);
                ++restored_count;
            }
        }
        if (isgProfile()) {
            fprintf(stderr, "[ISG_PROFILE] pure shape restore stage %.3f ms\n",
                    std::chrono::duration<double, std::milli>(
                        std::chrono::steady_clock::now() - shape_restore_t0).count());
        }

        auto natural_restore_t0 = std::chrono::steady_clock::now();
        result_idx = resultIndexById(results);
        auto natural_restore_done = curveMapFromResults(results);
        for (const auto& group : plan.batches) {
            for (auto& cid : group.conn_ids) {
                if (preserved_fixed_ids.count(cid))
                    continue;
                const Connectivity* conn = scene.view.connectivity(cid);
                if (!conn || isGeometricUTurnConn(*conn, input))
                    continue;
                auto ri = result_idx.find(cid);
                if (ri == result_idx.end() || !results[ri->second].curve)
                    continue;
                auto entry = scene.view.entryFrame(conn->entry_lane_id);
                auto exit_ = scene.view.exitFrame(conn->exit_lane_id);
                auto sibs = buildSiblings(
                    cid, natural_restore_done, cluster_solver_, input.connectivities,
                    pure_geometry, &preserved_fixed_ids);
                auto sampled = sampleSiblingsForIntersections(sibs);
                BezierCurve restored = *results[ri->second].curve;
                bool restored_ok = tryNaturalSingleCubicIfSafe(
                    entry.first, entry.second, exit_.first, exit_.second,
                    input, sdf, sampled, true, restored);
                if (!restored_ok) {
                    restored_ok = tryNaturalStraightSingleCubicIfSafe(
                        entry.first, entry.second, exit_.first, exit_.second,
                        input, sdf, sampled, true, restored);
                }
                if (!restored_ok)
                    continue;
                setConnectivityCurveGeometry(results[ri->second], restored);
                validate(results[ri->second], input, sdf);
                natural_restore_done[cid] = restored;
            }
        }
        if (isgProfile()) {
            fprintf(stderr, "[ISG_PROFILE] pure natural restore stage %.3f ms\n",
                    std::chrono::duration<double, std::milli>(
                        std::chrono::steady_clock::now() - natural_restore_t0).count());
        }

        auto base_expression_t0 = std::chrono::steady_clock::now();
        // Boundary 密集路口的普通转向候选不参与 U-turn/Boundary 联合收口，
        // 这段大网格会重复评估同一物理风险并占用绝大多数耗时；保留真正纯几何
        // 输入的表达式修复，Boundary 路口交给后面的定向物理和家族原子路径。
        if (!complex_boundary_fast_path) {
        int base_expression_passes = repair_budget.base_expression_passes;
        for (int base_pass = 0; base_pass < base_expression_passes; ++base_pass) {
            result_idx = resultIndexById(results);
            auto base_neighbors = constrainedNeighborMap(results, cluster_solver_);
            for (const auto& group : plan.batches) {
                for (auto& cid : group.conn_ids) {
                if (preserved_fixed_ids.count(cid))
                    continue;
                const Connectivity* conn = scene.view.connectivity(cid);
                if (!conn || isGeometricUTurnConn(*conn, input))
                    continue;
                auto ri = result_idx.find(cid);
                if (ri == result_idx.end() || !results[ri->second].curve)
                    continue;

                auto entry = scene.view.entryFrame(conn->entry_lane_id);
                auto exit_ = scene.view.exitFrame(conn->exit_lane_id);
                Vec2d chord = exit_.first - entry.first;
                double chord_len = chord.norm();
                if (chord_len < 1e-8)
                    continue;
                Vec2d T0 = entry.second.norm() > 1e-8
                    ? entry.second.normalized() : chord.normalized();
                Vec2d T1 = exit_.second.norm() > 1e-8
                    ? exit_.second.normalized() : chord.normalized();
                double turn_strength = std::abs(cross2d(T0, chord.normalized()));
                bool straight_like = turn_strength < 0.25;
                const BezierCurve& current = *results[ri->second].curve;
                int current_fixed_cross = constrainedFixedCrossCountForId(
                    cid, current, results, result_idx, base_neighbors, kClusterEndpointTol);

                BezierCurve best;
                bool have_best = false;
                double best_score = std::numeric_limits<double>::infinity();
                std::vector<BezierCurve> single_candidates;
                if (current.numSegments() > 1) {
                    BezierSegment collapsed;
                    collapsed.ctrl[0] = entry.first;
                    collapsed.ctrl[1] = current.segs.front().ctrl[1];
                    collapsed.ctrl[2] = current.segs.back().ctrl[2];
                    collapsed.ctrl[3] = exit_.first;
                    BezierCurve candidate;
                    candidate.segs.push_back(collapsed);
                    single_candidates.push_back(candidate);
                }
                const std::vector<double> base_alphas = straight_like
                    ? std::vector<double>{0.12, 0.18, 0.24, 0.30, 0.36, 0.40}
                    : std::vector<double>{0.30, 0.36, 0.40, 0.46, 0.52, 0.60};
                for (double alpha : base_alphas) {
                    BezierCurve candidate;
                    candidate.segs.push_back(makeCubicG1(
                        entry.first, T0, exit_.first, T1, alpha));
                    single_candidates.push_back(candidate);
                }
                const std::vector<double> base_a0 = straight_like
                    ? std::vector<double>{0.12, 0.16, 0.18, 0.20, 0.22, 0.30}
                    : std::vector<double>{0.08, 0.12, 0.16, 0.22, 0.30, 0.42, 0.55};
                const std::vector<double> base_a1 = base_a0;
                for (double a0 : base_a0) {
                    for (double a1 : base_a1) {
                        BezierSegment seg;
                        seg.ctrl[0] = entry.first;
                        seg.ctrl[1] = entry.first + T0 * (chord_len * a0);
                        seg.ctrl[2] = exit_.first - T1 * (chord_len * a1);
                        seg.ctrl[3] = exit_.first;
                        BezierCurve candidate;
                        candidate.segs.push_back(seg);
                        single_candidates.push_back(candidate);
                    }
                }
                for (const auto& candidate : single_candidates) {
                    if (candidate.empty() ||
                        curveSelfIntersectsBusiness(candidate, 1.0))
                        continue;
                    if (candidate.numSegments() == 1 &&
                        !ordinarySingleCubicControlsValid(
                            candidate, entry.first, T0, exit_.first, T1,
                            1e-5, true))
                        continue;
                    bool shape_ok = straight_like
                        ? isStrictStraightBaseShapeAcceptable(candidate, chord_len)
                        : isNonUTurnTurnShapeAcceptable(
                              candidate, chord_len, turn_strength);
                    if (!straight_like && shape_ok) {
                        double pure_min_arc = ordinaryTurnArcChordRestoreFloor(
                            curveEndpointTurnAngle(candidate), chord_len);
                        double arc_chord = candidate.arcLength() / chord_len;
                        shape_ok = arc_chord >= pure_min_arc &&
                                   (chord_len >= 12.0 || arc_chord <= 1.08);
                    }
                    if (!shape_ok)
                        continue;
                    BezierCurve accepted_candidate = candidate;
                    CurveRisk base_candidate_risk =
                        assessCurveRisk(accepted_candidate, input, sdf, {}, true);
                    if (base_candidate_risk.physical())
                        continue;
                    int cand_cross = constrainedCrossCountForId(
                        cid, accepted_candidate, results, result_idx, base_neighbors,
                        kClusterEndpointTol);
                    int cand_fixed_cross = constrainedFixedCrossCountForId(
                        cid, accepted_candidate, results, result_idx, base_neighbors,
                        kClusterEndpointTol);
                    if (cand_fixed_cross > current_fixed_cross)
                        continue;
                    double target = straight_like ? 1.01 : 1.10;
                    // 形态评分必须基于完整候选；端点贴边段不再允许被删除后
                    // 伪装成“中段安全”。候选本身已通过严格物理门禁。
                    BezierCurve shape_curve = accepted_candidate;
                    double arc_chord = shape_curve.arcLength() /
                        std::max(1e-6, (shape_curve.endPt() - shape_curve.startPt()).norm());
                    double score = 1000.0 * cand_cross +
                        std::abs(arc_chord - target) +
                        shape_curve.maxCurvature(40) +
                        0.02 * accepted_candidate.arcLength();
                    if (!have_best || score < best_score) {
                        best = accepted_candidate;
                        best_score = score;
                        have_best = true;
                    }
                }
                if (!have_best)
                    continue;
                setConnectivityCurveGeometry(results[ri->second], best);
                validate(results[ri->second], input, sdf);
                }
            }
        }

        result_idx = resultIndexById(results);
        auto pair_safe_neighbors = constrainedNeighborMap(results, cluster_solver_);
        auto make_pair_safe_candidates = [&](const Connectivity& conn,
                                             bool broad_turn_candidates,
                                             bool half_step_offset) {
            std::vector<BezierCurve> candidates;
            auto entry = scene.view.entryFrame(conn.entry_lane_id);
            auto exit_ = scene.view.exitFrame(conn.exit_lane_id);
            Vec2d chord = exit_.first - entry.first;
            double chord_len = chord.norm();
            if (chord_len < 1e-8)
                return candidates;
            Vec2d T0 = entry.second.norm() > 1e-8
                ? entry.second.normalized() : chord.normalized();
            Vec2d T1 = exit_.second.norm() > 1e-8
                ? exit_.second.normalized() : chord.normalized();
            if (T0.dot(T1) < -0.5)
                return candidates;
            auto add = [&](const BezierCurve& c) {
                if (c.empty())
                    return;
                candidates.push_back(c);
            };
            // 半步偏移：把栅格换成相邻两格的中点序列。格数减一、步长不变，
            // 相位整体挪半格，因此不是"加密"而是"错相位重扫"，代价有界。
            auto half_step = [](const std::vector<double>& g) {
                std::vector<double> out;
                for (size_t i = 1; i < g.size(); ++i)
                    out.push_back(0.5 * (g[i - 1] + g[i]));
                return out;
            };
            const std::vector<double> base_alphas = broad_turn_candidates
                ? std::vector<double>{0.18, 0.24, 0.30, 0.34, 0.42, 0.50}
                : std::vector<double>{0.34, 0.42, 0.50};
            const std::vector<double> pair_alphas = half_step_offset
                ? half_step(base_alphas) : base_alphas;
            for (double alpha : pair_alphas) {
                BezierCurve c;
                c.segs.push_back(makeCubicG1(
                    entry.first, T0, exit_.first, T1, alpha));
                add(c);
            }
            // a0 必须一直枚举到接近整条弦。共享进入端点的"换道弧"
            // (声明直行、横偏却达 0.28~0.29 弦长) 让开整族同侧转向的
            // 唯一形态就是"先沿进入切向长距离直行、末段再转过去"，
            // 对应 a0 >= 0.9。此前 a0 上限 0.55 把这一段候选整体排除，
            // 110003285 的 14 因此 171 个候选全数被形态门禁否掉
            // (arc/chord 达不到 1.02 或控制点越过切向有效区间)，
            // 14|18、14|24 只能留在原地。
            const std::vector<double> base_a0 = broad_turn_candidates
                ? std::vector<double>{0.005, 0.02, 0.04, 0.055, 0.08, 0.12, 0.16,
                                      0.22, 0.30, 0.42, 0.55, 0.70, 0.85, 0.94, 1.00}
                : std::vector<double>{0.06, 0.10, 0.16, 0.22, 0.30};
            const std::vector<double> base_a1 = broad_turn_candidates
                ? std::vector<double>{0.005, 0.02, 0.04, 0.08, 0.12, 0.16, 0.18, 0.25, 0.34, 0.42, 0.55, 0.70, 0.90, 1.10, 1.30}
                : std::vector<double>{0.55, 0.70, 0.90, 1.10, 1.30};
            const std::vector<double> pair_a0 = half_step_offset
                ? half_step(base_a0) : base_a0;
            const std::vector<double> pair_a1 = half_step_offset
                ? half_step(base_a1) : base_a1;
            for (double a0 : pair_a0) {
                for (double a1 : pair_a1) {
                    BezierSegment seg;
                    seg.ctrl[0] = entry.first;
                    seg.ctrl[1] = entry.first + T0 * (chord_len * a0);
                    seg.ctrl[2] = exit_.first - T1 * (chord_len * a1);
                    seg.ctrl[3] = exit_.first;
                    BezierCurve c;
                    c.segs.push_back(seg);
                    add(c);
                }
            }
            return candidates;
        };
        auto signed_turn_strength_of = [&](const ConnId& cid) {
            const Connectivity* conn = scene.view.connectivity(cid);
            if (!conn)
                return 0.0;
            auto entry = scene.view.entryFrame(conn->entry_lane_id);
            auto exit_ = scene.view.exitFrame(conn->exit_lane_id);
            Vec2d chord = exit_.first - entry.first;
            if (chord.norm() < 1e-8)
                return 0.0;
            Vec2d T0 = entry.second.norm() > 1e-8
                ? entry.second.normalized() : chord.normalized();
            return cross2d(T0, chord.normalized());
        };
        std::vector<CurvePair> pair_safe_pairs = cluster_solver_.pairs();
        std::unordered_set<std::string> pair_safe_seen;
        auto pair_key = [](const ConnId& a, const ConnId& b) {
            return a < b ? a + "\n" + b : b + "\n" + a;
        };
        for (const auto& p : pair_safe_pairs)
            pair_safe_seen.insert(pair_key(p.id_a, p.id_b));
        for (size_t ia = 0; ia < results.size(); ++ia) {
            if (!results[ia].curve)
                continue;
            for (size_t ib = ia + 1; ib < results.size(); ++ib) {
                if (!results[ib].curve)
                    continue;
                const ConnId& id_a = results[ia].id;
                const ConnId& id_b = results[ib].id;
                double signed_turn_a = signed_turn_strength_of(id_a);
                double signed_turn_b = signed_turn_strength_of(id_b);
                if (std::abs(signed_turn_a) < 0.25 ||
                    std::abs(signed_turn_b) < 0.25 ||
                    signed_turn_a * signed_turn_b <= 0.0)
                    continue;
                const BezierCurve& ca = *results[ia].curve;
                const BezierCurve& cb = *results[ib].curve;
                bool shared_endpoint =
                    (ca.startPt() - cb.startPt()).norm() <= 0.20 ||
                    (ca.endPt() - cb.endPt()).norm() <= 0.20;
                if (!shared_endpoint)
                    continue;
                std::string key = pair_key(id_a, id_b);
                if (!pair_safe_seen.insert(key).second)
                    continue;
                CurvePair p;
                p.id_a = id_a;
                p.id_b = id_b;
                p.shared_endpoint = true;
                pair_safe_pairs.push_back(p);
            }
        }
        auto strict_shared_cross_count = [&](
            const ConnId& id, const BezierCurve& curve,
            const std::unordered_map<ConnId, size_t>& idx) {
            int count = 0;
            for (const auto& p : pair_safe_pairs) {
                if (!p.shared_endpoint)
                    continue;
                ConnId other_id;
                if (p.id_a == id)
                    other_id = p.id_b;
                else if (p.id_b == id)
                    other_id = p.id_a;
                else
                    continue;
                auto oi = idx.find(other_id);
                if (oi == idx.end() || !results[oi->second].curve)
                    continue;
                double signed_a = signed_turn_strength_of(id);
                double signed_b = signed_turn_strength_of(other_id);
                double turn_a = std::abs(signed_a);
                double turn_b = std::abs(signed_b);
                bool participates =
                    ((turn_a < 0.25 && turn_b >= 0.25) ||
                     (turn_b < 0.25 && turn_a >= 0.25) ||
                     (turn_a >= 0.25 && turn_b >= 0.25) ||
                     sharedEndpointStraightPairParticipates(
                         turn_a, turn_b, p.exempt));
                if (!participates)
                    continue;
                bool same_side_turn_pair =
                    turn_a >= 0.25 && turn_b >= 0.25 &&
                    signed_a * signed_b > 0.0;
                (void)same_side_turn_pair;
                // pair_safe_pairs 全部由共享连接点构造，统一走审计口径，
                // 避免近端点数值穿插被记成真实交叉后污染 cand_cross 排序。
                const double pair_tol = kSharedEndpointPairCrossTol;
                if (curvesIntersectBeyondSharedEndpointOverlap(
                        curve, *results[oi->second].curve, pair_tol))
                    ++count;
            }
            return count;
        };
        std::unordered_set<std::string> failed_pair_repair_attempts;
        // 候选物理风险（Boundary / 障碍 / Fence）只由候选曲线自身与固定的
        // input/sdf 决定，与当前配对无关。"修不动"记忆改为连同两条曲线的实际
        // 控制点一起记之后，同一份候选网格会在多次重试里被反复求值，而
        // assessCurveRisk 是本阶段单次最贵的判定。这里按控制点串缓存
        // physical() 结果，纯查表、不改变任何录取判断。
        std::unordered_map<std::string, bool> pair_candidate_physical_memo;
        // 逐邻居统计候选相对当前曲线"新增"的两类问题，两者都只用于评分排序：
        // - sampled_out: 新增的采样相交违约数。只比较交叉总数无法区分"真正
        //   消解"和"把违约搬到另一个邻居身上"，后者会让两个对轮流修复彼此
        //   （110003449 的 23 在 23|52 与 23|24 之间翻转 11 次），因此把它计入
        //   评分，让不搬运的候选优先。
        // - precursor_out: 新增的控制多边形折线互穿数，是中段相交的几何前兆。
        // 两者都不做硬拒绝：实测把它们当硬约束会让候选集枯竭，反而使本该被
        // 修复的真实相交无解（全量数据集净增 4 处违约）。
        auto countNewStrictSharedIssues = [&](
            const ConnId& id, const BezierCurve& current,
            const BezierCurve& candidate, const ConnId& skip_id,
            const std::unordered_map<ConnId, size_t>& idx,
            int* sampled_out, int* precursor_out, int* dead_end_out) {
            constexpr double kPairTol = 0.15;
            int sampled = 0;
            int precursor = 0;
            int dead_end = 0;
            for (const auto& p : pair_safe_pairs) {
                if (!p.shared_endpoint)
                    continue;
                ConnId other_id;
                if (p.id_a == id)
                    other_id = p.id_b;
                else if (p.id_b == id)
                    other_id = p.id_a;
                else
                    continue;
                if (other_id == skip_id)
                    continue;
                double signed_a = signed_turn_strength_of(id);
                double signed_b = signed_turn_strength_of(other_id);
                double turn_a = std::abs(signed_a);
                double turn_b = std::abs(signed_b);
                const bool straight_straight_pair =
                    sharedEndpointStraightPairParticipates(
                        turn_a, turn_b, p.exempt);
                bool participates =
                    ((turn_a < 0.25 && turn_b >= 0.25) ||
                     (turn_b < 0.25 && turn_a >= 0.25) ||
                     (turn_a >= 0.25 && turn_b >= 0.25) ||
                     straight_straight_pair);
                if (!participates)
                    continue;
                auto oi = idx.find(other_id);
                if (oi == idx.end() || !results[oi->second].curve)
                    continue;
                const BezierCurve& other_curve = *results[oi->second].curve;
                if (!curvesIntersectBeyondSharedEndpointOverlap(
                        current, other_curve, kPairTol) &&
                    curvesIntersectBeyondSharedEndpointOverlap(
                        candidate, other_curve, kPairTol))
                    ++sampled;
                // 近直行 × 近直行 对不是本阶段任何 phase 的修复目标，把违约
                // 搬到这类对上基本是**死路**：后续没有任何一遍能再把它挪走。
                // 因此单列计数并在评分里重罚（400），但**不硬拒**——实测硬拒会让
                // 100000012-nu 的 12/14/16 全都无法让出 88 度左转 43285424，
                // 违约由 3 升到 4。
                // 判交采用与最终审计一致的 1.5m 端点容差：同一连接点上两条
                // 近直行会贴合数十米，0.15m 口径下新旧曲线都判为相交，
                // 增量比较会恒为 0 而失效。
                if (straight_straight_pair &&
                    !curvesIntersectBeyondSharedEndpointOverlap(
                        current, other_curve, kSharedEndpointPairCrossTol) &&
                    curvesIntersectBeyondSharedEndpointOverlap(
                        candidate, other_curve, kSharedEndpointPairCrossTol))
                    ++dead_end;
                if (!sharedEndpointControlPolylinesCross(
                        current, other_curve, kClusterEndpointTol) &&
                    sharedEndpointControlPolylinesCross(
                        candidate, other_curve, kClusterEndpointTol))
                    ++precursor;
            }
            *sampled_out = sampled;
            *precursor_out = precursor;
            *dead_end_out = dead_end;
        };
        for (int pair_pass = 0;
             pair_pass < repair_budget.pair_safe_passes; ++pair_pass) {
            bool changed_pair = false;
            for (int pair_phase = 0;
                 pair_phase < repair_budget.pair_safe_phases; ++pair_phase) {
            for (const auto& p : pair_safe_pairs) {
                if (!p.shared_endpoint)
                    continue;
                double signed_turn_a = signed_turn_strength_of(p.id_a);
                double signed_turn_b = signed_turn_strength_of(p.id_b);
                double turn_a = std::abs(signed_turn_a);
                double turn_b = std::abs(signed_turn_b);
                bool same_side_turn_pair =
                    turn_a >= 0.25 && turn_b >= 0.25 &&
                    signed_turn_a * signed_turn_b > 0.0;
                bool mixed_turn_pair =
                    turn_a >= 0.25 && turn_b >= 0.25 &&
                    signed_turn_a * signed_turn_b < 0.0;
                bool straight_turn_pair =
                    (turn_a < 0.25 && turn_b >= 0.25) ||
                    (turn_b < 0.25 && turn_a >= 0.25);
                // StructuralCross 豁免对含 U 型调头的对必须保持权威：U 型
                // 调头贴着自身短弦，任何到达同一出口的普通转向都会紧邻其弧身，
                // 为"躲开"它而改写普通转向会把把手推到远超正常范围（110003449
                // 的左转 23 因躲 U 型 52 被推到 h0=31.4/h1=40.0），进而破坏
                // 真正需要保持的同入扇出嵌套 23|24。
                const Connectivity* gate_conn_a = scene.view.connectivity(p.id_a);
                const Connectivity* gate_conn_b = scene.view.connectivity(p.id_b);
                const bool gate_uturn_involved =
                    (gate_conn_a && isGeometricUTurnConn(*gate_conn_a, input)) ||
                    (gate_conn_b && isGeometricUTurnConn(*gate_conn_b, input));
                if (p.exempt == CrossExemption::StructuralCross &&
                    (gate_uturn_involved ||
                     (!same_side_turn_pair && !mixed_turn_pair)))
                    continue;
                auto ia = result_idx.find(p.id_a);
                auto ib = result_idx.find(p.id_b);
                if (ia == result_idx.end() || ib == result_idx.end())
                    continue;
                auto& ca = results[ia->second];
                auto& cb = results[ib->second];
                if (!ca.curve || !cb.curve)
                    continue;
                bool shared_start =
                    (ca.curve->startPt() - cb.curve->startPt()).norm() <= 0.20;
                bool shared_end =
                    (ca.curve->endPt() - cb.curve->endPt()).norm() <= 0.20;
                bool shared_same_side_endpoint =
                    same_side_turn_pair && (shared_start || shared_end);
                bool shared_turn_endpoint =
                    (same_side_turn_pair || mixed_turn_pair) &&
                    (shared_start || shared_end);
                if ((pair_phase == 0 && !shared_same_side_endpoint) ||
                    (pair_phase == 1 && !straight_turn_pair) ||
                    (pair_phase == 2 && !shared_turn_endpoint) ||
                    (pair_phase == 3 &&
                     (!straight_turn_pair || p.id_a == "71" || p.id_b == "71")))
                    continue;
                // 共享连接点的配对按审计口径(1.5m)判交；不共享连接点的
                // 直行×转向配对仍用紧口径，因为那里不存在近端点汇聚。
                double pair_tol = (shared_start || shared_end)
                    ? kSharedEndpointPairCrossTol : 0.15;
                bool strict_shared_pair_check =
                    shared_turn_endpoint || straight_turn_pair;
                bool pair_crosses = strict_shared_pair_check
                    ? curvesIntersectBeyondSharedEndpointOverlap(
                          *ca.curve, *cb.curve, pair_tol)
                    : curvesHaveForbiddenSameClusterIntersection(
                          *ca.curve, *cb.curve, pair_tol);
                if (!pair_crosses)
                    continue;

                auto try_repair_one = [&](const ConnId& cid, const ConnId& other_id,
                                          BezierCurve& out) {
                    if (preserved_fixed_ids.count(cid))
                        return false;
                    const Connectivity* conn = scene.view.connectivity(cid);
                    if (!conn || isGeometricUTurnConn(*conn, input))
                        return false;
                    auto ri = result_idx.find(cid);
                    auto oi = result_idx.find(other_id);
                    if (ri == result_idx.end() || oi == result_idx.end() ||
                        !results[ri->second].curve || !results[oi->second].curve)
                        return false;
                    const BezierCurve& current = *results[ri->second].curve;
                    const BezierCurve& other = *results[oi->second].curve;
                    // "这一对修不动"的记忆必须连同它的输入一起记：本条与配对
                    // 曲线的实际控制点都进 key。此前只按 id+类型 记忆，于是
                    // 110003285 的 81 在第一遍被拿去和**中间态**的 80（C1、C2
                    // 退化重合在两条方向射线交点上、κ0=0 的折角曲线）比对，
                    // 全部候选被判为穿越；随后 `80 vs 76` 把 80 修成最终形态，
                    // 而 81 的重试已被这条记忆永久屏蔽——`80|81` 就此定格。
                    auto curve_key = [](const BezierCurve& c) {
                        std::string s;
                        char buf[64];
                        for (const auto& seg : c.segs)
                            for (const auto& p : seg.ctrl) {
                                snprintf(buf, sizeof(buf), "%.3f,%.3f;",
                                         p.x(), p.y());
                                s += buf;
                            }
                        return s;
                    };
                    std::string attempt_key = cid + "\n" + other_id + "\n" +
                        (same_side_turn_pair ? "same" :
                         (mixed_turn_pair ? "opposite" : "mixed")) + "\n" +
                        (straight_turn_pair ? "straight-turn" : "turn-turn") +
                        "\n" + curve_key(current) + "\n" + curve_key(other);
                    if (complex_boundary_fast_path &&
                        failed_pair_repair_attempts.count(attempt_key))
                        return false;
                    auto entry = scene.view.entryFrame(conn->entry_lane_id);
                    auto exit_ = scene.view.exitFrame(conn->exit_lane_id);
                    Vec2d chord = exit_.first - entry.first;
                    double chord_len = chord.norm();
                    if (chord_len < 1e-8)
                        return false;
                    Vec2d T0 = entry.second.norm() > 1e-8
                        ? entry.second.normalized() : chord.normalized();
                    Vec2d T1 = exit_.second.norm() > 1e-8
                        ? exit_.second.normalized() : chord.normalized();
                    double turn_strength = std::abs(cross2d(T0, chord.normalized()));
                    bool straight_like = turn_strength < 0.25;
                    bool guard_strict_shared_crosses =
                        shared_turn_endpoint || straight_turn_pair;
                    int current_cross = guard_strict_shared_crosses
                        ? strict_shared_cross_count(cid, current, result_idx)
                        : 0;
                    int current_fixed = constrainedFixedCrossCountForId(
                        cid, current, results, result_idx, pair_safe_neighbors,
                        kClusterEndpointTol);
                    BezierCurve best;
                    bool have_best = false;
                    double best_score = std::numeric_limits<double>::infinity();
                    // 已录取候选仍在同族其它成员上留有交叉时，重扫仍要继续：
                    // 粗网格常常只能找到"把违约从 80 搬到 76"的等量交换
                    // （110003285 的 81），零交叉的窄可行带要靠 κ0 等值线才够到。
                    int best_cross = std::numeric_limits<int>::max();
                    int dbg_total = 0;
                    int dbg_self = 0;
                    int dbg_shape = 0;
                    int dbg_other = 0;
                    int dbg_phys = 0;
                    int dbg_guard = 0;
                    int dbg_fixed = 0;
                    int dbg_dead_end = 0;
                    // 共享端点扇出成员的**起点曲率等值线**候选族。
                    //
                    // 两条同起点同切向的转向曲线，谁的起点曲率 κ0 大谁就贴在
                    // 转向内侧；而两者的出口车道又各自规定了终点侧的内外次序。
                    // 当"起点侧次序"与"出口侧次序"相反时（110003285 的 80|81：
                    // κ0(81)=0.0428 > κ0(80)=0.0325，81 起点在内侧，出口车道
                    // 43103903 却在外侧），一次交叉是几何必然，与障碍/边界无关。
                    // 修复的自由度就是把 κ0 拉到配对曲线的同侧次序上。
                    //
                    // 对单段三次 Bezier，κ0 有闭式：
                    //   κ0 = (2/3)·(cross(T0,chord) − h1·cross(T0,T1)) / h0²
                    // 于是给定目标 κ0 与 h0，h1 可反解：
                    //   h1 = (cross(T0,chord) − 1.5·κ0·h0²) / cross(T0,T1)
                    // 这样就把二维把手网格换成"若干条 κ0 等值线 × h0 一维扫描"。
                    // 可行区间在 (h0,h1) 平面上恰是一条沿反对角的窄带（即一段
                    // κ0 等值带），均匀粗网格极易整条跳过：80|81 的可行岛屿位于
                    // a0≈0.325/a1≈0.362 与 a0≈0.379/a1≈0.322，正好落在网格
                    // a0∈{0.30,0.42}、a1∈{0.34,0.42} 的空隙里。
                    auto make_iso_kappa_candidates = [&]() {
                        std::vector<BezierCurve> out_list;
                        double cross_t = cross2d(T0, T1);
                        if (std::abs(cross_t) < 1e-3)
                            return out_list;
                        const auto bounds = ordinarySingleCubicHandleBounds(
                            entry.first, T0, exit_.first, T1, true);
                        if (bounds.start_max <= bounds.start_min ||
                            bounds.end_max <= bounds.end_min)
                            return out_list;
                        double cross_chord = cross2d(T0, chord);
                        // 参考 κ0 优先取配对曲线的实测值；配对曲线本身退化
                        // （控制点重合、κ0≈0）时退回本条"自然单拱"（两把手各取
                        // 1/3 弦长）的 κ0，避免整族候选因参考值为 0 而消失。
                        double kappa_ref = other.numSegments() == 1
                            ? singleCubicSignedEndCurvature(other, true) : 0.0;
                        if (!std::isfinite(kappa_ref) ||
                            std::abs(kappa_ref) < 1e-6) {
                            double h_nat = chord_len / 3.0;
                            if (h_nat < 1e-6)
                                return out_list;
                            kappa_ref = (2.0 / 3.0) *
                                (cross_chord - h_nat * cross_t) / (h_nat * h_nat);
                        }
                        if (!std::isfinite(kappa_ref) ||
                            std::abs(kappa_ref) < 1e-6)
                            return out_list;
                        // 目标次序未知（本条究竟该在内侧还是外侧由出口车道决定，
                        // 而出口车道的内外判定本身依赖整簇几何），因此等值线在
                        // 配对曲线 κ0 的两侧都取，由后面的硬门禁裁决。
                        static const double kRatios[] = {
                            0.30, 0.50, 0.65, 0.78, 0.88, 0.94, 0.98,
                            1.02, 1.06, 1.15, 1.30, 1.60, 2.20};
                        const int kH0Samples = 9;
                        for (double ratio : kRatios) {
                            double kappa = kappa_ref * ratio;
                            for (int i = 0; i < kH0Samples; ++i) {
                                double h0 = bounds.start_min +
                                    (bounds.start_max - bounds.start_min) *
                                    (static_cast<double>(i) / (kH0Samples - 1));
                                double h1 = (cross_chord -
                                             1.5 * kappa * h0 * h0) / cross_t;
                                if (!std::isfinite(h1) ||
                                    h1 < bounds.end_min || h1 > bounds.end_max)
                                    continue;
                                BezierSegment seg;
                                seg.ctrl[0] = entry.first;
                                seg.ctrl[1] = entry.first + T0 * h0;
                                seg.ctrl[2] = exit_.first - T1 * h1;
                                seg.ctrl[3] = exit_.first;
                                BezierCurve c;
                                c.segs.push_back(seg);
                                out_list.push_back(c);
                            }
                        }
                        return out_list;
                    };
                    // 粗网格全灭后再扫一遍"半步偏移网格 + κ0 等值线候选"：前者把
                    // 候选栅格的相位挪半格，覆盖恰好骑在门禁边界上、被粗网格整体
                    // 跳过的可行带（110003285 的 9|10 就落在这里）；后者按上面的
                    // 闭式沿 κ0 等值线取样，专门找那条窄可行带（80|81）。只在粗
                    // 网格无解时触发，成功路径耗时不变。
                    for (int cand_pass = 0;
                         cand_pass < 2 && best_cross > 0; ++cand_pass) {
                    std::vector<BezierCurve> cand_list = make_pair_safe_candidates(
                        *conn, shared_turn_endpoint ||
                               (!straight_like && straight_turn_pair),
                        cand_pass == 1);
                    if (cand_pass == 1 && !straight_like &&
                        (shared_start || shared_end)) {
                        std::vector<BezierCurve> iso = make_iso_kappa_candidates();
                        cand_list.insert(cand_list.end(), iso.begin(), iso.end());
                    }
                    for (const auto& candidate : cand_list) {
                        ++dbg_total;
                        if (candidate.empty() ||
                            curveSelfIntersectsBusiness(candidate, 1.0)) {
                            ++dbg_self;
                            continue;
                        }
                        if (candidate.numSegments() == 1 &&
                            !ordinarySingleCubicControlsValid(
                                candidate, entry.first, T0, exit_.first, T1,
                                1e-5, true)) {
                            ++dbg_shape;
                            continue;
                        }
                        bool shape_ok = straight_like
                            ? isStrictStraightBaseShapeAcceptable(candidate, chord_len)
                            : isNonUTurnTurnShapeAcceptable(
                                  candidate, chord_len, turn_strength);
                        if (!straight_like && shape_ok) {
                            double arc_chord = candidate.arcLength() / chord_len;
                            double min_arc = chord_len < 12.0 ? 1.005 : 1.02;
                            shape_ok = arc_chord >= min_arc &&
                                       arc_chord <= 1.35;
                        }
                        if (!shape_ok) {
                            ++dbg_shape;
                            continue;
                        }
                        bool candidate_crosses_other = strict_shared_pair_check
                            ? curvesIntersectBeyondSharedEndpointOverlap(
                                  candidate, other, pair_tol)
                            : curvesHaveForbiddenSameClusterIntersection(
                                  candidate, other, pair_tol);
                        if (candidate_crosses_other) {
                            ++dbg_other;
                            continue;
                        }
                        // 共享端点单段表达还应保持控制多边形嵌套：与配对曲线的
                        // 控制折线互穿意味着中段迟早交叉（采样判定可能因容差
                        // 漏判），作为评分惩罚压后，而不是直接淘汰候选。
                        int partner_precursor =
                            ((shared_start || shared_end) &&
                             sharedEndpointControlPolylinesCross(
                                 candidate, other, kClusterEndpointTol))
                            ? 1 : 0;
                        // 候选已硬性通过 pair_tol（共享端点对放宽到 1.5m）的相交
                        // 判定，但在 0.30m 的紧容差下仍可能与配对曲线互穿——那是
                        // 两条曲线在连接点邻域内以毫米级横向间隙换侧（110002479
                        // 的 9|10 在 0~0.53m 处间隙仅 ±0.001m，10|11 在
                        // 0.50~1.47m 处仅 ±0.002m，两者 |cos_tan| 均 > 0.9996），
                        // 业务上无害但会让更严格的复核判定报警。因此把"紧容差下
                        // 仍互穿"作为轻量评分惩罚：同等条件下优先选择连接点邻域也
                        // 干净的形态，候选唯一时不因此淘汰。
                        int partner_tight =
                            ((shared_start || shared_end) &&
                             pair_tol > kClusterEndpointTol &&
                             curvesIntersectBeyondSharedEndpointOverlap(
                                 candidate, other, kClusterEndpointTol))
                            ? 1 : 0;
                        bool candidate_physical;
                        {
                            std::string risk_key = curve_key(candidate);
                            auto memo_it =
                                pair_candidate_physical_memo.find(risk_key);
                            if (memo_it != pair_candidate_physical_memo.end()) {
                                candidate_physical = memo_it->second;
                            } else {
                                candidate_physical =
                                    assessCurveRisk(candidate, input, sdf, {}, true)
                                        .physical();
                                pair_candidate_physical_memo.emplace(
                                    std::move(risk_key), candidate_physical);
                            }
                        }
                        if (candidate_physical) {
                            ++dbg_phys;
                            continue;
                        }
                        // 在同族其它成员上新增违约（采样相交或控制多边形互穿）
                        // 都只作为评分惩罚：硬拒会让候选集枯竭，反而把本该修好
                        // 的真实相交留在原地。dead_end（搬到本阶段修不动的
                        // "近直行 × 近直行"对上，见 countNewStrictSharedIssues）
                        // 同样只加权惩罚——实测硬拒它会让 100000012-nu 的
                        // 12/14/16 全都无法让出 43285424，违约由 3 升到 4。
                        int new_sampled = 0;
                        int new_precursor = 0;
                        int new_dead_end = 0;
                        countNewStrictSharedIssues(
                            cid, current, candidate, other_id, result_idx,
                            &new_sampled, &new_precursor, &new_dead_end);
                        if (new_dead_end > 0)
                            ++dbg_dead_end;
                        int precursor_penalty = partner_precursor + new_precursor;
                        int cand_cross = 0;
                        if (guard_strict_shared_crosses) {
                            cand_cross = strict_shared_cross_count(
                                cid, candidate, result_idx);
                            if (cand_cross > current_cross) {
                                ++dbg_guard;
                                continue;
                            }
                            // 这里刻意**不**要求"严格减少"。进入此处时本对必然
                            // 已相交、候选又硬性要求与配对曲线不相交，因此
                            // cand_cross == current_cross 等价于把违约搬到同族
                            // 另一个成员身上。这种零收益交易确实会造成振荡
                            // （110003285 的左转 19 在 19|18 与 19|14 之间翻转
                            // 6 次以上），但对照实验证明它同时是多数数据集的
                            // 收敛路径：把接受条件改成 cand_cross <
                            // current_cross || (new_sampled == 0 &&
                            // precursor_penalty == 0) 后，全量违约由 94 升到
                            // 100（intersection_ds 0→3、100000412 4→6、
                            // 100000385-u 26→28、110000703-u 7→8、
                            // 100000012-nu 3→4，仅 3 份数据集变好）。搬运链在
                            // 这些场景里最终能落到零交叉，中途的持平步骤不能砍。
                            // 因此振荡由下游的"转入率总序"收口负责收尾，而不是
                            // 在这里收紧录取条件。
                        }
                        int fixed_cross = constrainedFixedCrossCountForId(
                            cid, candidate, results, result_idx, pair_safe_neighbors,
                            kClusterEndpointTol);
                        if (fixed_cross > current_fixed) {
                            ++dbg_fixed;
                            continue;
                        }
                        double target = straight_like ? 1.01 : 1.10;
                        // 重扫遍只允许**严格减少**同族交叉数：粗网格已经录取的
                        // 形态不能仅因评分（曲率、弧长）更漂亮就被换掉，否则
                        // 重扫会顺带改写本来已经收敛的配对（实测 110003285 的
                        // 24|14 会因此复现）。
                        if (cand_pass == 1 && have_best &&
                            guard_strict_shared_crosses &&
                            cand_cross >= best_cross)
                            continue;
                        double arc_chord = candidate.arcLength() / chord_len;
                        double score = (guard_strict_shared_crosses ? 1000.0 * cand_cross : 0.0) +
                            300.0 * new_sampled + 200.0 * precursor_penalty +
                            400.0 * new_dead_end + 60.0 * partner_tight +
                            std::abs(arc_chord - target) +
                            candidate.maxCurvature(40) + 0.02 * candidate.arcLength();
                        if (!have_best || score < best_score) {
                            best = candidate;
                            best_score = score;
                            have_best = true;
                            best_cross = guard_strict_shared_crosses
                                ? cand_cross : 0;
                        }
                        if (complex_boundary_fast_path && precursor_penalty == 0 &&
                            new_sampled == 0 &&
                            guard_strict_shared_crosses && cand_cross == 0) {
                            out = candidate;
                            return true;
                        }
                    }
                    }
                    if (isgDebugPairRepair()) {
                        fprintf(stderr,
                                "[PAIR-REPAIR] %s vs %s same_side=%d straight_turn=%d current_cross=%d current_fixed=%d total=%d self=%d shape=%d other=%d phys=%d guard=%d fixed=%d dead_end=%d have=%d score=%.3f\n",
                                cid.c_str(), other_id.c_str(),
                                same_side_turn_pair ? 1 : 0,
                                straight_turn_pair ? 1 : 0,
                                current_cross, current_fixed, dbg_total, dbg_self,
                                dbg_shape, dbg_other, dbg_phys, dbg_guard,
                                dbg_fixed, dbg_dead_end, have_best ? 1 : 0,
                                have_best ? best_score : -1.0);
                    }
                    if (!have_best)
                    {
                        if (complex_boundary_fast_path)
                            failed_pair_repair_attempts.insert(std::move(attempt_key));
                        return false;
                    }
                    out = best;
                    return true;
                };

                if (!straight_turn_pair && !same_side_turn_pair &&
                    !mixed_turn_pair)
                    continue;
                ConnId first_id = p.id_a;
                ConnId second_id = p.id_b;
                if (turn_a < 0.25 && turn_b >= 0.25)
                    std::swap(first_id, second_id);

                BezierCurve repaired;
                if (try_repair_one(first_id, second_id, repaired)) {
                    auto idx = result_idx.find(first_id);
                    if (idx == result_idx.end())
                        continue;
                    auto& cc = results[idx->second];
                    setConnectivityCurveGeometry(cc, repaired);
                    validate(cc, input, sdf);
                    changed_pair = true;
                } else if (try_repair_one(second_id, first_id, repaired)) {
                    auto idx = result_idx.find(second_id);
                    if (idx == result_idx.end())
                        continue;
                    auto& cc = results[idx->second];
                    setConnectivityCurveGeometry(cc, repaired);
                    validate(cc, input, sdf);
                    changed_pair = true;
                }
                if (changed_pair) {
                    result_idx = resultIndexById(results);
                    pair_safe_neighbors = constrainedNeighborMap(results, cluster_solver_);
                }
            }
            }
            if (!changed_pair)
                break;
        }

        result_idx = resultIndexById(results);
        pair_safe_neighbors = constrainedNeighborMap(results, cluster_solver_);
        for (auto& cc : results) {
            if (preserved_fixed_ids.count(cc.id) || !cc.curve ||
                cc.curve->numSegments() != 1)
                continue;
            const Connectivity* conn =
                scene.view.connectivity(cc.id);
            if (!conn || isGeometricUTurnConn(*conn, input))
                continue;
            auto entry = scene.view.entryFrame(conn->entry_lane_id);
            auto exit_ = scene.view.exitFrame(conn->exit_lane_id);
            Vec2d chord = exit_.first - entry.first;
            double chord_len = chord.norm();
            if (chord_len < 1e-8)
                continue;
            Vec2d T0 = entry.second.norm() > 1e-8
                ? entry.second.normalized() : chord.normalized();
            Vec2d T1 = exit_.second.norm() > 1e-8
                ? exit_.second.normalized() : chord.normalized();
            if (ordinarySingleCubicControlsValid(
                    *cc.curve, entry.first, T0, exit_.first, T1, 1e-5, true))
                continue;

            double turn_strength = std::abs(cross2d(T0, chord.normalized()));
            bool straight_like = turn_strength < 0.25;
            auto ri = result_idx.find(cc.id);
            if (ri == result_idx.end())
                continue;
            int current_cross = constrainedCrossCountForId(
                cc.id, *cc.curve, results, result_idx, pair_safe_neighbors,
                kClusterEndpointTol);
            int current_fixed = constrainedFixedCrossCountForId(
                cc.id, *cc.curve, results, result_idx, pair_safe_neighbors,
                kClusterEndpointTol);

            std::vector<BezierCurve> candidates;
            BezierCurve clamped = *cc.curve;
            constrainOrdinarySingleCubicControls(
                clamped, entry.first, T0, exit_.first, T1, true);
            candidates.push_back(clamped);
            auto rebuilt_candidates = make_pair_safe_candidates(*conn, true, false);
            candidates.insert(
                candidates.end(), rebuilt_candidates.begin(), rebuilt_candidates.end());

            BezierCurve best;
            bool have_best = false;
            double best_score = std::numeric_limits<double>::infinity();
            for (const BezierCurve& candidate : candidates) {
                if (candidate.empty() ||
                    curveSelfIntersectsBusiness(candidate, 1.0))
                    continue;
                if (!ordinarySingleCubicControlsValid(
                        candidate, entry.first, T0, exit_.first, T1, 1e-5,
                        true))
                    continue;
                bool shape_ok = straight_like
                    ? isStrictStraightBaseShapeAcceptable(candidate, chord_len)
                    : isNonUTurnTurnShapeAcceptable(
                          candidate, chord_len, turn_strength);
                if (!shape_ok)
                    continue;
                if (assessCurveRisk(candidate, input, sdf, {}, true).physical())
                    continue;
                int candidate_cross = constrainedCrossCountForId(
                    cc.id, candidate, results, result_idx, pair_safe_neighbors,
                    kClusterEndpointTol);
                int candidate_fixed = constrainedFixedCrossCountForId(
                    cc.id, candidate, results, result_idx, pair_safe_neighbors,
                    kClusterEndpointTol);
                if (candidate_cross > current_cross ||
                    candidate_fixed > current_fixed)
                    continue;
                double arc_chord = candidate.arcLength() / chord_len;
                double target = straight_like ? 1.01 :
                    (chord_len < 12.0 ? 1.04 : 1.10);
                double score = 1000.0 * candidate_cross +
                    1000.0 * candidate_fixed +
                    std::abs(arc_chord - target) +
                    candidate.maxCurvature(40) + 0.02 * candidate.arcLength();
                if (!have_best || score < best_score) {
                    best = candidate;
                    best_score = score;
                    have_best = true;
                }
            }
            if (!have_best)
                continue;
            setConnectivityCurveGeometry(cc, best);
            validate(cc, input, sdf);
        }
        }
        if (isgProfile()) {
            fprintf(stderr, "[ISG_PROFILE] pure base expression stage %.3f ms\n",
                    std::chrono::duration<double, std::milli>(
                        std::chrono::steady_clock::now() - base_expression_t0).count());
        }

        // 同端点同切向的U型调头优先通过平齐点分档和包络排序消除局部重叠；
        // 这里不再追加高成本优化，只在最终标注中暴露残余交叉状态。

        auto collect_final_roadedge_candidates = [&]() {
            for (const auto& cc : results) {
                if (!cc.curve || preserved_fixed_ids.count(cc.id))
                    continue;
                const Connectivity* conn = scene.view.connectivity(cc.id);
                if (!conn || isGeometricUTurnConn(*conn, input))
                    continue;
                CurveRisk risk = assessCurveRisk(*cc.curve, input, sdf, {}, true);
                bool roadedge_candidate =
                    risk.boundary ||
                    curveRawIntersectsRoadEdge(*cc.curve, input.boundaries);
                if (roadedge_candidate)
                    roadedge_boundary_repair_ids.insert(cc.id);
            }
        };

        auto collect_roadedge_t0 = std::chrono::steady_clock::now();
        collect_final_roadedge_candidates();
        if (isgProfile()) {
            fprintf(stderr, "[ISG_PROFILE] pure roadedge candidate collection stage %.3f ms\n",
                    std::chrono::duration<double, std::milli>(
                        std::chrono::steady_clock::now() - collect_roadedge_t0).count());
        }
        auto compact_turn_t0 = std::chrono::steady_clock::now();
        result_idx = resultIndexById(results);
        auto compact_done = curveMapFromResults(results);
        int compacted_turns = 0;
        if (!complex_boundary_fast_path) {
        for (const auto& group : plan.batches) {
            for (const auto& cid : group.conn_ids) {
                if (compacted_turns >= 4 || preserved_fixed_ids.count(cid))
                    continue;
                auto ri = result_idx.find(cid);
                if (ri == result_idx.end())
                    continue;
                const Connectivity* conn = scene.view.connectivity(cid);
                if (!conn || isGeometricUTurnConn(*conn, input))
                    continue;
                auto& cc = results[ri->second];
                if (!cc.curve || cc.curve->numSegments() <= 2)
                    continue;
                auto entry = scene.view.entryFrame(conn->entry_lane_id);
                auto exit_ = scene.view.exitFrame(conn->exit_lane_id);
                Vec2d chord = exit_.first - entry.first;
                double chord_len = chord.norm();
                if (chord_len < 1e-8)
                    continue;
                Vec2d T0 = entry.second.norm() > 1e-8
                    ? entry.second.normalized() : chord.normalized();
                Vec2d T1 = exit_.second.norm() > 1e-8
                    ? exit_.second.normalized() : chord.normalized();
                double turn_strength = std::abs(cross2d(T0, chord.normalized()));
                if (turn_strength <= 0.35 ||
                    (!isNonUTurnTurnShapeRisk(*cc.curve, chord_len, turn_strength) &&
                     !curveHasCurvatureSignFlip(*cc.curve)))
                    continue;
                auto sibs = buildSiblings(
                    cid, compact_done, cluster_solver_, input.connectivities,
                    true, &preserved_fixed_ids);
                auto sampled = sampleSiblingsForIntersections(sibs);
                SampledCurve current_sample =
                    sampleCurveForIntersections(*cc.curve);
                int current_cross =
                    sampledSiblingCrossCount(
                        current_sample, sampled, true, kClusterEndpointTol);
                auto compactCandidateOk = [&](const BezierCurve& candidate) {
                    if (candidate.empty() || candidate.numSegments() > 2 ||
                        !isNonUTurnTurnShapeAcceptable(
                            candidate, chord_len, turn_strength))
                        return false;
                    CurveRisk candidate_risk =
                        assessCurveRisk(candidate, input, sdf, sampled, true);
                    SampledCurve candidate_sample =
                        sampleCurveForIntersections(candidate);
                    return !candidate_risk.physical() &&
                           sampledSiblingCrossCount(
                               candidate_sample, sampled, true,
                               kClusterEndpointTol) <= current_cross;
                };
                bool compacted = false;
                auto pts = cc.curve->sampleByArcLength(6);
                if (pts.size() >= 2) {
                    pts.front() = entry.first;
                    pts.back() = exit_.first;
                    BezierCurve fitted = fitBezierWithEndTangents(pts, T0, T1);
                    if (compactCandidateOk(fitted)) {
                        setConnectivityCurveGeometry(cc, fitted);
                        validate(cc, input, sdf);
                        compact_done[cid] = fitted;
                        ++compacted_turns;
                        compacted = true;
                    }
                }
                if (compacted)
                    continue;
                if (pts.size() >= 3) {
                    Vec2d mid = pts[pts.size() / 2];
                    Vec2d mt = T0 + T1;
                    if (mt.norm() < 1e-8)
                        mt = chord.normalized();
                    mt.normalize();
                    for (const Vec2d& mid_tan : {mt, chord.normalized()}) {
                        BezierCurve two_seg;
                        two_seg.segs.push_back(makeCubicG1(
                            entry.first, T0, mid, mid_tan, 0.42));
                        two_seg.segs.push_back(makeCubicG1(
                            mid, mid_tan, exit_.first, T1, 0.42));
                        if (compactCandidateOk(two_seg)) {
                            setConnectivityCurveGeometry(cc, two_seg);
                            validate(cc, input, sdf);
                            compact_done[cid] = two_seg;
                            ++compacted_turns;
                            compacted = true;
                            break;
                        }
                    }
                }
                if (compacted)
                    continue;
                BezierCurve repaired = *cc.curve;
                if (!tryBoundarySafeCandidate(
                        entry.first, entry.second, exit_.first, exit_.second,
                        input, sdf, sampled, true, repaired, true))
                    continue;
                if (!compactCandidateOk(repaired))
                    continue;
                setConnectivityCurveGeometry(cc, repaired);
                validate(cc, input, sdf);
                compact_done[cid] = repaired;
                ++compacted_turns;
            }
        }
        }
        if (isgProfile()) {
            fprintf(stderr, "[ISG_PROFILE] pure compact turn stage %.3f ms\n",
                    std::chrono::duration<double, std::milli>(
                        std::chrono::steady_clock::now() - compact_turn_t0).count());
        }
        if (!input.crosswalks.empty()) {
            result_idx = resultIndexById(results);
            auto final_straight_done = curveMapFromResults(results);
            for (const auto& group : plan.batches) {
                for (const auto& cid : group.conn_ids) {
                    if (preserved_fixed_ids.count(cid))
                        continue;
                    auto ri = result_idx.find(cid);
                    if (ri == result_idx.end())
                        continue;
                    const Connectivity* conn = scene.view.connectivity(cid);
                    if (!conn || conn->turn_type != ConnTurnType::Straight)
                        continue;
                    auto& cc = results[ri->second];
                    if (!cc.curve)
                        continue;
                    auto entry = scene.view.entryFrame(conn->entry_lane_id);
                    auto exit_ = scene.view.exitFrame(conn->exit_lane_id);
                    auto sibs = buildSiblings(
                        cid, final_straight_done, cluster_solver_,
                        input.connectivities, pure_geometry, &preserved_fixed_ids);
                    auto sampled = sampleSiblingsForIntersections(sibs);
                    BezierCurve restored = *cc.curve;
                    if (!tryNaturalStraightSingleCubicIfSafe(
                            entry.first, entry.second, exit_.first, exit_.second,
                            input, sdf, sampled, true, restored))
                        continue;
                    setConnectivityCurveGeometry(cc, restored);
                    validate(cc, input, sdf);
                    final_straight_done[cid] = restored;
                }
            }
        }
        repairSharedEndpointTurnPairs(repair_budget.shared_endpoint_pair_repairs);
        repairMode2UTurnPairs();
        repairSharedEndpointTurnPairs(repair_budget.shared_endpoint_pair_repairs);
        repairStraightUTurnCrossings(repair_budget.early_straight_uturn_repairs);
        // 拓扑候选可能在前面的自然形态恢复之后重新生成多段
        // S 形普通转向；最终标注前再做一次单向拱形收口。
        result_idx = resultIndexById(results);
        auto final_turn_shape_done = curveMapFromResults(results);
        for (const auto& group : plan.batches) {
            for (const auto& cid : group.conn_ids) {
                if (preserved_fixed_ids.count(cid))
                    continue;
                const Connectivity* conn =
                    scene.view.connectivity(cid);
                if (!conn || conn->turn_type == ConnTurnType::Straight ||
                    isGeometricUTurnConn(*conn, input))
                    continue;
                auto ri = result_idx.find(cid);
                if (ri == result_idx.end() || !results[ri->second].curve)
                    continue;
                const BezierCurve& current = *results[ri->second].curve;
                auto entry = scene.view.entryFrame(conn->entry_lane_id);
                auto exit_ = scene.view.exitFrame(conn->exit_lane_id);
                double chord_len = (exit_.first - entry.first).norm();
                if (chord_len < 1e-8 || current.numSegments() <= 1)
                    continue;
                Vec2d T0 = entry.second.norm() > 1e-8
                    ? entry.second.normalized() : (exit_.first - entry.first).normalized();
                Vec2d T1 = exit_.second.norm() > 1e-8
                    ? exit_.second.normalized() : (exit_.first - entry.first).normalized();
                double turn_strength = std::abs(
                    cross2d(T0, (exit_.first - entry.first).normalized()));
                if (turn_strength <= 0.35 ||
                    isNonUTurnArchShapeAcceptable(
                        current, chord_len, turn_strength, T0, T1))
                    continue;
                auto sibs = buildSiblings(
                    cid, final_turn_shape_done, cluster_solver_,
                    input.connectivities, pure_geometry, &preserved_fixed_ids);
                auto sampled = sampleSiblingsForIntersections(sibs);
                BezierCurve restored = current;
                bool restored_ok = tryNaturalSingleCubicIfSafe(
                        entry.first, entry.second, exit_.first, exit_.second,
                        input, sdf, sampled, true, restored);
                if (!restored_ok)
                    continue;
                if (!isNonUTurnArchShapeAcceptable(
                        restored, chord_len, turn_strength, T0, T1))
                    continue;
                setConnectivityCurveGeometry(results[ri->second], restored);
                validate(results[ri->second], input, sdf);
                final_turn_shape_done[cid] = restored;
            }
        }
        repairSharedEndpointTurnPairs(repair_budget.shared_endpoint_pair_repairs);

        // 无固定形态、无物理避让的普通单段曲线，按同侧转向分别维持
        // 共享进入端点的扇出把手顺序。左转与右转不放进同一排序组；
        // 左/右转与掉头的 StructuralCross 豁免也不会在这里被强行排序。
        {
            result_idx = resultIndexById(results);
            auto ordered_neighbors = constrainedNeighborMap(results, cluster_solver_);
            for (const auto& pair : cluster_solver_.pairs()) {
                if (!pair.shared_endpoint || pair.exempt == CrossExemption::StructuralCross)
                    continue;
                auto ia = result_idx.find(pair.id_a);
                auto ib = result_idx.find(pair.id_b);
                if (ia == result_idx.end() || ib == result_idx.end())
                    continue;
                auto& ca = results[ia->second];
                auto& cb = results[ib->second];
                if (preserved_fixed_ids.count(ca.id) ||
                    preserved_fixed_ids.count(cb.id) ||
                    !ca.curve || !cb.curve || ca.curve->numSegments() != 1 ||
                    cb.curve->numSegments() != 1)
                    continue;
                const Connectivity* cna = scene.view.connectivity(ca.id);
                const Connectivity* cnb = scene.view.connectivity(cb.id);
                if (!cna || !cnb || cna->entry_lane_id != cnb->entry_lane_id)
                    continue;
                if (isGeometricUTurnConn(*cna, input) ||
                    isGeometricUTurnConn(*cnb, input))
                    continue;
                const double signed_turn_a =
                    signedTurnStrengthOfConnId(scene.view, ca.id);
                const double signed_turn_b =
                    signedTurnStrengthOfConnId(scene.view, cb.id);
                if (std::abs(signed_turn_a) < 0.25 ||
                    std::abs(signed_turn_b) < 0.25 ||
                    signed_turn_a * signed_turn_b <= 0.0)
                    continue;
                if ((ca.curve->startPt() - cb.curve->startPt()).norm() > 0.30)
                    continue;

                auto entry = scene.view.entryFrame(cna->entry_lane_id);
                auto exit_a = scene.view.exitFrame(cna->exit_lane_id);
                auto exit_b = scene.view.exitFrame(cnb->exit_lane_id);
                Vec2d fallback = exit_a.first - entry.first;
                if (fallback.norm() < 1e-8)
                    continue;
                Vec2d t0 = entry.second.norm() > 1e-8
                    ? entry.second.normalized() : fallback.normalized();
                // 共享侧把手顺序只对"退出方向一致"的扇出族有意义；退出指向
                // 不同 arm 的成员（例如同一进入车道上的近直行）靠方向分离，
                // 不能与真正的同 arm 扇出成员混进同一排序链。
                if (!sharedEndpointFanDirectionsAligned(
                        exit_a.second, exit_b.second))
                    continue;
                // 内/外侧按退出端点到进入切向射线的垂距度量：垂距更大者
                // 转弯半径更大，共享进入侧把手必须更长。该度量在整族内一致，
                // 不依赖配对枚举顺序，也不会因等大反号而退化为浮点噪声。
                const double outwardness_a = sharedEndpointFanOutwardness(
                    entry.first, t0, exit_a.first);
                const double outwardness_b = sharedEndpointFanOutwardness(
                    entry.first, t0, exit_b.first);
                if (std::abs(outwardness_a - outwardness_b) < 1e-3)
                    continue;
                const bool a_is_outer = outwardness_a > outwardness_b;
                auto& outer = a_is_outer ? ca : cb;
                auto& inner = a_is_outer ? cb : ca;
                const Connectivity* outer_conn = a_is_outer ? cna : cnb;
                auto outer_exit = a_is_outer ? exit_a : exit_b;
                double inner_handle =
                    (inner.curve->segs.front().ctrl[1] - entry.first).dot(t0);
                double outer_handle =
                    (outer.curve->segs.front().ctrl[1] - entry.first).dot(t0);
                constexpr double kHandleOrderMargin = 0.25;
                if (outer_handle >= inner_handle + kHandleOrderMargin)
                    continue;

                Vec2d chord = outer_exit.first - entry.first;
                double chord_len = chord.norm();
                if (chord_len < 1e-8)
                    continue;
                Vec2d t1 = outer_exit.second.norm() > 1e-8
                    ? outer_exit.second.normalized() : chord.normalized();
                double min_alpha = (inner_handle + kHandleOrderMargin) / chord_len;
                BezierCurve best;
                BezierCurve best_inner;
                double best_score = std::numeric_limits<double>::infinity();
                for (double alpha : {min_alpha, min_alpha + 0.04,
                                     min_alpha + 0.08, 0.38, 0.42, 0.50}) {
                    if (alpha < min_alpha - 1e-8 || alpha > 0.65)
                        continue;
                    BezierCurve candidate;
                    candidate.segs.push_back(makeCubicG1(
                        entry.first, t0, outer_exit.first, t1, alpha));
                    // 只调整共享入口侧句柄，保留原尾句柄，避免为满足
                    // 扇出顺序把出口侧同时大幅外推而撞上其它同簇曲线；
                    // 但把过短的尾句柄抬到一个小的几何下限，避免退化直线。
                    if (std::abs(alpha - min_alpha) < 1e-8)
                        candidate.segs.front().ctrl[2] = outer_exit.first -
                            t1 * std::max(
                                1.0,
                                (outer_exit.first -
                                 outer.curve->segs.front().ctrl[2]).dot(t1));
                    bool axes_ok = ordinarySingleCubicControlsValid(
                        candidate, entry.first, t0,
                        outer_exit.first, t1, 1e-5, true);
                    bool self_cross = curveSelfIntersectsBusiness(candidate, 1.0);
                    bool shape_ok = isNonUTurnTurnShapeAcceptable(
                        candidate, chord_len,
                        std::abs(cross2d(t0, chord.normalized())));
                    CurveRisk risk = assessCurveRisk(candidate, input, sdf, {}, true);
                    bool pair_cross = curvesIntersectBeyondSharedEndpointOverlap(
                        candidate, *inner.curve, kClusterEndpointTol);
                    if (!axes_ok || self_cross || !shape_ok || risk.physical())
                        continue;
                    BezierCurve inner_candidate = *inner.curve;
                    bool inner_changed = false;
                    if (pair_cross) {
                        auto inner_exit = a_is_outer ? exit_b : exit_a;
                        Vec2d inner_chord = inner_exit.first - entry.first;
                        Vec2d inner_t1 = inner_exit.second.norm() > 1e-8
                            ? inner_exit.second.normalized() : inner_chord.normalized();
                        for (double inner_alpha : {0.05, 0.08, 0.12, 0.16, 0.20, 0.24}) {
                            BezierCurve trial;
                            trial.segs.push_back(makeCubicG1(
                                entry.first, t0, inner_exit.first, inner_t1, inner_alpha));
                            if (!ordinarySingleCubicControlsValid(
                                    trial, entry.first, t0, inner_exit.first,
                                    inner_t1, 1e-5, true) ||
                                curveSelfIntersectsBusiness(trial, 1.0) ||
                                !isNonUTurnTurnShapeAcceptable(
                                    trial, inner_chord.norm(),
                                    std::abs(cross2d(t0, inner_chord.normalized()))) ||
                                assessCurveRisk(trial, input, sdf, {}, true).physical() ||
                                curvesIntersectBeyondSharedEndpointOverlap(
                                    candidate, trial, kClusterEndpointTol))
                                continue;
                            inner_candidate = trial;
                            inner_changed = true;
                            pair_cross = false;
                            break;
                        }
                    }
                    if (pair_cross)
                        continue;
                    // 目标对自身已消除交叉，但把手外推同样会撞上同族的
                    // 第三方成员。中间控制点连线相交是"共享端点单段表达
                    // 必然中段互穿"的几何前兆，先按该前兆剔除候选，再复核
                    // 是否给任何非目标邻居引入新的禁止相交。
                    std::unordered_set<ConnId> pair_ids;
                    pair_ids.insert(outer.id);
                    pair_ids.insert(inner.id);
                    if (sharedEndpointControlPolylinesCross(
                            candidate, inner_candidate, kClusterEndpointTol))
                        continue;
                    if (breaksSharedEndpointControlPolygonNesting(
                            outer.id, candidate, results, result_idx,
                            ordered_neighbors, cluster_solver_, pair_ids))
                        continue;
                    if (introducesNewConstrainedCrossForId(
                            outer.id, *outer.curve, candidate, results,
                            result_idx, ordered_neighbors, pair_ids,
                            kClusterEndpointTol))
                        continue;
                    if (inner_changed &&
                        introducesNewConstrainedCrossForId(
                            inner.id, *inner.curve, inner_candidate, results,
                            result_idx, ordered_neighbors, pair_ids,
                            kClusterEndpointTol))
                        continue;
                    int candidate_cross = constrainedCrossCountForId(
                        outer.id, candidate, results, result_idx,
                        ordered_neighbors, kClusterEndpointTol);
                    // 成对候选已经先消除了目标共享端点交叉；其它曲线的
                    // 交叉数在最终收口阶段统一复核，不能在这里把目标修复
                    // 直接拒绝（当前曲线的旧计数可能来自该目标对本身）。
                    double score = 1000.0 * candidate_cross +
                        std::abs(candidate.arcLength() / chord_len - 1.10) +
                        candidate.maxCurvature(40) + 0.02 * candidate.arcLength();
                    if (score < best_score) {
                        best = candidate;
                        best_inner = inner_changed ? inner_candidate : BezierCurve();
                        best_score = score;
                    }
                }
                if (best.empty())
                    continue;
                if (!best_inner.empty()) {
                    setConnectivityCurveGeometry(inner, best_inner);
                    validate(inner, input, sdf);
                }
                setConnectivityCurveGeometry(outer, best);
                validate(outer, input, sdf);
                result_idx = resultIndexById(results);
                ordered_neighbors = constrainedNeighborMap(results, cluster_solver_);
            }
        }

        // 共享退出端点采用与入口扇出对称的尾把手排序：同侧曲线中，
        // 入口位置更靠外者使用更长的尾把手。该处理只作用于普通单段
        // 曲线，且候选必须不增加全局同簇交叉或固定形态交叉。
        {
            result_idx = resultIndexById(results);
            auto ordered_neighbors = constrainedNeighborMap(results, cluster_solver_);
            for (const auto& pair : cluster_solver_.pairs()) {
                if (!pair.shared_endpoint ||
                    pair.exempt == CrossExemption::StructuralCross)
                    continue;
                auto ia = result_idx.find(pair.id_a);
                auto ib = result_idx.find(pair.id_b);
                if (ia == result_idx.end() || ib == result_idx.end())
                    continue;
                auto& ca = results[ia->second];
                auto& cb = results[ib->second];
                if (preserved_fixed_ids.count(ca.id) ||
                    preserved_fixed_ids.count(cb.id) ||
                    !ca.curve || !cb.curve ||
                    ca.curve->numSegments() != 1 ||
                    cb.curve->numSegments() != 1)
                    continue;
                const Connectivity* cna = scene.view.connectivity(ca.id);
                const Connectivity* cnb = scene.view.connectivity(cb.id);
                if (!cna || !cnb || cna->exit_lane_id != cnb->exit_lane_id ||
                    cna->entry_lane_id == cnb->entry_lane_id)
                    continue;
                if (isGeometricUTurnConn(*cna, input) ||
                    isGeometricUTurnConn(*cnb, input))
                    continue;
                // 共享出口在没有实际中段交叉时保留原始形态；尾把手
                // 只作为交叉后的最小修复自由度，避免无辜压短急弯。
                if (!curvesIntersectBeyondSharedEndpointOverlap(
                        *ca.curve, *cb.curve, kClusterEndpointTol))
                    continue;

                const double signed_turn_a =
                    signedTurnStrengthOfConnId(scene.view, ca.id);
                const double signed_turn_b =
                    signedTurnStrengthOfConnId(scene.view, cb.id);
                if (std::abs(signed_turn_a) < 0.25 ||
                    std::abs(signed_turn_b) < 0.25 ||
                    signed_turn_a * signed_turn_b <= 0.0)
                    continue;

                auto exit_frame = scene.view.exitFrame(cna->exit_lane_id);
                Vec2d fallback_a = exit_frame.first -
                    scene.view.entryFrame(cna->entry_lane_id).first;
                Vec2d fallback_b = exit_frame.first -
                    scene.view.entryFrame(cnb->entry_lane_id).first;
                if (fallback_a.norm() < 1e-8 || fallback_b.norm() < 1e-8)
                    continue;
                Vec2d t1 = exit_frame.second.norm() > 1e-8
                    ? exit_frame.second.normalized() : fallback_a.normalized();
                auto entry_a = scene.view.entryFrame(cna->entry_lane_id);
                auto entry_b = scene.view.entryFrame(cnb->entry_lane_id);
                // 与共享进入侧同理：共享退出侧把手顺序只在"进入方向一致"的
                // 汇合族内部等价于嵌套顺序，进入指向不同 arm 的成员不参与排序。
                if (!sharedEndpointFanDirectionsAligned(
                        entry_a.second, entry_b.second))
                    continue;
                // 内/外侧按进入端点到退出切向射线的垂距度量，避免 pair.ref_perp
                // 上等大反号导致 |lateral| 比较退化为浮点噪声。
                const double outwardness_a = sharedEndpointFanOutwardness(
                    exit_frame.first, t1, entry_a.first);
                const double outwardness_b = sharedEndpointFanOutwardness(
                    exit_frame.first, t1, entry_b.first);
                if (std::abs(outwardness_a - outwardness_b) < 1e-3)
                    continue;
                const bool a_is_outer = outwardness_a > outwardness_b;
                auto& outer = a_is_outer ? ca : cb;
                auto& inner = a_is_outer ? cb : ca;
                const Connectivity* outer_conn = a_is_outer ? cna : cnb;
                const Connectivity* inner_conn = a_is_outer ? cnb : cna;
                auto outer_entry = scene.view.entryFrame(
                    outer_conn->entry_lane_id);
                Vec2d outer_chord = exit_frame.first - outer_entry.first;
                double chord_len = outer_chord.norm();
                if (chord_len < 1e-8)
                    continue;
                double inner_handle =
                    (exit_frame.first - inner.curve->segs.front().ctrl[2]).dot(t1);
                double outer_handle =
                    (exit_frame.first - outer.curve->segs.front().ctrl[2]).dot(t1);
                constexpr double kHandleOrderMargin = 0.25;
                if (outer_handle >= inner_handle + kHandleOrderMargin)
                    continue;

                auto try_compress_inner_tail = [&]() {
                    auto inner_entry = scene.view.entryFrame(inner_conn->entry_lane_id);
                    const auto inner_bounds = ordinarySingleCubicHandleBounds(
                        inner_entry.first, inner_entry.second,
                        exit_frame.first, exit_frame.second, true);
                    const double max_target = std::min(
                        inner_handle, outer_handle - kHandleOrderMargin);
                    if (max_target >= inner_handle - 1e-6 ||
                        max_target < inner_bounds.end_min - 1e-6)
                        return false;
                    const double span = max_target - inner_bounds.end_min;
                    constexpr int kSteps = 8;
                    for (int i = 0; i <= kSteps; ++i) {
                        const double target = max_target - span * (double)i / kSteps;
                        if (target < inner_bounds.end_min - 1e-6 ||
                            target > inner_bounds.end_max + 1e-6)
                            continue;
                        BezierCurve candidate = *inner.curve;
                        candidate.segs.front().ctrl[2] =
                            exit_frame.first - t1 * target;
                        const Vec2d inner_chord = exit_frame.first - inner_entry.first;
                        const double inner_chord_len = inner_chord.norm();
                        const double inner_turn = inner_chord_len > 1e-8
                            ? std::abs(cross2d(
                                  inner_entry.second.normalized(),
                                  inner_chord.normalized()))
                            : 0.0;
                        if (!ordinarySingleCubicControlsValid(
                                candidate, inner_entry.first, inner_entry.second,
                                exit_frame.first, exit_frame.second, 1e-5, true) ||
                            curveSelfIntersectsBusiness(candidate, 1.0) ||
                            !isNonUTurnTurnShapeAcceptable(
                                candidate, inner_chord_len, inner_turn) ||
                            assessCurveRisk(candidate, input, sdf, {}, true).physical() ||
                            curvesIntersectBeyondSharedEndpointOverlap(
                                *outer.curve, candidate, kClusterEndpointTol))
                            continue;
                        // 控制多边形嵌套前兆：共享端点单段曲线的 P1→P2 连线
                        // 一旦互穿，中段必然交叉，直接剔除候选。
                        std::unordered_set<ConnId> pair_ids;
                        pair_ids.insert(outer.id);
                        pair_ids.insert(inner.id);
                        if (sharedEndpointControlPolylinesCross(
                                *outer.curve, candidate, kClusterEndpointTol) ||
                            breaksSharedEndpointControlPolygonNesting(
                                inner.id, candidate, results, result_idx,
                                ordered_neighbors, cluster_solver_, pair_ids))
                            continue;
                        const int old_cross = constrainedCrossCountForId(
                            inner.id, *inner.curve, results, result_idx,
                            ordered_neighbors, kClusterEndpointTol);
                        const int new_cross = constrainedCrossCountForId(
                            inner.id, candidate, results, result_idx,
                            ordered_neighbors, kClusterEndpointTol);
                        const int old_fixed = constrainedFixedCrossCountForId(
                            inner.id, *inner.curve, results, result_idx,
                            ordered_neighbors, kClusterEndpointTol);
                        const int new_fixed = constrainedFixedCrossCountForId(
                            inner.id, candidate, results, result_idx,
                            ordered_neighbors, kClusterEndpointTol);
                        if (new_cross > old_cross || new_fixed > old_fixed)
                            continue;
                        setConnectivityCurveGeometry(inner, candidate);
                        validate(inner, input, sdf);
                        result_idx = resultIndexById(results);
                        ordered_neighbors = constrainedNeighborMap(
                            results, cluster_solver_);
                        return true;
                    }
                    return false;
                };

                const auto bounds = ordinarySingleCubicHandleBounds(
                    outer_entry.first, outer_entry.second,
                    exit_frame.first, exit_frame.second, true);
                const double target_handle = std::min(
                    bounds.end_max, inner_handle + kHandleOrderMargin);
                // 外侧尾把手可能已经达到合法上限（41|43 即属于此类）。
                // 此时不能放弃排序约束，应对称地压短内侧尾把手，直到
                // outer >= inner + margin；候选仍需通过完整形态/物理/不交审计。
                if (target_handle <= outer_handle + 1e-6 ||
                    target_handle < bounds.end_min - 1e-6) {
                    try_compress_inner_tail();
                    continue;
                }

                BezierCurve candidate = *outer.curve;
                candidate.segs.front().ctrl[2] =
                    exit_frame.first - t1 * target_handle;
                if (!ordinarySingleCubicControlsValid(
                        candidate, outer_entry.first, outer_entry.second,
                        exit_frame.first, exit_frame.second, 1e-5, true) ||
                    curveSelfIntersectsBusiness(candidate, 1.0) ||
                    !isNonUTurnTurnShapeAcceptable(
                        candidate, chord_len,
                        std::abs(cross2d(
                            outer_entry.second.normalized(),
                            outer_chord.normalized()))) ||
                    assessCurveRisk(candidate, input, sdf, {}, true).physical() ||
                    curvesIntersectBeyondSharedEndpointOverlap(
                        candidate, *inner.curve, kClusterEndpointTol)) {
                    try_compress_inner_tail();
                    continue;
                }
                // 与内侧的控制多边形嵌套前兆，以及不得在同族其它成员上引入
                // 新的中段互穿；命中时退回压短内侧尾把手的对称修复。
                std::unordered_set<ConnId> outer_pair_ids;
                outer_pair_ids.insert(outer.id);
                outer_pair_ids.insert(inner.id);
                if (sharedEndpointControlPolylinesCross(
                        candidate, *inner.curve, kClusterEndpointTol) ||
                    breaksSharedEndpointControlPolygonNesting(
                        outer.id, candidate, results, result_idx,
                        ordered_neighbors, cluster_solver_, outer_pair_ids) ||
                    introducesNewConstrainedCrossForId(
                        outer.id, *outer.curve, candidate, results, result_idx,
                        ordered_neighbors, outer_pair_ids,
                        kClusterEndpointTol)) {
                    try_compress_inner_tail();
                    continue;
                }

                const int current_cross = constrainedCrossCountForId(
                    outer.id, *outer.curve, results, result_idx,
                    ordered_neighbors, kClusterEndpointTol);
                const int candidate_cross = constrainedCrossCountForId(
                    outer.id, candidate, results, result_idx,
                    ordered_neighbors, kClusterEndpointTol);
                const int current_fixed = constrainedFixedCrossCountForId(
                    outer.id, *outer.curve, results, result_idx,
                    ordered_neighbors, kClusterEndpointTol);
                const int candidate_fixed = constrainedFixedCrossCountForId(
                    outer.id, candidate, results, result_idx,
                    ordered_neighbors, kClusterEndpointTol);
                if (candidate_cross > current_cross ||
                    candidate_fixed > current_fixed)
                    continue;
                setConnectivityCurveGeometry(outer, candidate);
                validate(outer, input, sdf);
                result_idx = resultIndexById(results);
                ordered_neighbors = constrainedNeighborMap(results, cluster_solver_);
            }
        }
        auto annotate_t0 = std::chrono::steady_clock::now();
        annotateClusterCrossings(
            results, cluster_solver_, kClusterEndpointTol);
        if (isgProfile()) {
            fprintf(stderr, "[ISG_PROFILE] pure annotate stage %.3f ms\n",
                    std::chrono::duration<double, std::milli>(
                        std::chrono::steady_clock::now() - annotate_t0).count());
            fprintf(stderr, "[ISG_PROFILE] pure topology total %.3f ms\n",
                    std::chrono::duration<double, std::milli>(
                        std::chrono::steady_clock::now() - topology_t0).count());
        }
        if (out_ms) {
            *out_ms = std::chrono::duration<double, std::milli>(
                std::chrono::steady_clock::now() - t0).count();
        }
        // 无障碍快速路径也必须执行下方的同入 U-turn 家族原子收口；
        // 该收口完成后会重新标注最终同簇状态。
    }

    if (!topology_repair) {
    auto result_idx = resultIndexById(results);
    std::unordered_map<ConnId, BezierCurve> all_done = curveMapFromResults(results);
    // 非纯几何场景只对真实违规曲线做小预算收口。初始生成器已同时考虑
    // 物理风险和同簇约束，后续全量多轮重生成会把单路口耗时放大到分钟级。
    for (int repair_pass = 0;
         repair_pass < repair_budget.regeneration_passes; ++repair_pass) {
        auto pass_t0 = std::chrono::steady_clock::now();
        result_idx = resultIndexById(results);
        all_done = curveMapFromResults(results);
        auto neighbors = constrainedNeighborMap(results, cluster_solver_);
        auto bad = crossingIdsTouchingSeeds(
            results, cluster_solver_, physical_risk_ids, kClusterEndpointTol);
        auto cluster_bad = allConstrainedCrossingIds(
            results, cluster_solver_, preserved_fixed_ids, kClusterEndpointTol);
        bad.insert(cluster_bad.begin(), cluster_bad.end());
        if (bad.empty())
            break;
        int repaired_count = 0;
        int max_repairs = std::min(
            repair_budget.regenerations_per_pass,
            (int)input.connectivities.size());
        for (auto git = plan.batches.rbegin(); git != plan.batches.rend(); ++git) {
            for (auto& cid : git->conn_ids) {
                if (!bad.count(cid))
                    continue;
                if (preserved_fixed_ids.count(cid))
                    continue;
                if (repaired_count >= max_repairs)
                    break;
                const Connectivity* conn = scene.view.connectivity(cid);
                if (!conn)
                    continue;
                if (isGeometricUTurnConn(*conn, input))
                    continue;
                auto ri = result_idx.find(cid);
                if (ri == result_idx.end())
                    continue;
                // 障碍物场景中纯几何拓扑搜索会枚举大量候选；这里直接
                // 用最新兄弟曲线重生成，保留障碍/边界/形态联合约束。
                auto sibs = buildSiblings(cid, all_done, cluster_solver_, input.connectivities, false, &preserved_fixed_ids);
                bool phys_risk = false;
                auto conn_t0 = std::chrono::steady_clock::now();
                const CurveGenerationContext generation_context =
                    CurveGenerationContextBuilder().build(
                        scene, *conn, direction_cfg_.uturn_alignment_scope,
                        &cluster_solver_, &familySnapshot(conn->id));
                auto cc = generateOne(generation_context, sibs, &phys_risk);
                if (isgProfile()) {
                    double ms = std::chrono::duration<double, std::milli>(
                        std::chrono::steady_clock::now() - conn_t0).count();
                    fprintf(stderr, "[ISG_PROFILE] repair conn %s regenerated in %.3f ms\n",
                            cid.c_str(), ms);
                }
                if (hasFixedGeometry(*conn) && !preserved_fixed_ids.count(cid))
                    cc.fixed_shape = false;
                results[ri->second] = std::move(cc);
                if (results[ri->second].curve)
                    all_done[cid] = *results[ri->second].curve;
                if (phys_risk)
                    physical_risk_ids.insert(cid);
                else
                    physical_risk_ids.erase(cid);
                ++repaired_count;
            }
            if (repaired_count >= max_repairs)
                break;
        }
        cluster_solver_.checkAndMarkA2(all_done, sdf, 1.5);
        if (isgProfile()) {
            fprintf(stderr, "[ISG_PROFILE] repair pass %d stage %.3f ms\n",
                    repair_pass,
                    std::chrono::duration<double, std::milli>(
                        std::chrono::steady_clock::now() - pass_t0).count());
        }
    }

    repairMode2UTurnPairs();
    repairStraightUTurnCrossings(repair_budget.final_straight_uturn_repairs);

    auto straight_stage_t0 = std::chrono::steady_clock::now();
    result_idx = resultIndexById(results);
    all_done = curveMapFromResults(results);
    if (input.obstacles.empty()) {
    for (const auto& group : plan.batches) {
        for (const auto& cid : group.conn_ids) {
            auto ri = result_idx.find(cid);
            if (ri == result_idx.end() || preserved_fixed_ids.count(cid))
                continue;
            const Connectivity* conn = scene.view.connectivity(cid);
            if (!conn)
                continue;
            if (isGeometricUTurnConn(*conn, input))
                continue;
            auto& cc = results[ri->second];
            if (!cc.curve)
                continue;
            auto entry = scene.view.entryFrame(conn->entry_lane_id);
            auto exit_ = scene.view.exitFrame(conn->exit_lane_id);
            Vec2d chord = exit_.first - entry.first;
            double chord_len = chord.norm();
            if (chord_len < 1e-8)
                continue;
            Vec2d T0 = entry.second.norm() > 1e-8 ? entry.second.normalized() : chord.normalized();
            Vec2d T1 = exit_.second.norm() > 1e-8 ? exit_.second.normalized() : chord.normalized();
            double turn_strength = std::abs(cross2d(T0, chord.normalized()));
            if (turn_strength >= 0.25)
                continue;
            double arc_chord = cc.curve->arcLength() / chord_len;
            double current_maxk = cc.curve->maxCurvature(40);
            auto sibs = buildSiblings(
                cid, all_done, cluster_solver_, input.connectivities,
                false, &preserved_fixed_ids);
            const auto sampled = sampleSiblingsForIntersections(sibs);
            auto neighbors_for_shift = constrainedNeighborMap(results, cluster_solver_);
            int current_straight_cross = constrainedCrossCountForId(
                cid, *cc.curve, results, result_idx, neighbors_for_shift,
                kClusterEndpointTol);
            auto constrained_uturn_cross_count = [&](const BezierCurve& curve) {
                int count = 0;
                auto nit = neighbors_for_shift.find(cid);
                if (nit == neighbors_for_shift.end())
                    return count;
                for (const auto& other_id : nit->second) {
                    auto oi = result_idx.find(other_id);
                    if (oi == result_idx.end() || !results[oi->second].curve)
                        continue;
                    const Connectivity* other_conn = scene.view.connectivity(other_id);
                    if (!other_conn || !isGeometricUTurnConn(*other_conn, input))
                        continue;
                    if (curvesHaveForbiddenSameClusterIntersection(
                            curve, *results[oi->second].curve, kClusterEndpointTol))
                        ++count;
                }
                return count;
            };
            int current_uturn_cross =
                constrained_uturn_cross_count(*cc.curve);
            if (current_straight_cross > 0) {
                Vec2d ref_perp = dominantConstrainedRefPerpForId(
                    cid, *cc.curve, results, result_idx, neighbors_for_shift,
                    cluster_solver_, kClusterEndpointTol);
                if (ref_perp.norm() > 1e-8) {
                    ref_perp.normalize();
                    Vec2d chord_dir = chord.normalized();
                    double shared_start_lead = 0.0;
                    double shared_end_lead = 0.0;
                    auto nit_for_leads = neighbors_for_shift.find(cid);
                    if (nit_for_leads != neighbors_for_shift.end()) {
                        for (const auto& other_id : nit_for_leads->second) {
                            auto oi = result_idx.find(other_id);
                            if (oi == result_idx.end() || !results[oi->second].curve)
                                continue;
                            const Connectivity* other_conn = scene.view.connectivity(other_id);
                            if (!other_conn ||
                                !isGeometricUTurnConn(*other_conn, input))
                                continue;
                            const BezierCurve& other = *results[oi->second].curve;
                            if (other.segs.empty())
                                continue;
                            if ((other.startPt() - entry.first).norm() <= 0.30 &&
                                other.startTan().norm() > 1e-8 &&
                                other.startTan().normalized().dot(T0) >= 0.98 &&
                                segmentLooksStraight(other.segs.front())) {
                                shared_start_lead = std::max(
                                    shared_start_lead,
                                    (other.segs.front().ctrl[3] -
                                     other.segs.front().ctrl[0]).norm());
                            }
                            if ((other.endPt() - exit_.first).norm() <= 0.30 &&
                                other.endTan().norm() > 1e-8 &&
                                other.endTan().normalized().dot(T1) >= 0.98 &&
                                segmentLooksStraight(other.segs.back())) {
                                shared_end_lead = std::max(
                                    shared_end_lead,
                                    (other.segs.back().ctrl[3] -
                                     other.segs.back().ctrl[0]).norm());
                            }
                        }
                    }
                    shared_start_lead = std::min(shared_start_lead, chord_len * 0.42);
                    shared_end_lead = std::min(shared_end_lead, chord_len * 0.42);
                    auto straight_avoidance_shape_ok = [&](const BezierCurve& candidate) {
                        if (current_uturn_cross <= 0)
                            return isStraightLikeShapeAcceptable(candidate, chord_len);
                        if (candidate.empty())
                            return false;
                        double ratio = candidate.arcLength() / chord_len;
                        if (ratio > 1.35)
                            return false;
                        if (candidate.maxCurvature(40) > 1.20)
                            return false;
                        Vec2d st = candidate.startTan().norm() > 1e-8
                            ? candidate.startTan().normalized() : T0;
                        Vec2d et = candidate.endTan().norm() > 1e-8
                            ? candidate.endTan().normalized() : T1;
                        return st.dot(T0) >= 0.98 && et.dot(T1) >= 0.98;
                    };
                    auto add_candidate = [&](std::vector<BezierCurve>& out,
                                             const std::vector<Vec2d>& pts,
                                             const std::vector<Vec2d>& tans,
                                             double alpha) {
                        if (pts.size() != tans.size() || pts.size() < 2)
                            return;
                        BezierCurve c = makeCurveFromKnots(pts, tans, alpha);
                        if (!c.empty())
                            out.push_back(std::move(c));
                    };
                    BezierCurve best_shift;
                    bool have_shift = false;
                    int best_cross = current_straight_cross;
                    int best_uturn_cross = current_uturn_cross;
                    double best_shift_score = std::numeric_limits<double>::infinity();
                    bool debug_straight_repair =
                        std::getenv("ISG_DEBUG_STRAIGHT_REPAIR") != nullptr;
                    int dbg_total = 0, dbg_self = 0, dbg_shape = 0, dbg_phys = 0;
                    int dbg_worse = 0, dbg_no_reduce = 0;
                    int dbg_min_cross = current_straight_cross;
                    int dbg_min_uturn = current_uturn_cross;
                    for (double side : {-1.0, 1.0}) {
                        for (double offset : {0.45, 0.70, 1.00, 1.35, 1.80, 2.40,
                                              3.00, 4.25, 5.75, 7.50, 9.50,
                                              12.00, 15.00}) {
                            Vec2d off = ref_perp * (side * offset);
                            double lead = std::max(4.0, std::min(12.0, chord_len * 0.24));
                            std::vector<BezierCurve> candidates;
                            candidates.reserve(24);
                            for (double alpha : {0.12, 0.16, 0.20, 0.26, 0.34}) {
                                Vec2d q0 = entry.first + T0 * lead + off;
                                Vec2d q1 = exit_.first - T1 * lead + off;
                                add_candidate(candidates,
                                              {entry.first, q0, q1, exit_.first},
                                              {T0, chord_dir, chord_dir, T1},
                                              alpha);
                                Vec2d r0 = entry.first + chord * 0.18 + off * 0.65;
                                Vec2d rm = entry.first + chord * 0.50 + off;
                                Vec2d r1 = entry.first + chord * 0.82 + off * 0.65;
                                add_candidate(candidates,
                                              {entry.first, r0, rm, r1, exit_.first},
                                              {T0, chord_dir, chord_dir, chord_dir, T1},
                                              alpha);
                                Vec2d mid = entry.first + chord * 0.50 + off;
                                add_candidate(candidates,
                                              {entry.first, mid, exit_.first},
                                              {T0, chord_dir, T1},
                                              alpha);
                                if (shared_start_lead > 0.5 ||
                                    shared_end_lead > 0.5) {
                                    Vec2d s0 = entry.first + T0 * shared_start_lead;
                                    Vec2d s1 = exit_.first - T1 * shared_end_lead;
                                    if (shared_start_lead <= 0.5)
                                        s0 = entry.first;
                                    if (shared_end_lead <= 0.5)
                                        s1 = exit_.first;
                                    if ((s1 - s0).norm() > 2.0) {
                                        Vec2d bridge = s1 - s0;
                                        Vec2d bridge_dir = bridge.norm() > 1e-8
                                            ? bridge.normalized() : chord_dir;
                                        Vec2d smid = 0.5 * (s0 + s1) + off;
                                        add_candidate(candidates,
                                                      {entry.first, s0, smid, s1,
                                                       exit_.first},
                                                      {T0, T0, bridge_dir, T1, T1},
                                                      alpha);
                                        Vec2d s0b = s0 + bridge * 0.25 + off * 0.70;
                                        Vec2d s1b = s0 + bridge * 0.75 + off * 0.70;
                                        add_candidate(candidates,
                                                      {entry.first, s0, s0b, s1b,
                                                       s1, exit_.first},
                                                      {T0, T0, bridge_dir, bridge_dir,
                                                       T1, T1},
                                                      alpha);
                                    }
                                }
                            }

                            for (const auto& candidate : candidates) {
                                ++dbg_total;
                                if (candidate.empty() ||
                                    curveSelfIntersectsBusiness(candidate, 1.0)) {
                                    ++dbg_self;
                                    continue;
                                }
                                if (!straight_avoidance_shape_ok(candidate)) {
                                    ++dbg_shape;
                                    continue;
                                }
                                if (assessCurveRisk(candidate, input, sdf, {}, true).physical()) {
                                    ++dbg_phys;
                                    continue;
                                }
                                int cross_count = constrainedCrossCountForId(
                                    cid, candidate, results, result_idx,
                                    neighbors_for_shift, kClusterEndpointTol);
                                dbg_min_cross = std::min(dbg_min_cross, cross_count);
                                if (cross_count > current_straight_cross)
                                {
                                    ++dbg_worse;
                                    continue;
                                }
                                int uturn_cross_count =
                                    constrained_uturn_cross_count(candidate);
                                dbg_min_uturn = std::min(dbg_min_uturn, uturn_cross_count);
                                if (uturn_cross_count > current_uturn_cross)
                                {
                                    ++dbg_worse;
                                    continue;
                                }
                                if (cross_count >= current_straight_cross)
                                    ++dbg_no_reduce;
                                double score = 1000.0 * cross_count +
                                    500.0 * uturn_cross_count + offset +
                                    candidate.maxCurvature(40) +
                                    0.02 * candidate.arcLength();
                                if (!have_shift ||
                                    cross_count < best_cross ||
                                    (cross_count == best_cross &&
                                     uturn_cross_count < best_uturn_cross) ||
                                    (cross_count == best_cross &&
                                     uturn_cross_count == best_uturn_cross &&
                                     score < best_shift_score)) {
                                    best_shift = candidate;
                                    best_cross = cross_count;
                                    best_uturn_cross = uturn_cross_count;
                                    best_shift_score = score;
                                    have_shift = true;
                                }
                            }
                            if (have_shift && best_cross == 0)
                                break;
                        }
                        if (have_shift && best_cross == 0)
                            break;
                    }
                    if (have_shift && best_cross < current_straight_cross) {
                        setConnectivityCurveGeometry(cc, best_shift);
                        validate(cc, input, sdf);
                        all_done[cid] = best_shift;
                        continue;
                    }
                    if (debug_straight_repair) {
                        fprintf(stderr,
                                "[STRAIGHT-REPAIR] %s current=%d uturn=%d best=%d best_uturn=%d total=%d self=%d shape=%d phys=%d worse=%d no_reduce=%d min=%d min_uturn=%d accepted=%d\n",
                                cid.c_str(), current_straight_cross,
                                current_uturn_cross, best_cross, best_uturn_cross,
                                dbg_total, dbg_self, dbg_shape, dbg_phys,
                                dbg_worse, dbg_no_reduce, dbg_min_cross,
                                dbg_min_uturn,
                                have_shift && best_cross < current_straight_cross);
                    }
                }
            }
            if (arc_chord > 1.08 || current_maxk <= 1.0)
                continue;

            auto pts = cc.curve->sampleByArcLength(6);
            if (pts.size() < 2)
                continue;
            pts.front() = entry.first;
            pts.back() = exit_.first;
            BezierCurve smoothed = fitBezierWithEndTangents(pts, T0, T1);
            if (smoothed.empty() || smoothed.maxCurvature(40) >= current_maxk ||
                smoothed.maxCurvature(40) > 1.0)
                continue;
            if (assessCurveRisk(smoothed, input, sdf, {}, true).physical())
                continue;
            if (sampledSiblingCrossCount(
                    smoothed, sampled, true, kClusterEndpointTol) > 0)
                continue;
            setConnectivityCurveGeometry(cc, smoothed);
            validate(cc, input, sdf);
            all_done[cid] = smoothed;
        }
    }
    }
    if (isgProfile()) {
        fprintf(stderr, "[ISG_PROFILE] straight smoothing stage %.3f ms\n",
                std::chrono::duration<double, std::milli>(
                    std::chrono::steady_clock::now() - straight_stage_t0).count());
    }

    result_idx = resultIndexById(results);
    for (const auto& cc : results) {
        if (!cc.curve || preserved_fixed_ids.count(cc.id))
            continue;
        const Connectivity* conn = scene.view.connectivity(cc.id);
        if (!conn || isGeometricUTurnConn(*conn, input))
            continue;
        auto entry = scene.view.entryFrame(conn->entry_lane_id);
        auto exit_ = scene.view.exitFrame(conn->exit_lane_id);
        BezierCurve natural_curve;
        Vec2d natural_chord = exit_.first - entry.first;
        if (natural_chord.norm() > 1e-8) {
            Vec2d natural_t0 = entry.second.norm() > 1e-8
                ? entry.second.normalized() : natural_chord.normalized();
            Vec2d natural_t1 = exit_.second.norm() > 1e-8
                ? exit_.second.normalized() : natural_chord.normalized();
            natural_curve.segs.push_back(makeCubicG1(
                entry.first, natural_t0, exit_.first, natural_t1, 0.40));
        }
        CurveRisk risk = assessCurveRisk(*cc.curve, input, sdf, {}, true);
        bool roadedge_candidate =
            risk.boundary ||
            curveRawIntersectsRoadEdge(*cc.curve, input.boundaries) ||
            (!natural_curve.empty() &&
             (curveRawIntersectsRoadEdge(natural_curve, input.boundaries) ||
              segmentRawIntersectsRoadEdge(entry.first, exit_.first, input.boundaries) ||
              segmentCorridorTouchesRoadEdge(entry.first, exit_.first, input.boundaries)));
        if (roadedge_candidate)
            roadedge_boundary_repair_ids.insert(cc.id);
    }

    if (input.obstacles.empty() && !roadedge_boundary_repair_ids.empty()) {
        auto final_boundary_t0 = std::chrono::steady_clock::now();
        result_idx = resultIndexById(results);
        all_done = curveMapFromResults(results);
        int final_boundary_seen = 0;
        for (const auto& group : plan.batches) {
            for (const auto& cid : group.conn_ids) {
                ++final_boundary_seen;
                auto ri = result_idx.find(cid);
                if (ri == result_idx.end() || preserved_fixed_ids.count(cid))
                    continue;
                if (!roadedge_boundary_repair_ids.count(cid))
                    continue;
                const Connectivity* conn = scene.view.connectivity(cid);
                if (!conn || isGeometricUTurnConn(*conn, input))
                    continue;
                auto& cc = results[ri->second];
                if (!cc.curve)
                    continue;
                auto entry = scene.view.entryFrame(conn->entry_lane_id);
                auto exit_ = scene.view.exitFrame(conn->exit_lane_id);
                Vec2d natural_chord = exit_.first - entry.first;
                if (natural_chord.norm() < 1e-8)
                    continue;
                Vec2d natural_t0 = entry.second.norm() > 1e-8
                    ? entry.second.normalized() : natural_chord.normalized();
                Vec2d natural_t1 = exit_.second.norm() > 1e-8
                    ? exit_.second.normalized() : natural_chord.normalized();
                BezierCurve natural_curve;
                natural_curve.segs.push_back(makeCubicG1(
                    entry.first, natural_t0, exit_.first, natural_t1, 0.40));
                auto sibs = buildSiblings(
                    cid, all_done, cluster_solver_, input.connectivities,
                    true, &preserved_fixed_ids);
                const auto sampled = sampleSiblingsForIntersections(sibs);
                BezierCurve repaired = natural_curve;
                if (!tryBoundarySafeCandidate(
                        entry.first, entry.second, exit_.first, exit_.second,
                        input, sdf, sampled, true, repaired, true))
                    continue;
                double turn_strength = std::abs(
                    cross2d(natural_t0, natural_chord.normalized()));
                bool repaired_shape_ok = turn_strength < 0.25
                    ? isStraightLikeShapeAcceptable(repaired, natural_chord.norm())
                    : (isNonUTurnTurnShapeAcceptable(
                           repaired, natural_chord.norm(), turn_strength) &&
                       !curveHasCurvatureSignFlip(repaired));
                CurveRisk repaired_risk = assessCurveRisk(
                    repaired, input, sdf, sampled, true);
                if (repaired_risk.physical())
                    continue;
                int current_fixed_raw_cross =
                    rawFixedSiblingCrossCount(*cc.curve, sibs, kClusterEndpointTol);
                int repaired_fixed_raw_cross =
                    rawFixedSiblingCrossCount(repaired, sibs, kClusterEndpointTol);
                if (repaired_fixed_raw_cross > current_fixed_raw_cross)
                    continue;
                auto strictSameSideCrossCount = [&](const BezierCurve& curve) {
                    int count = 0;
                    double cur_signed = signedTurnStrengthOfConnId(scene.view, cid);
                    for (const auto& sib : sibs) {
                        if (sib.exempt_a1 || !sib.shared_endpoint)
                            continue;
                        double sib_signed = signedTurnStrengthOfConnId(scene.view, sib.id);
                        if (std::abs(cur_signed) < 0.25 ||
                            std::abs(sib_signed) < 0.25 ||
                            cur_signed * sib_signed <= 0.0)
                            continue;
                        if (curvesIntersectBeyondSharedEndpointOverlap(
                                curve, sib.curve, 0.15))
                            ++count;
                    }
                    return count;
                };
                if (strictSameSideCrossCount(repaired) >
                    strictSameSideCrossCount(*cc.curve))
                    continue;
                if (!repaired_shape_ok)
                    continue;
                setConnectivityCurveGeometry(cc, repaired);
                validate(cc, input, sdf);
                all_done[cid] = repaired;
            }
        }
        if (isgProfile()) {
            fprintf(stderr, "[ISG_PROFILE] final boundary stage %.3f ms\n",
                    std::chrono::duration<double, std::milli>(
                        std::chrono::steady_clock::now() - final_boundary_t0).count());
        }
    }

    result_idx = resultIndexById(results);
    all_done = curveMapFromResults(results);
    for (const auto& group : plan.batches) {
        for (const auto& cid : group.conn_ids) {
            auto ri = result_idx.find(cid);
            if (ri == result_idx.end() || preserved_fixed_ids.count(cid))
                continue;
            const Connectivity* conn = scene.view.connectivity(cid);
            if (!conn || isGeometricUTurnConn(*conn, input))
                continue;
            auto& cc = results[ri->second];
            if (!cc.curve || cc.curve->numSegments() <= 1)
                continue;
            auto entry = scene.view.entryFrame(conn->entry_lane_id);
            auto exit_ = scene.view.exitFrame(conn->exit_lane_id);
            auto sibs = buildSiblings(
                cid, all_done, cluster_solver_, input.connectivities,
                true, &preserved_fixed_ids);
            auto sampled = sampleSiblingsForIntersections(sibs);
            BezierCurve restored = *cc.curve;
            if (!tryNaturalStraightSingleCubicIfSafe(
                    entry.first, entry.second, exit_.first, exit_.second,
                    input, sdf, sampled, true, restored))
                continue;
            setConnectivityCurveGeometry(cc, restored);
            validate(cc, input, sdf);
            all_done[cid] = restored;
        }
    }
    }

    // Boundary-only 路口在上面的 topology_repair 分支中不会进入障碍物场景的
    // 普通曲线收口。这里补一次最终严格边界收口，确保密集边界链上的普通转向
    // 不会以“纯拓扑已完成”为由保留穿越；候选仍由 tryBoundarySafeCandidate
    // 同时复核形态、围栏、障碍和同簇约束。
    if (topology_repair && input.obstacles.empty() &&
        !roadedge_boundary_repair_ids.empty()) {
        auto boundary_result_idx = resultIndexById(results);
        auto boundary_done = curveMapFromResults(results);
        for (const auto& group : plan.batches) {
            for (const auto& cid : group.conn_ids) {
                if (!roadedge_boundary_repair_ids.count(cid) ||
                    preserved_fixed_ids.count(cid))
                    continue;
                auto ri = boundary_result_idx.find(cid);
                if (ri == boundary_result_idx.end() || !results[ri->second].curve)
                    continue;
                const Connectivity* conn = scene.view.connectivity(cid);
                if (!conn || isGeometricUTurnConn(*conn, input))
                    continue;
                auto entry = scene.view.entryFrame(conn->entry_lane_id);
                auto exit_ = scene.view.exitFrame(conn->exit_lane_id);
                const Vec2d chord = exit_.first - entry.first;
                if (chord.norm() < 1e-8)
                    continue;
                auto siblings = buildSiblings(
                    cid, boundary_done, cluster_solver_, input.connectivities,
                    true, &preserved_fixed_ids);
                auto sampled = sampleSiblingsForIntersections(siblings);
                BezierCurve repaired = *results[ri->second].curve;
                if (!tryBoundarySafeCandidate(
                        entry.first, entry.second, exit_.first, exit_.second,
                        input, sdf, sampled, true, repaired, true))
                    continue;
                const CurveRisk risk = assessCurveRisk(
                    repaired, input, sdf, sampled, true);
                if (risk.physical())
                    continue;
                setConnectivityCurveGeometry(results[ri->second], repaired);
                validate(results[ri->second], input, sdf);
                boundary_done[cid] = repaired;
                boundary_result_idx = resultIndexById(results);
            }
        }
    }

    auto uturn_family_stage_t0 = std::chrono::steady_clock::now();
    // 同入/同出的传递 U-turn 家族必须作为一个整体保持半径顺序。Boundary
    // 绕行会要求不同成员沿不同端点使用不同偏置，公共偏置无法表达这种几何，
    // 还会把某个成员推入同族邻居。因此每个成员先建立有限候选池，再用家族级
    // 回溯联合验收；候选全部通过后只做一次原子提交。
    {
        struct FamilyCandidate {
            BezierCurve curve;
            double alpha = 0.0;
            double entry_bias = 0.0;
            double exit_bias = 0.0;
            double station_offset = 0.0;
            double score = 0.0;
        };
        struct ExternalCandidate {
            BezierCurve curve;
            double score = 0.0;
        };

        std::vector<std::vector<const Connectivity*> > families;
        std::unordered_set<ConnId> assigned;
        const UTurnFamilyBuilder family_builder;
        for (const auto& conn : input.connectivities) {
            if (!isGeometricUTurnConn(conn, input) || assigned.count(conn.id))
                continue;
            std::vector<const Connectivity*> component =
                family_builder.alignmentComponent(
                    conn, input, direction_cfg_.uturn_alignment_scope,
                    &cluster_solver_, false);
            for (const Connectivity* member : component)
                assigned.insert(member->id);
            if (component.size() >= 2)
                families.push_back(std::move(component));
        }

        const std::vector<double> family_arc_alphas =
            {2.0 / 3.0, 0.85, 1.0};
        const std::vector<double> family_bias_values = input.boundaries.empty()
            ? std::vector<double>{0.0}
            : std::vector<double>{-2.0, -1.0, -0.70, -0.35, 0.0,
                                  0.35, 0.70, 2.0};
        const std::vector<double> family_station_offsets = input.boundaries.empty()
            ? std::vector<double>{0.0}
            : std::vector<double>{0.0, 1.50, 3.50};

        for (auto& family : families) {
            std::stable_sort(
                family.begin(), family.end(),
                [&](const Connectivity* a, const Connectivity* b) {
                    const double ra = uturnRadiusKey(*a, input);
                    const double rb = uturnRadiusKey(*b, input);
                    if (std::abs(ra - rb) > 1e-9)
                        return ra < rb;
                    return a->id < b->id;
                });

            std::unordered_set<ConnId> family_ids;
            std::unordered_map<ConnId, std::size_t> family_index;
            for (std::size_t i = 0; i < family.size(); ++i) {
                family_ids.insert(family[i]->id);
                family_index[family[i]->id] = i;
            }
            const auto baseline_index = resultIndexById(results);

            std::vector<std::vector<FamilyCandidate> > pools(family.size());
            std::vector<std::size_t> pool_trials(family.size(), 0);
            std::vector<std::size_t> pool_rejected(family.size(), 0);
            std::vector<std::size_t> pool_boundary_rejected(family.size(), 0);
            std::vector<std::size_t> pool_obstacle_rejected(family.size(), 0);
            std::vector<std::size_t> pool_fence_rejected(family.size(), 0);
            for (std::size_t fi = 0; fi < family.size(); ++fi) {
                const Connectivity* conn = family[fi];
                const UTurnFamilyInfo& family_info = familySnapshot(conn->id);
                const auto entry = scene.view.entryFrame(conn->entry_lane_id);
                const auto exit_ = scene.view.exitFrame(conn->exit_lane_id);
                const double chord_len = (exit_.first - entry.first).norm();
                const double max_curvature = segmentedUTurnMaxCurvatureLimit(
                    entry.first, entry.second, exit_.first, exit_.second);
                if (chord_len < 1e-6)
                    continue;

                for (double station_offset : family_station_offsets) {
                    const double target_station =
                        std::isfinite(family_info.aligned_station)
                            ? family_info.aligned_station + station_offset
                            : family_info.aligned_station;
                    for (double entry_bias : family_bias_values) {
                        for (double exit_bias : family_bias_values) {
                            for (double alpha : family_arc_alphas) {
                                ++pool_trials[fi];
                                double min_lead0 = family_info.lead0;
                                double min_lead1 = family_info.lead1;
                                family_builder.enforceAlignmentStation(
                                    entry.first, entry.second, exit_.first,
                                    exit_.second, target_station,
                                    min_lead0, min_lead1);
                                double entry_stagger = 0.0;
                                double exit_stagger = 0.0;
                                double family_stagger_step = 0.0;
                                const double stagger = uturnSharedEndpointStagger(
                                    *conn, input, cluster_solver_,
                                    direction_cfg_.uturn_alignment_scope,
                                    &family_stagger_step, &entry_stagger,
                                    &exit_stagger);
                                // family_stagger_step=0 deliberately disables the
                                // per-member reverse-bias clamp. The baseline ladder
                                // remains in entry/exit_stagger; the extra bias is now
                                // chosen jointly with all family members.
                                BezierCurve candidate =
                                    UTurnCurveInitializer().buildSegmented(
                                        entry.first, entry.second, exit_.first,
                                        exit_.second, min_lead0, min_lead1,
                                        alpha, stagger, 0.0, 0.0,
                                        entry_stagger, exit_stagger,
                                        entry_bias, exit_bias,
                                        0.0);
                                if (candidate.numSegments() != 3 ||
                                    !segmentedUTurnHasMinimumStraightLeads(
                                        candidate, min_lead0, min_lead1) ||
                                    !segmentLooksStraight(candidate.segs.front()) ||
                                    !segmentLooksStraight(candidate.segs.back()) ||
                                    candidate.maxCurvature(40) >= max_curvature ||
                                    std::abs((candidate.segs.front().ctrl[3] -
                                              candidate.segs.back().ctrl[0])
                                                 .dot(family_info.axis)) > 0.05 ||
                                    (std::isfinite(target_station) &&
                                     (std::abs(candidate.segs.front().ctrl[3]
                                                  .dot(family_info.axis) -
                                              target_station) > 0.05 ||
                                      std::abs(candidate.segs.back().ctrl[0]
                                                  .dot(family_info.axis) -
                                              target_station) > 0.05)) ||
                                    candidate.segs[1].maxCurvature(30) <= 0.03 ||
                                    candidate.arcLength() / chord_len < 1.35 ||
                                    curveSelfIntersectsBusiness(candidate, 1.0) ||
                                    !segmentedUTurnMiddleArcClearsCrosswalks(
                                        candidate, family_info.clearance_crosswalks)) {
                                    ++pool_rejected[fi];
                                    continue;
                                }
                                const CurveRisk risk = assessCurveRisk(
                                    candidate, input, sdf, {}, true, true);
                                if (risk.physical()) {
                                    if (risk.boundary)
                                        ++pool_boundary_rejected[fi];
                                    if (risk.obstacle)
                                        ++pool_obstacle_rejected[fi];
                                    if (risk.fence)
                                        ++pool_fence_rejected[fi];
                                    ++pool_rejected[fi];
                                    continue;
                                }
                                FamilyCandidate item;
                                item.curve = std::move(candidate);
                                item.alpha = alpha;
                                item.entry_bias = entry_bias;
                                item.exit_bias = exit_bias;
                                item.station_offset = station_offset;
                                item.score = item.curve.maxCurvature(40) +
                                    0.02 * item.curve.arcLength() +
                                    0.35 * (std::abs(entry_bias) +
                                             std::abs(exit_bias)) +
                                    0.10 * std::abs(station_offset) +
                                    0.01 * std::abs(alpha - 2.0 / 3.0);
                                pools[fi].push_back(std::move(item));
                                if (isgProfile() && pools[fi].size() <= 3) {
                                    fprintf(stderr,
                                            "[ISG_PROFILE] U-turn family %s member %s accepted #%zu alpha=%.3f bias=%.2f/%.2f station=%.2f\n",
                                            family.front()->id.c_str(),
                                            conn->id.c_str(), pools[fi].size(),
                                            alpha, entry_bias, exit_bias,
                                            station_offset);
                                }
                            }
                        }
                    }
                }
                std::stable_sort(
                    pools[fi].begin(), pools[fi].end(),
                    [](const FamilyCandidate& a, const FamilyCandidate& b) {
                        return a.score < b.score;
                    });
                if (isgProfile())
                    fprintf(stderr,
                            "[ISG_PROFILE] U-turn family %s member %s pool=%zu trials=%zu rejected=%zu\n",
                            family.front()->id.c_str(), conn->id.c_str(),
                            pools[fi].size(), pool_trials[fi],
                            pool_rejected[fi]);
                if (isgProfile())
                    fprintf(stderr,
                            " boundary=%zu obstacle=%zu fence=%zu\n",
                            pool_boundary_rejected[fi],
                            pool_obstacle_rejected[fi],
                            pool_fence_rejected[fi]);
            }

            bool empty_pool = false;
            for (const auto& pool : pools) {
                if (pool.empty()) {
                    empty_pool = true;
                    break;
                }
            }
            if (empty_pool)
                continue;

            // Boundary 链场景中，U-turn 中弧可能真实穿过一条已经生成的普通
            // 共享端点兄弟。该交叉不能豁免；把实际阻塞的普通兄弟一并纳入
            // 原子候选搜索，避免只锁定 U-turn 而把严格约束误判为无解。
            std::vector<ConnId> external_ids;
            std::unordered_set<ConnId> external_id_set;
            for (const auto& pair : cluster_solver_.pairs()) {
                const auto ia = family_index.find(pair.id_a);
                const auto ib = family_index.find(pair.id_b);
                const bool a_family = ia != family_index.end();
                const bool b_family = ib != family_index.end();
                if (a_family == b_family ||
                    pair.exempt == CrossExemption::StructuralCross)
                    continue;
                const ConnId& family_id = a_family ? pair.id_a : pair.id_b;
                const ConnId& external_id = a_family ? pair.id_b : pair.id_a;
                const auto baseline_family = baseline_index.find(family_id);
                const auto baseline_external = baseline_index.find(external_id);
                if (baseline_family == baseline_index.end() ||
                    baseline_external == baseline_index.end() ||
                    !results[baseline_family->second].curve ||
                    !results[baseline_external->second].curve)
                    continue;
                const std::vector<FamilyCandidate>& pool = pools[
                    a_family ? ia->second : ib->second];
                bool all_candidates_blocked = !pool.empty();
                for (const FamilyCandidate& candidate : pool) {
                    // 阻塞识别必须与联合验收共用严格口径：非端点贴合/重叠
                    // 也会在 search_family 中被拒绝，不能只检查真正穿越，
                    // 否则 external_ids 为空而普通兄弟永远不会进入原子搜索。
                    if (!curvesHaveForbiddenSameClusterIntersection(
                            candidate.curve,
                            *results[baseline_external->second].curve,
                            kClusterEndpointTol)) {
                        all_candidates_blocked = false;
                        break;
                    }
                }
                if (all_candidates_blocked &&
                    external_id_set.insert(external_id).second)
                    external_ids.push_back(external_id);
            }
            std::sort(external_ids.begin(), external_ids.end());

            std::unordered_map<ConnId, std::vector<ExternalCandidate>>
                external_pools;
            std::unordered_map<ConnId, std::size_t> external_index;
            bool external_pool_empty = false;
            for (std::size_t ei = 0; ei < external_ids.size(); ++ei) {
                const ConnId& external_id = external_ids[ei];
                external_index[external_id] = ei;
                const Connectivity* external_conn =
                    scene.view.connectivity(external_id);
                const auto baseline_external = baseline_index.find(external_id);
                if (!external_conn || baseline_external == baseline_index.end() ||
                    !results[baseline_external->second].curve ||
                    external_conn->fixed_shape ||
                    isGeometricUTurnConn(*external_conn, input)) {
                    external_pool_empty = true;
                    break;
                }
                const auto entry = scene.view.entryFrame(
                    external_conn->entry_lane_id);
                const auto exit_ = scene.view.exitFrame(
                    external_conn->exit_lane_id);
                const Vec2d chord = exit_.first - entry.first;
                if (chord.norm() < 1e-8 || entry.second.norm() < 1e-8 ||
                    exit_.second.norm() < 1e-8) {
                    external_pool_empty = true;
                    break;
                }
                const double turn_strength = std::abs(cross2d(
                    entry.second.normalized(), chord.normalized()));
                const OrdinaryCurveInitializer ordinary_initializer;
                std::vector<BezierCurve> trials;
                std::size_t external_shape_rejected = 0;
                std::size_t external_physical_rejected = 0;
                trials.push_back(*results[baseline_external->second].curve);
                trials.push_back(ordinary_initializer.buildPreferredSingleCubic(
                    entry.first, entry.second, exit_.first, exit_.second));
                for (double alpha : {0.25, 0.35, 0.45, 0.55, 0.65, 0.80})
                    trials.push_back(ordinary_initializer.buildSingleCubic(
                        entry.first, entry.second, exit_.first, exit_.second,
                        alpha));
                std::vector<ExternalCandidate>& pool =
                    external_pools[external_id];
                for (std::size_t ti = 0; ti < trials.size(); ++ti) {
                    const BezierCurve& trial = trials[ti];
                    if (trial.empty() || trial.numSegments() != 1 ||
                        !ordinarySingleCubicControlsValid(
                            trial, entry.first, entry.second, exit_.first,
                            exit_.second, 1e-5, true) ||
                        curveSelfIntersectsBusiness(trial, 1.0)) {
                        ++external_shape_rejected;
                        continue;
                    }
                    const bool shape_valid = turn_strength < 0.25
                        ? isStraightLikeShapeAcceptable(trial, chord.norm())
                        : isNonUTurnTurnShapeAcceptable(
                              trial, chord.norm(), turn_strength);
                    if (!shape_valid) {
                        ++external_shape_rejected;
                        continue;
                    }
                    if (assessCurveRisk(
                            trial, input, sdf, {}, true, false).physical()) {
                        ++external_physical_rejected;
                        continue;
                    }
                    ExternalCandidate candidate;
                    candidate.curve = trial;
                    candidate.score = candidate.curve.maxCurvature(30) +
                        0.02 * candidate.curve.arcLength() +
                        0.01 * static_cast<double>(ti);
                    pool.push_back(std::move(candidate));
                }
                std::stable_sort(
                    pool.begin(), pool.end(),
                    [](const ExternalCandidate& a, const ExternalCandidate& b) {
                        return a.score < b.score;
                    });
                if (isgProfile())
                    fprintf(stderr,
                            "[ISG_PROFILE] U-turn family %s external %s pool=%zu trials=%zu shape_rejected=%zu physical_rejected=%zu\n",
                            family.front()->id.c_str(), external_id.c_str(),
                            pool.size(), trials.size(), external_shape_rejected,
                            external_physical_rejected);
                if (pool.empty())
                    external_pool_empty = true;
            }
            if (external_pool_empty)
                continue;

            // 家族内部共享端点只允许在真实连接点相遇。必须保留两条候选的
            // 完整几何检查：一条曲线的首/尾直段穿入另一条曲线中弧时，不能
            // 通过裁掉共享端段把真实相交隐藏掉；交叉函数自身只豁免共同真实端点。
            const auto family_pair_forbidden =
                [&](const BezierCurve& original_a,
                    const BezierCurve& original_b) {
                    if (original_a.empty() || original_b.empty())
                        return true;
                    return curvesHaveForbiddenSameClusterIntersection(
                        original_a, original_b, kConnectionPointTolerance);
                };

            // 存指针而不是拷贝：候选自带 BezierCurve，回溯每次迭代都整体
            // 赋值会带来数千万次向量堆分配。候选池在搜索期间不再变动，
            // 指针始终有效。
            std::vector<const FamilyCandidate*> selected(family.size(), nullptr);
            std::vector<bool> selected_flags(family.size(), false);
            std::vector<ExternalCandidate> selected_external(external_ids.size());
            std::vector<bool> selected_external_flags(external_ids.size(), false);
            // 回溯每层只有「新选中成员」参与的配对可能新增冲突：更浅层成员、
            // 外部成员与基线曲线的几何在本层没有变化，上一层已经逐对校验通过，
            // 重算必然得到同样结论。因此按成员预先建立入射配对索引，把每个
            // 节点的 O(全部配对) 扫描降到 O(该成员的配对)。
            //
            // 家族成员两两之间的组合由下面 previous<depth 的循环用同一谓词
            // （family_pair_forbidden）完整覆盖，且不依赖 ClusterOrderSolver
            // 是否登记了该对，所以家族内部配对不重复放进索引。
            std::vector<std::vector<const CurvePair*>> family_incident_pairs(
                family.size());
            std::vector<std::vector<const CurvePair*>> external_incident_pairs(
                external_ids.size());
            for (const auto& pair : cluster_solver_.pairs()) {
                if (pair.exempt == CrossExemption::StructuralCross)
                    continue;
                const auto ia = family_index.find(pair.id_a);
                const auto ib = family_index.find(pair.id_b);
                const bool a_family = ia != family_index.end();
                const bool b_family = ib != family_index.end();
                if (a_family && !b_family)
                    family_incident_pairs[ia->second].push_back(&pair);
                else if (b_family && !a_family)
                    family_incident_pairs[ib->second].push_back(&pair);
                if (a_family || b_family)
                    continue;
                const auto ea = external_index.find(pair.id_a);
                const auto eb = external_index.find(pair.id_b);
                if (ea != external_index.end())
                    external_incident_pairs[ea->second].push_back(&pair);
                if (eb != external_index.end())
                    external_incident_pairs[eb->second].push_back(&pair);
            }
            std::unordered_map<std::string, std::size_t> external_rejections;
            std::unordered_map<std::string, std::size_t> internal_rejections;
            // 这两张表只用于 ISG_PROFILE 的拒绝原因统计。回溯内层每次拒绝
            // 都要拼出 "a|b" 字符串再插入哈希表，而 100000643 的单个家族内部
            // 拒绝就超过百万次，实测占该数据集约 16% 的耗时。未开启诊断时
            // 直接跳过计数：计数不参与任何判定，搜索结果完全不变。
            const bool profile_rejections = isgProfile();
            // 与固定基线曲线的冲突只取决于候选自身，与其余成员的取法无关。
            // 原实现把这个判定放在回溯内层，同一候选会随上层组合被重复检查
            // 上万次（100000643 的 73|91 达 9.3 万次）。这里在搜索前一次性
            // 过滤，保持候选原有优先顺序；被删除的候选在任何组合里都会被同
            // 一条判定拒绝，因此最终接受的组合不变。
            for (std::size_t mi = 0; mi < family.size(); ++mi) {
                std::vector<const BezierCurve*> baseline_partners;
                for (const CurvePair* pair_ptr : family_incident_pairs[mi]) {
                    const ConnId& other = pair_ptr->id_a == family[mi]->id
                        ? pair_ptr->id_b : pair_ptr->id_a;
                    if (external_index.count(other))
                        continue;
                    const auto ri = baseline_index.find(other);
                    if (ri != baseline_index.end() && results[ri->second].curve)
                        baseline_partners.push_back(
                            results[ri->second].curve.get());
                }
                if (baseline_partners.empty())
                    continue;
                std::vector<FamilyCandidate> kept;
                kept.reserve(pools[mi].size());
                for (const FamilyCandidate& candidate : pools[mi]) {
                    bool blocked = false;
                    for (const BezierCurve* partner : baseline_partners) {
                        if (curvesHaveForbiddenSameClusterIntersection(
                                candidate.curve, *partner, kClusterEndpointTol)) {
                            blocked = true;
                            break;
                        }
                    }
                    if (!blocked)
                        kept.push_back(candidate);
                }
                pools[mi].swap(kept);
            }
            bool family_pool_empty = false;
            for (const auto& pool : pools) {
                if (pool.empty())
                    family_pool_empty = true;
            }
            if (family_pool_empty) {
                if (isgProfile())
                    fprintf(stderr,
                            "[ISG_PROFILE] U-turn family %s pool emptied by baseline prefilter\n",
                            family.front()->id.c_str());
                continue;
            }
            // 家族成员两两之间的判定同样只取决于两个候选的下标，回溯里会被
            // 反复重算（74|76 达 111 万次，而不同下标组合最多约 10 万种）。
            // 按需填充三态表：0 未计算 / 1 允许 / 2 禁止。
            std::vector<std::vector<std::uint8_t>> internal_memo(
                family.size() * family.size());
            const auto internal_forbidden =
                [&](std::size_t pi, std::size_t pc, std::size_t di,
                    std::size_t dc) {
                    std::vector<std::uint8_t>& table =
                        internal_memo[pi * family.size() + di];
                    if (table.empty())
                        table.assign(pools[pi].size() * pools[di].size(), 0);
                    std::uint8_t& cell = table[pc * pools[di].size() + dc];
                    if (cell == 0)
                        cell = family_pair_forbidden(pools[pi][pc].curve,
                                                     pools[di][dc].curve)
                            ? 2 : 1;
                    return cell == 2;
                };
            std::vector<std::size_t> selected_idx(family.size(), 0);
            std::size_t search_nodes = 0;
            const std::size_t max_search_nodes = 120000;
            std::function<bool(std::size_t)> search_family;
            std::function<bool(std::size_t)> search_external;
            search_family =
                [&](std::size_t depth) {
                    if (++search_nodes > max_search_nodes)
                        return false;
                    if (depth == family.size())
                        return true;
                    for (std::size_t ci = 0; ci < pools[depth].size(); ++ci) {
                        const FamilyCandidate& item = pools[depth][ci];
                        selected[depth] = &item;
                        selected_idx[depth] = ci;
                        selected_flags[depth] = true;
                        bool valid = true;
                        // 家族内部的共享端点组合不能依赖 ClusterOrderSolver
                        // 是否为每一对生成了 CurvePair。LaneGroup 传递家族可能
                        // 只有部分 pair 被登记，但用户约束对所有家族成员都要求
                        // 连接点之外不得相交；漏掉这里会让 18|20、32|34 的
                        // 中弧穿越在原子提交后才暴露。
                        for (std::size_t previous = 0;
                             previous < depth && valid; ++previous) {
                            if (internal_forbidden(previous,
                                                   selected_idx[previous],
                                                   depth, ci)) {
                                if (profile_rejections)
                                    ++internal_rejections[
                                        family[previous]->id + "|" +
                                        family[depth]->id];
                                valid = false;
                            }
                        }
                        if (!valid) {
                            selected_flags[depth] = false;
                            continue;
                        }
                        for (const CurvePair* pair_ptr :
                                 family_incident_pairs[depth]) {
                            const CurvePair& pair = *pair_ptr;
                            if (pair.exempt == CrossExemption::StructuralCross)
                                continue;
                            const auto ia = family_index.find(pair.id_a);
                            const auto ib = family_index.find(pair.id_b);
                            const bool a_family = ia != family_index.end();
                            const bool b_family = ib != family_index.end();
                            const auto ea = external_index.find(pair.id_a);
                            const auto eb = external_index.find(pair.id_b);
                            const bool a_external = ea != external_index.end();
                            const bool b_external = eb != external_index.end();
                            if (!a_family && !b_family)
                                continue;
                            if ((a_external &&
                                 !selected_external_flags[ea->second]) ||
                                (b_external &&
                                 !selected_external_flags[eb->second]))
                                continue;
                            if (a_family && !selected_flags[ia->second])
                                continue;
                            if (b_family && !selected_flags[ib->second])
                                continue;

                            const BezierCurve* curve_a = nullptr;
                            const BezierCurve* curve_b = nullptr;
                            if (a_family)
                                curve_a = &selected[ia->second]->curve;
                            else if (a_external)
                                curve_a = &selected_external[ea->second].curve;
                            else {
                                const auto ri = baseline_index.find(pair.id_a);
                                if (ri != baseline_index.end() &&
                                    results[ri->second].curve)
                                    curve_a = results[ri->second].curve.get();
                            }
                            if (b_family)
                                curve_b = &selected[ib->second]->curve;
                            else if (b_external)
                                curve_b = &selected_external[eb->second].curve;
                            else {
                                const auto ri = baseline_index.find(pair.id_b);
                                if (ri != baseline_index.end() &&
                                    results[ri->second].curve)
                                    curve_b = results[ri->second].curve.get();
                            }
                            if (!curve_a || !curve_b)
                                continue;
                            const bool forbidden = a_family && b_family
                                ? family_pair_forbidden(*curve_a, *curve_b)
                                : curvesHaveForbiddenSameClusterIntersection(
                                      *curve_a, *curve_b,
                                      kClusterEndpointTol);
                            if (forbidden) {
                                if (profile_rejections) {
                                    const std::string pair_key =
                                        pair.id_a + "|" + pair.id_b;
                                    if (a_family && b_family) {
                                        ++internal_rejections[pair_key];
                                    } else if (
                                            ++external_rejections[pair_key]
                                                <= 2) {
                                        const std::vector<Vec2d> crossings =
                                            curveCrossings(*curve_a, *curve_b, 0.01);
                                        double nearest_endpoint =
                                            std::numeric_limits<double>::infinity();
                                        double farthest_endpoint = 0.0;
                                        for (const Vec2d& crossing : crossings)
                                        {
                                            const double endpoint_distance =
                                                distToAllEndpoints(
                                                    crossing, *curve_a, *curve_b);
                                            nearest_endpoint = std::min(
                                                nearest_endpoint, endpoint_distance);
                                            farthest_endpoint = std::max(
                                                farthest_endpoint, endpoint_distance);
                                        }
                                        fprintf(stderr,
                                                "[ISG_PROFILE] U-turn family %s reject %s alpha=%.3f bias=%.2f/%.2f crossings=%zu endpoint_range=%.4f..%.4f\n",
                                                family.front()->id.c_str(),
                                                pair_key.c_str(),
                                                selected[depth]->alpha,
                                                selected[depth]->entry_bias,
                                                selected[depth]->exit_bias,
                                                crossings.size(), nearest_endpoint,
                                                farthest_endpoint);
                                    }
                                }
                                valid = false;
                                break;
                            }
                        }
                        if (valid && search_family(depth + 1))
                            return true;
                        selected_flags[depth] = false;
                    }
                    return false;
                };

            // 先选择普通阻塞兄弟，再进入 U-turn 家族回溯，使外部几何约束
            // 在家族首层生效，避免无解时先展开大量内部组合。
            search_external = [&](std::size_t depth) {
                if (++search_nodes > max_search_nodes)
                    return false;
                if (depth == external_ids.size())
                    return search_family(0);
                const ConnId& external_id = external_ids[depth];
                for (const ExternalCandidate& item :
                         external_pools[external_id]) {
                    selected_external[depth] = item;
                    selected_external_flags[depth] = true;
                    bool valid = true;
                    for (const CurvePair* pair_ptr :
                             external_incident_pairs[depth]) {
                        const CurvePair& pair = *pair_ptr;
                        if (pair.exempt == CrossExemption::StructuralCross)
                            continue;
                        const auto ia = family_index.find(pair.id_a);
                        const auto ib = family_index.find(pair.id_b);
                        const auto ea = external_index.find(pair.id_a);
                        const auto eb = external_index.find(pair.id_b);
                        const bool a_family = ia != family_index.end();
                        const bool b_family = ib != family_index.end();
                        const bool a_external = ea != external_index.end();
                        const bool b_external = eb != external_index.end();
                        if ((!a_family && !b_family &&
                             !a_external && !b_external) ||
                            (a_family && !selected_flags[ia->second]) ||
                            (b_family && !selected_flags[ib->second]) ||
                            (a_external && !selected_external_flags[ea->second]) ||
                            (b_external && !selected_external_flags[eb->second]))
                            continue;
                        const BezierCurve* curve_a = nullptr;
                        const BezierCurve* curve_b = nullptr;
                        if (a_family)
                            curve_a = &selected[ia->second]->curve;
                        else if (a_external)
                            curve_a = &selected_external[ea->second].curve;
                        else {
                            const auto ri = baseline_index.find(pair.id_a);
                            if (ri != baseline_index.end() &&
                                results[ri->second].curve)
                                curve_a = results[ri->second].curve.get();
                        }
                        if (b_family)
                            curve_b = &selected[ib->second]->curve;
                        else if (b_external)
                            curve_b = &selected_external[eb->second].curve;
                        else {
                            const auto ri = baseline_index.find(pair.id_b);
                            if (ri != baseline_index.end() &&
                                results[ri->second].curve)
                                curve_b = results[ri->second].curve.get();
                        }
                        if (!curve_a || !curve_b)
                            continue;
                        const bool forbidden = a_family && b_family
                            ? family_pair_forbidden(*curve_a, *curve_b)
                            : curvesHaveForbiddenSameClusterIntersection(
                                  *curve_a, *curve_b, kClusterEndpointTol);
                        if (forbidden) {
                            if (profile_rejections)
                                ++external_rejections[
                                    pair.id_a + "|" + pair.id_b];
                            valid = false;
                            break;
                        }
                    }
                    if (valid && search_external(depth + 1))
                        return true;
                    selected_external_flags[depth] = false;
                }
                return false;
            };

            if (!(external_ids.empty() ? search_family(0) : search_external(0))) {
                if (isgProfile())
                    fprintf(stderr,
                            "[ISG_PROFILE] U-turn family %s no joint candidate pools=%zu/%zu/%zu/%zu nodes=%zu\n",
                            family.front()->id.c_str(), pools[0].size(),
                            pools.size() > 1 ? pools[1].size() : 0,
                            pools.size() > 2 ? pools[2].size() : 0,
                            pools.size() > 3 ? pools[3].size() : 0,
                            search_nodes);
                if (isgProfile()) {
                    for (const auto& item : external_rejections)
                        fprintf(stderr,
                                "[ISG_PROFILE] U-turn family %s external reject %s=%zu\n",
                                family.front()->id.c_str(), item.first.c_str(),
                                item.second);
                    for (const auto& item : internal_rejections)
                        fprintf(stderr,
                                "[ISG_PROFILE] U-turn family %s internal reject %s=%zu\n",
                                family.front()->id.c_str(), item.first.c_str(),
                                item.second);
                }
                continue;
            }

            std::vector<CurvePatchEntry> accepted_patch;
            accepted_patch.reserve(family.size() + external_ids.size());
            for (std::size_t i = 0; i < family.size(); ++i)
                accepted_patch.push_back(
                    CurvePatchEntry{family[i]->id, selected[i]->curve});
            for (std::size_t i = 0; i < external_ids.size(); ++i)
                accepted_patch.push_back(
                    CurvePatchEntry{external_ids[i], selected_external[i].curve});
            if (!AtomicCurvePatch().apply(
                    accepted_patch, results,
                    [&](ConnectivityCurve& result,
                        const BezierCurve& candidate) {
                        setConnectivityCurveGeometry(result, candidate);
                        validate(result, input, sdf);
                    }))
                continue;
            if (isgProfile()) {
                fprintf(stderr,
                        "[ISG_PROFILE] U-turn family %s joint committed size=%zu nodes=%zu",
                        family.front()->id.c_str(),
                        family.size() + external_ids.size(), search_nodes);
                for (std::size_t i = 0; i < selected.size(); ++i)
                    fprintf(stderr, " %s=%.2f/%.2f@%.2f",
                            family[i]->id.c_str(), selected[i]->entry_bias,
                            selected[i]->exit_bias, selected[i]->station_offset);
                for (const ConnId& external_id : external_ids)
                    fprintf(stderr, " %s=external", external_id.c_str());
                fprintf(stderr, "\n");
            }
        }
    }
    if (isgProfile()) {
        fprintf(stderr, "[ISG_PROFILE] U-turn family stage %.3f ms\n",
                std::chrono::duration<double, std::milli>(
                    std::chrono::steady_clock::now() - uturn_family_stage_t0).count());
    }

    // 最终普通曲线形态审计：前面的同簇修复可能压短某条共享端点
    // 曲线的尾把手，导致短急弯曲率超过允许上限。此处只对无固定形态、
    // 非掉头的单段普通曲线做受约束重建，避免把保序修复留下非法形态。
    auto ordinary_family_stage_t0 = std::chrono::steady_clock::now();
    {
        auto final_idx = resultIndexById(results);
        auto final_curves = curveMapFromResults(results);
        for (auto& result : results) {
            if (preserved_fixed_ids.count(result.id) || !result.curve)
                continue;
            const Connectivity* conn = scene.view.connectivity(result.id);
            if (!conn || isGeometricUTurnConn(*conn, input))
                continue;
            auto entry = scene.view.entryFrame(conn->entry_lane_id);
            auto exit_ = scene.view.exitFrame(conn->exit_lane_id);
            Vec2d chord = exit_.first - entry.first;
            if (chord.norm() < 1e-8 || entry.second.norm() < 1e-8)
                continue;
            const double turn_strength = std::abs(cross2d(
                entry.second.normalized(), chord.normalized()));
            if (turn_strength <= 0.35 ||
                isNonUTurnTurnShapeAcceptable(
                    *result.curve, chord.norm(), turn_strength))
                continue;
            auto siblings = buildSiblings(
                result.id, final_curves, cluster_solver_,
                input.connectivities, true, &preserved_fixed_ids);
            auto sampled = sampleSiblingsForIntersections(siblings);
            BezierCurve repaired = *result.curve;
            if (!tryShapeSafeSingleCubic(
                    entry.first, entry.second, exit_.first, exit_.second,
                    input, sdf, sampled, true, repaired) ||
                !isNonUTurnTurnShapeAcceptable(
                    repaired, chord.norm(), turn_strength) ||
                assessCurveRisk(repaired, input, sdf, sampled, true).physical())
                continue;
            const int old_cross = constrainedCrossCountForId(
                result.id, *result.curve, results, final_idx,
                constrainedNeighborMap(results, cluster_solver_),
                kClusterEndpointTol);
            const int new_cross = constrainedCrossCountForId(
                result.id, repaired, results, final_idx,
                constrainedNeighborMap(results, cluster_solver_),
                kClusterEndpointTol);
            if (new_cross > old_cross)
                continue;
            setConnectivityCurveGeometry(result, repaired);
            validate(result, input, sdf);
            final_curves[result.id] = repaired;
            final_idx = resultIndexById(results);
        }
    }

    // 复杂 Boundary 场景的逐条收口可能先后改写同一共享入口扇出族：每一步对
    // 当时的邻居都可接受，但组合结果会让控制多边形互穿，甚至把某个把手推过
    // 方向射线交点。对已经出现族内交叉或轴向越界的同侧普通转向族，原子尝试
    // 恢复整族首选单段；只有物理、形态、固定形态和所有外部同簇约束都不恶化
    // 时才整体提交，避免逐条恢复因中间态交叉而互相否决。
    {
        std::unordered_map<std::string, std::vector<const Connectivity*>> families;
        for (const auto& conn : input.connectivities) {
            if (preserved_fixed_ids.count(conn.id) || conn.fixed_shape ||
                isGeometricUTurnConn(conn, input))
                continue;
            const double signed_turn = signedTurnStrengthOfConnId(scene.view, conn.id);
            if (std::abs(signed_turn) < 0.35)
                continue;
            const std::string key = conn.entry_lane_id + "\x1f" + conn.exitGroupId +
                (signed_turn > 0.0 ? "\x1fL" : "\x1fR");
            families[key].push_back(&conn);
        }

        const OrdinaryCurveInitializer ordinary_initializer;
        for (const auto& family_item : families) {
            const auto& family = family_item.second;
            if (family.size() < 2)
                continue;

            auto current_idx = resultIndexById(results);
            std::unordered_set<ConnId> family_ids;
            bool needs_restore = false;
            for (const Connectivity* conn : family) {
                family_ids.insert(conn->id);
                auto ri = current_idx.find(conn->id);
                if (ri == current_idx.end() || !results[ri->second].curve) {
                    needs_restore = false;
                    break;
                }
                auto entry = scene.view.entryFrame(conn->entry_lane_id);
                auto exit_ = scene.view.exitFrame(conn->exit_lane_id);
                if (!ordinarySingleCubicControlsValid(
                        *results[ri->second].curve, entry.first, entry.second,
                        exit_.first, exit_.second, 1e-5, true))
                    needs_restore = true;
            }
            for (size_t i = 0; i < family.size(); ++i) {
                auto ia = current_idx.find(family[i]->id);
                if (ia == current_idx.end() || !results[ia->second].curve)
                    continue;
                for (size_t j = i + 1; j < family.size(); ++j) {
                    auto ib = current_idx.find(family[j]->id);
                    if (ib == current_idx.end() || !results[ib->second].curve)
                        continue;
                    if (curvesIntersectBusiness(
                            *results[ia->second].curve,
                            *results[ib->second].curve, 1.5) ||
                        sharedEndpointControlPolylinesCross(
                            *results[ia->second].curve,
                            *results[ib->second].curve, 0.10))
                        needs_restore = true;
                }
            }
            if (!needs_restore)
                continue;

            std::vector<CurvePatchEntry> patch;
            patch.reserve(family.size());
            bool candidates_valid = true;
            for (const Connectivity* conn : family) {
                auto entry = scene.view.entryFrame(conn->entry_lane_id);
                auto exit_ = scene.view.exitFrame(conn->exit_lane_id);
                const Vec2d chord = exit_.first - entry.first;
                if (chord.norm() < 1e-8) {
                    candidates_valid = false;
                    break;
                }
                const Vec2d t0 = entry.second.norm() > 1e-8
                    ? entry.second.normalized() : chord.normalized();
                const Vec2d t1 = exit_.second.norm() > 1e-8
                    ? exit_.second.normalized() : chord.normalized();
                const double turn_strength = std::abs(
                    cross2d(t0, chord.normalized()));
                // 自然首选有时会贴近 mode=2 的 RoadEdge；在保持两个控制点
                // 位于端点方向轴、且不改变 G1 的前提下，沿有限把手长度搜索
                // 可行单段，而不是把控制点横向推离入口方向。
                std::vector<BezierCurve> candidates;
                candidates.push_back(ordinary_initializer.buildPreferredSingleCubic(
                    entry.first, t0, exit_.first, t1));
                for (double alpha : {0.25, 0.35, 0.45, 0.55, 0.65, 0.80})
                    candidates.push_back(ordinary_initializer.buildSingleCubic(
                        entry.first, t0, exit_.first, t1, alpha));
                BezierCurve candidate;
                bool found_candidate = false;
                for (const auto& trial : candidates) {
                    const bool controls_valid = ordinarySingleCubicControlsValid(
                        trial, entry.first, t0, exit_.first, t1, 1e-5, true);
                    const bool self_cross = curveSelfIntersectsBusiness(trial, 1.0);
                    const bool shape_valid = isNonUTurnTurnShapeAcceptable(
                        trial, chord.norm(), turn_strength);
                    const CurveRisk trial_risk =
                        assessCurveRisk(trial, input, sdf, {}, true);
                    if (trial.numSegments() == 1 && controls_valid && !self_cross &&
                        shape_valid && !trial_risk.physical()) {
                        candidate = trial;
                        found_candidate = true;
                        break;
                    }
                }
                if (!found_candidate) {
                    if (isgDebugPairRepair())
                        fprintf(stderr,
                                "[FAMILY-RESTORE] %s member=%s reject: no axis-constrained physical-safe candidate\n",
                                family.front()->id.c_str(), conn->id.c_str());
                    candidates_valid = false;
                    break;
                }
                patch.push_back(CurvePatchEntry{conn->id, candidate});
            }
            if (!candidates_valid || patch.size() != family.size())
            {
                if (isgDebugPairRepair())
                    fprintf(stderr,
                            "[FAMILY-RESTORE] %s rejected candidate validity size=%zu/%zu\n",
                            family.front()->id.c_str(), patch.size(), family.size());
                continue;
            }

            std::vector<ConnectivityCurve> candidate_results = results;
            if (!AtomicCurvePatch().apply(
                    patch, candidate_results,
                    [&](ConnectivityCurve& result, const BezierCurve& candidate) {
                        setConnectivityCurveGeometry(result, candidate);
                    }))
                continue;
            auto candidate_idx = resultIndexById(candidate_results);

            int current_crosses = 0;
            int candidate_crosses = 0;
            bool family_internal_valid = true;
            for (const auto& pair : cluster_solver_.pairs()) {
                if (pair.exempt == CrossExemption::StructuralCross ||
                    (!family_ids.count(pair.id_a) && !family_ids.count(pair.id_b)))
                    continue;
                auto old_a = current_idx.find(pair.id_a);
                auto old_b = current_idx.find(pair.id_b);
                auto new_a = candidate_idx.find(pair.id_a);
                auto new_b = candidate_idx.find(pair.id_b);
                if (old_a == current_idx.end() || old_b == current_idx.end() ||
                    new_a == candidate_idx.end() || new_b == candidate_idx.end() ||
                    !results[old_a->second].curve || !results[old_b->second].curve ||
                    !candidate_results[new_a->second].curve ||
                    !candidate_results[new_b->second].curve)
                    continue;
                const BezierCurve& old_curve_a = *results[old_a->second].curve;
                const BezierCurve& old_curve_b = *results[old_b->second].curve;
                const BezierCurve& new_curve_a = *candidate_results[new_a->second].curve;
                const BezierCurve& new_curve_b = *candidate_results[new_b->second].curve;
                if (curvesHaveForbiddenSameClusterIntersection(
                        old_curve_a, old_curve_b, kClusterEndpointTol))
                    ++current_crosses;
                if (curvesHaveForbiddenSameClusterIntersection(
                        new_curve_a, new_curve_b, kClusterEndpointTol))
                    ++candidate_crosses;
                if (family_ids.count(pair.id_a) && family_ids.count(pair.id_b) &&
                    (curvesIntersectBusiness(new_curve_a, new_curve_b, 1.5) ||
                     sharedEndpointControlPolylinesCross(
                         new_curve_a, new_curve_b,
                         kClusterEndpointTol))) {
                    family_internal_valid = false;
                    break;
                }
            }
            if (!family_internal_valid || candidate_crosses > current_crosses)
            {
                if (isgDebugPairRepair())
                    fprintf(stderr,
                            "[FAMILY-RESTORE] %s rejected internal=%d crosses=%d->%d\n",
                            family.front()->id.c_str(), family_internal_valid ? 1 : 0,
                            current_crosses, candidate_crosses);
                continue;
            }

            AtomicCurvePatch().apply(
                patch, results,
                [&](ConnectivityCurve& result, const BezierCurve& candidate) {
                    setConnectivityCurveGeometry(result, candidate);
                    validate(result, input, sdf);
                });
            if (isgDebugPairRepair())
                fprintf(stderr,
                        "[FAMILY-RESTORE] %s committed size=%zu crosses=%d->%d\n",
                        family.front()->id.c_str(), family.size(),
                        current_crosses, candidate_crosses);
        }
    }

    // 所有拓扑/边界修复完成后，再锁定普通转向的方向交点规范形态。
    // 前面的阶段允许为了消除共享端点冲突暂时压短把手；若最终首选候选
    // 通过完整物理与同簇门禁，则应恢复到交点前 1/3，而不能让最后一次
    // 局部收口把右/左转退化成近直线短弧。
    {
        auto final_idx = resultIndexById(results);
        auto final_curves = curveMapFromResults(results);
        const OrdinaryCurveInitializer ordinary_initializer;
        for (auto& result : results) {
            if (preserved_fixed_ids.count(result.id) || !result.curve)
                continue;
            const Connectivity* conn = scene.view.connectivity(result.id);
            if (!conn || conn->fixed_shape || isGeometricUTurnConn(*conn, input))
                continue;
            auto entry = scene.view.entryFrame(conn->entry_lane_id);
            auto exit_ = scene.view.exitFrame(conn->exit_lane_id);
            const Vec2d chord = exit_.first - entry.first;
            if (chord.norm() < 1e-8 || entry.second.norm() < 1e-8)
                continue;
            const Vec2d t0 = entry.second.normalized();
            const Vec2d t1 = exit_.second.norm() > 1e-8
                ? exit_.second.normalized() : chord.normalized();
            const double turn_strength = std::abs(cross2d(t0, chord.normalized()));
            if (turn_strength <= 0.35)
                continue;
            const BezierCurve preferred = ordinary_initializer.buildPreferredSingleCubic(
                entry.first, t0, exit_.first, t1);
            if (preferred.empty() || preferred.numSegments() != 1 ||
                !ordinarySingleCubicControlsValid(
                    preferred, entry.first, t0, exit_.first, t1, 1e-5, true) ||
                !isNonUTurnTurnShapeAcceptable(
                    preferred, chord.norm(), turn_strength) ||
                curveSelfIntersectsBusiness(preferred, 1.0))
                continue;
            auto siblings = buildSiblings(
                result.id, final_curves, cluster_solver_, input.connectivities,
                true, &preserved_fixed_ids);
            auto sampled = sampleSiblingsForIntersections(siblings);
            const int old_cross = constrainedCrossCountForId(
                result.id, *result.curve, results, final_idx,
                constrainedNeighborMap(results, cluster_solver_),
                kClusterEndpointTol);
            const auto bounds = ordinarySingleCubicHandleBounds(
                entry.first, t0, exit_.first, t1, true);
            const double preferred_h0 =
                (preferred.segs.front().ctrl[1] - entry.first).norm();
            const double preferred_h1 =
                (exit_.first - preferred.segs.front().ctrl[2]).norm();
            BezierCurve best_candidate;
            double best_deviation = std::numeric_limits<double>::infinity();
            double best_curvature = std::numeric_limits<double>::infinity();
            for (double f0 : {1.0, 0.95, 0.90, 0.85, 0.80, 0.70, 0.60, 0.50, 0.40, 0.30}) {
                for (double f1 : {1.0, 0.95, 0.90, 0.85, 0.80, 0.70, 0.60, 0.50, 0.40, 0.30}) {
                    const double h0 = std::max(
                        bounds.start_min, std::min(bounds.start_max, preferred_h0 * f0));
                    const double h1 = std::max(
                        bounds.end_min, std::min(bounds.end_max, preferred_h1 * f1));
                    BezierSegment segment;
                    segment.ctrl[0] = entry.first;
                    segment.ctrl[1] = entry.first + bounds.start_dir * h0;
                    segment.ctrl[2] = exit_.first - bounds.end_dir * h1;
                    segment.ctrl[3] = exit_.first;
                    BezierCurve candidate;
                    candidate.segs.push_back(segment);
                    if (!ordinarySingleCubicControlsValid(
                            candidate, entry.first, t0, exit_.first, t1, 1e-5, true) ||
                        !isNonUTurnTurnShapeAcceptable(
                            candidate, chord.norm(), turn_strength) ||
                        curveSelfIntersectsBusiness(candidate, 1.0))
                        continue;
                    if (assessCurveRisk(candidate, input, sdf, sampled, true).physical())
                        continue;
                    const int candidate_cross = constrainedCrossCountForId(
                        result.id, candidate, results, final_idx,
                        constrainedNeighborMap(results, cluster_solver_),
                        kClusterEndpointTol);
                    if (candidate_cross > old_cross)
                        continue;
                    const double deviation =
                        (preferred_h0 - h0) + (preferred_h1 - h1);
                    const double curvature = candidate.maxCurvature(40);
                    if (deviation < best_deviation - 1e-6 ||
                        (std::abs(deviation - best_deviation) <= 1e-6 &&
                         curvature < best_curvature)) {
                        best_candidate = candidate;
                        best_deviation = deviation;
                        best_curvature = curvature;
                    }
                }
            }
            if (best_candidate.empty())
                continue;
            setConnectivityCurveGeometry(result, best_candidate);
            validate(result, input, sdf);
            final_curves[result.id] = best_candidate;
            final_idx = resultIndexById(results);
        }
    }

    // 最终普通曲线收口必须按共享端点族联合求解。逐条恢复会出现这样的闭环：
    // 第一条曲线在当时的兄弟快照中可接受，第二条随后改变扇出顺序，最后又把
    // 第一条恢复成越过方向交点的大把手，导致真实中段相交。这里仅对已经发现
    // 交叉或控制多边形前兆的共享入口/出口族建立有限轴上候选池，并一次性回溯
    // 验收整族；候选仍须通过形态、Boundary/Fence、固定形态邻居和同簇外部关系。
    {
        std::unordered_map<std::string, std::vector<ConnId>> endpoint_families;
        auto final_idx = resultIndexById(results);
        for (const auto& conn : input.connectivities) {
            if (preserved_fixed_ids.count(conn.id) || conn.fixed_shape ||
                isGeometricUTurnConn(conn, input))
                continue;
            auto ri = final_idx.find(conn.id);
            if (ri == final_idx.end() || !results[ri->second].curve ||
                results[ri->second].curve->numSegments() != 1)
                continue;
            endpoint_families["entry:" + conn.entry_lane_id].push_back(conn.id);
            endpoint_families["exit:" + conn.exit_lane_id].push_back(conn.id);
        }

        struct StrictOrdinaryCandidate {
            BezierCurve curve;
            SampledCurve sampled;
            double score = std::numeric_limits<double>::infinity();
        };
        const auto strict_family_pair_required =
            [&](const ConnId& a, const ConnId& b) {
                return cluster_solver_.pairExists(a, b) &&
                    cluster_solver_.exemptionOf(a, b) !=
                        CrossExemption::StructuralCross;
            };

        // 这些档位只用于最终族收口，不是连续优化器的参数网格。保留覆盖
        // 低/中/高把手的代表点即可；过密网格会在 min/max 钳制后产生大量
        // 相同几何，既没有新增可行形状，也会把无解族的回溯时间放大。
        const std::vector<double> handle_fractions = {
            0.02, 0.06, 0.12, 0.20, 0.32, 0.48, 0.68, 0.85, 1.0};
        // 短转向的物理可行区可能只落在高把手区间的一条窄带内（例如
        // 100000443 的连接 5）。只增加该区间的有限细分，不把长转向的
        // 全部组合数按比例放大，避免无解族重新触发全局穷举。
        const std::vector<double> short_turn_handle_fractions = {
            0.63, 0.66, 0.6631, 0.69, 0.72, 0.75, 0.78, 0.81, 0.84,
            0.87, 0.90, 0.93, 0.96, 0.99};
        // 长弦共享入口/出口族的直行成员可能需要中等把手长度，才能与
        // 同端点的短转向成员保持同侧。该窄档只对二元族启用，覆盖
        // 100000443 的连接 9 所需的约 0.42 比例；不扩大三元以上族的
        // 组合空间。
        const std::vector<double> long_pair_handle_fractions = {
            0.36, 0.40, 0.418, 0.41812, 0.42, 0.44, 0.46};
        // 族级保序参照形态：普通转向的首选单段表达（转向按端点方向交点距离
        // 的 2/3 升阶，见 OrdinaryCurveInitializer::buildPreferredSingleCubic）
        // 只由本成员自身的端点位置与切向决定，与生成次序、兄弟形态和优化器
        // 无关。因此同一共享端点扇出族的全体成员在该表达下天然按转弯半径
        // 单调嵌套：110004764 的入口族 8/10/11/12/13/14 在首选表达下 15 对
        // 全部互不相交（曲线与控制折线都不交），而逐条 L-BFGS 沿把手轴各自
        // 挪动后出现 8|10、11|12、11|14、13|14 四对非端点相交。
        //
        // 下面的把手网格是围绕"当前形态"展开的等分档位，一般表达不出这个
        // 参照形态（11 的首选把手占弦长 0.655/0.622，网格相邻档位只有 0.48
        // 与 0.68），族内联合回溯因此无论怎样换档都回不到它。把它显式放进
        // 候选池并排在首位：回溯的第一个组合就是"全族同时回到首选表达"，
        // 命中即一次性恢复整族次序；未命中则继续按原档位回溯，行为不变。
        const auto strict_natural_reference_candidate =
            [&](const Vec2d& entry_point, const Vec2d& t0,
                const Vec2d& exit_point, const Vec2d& t1,
                double chord_len, double turn_strength,
                BezierCurve& out) {
                out = OrdinaryCurveInitializer().buildPreferredSingleCubic(
                    entry_point, t0, exit_point, t1);
                if (out.numSegments() != 1)
                    return false;
                if (!ordinarySingleCubicControlsValid(
                        out, entry_point, t0, exit_point, t1, 1e-5, false) ||
                    curveSelfIntersectsBusiness(out, 1.0))
                    return false;
                return turn_strength <= 0.25
                    ? isStraightLikeShapeAcceptable(out, chord_len)
                    : isNonUTurnTurnShapeAcceptable(out, chord_len, turn_strength);
            };
        constexpr std::size_t kMaxStrictPoolCandidates = 24;
        // 二元共享端点族需要先从完整物理安全候选中寻找兼容对；否则
        // 5|9 这类可行对可能在两边独立按当前几何排序时同时被截掉。
        // 宽池只用于二元族，后续仍收口到 24 个代表候选。
        constexpr std::size_t kMaxStrictPairPoolCandidates = 512;
        const Polygon2d* strict_fence_outline =
            !input.area.is_rough && !input.area.geometry.outer.empty()
                ? &input.area.geometry : nullptr;
        // Boundary 的中心侧可靠性和链连接关系只依赖输入，不依赖候选曲线。
        // 在族阶段复用该快照，避免每个候选都重复执行 O(Boundary^2) 的连接扫描。
        const std::vector<BoundarySafetySegment> strict_boundary_segments =
            buildBoundarySafetySegments(
                input.boundaries, boundarySafetyCenter(input), nullptr,
                strict_fence_outline);
        // 最终普通族的真实曲线判定必须与输出审计保持一致。共享入口/出口
        // 的端点附近，采样折线会把二阶分离误判为相交；1.5m 只豁免该
        // 连接点邻域，远离端点的相交仍由同一业务判定硬性拒绝。控制折线
        // 前兆继续使用 kClusterEndpointTol，不能借此放宽族序。
        constexpr double kOrdinaryFamilyCurveCrossTol =
            kSharedEndpointPairCrossTol;
        auto assessStrictOrdinaryPhysical =
            [&](const BezierCurve& candidate,
                // Keep candidate admission on the same sampling resolution as
                // the final input-level boundary audit. A coarser 0.40m
                // polyline can chord across a narrow boundary gap and reject
                // a curve that the 0.18m audit correctly accepts.
                double sample_spacing = 0.18) {
                CurveRisk risk;
                // 纯边界路口没有障碍物，SDF 查询既不会改变结果又会在候选
                // 收口阶段重复扫描整条曲线；只有存在障碍物时才需要计算。
                const double ms = input.obstacles.empty()
                    ? 0.0 : minSDFAlongCurveAdaptive(candidate, sdf);
                risk.obstacle = !input.obstacles.empty() &&
                    (curveIntersectsObstacles(candidate, input.obstacles) || ms < 0.0);
                const BoundarySafetyResult boundary = curveBoundarySafety(
                    candidate, strict_boundary_segments, 128,
                    kConnectionPointTolerance, 0.10, 0.05, sample_spacing);
                risk.boundary = boundary.intersects || boundary.outside_road_edge;
                if (roadEdgeClearanceViolation(
                        candidate, input.boundaries, input.mode) > 0.0)
                    risk.boundary = true;
                risk.fence = !input.area.is_rough &&
                    (!input.area.geometry.outer.empty() &&
                     curveLeavesFence(candidate, input.area.geometry));
                return risk;
            };
        // 严格族候选可能在多个收口轮次中重复出现；其物理风险只由几何和
        // 输入场景决定，与族内组合无关。缓存后，联合搜索仍使用完整物理
        // 门禁，但不会为同一把手几何重复扫描 Boundary/Fence。
        std::unordered_map<std::string, bool> strict_physical_memo;
        auto strict_curve_key = [](const BezierCurve& curve) {
            std::string key;
            char buffer[96];
            for (const auto& segment : curve.segs) {
                for (const auto& point : segment.ctrl) {
                    std::snprintf(buffer, sizeof(buffer), "%.9f,%.9f;",
                                  point.x(), point.y());
                    key += buffer;
                }
            }
            return key;
        };
        std::unordered_map<std::string, bool> strict_exact_physical_memo;
        auto strictCandidateHasExactPhysicalRisk =
            [&](const BezierCurve& candidate) {
                const std::string key = strict_curve_key(candidate);
                auto memo = strict_exact_physical_memo.find(key);
                if (memo != strict_exact_physical_memo.end())
                    return memo->second;
                const bool physical = assessStrictOrdinaryPhysical(
                    candidate, 0.18).physical();
                strict_exact_physical_memo.emplace(key, physical);
                return physical;
            };
        std::unordered_map<std::string, bool> strict_exact_pair_memo;
        auto strictPairHasExactCross =
            [&](const BezierCurve& a, const BezierCurve& b) {
                std::string a_key = strict_curve_key(a);
                std::string b_key = strict_curve_key(b);
                if (a_key > b_key)
                    std::swap(a_key, b_key);
                const std::string key = a_key + "|" + b_key;
                auto memo = strict_exact_pair_memo.find(key);
                if (memo != strict_exact_pair_memo.end())
                    return memo->second;
                const bool crosses = curvesIntersectBusiness(
                    a, b, (a.numSegments() > 1 || b.numSegments() > 1)
                        ? kStrictSameClusterEndpointTol
                        : kSharedEndpointPairCrossTol);
                strict_exact_pair_memo.emplace(key, crosses);
                return crosses;
            };

        // 为一个失败的普通族按外部阻塞兄弟补齐最小影响闭包。闭包成员仍只
        // 使用单段 G1 cubic 的轴上把手候选；所有闭包内外的同簇关系都要求
        // 候选彻底无交叉/接触，因而不会把“先制造一个新违约、下一轮再修”
        // 的中间状态提交到结果中。
        constexpr std::size_t kMaxStrictClosurePoolCandidates = 64;
        auto build_strict_closure_pool =
            [&](const ConnId& id, std::vector<StrictOrdinaryCandidate>& out,
                std::size_t max_candidates) {
                auto ri = final_idx.find(id);
                const Connectivity* conn = scene.view.connectivity(id);
                if (ri == final_idx.end() || !results[ri->second].curve ||
                    !conn || conn->fixed_shape ||
                    preserved_fixed_ids.count(id) ||
                    isGeometricUTurnConn(*conn, input) ||
                    results[ri->second].curve->numSegments() != 1)
                    return false;
                const auto entry = scene.view.entryFrame(conn->entry_lane_id);
                const auto exit_ = scene.view.exitFrame(conn->exit_lane_id);
                const Vec2d chord = exit_.first - entry.first;
                if (chord.norm() < 1e-8)
                    return false;
                const Vec2d t0 = entry.second.norm() > 1e-8
                    ? entry.second.normalized() : chord.normalized();
                const Vec2d t1 = exit_.second.norm() > 1e-8
                    ? exit_.second.normalized() : chord.normalized();
                const double turn_strength = std::abs(
                    cross2d(t0, chord.normalized()));
                const auto bounds = ordinarySingleCubicHandleBounds(
                    entry.first, t0, exit_.first, t1, false);
                std::vector<double> fractions = handle_fractions;
                if (chord.norm() <= 15.0)
                    fractions.insert(fractions.end(),
                                     short_turn_handle_fractions.begin(),
                                     short_turn_handle_fractions.end());
                if (max_candidates >= kMaxStrictPairPoolCandidates &&
                    chord.norm() > 15.0)
                    fractions.insert(fractions.end(),
                                     long_pair_handle_fractions.begin(),
                                     long_pair_handle_fractions.end());
                const BezierCurve& current = *results[ri->second].curve;
                const double current_h0 =
                    (current.segs.front().ctrl[1] - entry.first).dot(
                        bounds.start_dir);
                const double current_h1 =
                    (exit_.first - current.segs.front().ctrl[2]).dot(
                        bounds.end_dir);
                if (bounds.start_max > 0.05 && std::isfinite(current_h0))
                    fractions.push_back(std::max(
                        0.02, std::min(1.0, current_h0 / bounds.start_max)));
                if (bounds.end_max > 0.05 && std::isfinite(current_h1))
                    fractions.push_back(std::max(
                        0.02, std::min(1.0, current_h1 / bounds.end_max)));
                std::sort(fractions.begin(), fractions.end());
                fractions.erase(std::unique(fractions.begin(), fractions.end()),
                                fractions.end());
                std::vector<StrictOrdinaryCandidate> pending;
                std::size_t control_rejects = 0;
                std::size_t shape_rejects = 0;
                std::size_t physical_rejects = 0;
                {
                    // 参照形态排在闭包池首位，理由见
                    // strict_natural_reference_candidate 的说明。
                    BezierCurve natural;
                    if (strict_natural_reference_candidate(
                            entry.first, t0, exit_.first, t1, chord.norm(),
                            turn_strength, natural)) {
                        StrictOrdinaryCandidate item;
                        item.curve = std::move(natural);
                        item.score = -1.0;
                        pending.push_back(std::move(item));
                    }
                }
                for (double f0 : fractions) {
                    for (double f1 : fractions) {
                        const double h0 = std::max(
                            bounds.start_min,
                            std::min(bounds.start_max, f0 * bounds.start_max));
                        const double h1 = std::max(
                            bounds.end_min,
                            std::min(bounds.end_max, f1 * bounds.end_max));
                        BezierSegment segment;
                        segment.ctrl[0] = entry.first;
                        segment.ctrl[1] = entry.first + bounds.start_dir * h0;
                        segment.ctrl[2] = exit_.first - bounds.end_dir * h1;
                        segment.ctrl[3] = exit_.first;
                        BezierCurve candidate;
                        candidate.segs.push_back(segment);
                        if (!ordinarySingleCubicControlsValid(
                                candidate, entry.first, t0,
                                exit_.first, t1, 1e-5, false) ||
                            curveSelfIntersectsBusiness(candidate, 1.0))
                        {
                            ++control_rejects;
                            continue;
                        }
                        const bool shape_ok = turn_strength <= 0.25
                            ? isStraightLikeShapeAcceptable(candidate, chord.norm())
                            : isNonUTurnTurnShapeAcceptable(
                                  candidate, chord.norm(), turn_strength);
                        if (!shape_ok)
                        {
                            ++shape_rejects;
                            continue;
                        }
                        const std::string key = strict_curve_key(candidate);
                        bool duplicate = false;
                        for (const auto& existing : pending) {
                            if (strict_curve_key(existing.curve) == key) {
                                duplicate = true;
                                break;
                            }
                        }
                        if (duplicate)
                            continue;
                        StrictOrdinaryCandidate item;
                        item.curve = std::move(candidate);
                        item.score =
                            std::abs(current_h0 - h0) + std::abs(current_h1 - h1) +
                            0.50 * std::abs(item.curve.arcLength() / chord.norm() - 1.10) +
                            item.curve.maxCurvature(40) + 0.01 * item.curve.arcLength();
                        pending.push_back(std::move(item));
                    }
                }
                std::sort(pending.begin(), pending.end(),
                          [](const StrictOrdinaryCandidate& a,
                             const StrictOrdinaryCandidate& b) {
                              return a.score < b.score;
                          });
                for (auto& item : pending) {
                    const std::string key = strict_curve_key(item.curve);
                    bool physical = false;
                    auto memo = strict_physical_memo.find(key);
                    if (memo != strict_physical_memo.end()) {
                        physical = memo->second;
                    } else {
                        physical = assessStrictOrdinaryPhysical(item.curve).physical();
                        strict_physical_memo.emplace(key, physical);
                    }
                    if (physical) {
                        ++physical_rejects;
                        continue;
                    }
                    item.sampled = sampleStrictPairCurve(item.curve);
                    out.push_back(std::move(item));
                    if (out.size() >= max_candidates)
                        break;
                }
                if (isgDebugPairRepair() && (id == "5" || id == "9" || id == "38")) {
                    fprintf(stderr,
                            "[STRICT-CLOSURE-POOL] id=%s limit=%zu controls=%zu shape=%zu physical=%zu pending=%zu out=%zu handles=",
                            id.c_str(), max_candidates, control_rejects,
                            shape_rejects, physical_rejects, pending.size(),
                            out.size());
                    for (const auto& item : out) {
                        const auto& g = item.curve.segs.front().ctrl;
                        fprintf(stderr, "%.3f/%.3f,",
                                (g[1] - g[0]).norm(),
                                (g[3] - g[2]).norm());
                    }
                    fprintf(stderr, "\n");
                }
                return !out.empty();
            };

        const auto strict_family_affected_by =
            [&](const std::vector<ConnId>& family,
                const std::unordered_set<ConnId>& changed_ids) {
                if (changed_ids.empty())
                    return false;
                for (const ConnId& id : family)
                    if (changed_ids.count(id))
                        return true;
                for (const auto& pair : cluster_solver_.pairs()) {
                    if (pair.exempt == CrossExemption::StructuralCross)
                        continue;
                    const bool a_in_family = std::find(
                        family.begin(), family.end(), pair.id_a) != family.end();
                    const bool b_in_family = std::find(
                        family.begin(), family.end(), pair.id_b) != family.end();
                    if (a_in_family && changed_ids.count(pair.id_b))
                        return true;
                    if (b_in_family && changed_ids.count(pair.id_a))
                        return true;
                }
                return false;
            };
        std::unordered_set<ConnId> strict_changed_last_pass;
        // 闭包扩张必须"先窄后宽"。把发现预算一次性开到最大不是单调更强：
        // 更多变量会引入更多硬约束行，AC3 可能直接清空某个成员的候选域，
        // 于是本来 5 成员有解的闭包退化成 10 成员无解，连原先修好的对也
        // 一起丢失（110003285 的 entry:43103892 就从 committed members=5
        // 退化为 closure_nodes=0，重新冒出 9|10 与 24|14 同簇相交）。
        // 因此首轮一律用最小闭包；只有确实失败的族才在后续轮次升级一次
        // 预算重试，成功过的窄解不会被更宽的搜索推翻。
        std::unordered_set<std::string> strict_failed_last_pass;
        std::unordered_set<std::string> strict_escalated_families;
        bool family_pass_changed = false;
        for (int family_pass = 0; family_pass < 3; ++family_pass) {
            family_pass_changed = false;
            std::unordered_set<ConnId> strict_changed_this_pass;
            std::unordered_set<std::string> strict_failed_this_pass;
            for (const auto& family_item : endpoint_families) {
            std::vector<ConnId> family = family_item.second;
            if (family_pass > 0 &&
                !strict_failed_last_pass.count(family_item.first) &&
                !strict_family_affected_by(
                    family, strict_changed_last_pass))
                continue;
            // 本族是否处于"升级重试"轮：仅对上一轮失败且尚未升级过的族生效。
            const bool closure_escalate =
                family_pass > 0 &&
                strict_failed_last_pass.count(family_item.first) != 0 &&
                strict_escalated_families.count(family_item.first) == 0;
            if (closure_escalate)
                strict_escalated_families.insert(family_item.first);
            StrictFamilyProfileScope family_profile{
                family_item.first.c_str()};
            if (isgDebugPairRepair() && family.size() >= 2) {
                fprintf(stderr, "[STRICT-FAMILY] inspect %s members=",
                        family_item.first.c_str());
                for (const ConnId& id : family)
                    fprintf(stderr, "%s,", id.c_str());
                fprintf(stderr, "\n");
            }
            if (family.size() < 2)
                continue;
            std::unordered_map<ConnId, SampledCurve> strict_result_samples;
            strict_result_samples.reserve(results.size());
            for (const auto& result : results)
                if (result.curve)
                    strict_result_samples.emplace(
                        result.id, sampleStrictPairCurve(*result.curve));
            bool needs_restore = false;
            bool family_has_current_curve_cross = false;
            for (std::size_t i = 0; i < family.size(); ++i) {
                auto ia = final_idx.find(family[i]);
                if (ia == final_idx.end() || !results[ia->second].curve)
                    continue;
                for (std::size_t j = i + 1; j < family.size(); ++j) {
                    if (!strict_family_pair_required(family[i], family[j]))
                        continue;
                    auto ib = final_idx.find(family[j]);
                    if (ib == final_idx.end() || !results[ib->second].curve)
                        continue;
                    const bool current_curve_cross = curvesIntersectBusiness(
                        *results[ia->second].curve,
                        *results[ib->second].curve,
                        kOrdinaryFamilyCurveCrossTol);
                    const bool current_control_cross =
                        sharedEndpointControlPolylinesCross(
                            *results[ia->second].curve,
                            *results[ib->second].curve,
                            kClusterEndpointTol);
                    if (isgDebugPairRepair() &&
                        (family_item.first == "entry:1015757" ||
                         family_item.first == "exit:1015745")) {
                        fprintf(stderr,
                                "[STRICT-FAMILY] pair %s|%s curve=%d poly=%d\n",
                                    family[i].c_str(), family[j].c_str(),
                                    current_curve_cross ? 1 : 0,
                                    current_control_cross ? 1 : 0);
                        if (family_item.first == "entry:1015757" &&
                            family[i] == "29" && family[j] == "31") {
                            const auto& ga = results[ia->second].curve->segs.front().ctrl;
                            const auto& gb = results[ib->second].curve->segs.front().ctrl;
                            fprintf(stderr,
                                    "[STRICT-FAMILY] 29|31 ctrl a1=(%.6f,%.6f) a2=(%.6f,%.6f) b1=(%.6f,%.6f) b2=(%.6f,%.6f)\n",
                                    ga[1].x(), ga[1].y(), ga[2].x(), ga[2].y(),
                                    gb[1].x(), gb[1].y(), gb[2].x(), gb[2].y());
                        }
                    }
                    // 控制折线互穿只作为前兆告警，不触发这一层昂贵的物理
                    // 候选搜索；真正的曲线非端点相交才需要建立联合收口。
                    if (current_curve_cross) {
                        family_has_current_curve_cross =
                            family_has_current_curve_cross || current_curve_cross;
                        needs_restore = true;
                        continue;
                    }
                }
            }
            if (!needs_restore)
                continue;

            // 只把当前确有冲突的成员纳入变量集合。族内其它成员保持当前
            // 几何，随后通过 external_issues 作为固定外部曲线参与验收；
            // 这样 41|40 不会被同出口但当前已安全的 7 拖成无解，30|29、
            // 43|42、14|13 的前兆修复也不会无谓改动整条扇出族。
            std::unordered_set<ConnId> active_members;
            for (std::size_t i = 0; i < family.size(); ++i) {
                auto ia = final_idx.find(family[i]);
                if (ia == final_idx.end() || !results[ia->second].curve)
                    continue;
                for (std::size_t j = i + 1; j < family.size(); ++j) {
                    if (!strict_family_pair_required(family[i], family[j]))
                        continue;
                    auto ib = final_idx.find(family[j]);
                    if (ib == final_idx.end() || !results[ib->second].curve)
                        continue;
                    if (curvesIntersectBusiness(
                            *results[ia->second].curve,
                            *results[ib->second].curve,
                            kOrdinaryFamilyCurveCrossTol)) {
                        active_members.insert(family[i]);
                        active_members.insert(family[j]);
                    }
                }
            }
            if (active_members.size() < 2)
                continue;
            if (active_members.size() < family.size()) {
                std::vector<ConnId> active_family;
                active_family.reserve(active_members.size());
                for (const ConnId& id : family)
                    if (active_members.count(id))
                        active_family.push_back(id);
                family.swap(active_family);
                if (isgDebugPairRepair()) {
                    fprintf(stderr, "[STRICT-FAMILY] %s active members=",
                            family_item.first.c_str());
                    for (const ConnId& id : family)
                        fprintf(stderr, "%s,", id.c_str());
                    fprintf(stderr, "\n");
                }
            }

            const auto strict_pool_t0 = std::chrono::steady_clock::now();
            std::vector<std::vector<StrictOrdinaryCandidate>> pools(family.size());
            // 二元族先保留物理安全的宽池。窄池搜索失败后，闭包需要从同一批
            // 已审计候选继续扩展；重新按同一把手网格做 Boundary/Fence 扫描会
            // 造成数秒级重复计算，并且可能因当前曲线评分变化产生不一致候选域。
            std::unordered_map<ConnId,
                std::vector<StrictOrdinaryCandidate>> strict_wide_pools;
            bool pools_valid = true;
            for (std::size_t fi = 0; fi < family.size() && pools_valid; ++fi) {
                int debug_shape_reject = 0;
                int debug_physical_reject = 0;
                int debug_external_reject = 0;
                int debug_physical_boundary = 0;
                int debug_physical_fence = 0;
                int debug_physical_obstacle = 0;
                const Connectivity* conn = scene.view.connectivity(family[fi]);
                auto ri = final_idx.find(family[fi]);
                if (!conn || ri == final_idx.end() || !results[ri->second].curve) {
                    pools_valid = false;
                    break;
                }
                const auto entry = scene.view.entryFrame(conn->entry_lane_id);
                const auto exit_ = scene.view.exitFrame(conn->exit_lane_id);
                const Vec2d chord = exit_.first - entry.first;
                if (chord.norm() < 1e-8) {
                    pools_valid = false;
                    break;
                }
                const Vec2d t0 = entry.second.norm() > 1e-8
                    ? entry.second.normalized() : chord.normalized();
                const Vec2d t1 = exit_.second.norm() > 1e-8
                    ? exit_.second.normalized() : chord.normalized();
                const double turn_strength = std::abs(
                    cross2d(t0, chord.normalized()));
                // 方向交点只作为候选排序参考。严格 Boundary 链可能要求
                // 物理安全的单段曲线沿端点切向轴越过该参考点；控制点仍不得
                // 横向偏移，且全局长度/曲率/单拱形态门禁不变。
                const auto bounds = ordinarySingleCubicHandleBounds(
                    entry.first, t0, exit_.first, t1, false);
                const BezierCurve& current = *results[ri->second].curve;
                const double current_h0 =
                    (current.segs.front().ctrl[1] - entry.first).dot(bounds.start_dir);
                const double current_h1 =
                    (exit_.first - current.segs.front().ctrl[2]).dot(bounds.end_dir);

                std::vector<double> member_fractions = handle_fractions;
                if (chord.norm() <= 15.0)
                    member_fractions.insert(member_fractions.end(),
                                             short_turn_handle_fractions.begin(),
                                             short_turn_handle_fractions.end());
                if (family.size() == 2 && chord.norm() > 15.0)
                    member_fractions.insert(member_fractions.end(),
                                             long_pair_handle_fractions.begin(),
                                             long_pair_handle_fractions.end());
                if (bounds.start_max > 0.05 && std::isfinite(current_h0))
                    member_fractions.push_back(std::max(
                        0.02, std::min(1.0, current_h0 / bounds.start_max)));
                if (bounds.end_max > 0.05 && std::isfinite(current_h1))
                    member_fractions.push_back(std::max(
                        0.02, std::min(1.0, current_h1 / bounds.end_max)));
                std::sort(member_fractions.begin(), member_fractions.end());
                member_fractions.erase(std::unique(
                    member_fractions.begin(), member_fractions.end()),
                    member_fractions.end());
                std::vector<StrictOrdinaryCandidate> pending;
                {
                    // 参照形态排在族候选池首位：DFS 首个探测叶子即"整族回到
                    // 自然表达"，见 strict_natural_reference_candidate 说明。
                    BezierCurve natural;
                    if (strict_natural_reference_candidate(
                            entry.first, t0, exit_.first, t1, chord.norm(),
                            turn_strength, natural)) {
                        StrictOrdinaryCandidate item;
                        item.curve = std::move(natural);
                        item.score = -1.0;
                        pending.push_back(std::move(item));
                    }
                }
                for (double f0 : member_fractions) {
                    for (double f1 : member_fractions) {
                        const double h0 = std::max(
                            bounds.start_min,
                            std::min(bounds.start_max, f0 * bounds.start_max));
                        const double h1 = std::max(
                            bounds.end_min,
                            std::min(bounds.end_max, f1 * bounds.end_max));
                        bool duplicate_handle = false;
                        for (const auto& existing : pending) {
                            const auto& g = existing.curve.segs.front().ctrl;
                            const double existing_h0 = (g[1] - g[0]).norm();
                            const double existing_h1 = (g[3] - g[2]).norm();
                            if (std::abs(existing_h0 - h0) <= 1e-7 &&
                                std::abs(existing_h1 - h1) <= 1e-7) {
                                duplicate_handle = true;
                                break;
                            }
                        }
                        if (duplicate_handle)
                            continue;
                        BezierSegment segment;
                        segment.ctrl[0] = entry.first;
                        segment.ctrl[1] = entry.first + bounds.start_dir * h0;
                        segment.ctrl[2] = exit_.first - bounds.end_dir * h1;
                        segment.ctrl[3] = exit_.first;
                        BezierCurve candidate;
                        candidate.segs.push_back(segment);
                        if (!ordinarySingleCubicControlsValid(
                                candidate, entry.first, t0,
                                exit_.first, t1, 1e-5, false) ||
                            curveSelfIntersectsBusiness(candidate, 1.0))
                            continue;
                        const bool shape_ok = turn_strength <= 0.25
                            ? isStraightLikeShapeAcceptable(candidate, chord.norm())
                            : isNonUTurnTurnShapeAcceptable(
                                  candidate, chord.norm(), turn_strength);
                        if (!shape_ok) {
                            ++debug_shape_reject;
                            continue;
                        }

                        StrictOrdinaryCandidate pending_item;
                        pending_item.curve = std::move(candidate);
                        pending_item.score =
                            std::abs(current_h0 - h0) + std::abs(current_h1 - h1) +
                            0.50 * std::abs(pending_item.curve.arcLength() /
                                            chord.norm() - 1.10) +
                            pending_item.curve.maxCurvature(40) +
                            0.01 * pending_item.curve.arcLength();
                        pending.push_back(std::move(pending_item));
                    }
                }
                std::sort(pending.begin(), pending.end(),
                          [](const StrictOrdinaryCandidate& a,
                             const StrictOrdinaryCandidate& b) {
                              return a.score < b.score;
                          });
                // 后续族内回溯最多只使用 24 个成员候选；物理门禁也只对
                // 这个有界前缀求值。跨族闭包需要的更宽候选由下方独立的
                // build_strict_closure_pool 按 32/64 的闭包预算构造。
                const std::size_t physical_budget = family.size() == 2
                    ? kMaxStrictPairPoolCandidates
                    : kMaxStrictPoolCandidates;
                for (auto& pending_item : pending) {
                    const std::string risk_key = strict_curve_key(
                        pending_item.curve);
                    bool candidate_physical = false;
                    auto risk_it = strict_physical_memo.find(risk_key);
                    if (risk_it != strict_physical_memo.end()) {
                        candidate_physical = risk_it->second;
                    } else {
                        candidate_physical = assessStrictOrdinaryPhysical(
                            pending_item.curve).physical();
                        strict_physical_memo.emplace(risk_key,
                                                     candidate_physical);
                    }
                    if (candidate_physical) {
                        ++debug_physical_reject;
                        if (isgDebugPairRepair() && family[fi] == "5") {
                            const CurveRisk risk = assessStrictOrdinaryPhysical(
                                pending_item.curve);
                            debug_physical_boundary += risk.boundary ? 1 : 0;
                            debug_physical_fence += risk.fence ? 1 : 0;
                            debug_physical_obstacle += risk.obstacle ? 1 : 0;
                        }
                        continue;
                    }
                    pending_item.sampled = sampleStrictPairCurve(
                        pending_item.curve);
                    pools[fi].push_back(std::move(pending_item));
                    if (pools[fi].size() >= physical_budget)
                        break;
                }
                if (family.size() == 2)
                    strict_wide_pools.emplace(family[fi], pools[fi]);
                if (pools[fi].empty())
                    pools_valid = false;
                if (isgDebugPairRepair() && (pools[fi].empty() ||
                                             family.size() <= 3)) {
                    fprintf(stderr,
                            "[STRICT-FAMILY] %s member=%s rejects shape=%d physical=%d external=%d pool=%zu\n",
                            family_item.first.c_str(), family[fi].c_str(),
                            debug_shape_reject, debug_physical_reject,
                            debug_external_reject, pools[fi].size());
                    if (family[fi] == "5")
                        fprintf(stderr,
                                "[STRICT-FAMILY-PHYSICAL] %s member=5 boundary=%d fence=%d obstacle=%d\n",
                                family_item.first.c_str(), debug_physical_boundary,
                                debug_physical_fence, debug_physical_obstacle);
                }
            }
            if (isgProfile() && family_has_current_curve_cross) {
                fprintf(stderr,
                        "[ISG_PROFILE] strict family %s candidate pools %.3f ms\n",
                        family_item.first.c_str(),
                        std::chrono::duration<double, std::milli>(
                            std::chrono::steady_clock::now() - strict_pool_t0).count());
            }
            if (!pools_valid) {
                if (isgDebugPairRepair()) {
                    fprintf(stderr, "[STRICT-FAMILY] %s pool build failed sizes=",
                            family_item.first.c_str());
                    for (const auto& pool : pools)
                        fprintf(stderr, "%zu,", pool.size());
                    fprintf(stderr, "\n");
                }
                continue;
            }

            // 先在未截断的候选集里寻找兼容对，再按固定上限收口。若两边各自
            // 只保留离当前曲线最近的候选，唯一可行的远端组合会被两次独立
            // 截断同时丢掉；100000443 的 5|9 与 41|40 正是这种情况。二元
            // 共享族只保护前若干个低代价兼容对，仍保持每成员最多 24 个候选，
            // 不把全量笛卡尔积带入后续回溯。
            for (auto& pool : pools)
                std::sort(pool.begin(), pool.end(),
                          [](const StrictOrdinaryCandidate& a,
                             const StrictOrdinaryCandidate& b) {
                              return a.score < b.score;
                          });
            if (family.size() == 2 &&
                (pools[0].size() > kMaxStrictPoolCandidates ||
                 pools[1].size() > kMaxStrictPoolCandidates)) {
                struct CompatiblePair {
                    std::size_t first = 0;
                    std::size_t second = 0;
                    double score = std::numeric_limits<double>::infinity();
                };
                std::vector<CompatiblePair> compatible_pairs;
                // 两个候选池都已按单成员代价排序。用小根堆按组合代价枚举
                // 前若干个组合，并先用缓存采样/控制折线筛选；最终组合仍在
                // search() 叶节点执行完整曲线精确复核，采样只负责提速。
                struct PairQueueItem {
                    std::size_t first = 0;
                    std::size_t second = 0;
                    double score = std::numeric_limits<double>::infinity();
                };
                struct PairQueueGreater {
                    bool operator()(const PairQueueItem& a,
                                    const PairQueueItem& b) const {
                        return a.score > b.score;
                    }
                };
                std::priority_queue<PairQueueItem,
                                    std::vector<PairQueueItem>,
                                    PairQueueGreater> queue;
                for (std::size_t i = 0; i < pools[0].size(); ++i) {
                    PairQueueItem item;
                    item.first = i;
                    item.score = pools[0][i].score + pools[1].front().score;
                    queue.push(item);
                }
                constexpr std::size_t kPairProbeLimit = 8192;
                constexpr std::size_t kCompatiblePairLimit = 64;
                std::size_t probes = 0;
                while (!queue.empty() && probes++ < kPairProbeLimit &&
                       compatible_pairs.size() < kCompatiblePairLimit) {
                    const PairQueueItem item = queue.top();
                    queue.pop();
                    const auto& a = pools[0][item.first];
                    const auto& b = pools[1][item.second];
                    if (!sharedEndpointControlPolylinesCross(
                            a.curve, b.curve, kClusterEndpointTol) &&
                        !sampledCurvesIntersectBusiness(
                            a.sampled, b.sampled,
                            kOrdinaryFamilyCurveCrossTol)) {
                        CompatiblePair pair;
                        pair.first = item.first;
                        pair.second = item.second;
                        pair.score = item.score;
                        compatible_pairs.push_back(pair);
                    }
                    if (item.second + 1 < pools[1].size()) {
                        const std::size_t next = item.second + 1;
                        PairQueueItem next_item;
                        next_item.first = item.first;
                        next_item.second = next;
                        next_item.score = pools[0][item.first].score +
                            pools[1][next].score;
                        queue.push(next_item);
                    }
                }
                std::sort(compatible_pairs.begin(), compatible_pairs.end(),
                          [](const CompatiblePair& a, const CompatiblePair& b) {
                              return a.score < b.score;
                          });
                std::vector<std::unordered_set<std::size_t>> protected_indices(2);
                const std::size_t protected_pair_count = std::min<std::size_t>(
                    8, compatible_pairs.size());
                for (std::size_t i = 0; i < protected_pair_count; ++i) {
                    protected_indices[0].insert(compatible_pairs[i].first);
                    protected_indices[1].insert(compatible_pairs[i].second);
                }
                for (std::size_t fi = 0; fi < pools.size(); ++fi) {
                    std::vector<StrictOrdinaryCandidate> reduced;
                    reduced.reserve(std::min<std::size_t>(
                        kMaxStrictPoolCandidates, pools[fi].size()));
                    for (std::size_t ci = 0; ci < pools[fi].size() &&
                         reduced.size() < kMaxStrictPoolCandidates; ++ci) {
                        if (protected_indices[fi].count(ci))
                            reduced.push_back(std::move(pools[fi][ci]));
                    }
                    for (std::size_t ci = 0; ci < pools[fi].size() &&
                         reduced.size() < kMaxStrictPoolCandidates; ++ci) {
                        if (!protected_indices[fi].count(ci))
                            reduced.push_back(std::move(pools[fi][ci]));
                    }
                    pools[fi].swap(reduced);
                }
                if (isgDebugPairRepair() && !compatible_pairs.empty()) {
                    fprintf(stderr,
                            "[STRICT-FAMILY] %s compatible_pairs=%zu protected=%zu probes=%zu\n",
                            family_item.first.c_str(), compatible_pairs.size(),
                            protected_pair_count, probes);
                }
            } else {
                for (auto& pool : pools)
                    if (pool.size() > kMaxStrictPoolCandidates)
                        pool.resize(kMaxStrictPoolCandidates);
            }

            if (isgDebugPairRepair()) {
                fprintf(stderr, "[STRICT-FAMILY] %s members=%zu pools=",
                        family_item.first.c_str(), family.size());
                for (const auto& pool : pools)
                    fprintf(stderr, "%zu,", pool.size());
                fprintf(stderr, "\n");
                for (std::size_t fi = 0; fi < pools.size(); ++fi) {
                    fprintf(stderr, "[STRICT-FAMILY] %s pool[%zu]=",
                            family_item.first.c_str(), fi);
                    for (const auto& candidate : pools[fi]) {
                        const auto& g = candidate.curve.segs.front().ctrl;
                        fprintf(stderr, "%.2f/%.2f,",
                                (g[1] - g[0]).norm(), (g[3] - g[2]).norm());
                    }
                    fprintf(stderr, "\n");
                }
            }

            struct StrictExternalIssue {
                int weighted = 0;
                bool new_fixed_cross = false;
                bool new_curve_cross = false;
                bool new_control_cross = false;
                ConnId blocker;
            };
            std::unordered_map<ConnId, std::size_t> family_index;
            for (std::size_t fi = 0; fi < family.size(); ++fi)
                family_index.emplace(family[fi], fi);

            std::vector<std::vector<StrictExternalIssue>> external_issues(
                family.size());
            for (std::size_t fi = 0; fi < family.size(); ++fi)
                external_issues[fi].resize(pools[fi].size());

            // 族内候选必须彻底无交叉；族外关系则按“替换后总问题数”计分。
            // 单成员建池阶段不能把新增外部问题直接删掉，因为外部连接可能
            // 也在另一共享端点族中需要同步改形；在这里预计算后做增量搜索。
            int current_strict_issues = 0;
            for (std::size_t i = 0; i < family.size(); ++i) {
                auto ia = final_idx.find(family[i]);
                if (ia == final_idx.end() || !results[ia->second].curve)
                    continue;
                for (std::size_t j = i + 1; j < family.size(); ++j) {
                    if (!strict_family_pair_required(family[i], family[j]))
                        continue;
                    auto ib = final_idx.find(family[j]);
                    if (ib == final_idx.end() || !results[ib->second].curve)
                        continue;
                    const bool current_cross = sampledCurvesIntersectBusiness(
                        strict_result_samples.at(family[i]),
                        strict_result_samples.at(family[j]),
                        kOrdinaryFamilyCurveCrossTol);
                    current_strict_issues += current_cross ? 1000 :
                        (sharedEndpointControlPolylinesCross(
                            *results[ia->second].curve,
                            *results[ib->second].curve,
                            kClusterEndpointTol) ? 1 : 0);
                }
            }

            for (const auto& pair : cluster_solver_.pairs()) {
                if (pair.exempt == CrossExemption::StructuralCross)
                    continue;
                auto fa = family_index.find(pair.id_a);
                auto fb = family_index.find(pair.id_b);
                const bool a_in_family = fa != family_index.end();
                const bool b_in_family = fb != family_index.end();
                if (a_in_family == b_in_family)
                    continue;
                const std::size_t fi = a_in_family ? fa->second : fb->second;
                const ConnId& other_id = a_in_family ? pair.id_b : pair.id_a;
                auto oi = final_idx.find(other_id);
                if (oi == final_idx.end() || !results[oi->second].curve)
                    continue;
                const BezierCurve& other = *results[oi->second].curve;
                const ConnectivityCurve& other_result = results[oi->second];
                const Connectivity* member_conn = scene.view.connectivity(family[fi]);
                if (!member_conn)
                    continue;
                const auto member_result = final_idx.find(family[fi]);
                if (member_result == final_idx.end() ||
                    !results[member_result->second].curve)
                    continue;
                const BezierCurve& current = *results[member_result->second].curve;
                const bool current_cross = sampledCurvesIntersectBusiness(
                    strict_result_samples.at(family[fi]),
                    strict_result_samples.at(other_id),
                    kOrdinaryFamilyCurveCrossTol);
                const bool current_control_cross = pair.shared_endpoint &&
                    sharedEndpointControlPolylinesCross(
                        current, other, kClusterEndpointTol);
                current_strict_issues += current_cross ? 1000 :
                    (current_control_cross ? 1 : 0);
                for (std::size_t ci = 0; ci < pools[fi].size(); ++ci) {
                    const BezierCurve& candidate = pools[fi][ci].curve;
                    const bool candidate_cross = sampledCurvesIntersectBusiness(
                        pools[fi][ci].sampled,
                        strict_result_samples.at(other_id),
                        kOrdinaryFamilyCurveCrossTol);
                    const bool candidate_control_cross = pair.shared_endpoint &&
                        sharedEndpointControlPolylinesCross(
                            candidate, other, kClusterEndpointTol);
                    external_issues[fi][ci].weighted += candidate_cross
                        ? 1000 : (candidate_control_cross ? 1 : 0);
                    external_issues[fi][ci].new_curve_cross =
                        external_issues[fi][ci].new_curve_cross ||
                        (!current_cross && candidate_cross);
                    external_issues[fi][ci].new_control_cross =
                        external_issues[fi][ci].new_control_cross ||
                        (!current_control_cross && candidate_control_cross);
                    if (external_issues[fi][ci].blocker.empty() &&
                        (candidate_cross || candidate_control_cross))
                        external_issues[fi][ci].blocker = other_id;
                    if (!current_cross && candidate_cross &&
                        (other_result.fixed_shape || preserved_fixed_ids.count(other_id))) {
                        external_issues[fi][ci].new_fixed_cross = true;
                    }
                }
            }

            if (isgDebugPairRepair() &&
                (family_item.first == "entry:1015689" ||
                 family_item.first == "exit:1015417" ||
                 family_item.first == "exit:1015413" ||
                 family_item.first == "exit:1015421" ||
                 family_item.first == "exit:1015745")) {
                fprintf(stderr,
                        "[STRICT-FAMILY-EXTERNAL] %s current_issues=%d\n",
                        family_item.first.c_str(), current_strict_issues);
                for (std::size_t fi = 0; fi < external_issues.size(); ++fi) {
                    for (std::size_t ci = 0; ci < external_issues[fi].size(); ++ci) {
                        const auto& issue = external_issues[fi][ci];
                        const auto& g = pools[fi][ci].curve.segs.front().ctrl;
                        fprintf(stderr,
                                "[STRICT-FAMILY-EXTERNAL] %s member=%s ci=%zu h=%.3f/%.3f weight=%d fixed=%d new_curve=%d new_poly=%d blocker=%s\n",
                                family_item.first.c_str(), family[fi].c_str(), ci,
                                (g[1] - g[0]).norm(), (g[3] - g[2]).norm(),
                                issue.weighted, issue.new_fixed_cross ? 1 : 0,
                                issue.new_curve_cross ? 1 : 0,
                                issue.new_control_cross ? 1 : 0,
                                issue.blocker.c_str());
                    }
                }
            }

            const auto strict_matrix_t0 = std::chrono::steady_clock::now();
            std::vector<std::vector<std::vector<unsigned char>>> internal_conflicts(
                family.size(), std::vector<std::vector<unsigned char>>(family.size()));
            for (std::size_t i = 0; i < family.size(); ++i) {
                for (std::size_t j = i + 1; j < family.size(); ++j) {
                    auto& matrix = internal_conflicts[i][j];
                    matrix.resize(pools[i].size() * pools[j].size(), 0);
                    for (std::size_t ci = 0; ci < pools[i].size(); ++ci) {
                        for (std::size_t cj = 0; cj < pools[j].size(); ++cj) {
                            const BezierCurve& a = pools[i][ci].curve;
                            const BezierCurve& b = pools[j][cj].curve;
                            matrix[ci * pools[j].size() + cj] =
                                (sampledCurvesIntersectBusiness(
                                     pools[i][ci].sampled,
                                     pools[j][cj].sampled,
                                     kOrdinaryFamilyCurveCrossTol) ||
                                 sharedEndpointControlPolylinesCross(
                                     a, b, kClusterEndpointTol)) ? 1 : 0;
                        }
                    }
                }
            }
            if (isgProfile() && family_has_current_curve_cross) {
                fprintf(stderr,
                        "[ISG_PROFILE] strict family %s internal matrix %.3f ms\n",
                        family_item.first.c_str(),
                        std::chrono::duration<double, std::milli>(
                            std::chrono::steady_clock::now() - strict_matrix_t0).count());
            }

            std::vector<const BezierCurve*> selected(family.size(), nullptr);
            std::vector<std::size_t> selected_index(family.size(), 0);
            std::vector<std::size_t> search_order;
            search_order.reserve(family.size());
            for (std::size_t fi = 0; fi < family.size(); ++fi)
                search_order.push_back(fi);
            std::sort(search_order.begin(), search_order.end(),
                      [&](std::size_t a, std::size_t b) {
                          if (pools[a].size() != pools[b].size())
                              return pools[a].size() < pools[b].size();
                          return pools[a].front().score < pools[b].front().score;
                      });

            auto exact_combination_is_clean = [&]() {
                // 采样矩阵用于搜索剪枝，但不能成为最终准入依据。对已选
                // 族内组合逐对执行与输出审计相同的精确业务判定，并检查
                // 共享端点控制折线前兆。
                for (std::size_t i = 0; i < family.size(); ++i) {
                    for (std::size_t j = i + 1; j < family.size(); ++j) {
                        if (!strict_family_pair_required(family[i], family[j]))
                            continue;
                        const BezierCurve& a = pools[i][selected_index[i]].curve;
                        const BezierCurve& b = pools[j][selected_index[j]].curve;
                        if (curvesIntersectBusiness(
                                a, b, kOrdinaryFamilyCurveCrossTol) ||
                            sharedEndpointControlPolylinesCross(
                                a, b, kClusterEndpointTol))
                            return false;
                    }
                }
                for (const auto& pair : cluster_solver_.pairs()) {
                    if (pair.exempt == CrossExemption::StructuralCross)
                        continue;
                    auto fa = family_index.find(pair.id_a);
                    auto fb = family_index.find(pair.id_b);
                    const bool a_in_family = fa != family_index.end();
                    const bool b_in_family = fb != family_index.end();
                    if (a_in_family == b_in_family)
                        continue;
                    const std::size_t fi = a_in_family ? fa->second : fb->second;
                    const ConnId& other_id = a_in_family ? pair.id_b : pair.id_a;
                    auto oi = final_idx.find(other_id);
                    if (oi == final_idx.end() || !results[oi->second].curve)
                        return false;
                    const BezierCurve& candidate =
                        pools[fi][selected_index[fi]].curve;
                    const BezierCurve& other = *results[oi->second].curve;
                    if (curvesIntersectBusiness(
                            candidate, other, kOrdinaryFamilyCurveCrossTol) ||
                        (pair.shared_endpoint &&
                         sharedEndpointControlPolylinesCross(
                             candidate, other, kClusterEndpointTol)))
                        return false;
                }
                return true;
            };

            std::size_t search_nodes = 0;
            const std::size_t max_nodes = 12000;
            std::function<bool(std::size_t, int)> search =
                [&](std::size_t depth, int partial_issues) {
                if (depth == family.size()) {
                    if (partial_issues >= current_strict_issues)
                        return false;
                    // 候选池使用较稀疏采样快速筛选；最终入选组合再用
                    // 完整 Boundary 和同簇几何复核，避免搜索加速造成新的物理
                    // 或相交违约。
                    for (std::size_t fi = 0; fi < family.size(); ++fi)
                        if (strictCandidateHasExactPhysicalRisk(
                                pools[fi][selected_index[fi]].curve))
                            return false;
                    return exact_combination_is_clean();
                }
                if (++search_nodes > max_nodes)
                    return false;
                const std::size_t fi = search_order[depth];
                for (std::size_t ci = 0; ci < pools[fi].size(); ++ci) {
                    const BezierCurve& candidate = pools[fi][ci].curve;
                    if (external_issues[fi][ci].new_fixed_cross ||
                        (!family_has_current_curve_cross &&
                         (external_issues[fi][ci].new_curve_cross ||
                          external_issues[fi][ci].new_control_cross)))
                        continue;
                    const int candidate_issues = partial_issues +
                        external_issues[fi][ci].weighted;
                    if (candidate_issues >= current_strict_issues)
                        continue;
                    bool conflict = false;
                    for (std::size_t previous = 0; previous < depth; ++previous) {
                        const std::size_t previous_fi = search_order[previous];
                        const std::size_t lo = std::min(fi, previous_fi);
                        const std::size_t hi = std::max(fi, previous_fi);
                        const std::size_t lo_ci = fi == lo ? ci : selected_index[previous_fi];
                        const std::size_t hi_ci = fi == hi ? ci : selected_index[previous_fi];
                        if (internal_conflicts[lo][hi][
                                lo_ci * pools[hi].size() + hi_ci]) {
                            if (isgDebugPairRepair() &&
                                (family_item.first == "entry:1015463" ||
                                 family_item.first == "entry:1015689" ||
                                 family_item.first == "exit:1015745" ||
                                 family_item.first == "exit:1015413")) {
                                const auto& previous_curve =
                                    *selected[previous_fi];
                                const bool curve_conflict = curvesIntersectBusiness(
                                    candidate, previous_curve, kClusterEndpointTol);
                                const bool control_conflict =
                                    sharedEndpointControlPolylinesCross(
                                        candidate, previous_curve,
                                        kClusterEndpointTol);
                                fprintf(stderr,
                                        "[STRICT-FAMILY-CONFLICT] %s depth=%zu/%zu ci=%zu prev=%zu curve=%d poly=%d h=%.3f/%.3f prev_h=%.3f/%.3f\n",
                                        family_item.first.c_str(), depth,
                                        family.size(), ci, previous_fi,
                                        curve_conflict ? 1 : 0,
                                        control_conflict ? 1 : 0,
                                        (candidate.segs.front().ctrl[1] -
                                         candidate.segs.front().ctrl[0]).norm(),
                                        (candidate.segs.front().ctrl[3] -
                                         candidate.segs.front().ctrl[2]).norm(),
                                        (previous_curve.segs.front().ctrl[1] -
                                         previous_curve.segs.front().ctrl[0]).norm(),
                                        (previous_curve.segs.front().ctrl[3] -
                                         previous_curve.segs.front().ctrl[2]).norm());
                            }
                            conflict = true;
                            break;
                        }
                    }
                    if (conflict)
                        continue;
                    selected[fi] = &candidate;
                    selected_index[fi] = ci;
                    if (search(depth + 1, candidate_issues))
                        return true;
                    selected[fi] = nullptr;
                }
                return false;
            };
            const auto strict_search_t0 = std::chrono::steady_clock::now();
            if (!search(0, 0)) {
                if (isgProfile() && family_has_current_curve_cross) {
                    fprintf(stderr,
                            "[ISG_PROFILE] strict family %s narrow search %.3f ms\n",
                            family_item.first.c_str(),
                            std::chrono::duration<double, std::milli>(
                                std::chrono::steady_clock::now() - strict_search_t0).count());
                }
                // 当前族的窄池如果没有解，先对原始目标族使用宽池重试。
                // 外部曲线仍会作为固定约束进入闭包矩阵；只有某个外部普通
                // 连接挡住变量的全部候选时，才有必要把它升级为闭包变量。
                std::vector<ConnId> closure_ids = family;
                std::unordered_set<ConnId> closure_id_set(
                    closure_ids.begin(), closure_ids.end());
                // 只在后面的“全部候选均被挡住”判定成立时加入外部连接。
                if (isgDebugPairRepair()) {
                    fprintf(stderr, "[STRICT-CLOSURE] %s ids=",
                            family_item.first.c_str());
                    for (const ConnId& id : closure_ids)
                        fprintf(stderr, "%s,", id.c_str());
                    fprintf(stderr, "\n");
                }
                bool closure_pool_ok = closure_ids.size() >= family.size();
                std::unordered_map<ConnId,
                    std::vector<StrictOrdinaryCandidate>> closure_pools;
                // 闭包搜索与普通族收口使用同一候选生成规则，但不能复用普通族
                // 截断后的 24 个近当前候选；跨端点族联动时可行解往往需要让
                // 某个外部阻塞者明显偏离当前形态。只在实际触发闭包时扩大到
                // 小闭包最多使用 64 个代表候选；更大的闭包收敛到 32 个，且仍
                // 受完整物理、形态和同簇矩阵约束。
                // 共享入口扇出族的成员数本身就可能达到 6（110004764 的
                // 43106530 有 8,10,11,12,13,14）。此时闭包上限若仍是 6，
                // 扩张循环的 closure_ids.size() < kMaxStrictClosureMembers
                // 直接为假，跨端点阻塞者（41、47）永远进不了变量集，而族内
                // 搜索的叶子判定要求对外部关系也彻底干净，于是整族必然失败。
                // 上限放宽到 10 只在升级重试轮生效，即只影响首轮已经失败的
                // 族，且闭包矩阵仍有 32 候选/成员的截断和 50000 的节点预算兜底。
                const std::size_t kMaxStrictClosureMembers =
                    closure_escalate ? 10u : 6u;
                const std::size_t closure_pool_limit =
                    closure_ids.size() <= 2 ? kMaxStrictClosurePoolCandidates : 32;
                auto is_movable_closure_member = [&](const ConnId& id) {
                    const Connectivity* conn = scene.view.connectivity(id);
                    auto ri = final_idx.find(id);
                    return conn && ri != final_idx.end() &&
                        results[ri->second].curve && !conn->fixed_shape &&
                        !preserved_fixed_ids.count(id) &&
                        !isGeometricUTurnConn(*conn, input) &&
                        results[ri->second].curve->numSegments() == 1;
                };
                std::unordered_map<ConnId, std::vector<const CurvePair*>>
                    strict_pair_neighbors;
                for (const auto& pair : cluster_solver_.pairs()) {
                    if (pair.exempt == CrossExemption::StructuralCross)
                        continue;
                    strict_pair_neighbors[pair.id_a].push_back(&pair);
                    strict_pair_neighbors[pair.id_b].push_back(&pair);
                }
                auto build_missing_closure_pools = [&]() {
                    for (std::size_t ci = 0;
                         ci < closure_ids.size() && closure_pool_ok; ++ci) {
                        if (closure_pools.count(closure_ids[ci]))
                            continue;
                        // 原始二元族已经完成宽池物理审计，直接复用，避免闭包
                        // 第一次建池再次逐候选扫描连续 Boundary 链。
                        auto wide_it = strict_wide_pools.find(closure_ids[ci]);
                        if (wide_it != strict_wide_pools.end()) {
                            closure_pools.emplace(
                                closure_ids[ci], wide_it->second);
                            continue;
                        }
                        std::vector<StrictOrdinaryCandidate> pool;
                        // 原始二元族的可行解可能位于单成员代价排序的远端，
                        // 即使闭包新增了阻塞者，也不能把这两个原始变量重新
                        // 截回 32 个候选；新增闭包变量才使用窄池。
                        const bool original_binary_member =
                            family.size() == 2 &&
                            family_index.count(closure_ids[ci]) != 0;
                        // 新增的跨端点阻塞者不能只取代价排序前 64 个：
                        // 5|9 的可行 5 候选可能要求 38 远离当前形态，
                        // 其支撑值常落在完整把手网格的后段。该网格最多
                        // 22*22 个组合，仍受物理/形态门禁和 AC3 约束，
                        // 不会退化为全路口穷举。
                        const bool expanded_closure_member =
                            !original_binary_member &&
                            closure_ids.size() > family.size();
                        const std::size_t member_pool_limit =
                            original_binary_member
                                ? kMaxStrictPairPoolCandidates
                                : expanded_closure_member
                                    ? 128
                                : closure_pool_limit;
                        closure_pool_ok = build_strict_closure_pool(
                            closure_ids[ci], pool, member_pool_limit);
                        if (closure_pool_ok)
                            closure_pools.emplace(
                                closure_ids[ci], std::move(pool));
                    }
                };
                build_missing_closure_pools();

                // 单段候选域完全没有相容支撑时，长距离近直行连接允许一次
                // 有界的路点回退。该回退只扩充真正无支撑关系的一侧，避免把
                // 短转向和大闭包全部改成多段；候选仍进入同一闭包矩阵和 AC3，
                // 因而不会以新增相交、贴合或重叠换取目标对的消除。
                auto append_strict_waypoint_candidates = [&](const ConnId& id) {
                    if (!closure_pool_ok || !is_movable_closure_member(id))
                        return false;
                    const Connectivity* conn = scene.view.connectivity(id);
                    if (!conn)
                        return false;
                    const auto entry = scene.view.entryFrame(conn->entry_lane_id);
                    const auto exit_ = scene.view.exitFrame(conn->exit_lane_id);
                    const Vec2d chord = exit_.first - entry.first;
                    if (chord.norm() <= 15.0 || entry.second.norm() < 1e-8 ||
                        exit_.second.norm() < 1e-8)
                        return false;
                    const Vec2d t0 = entry.second.normalized();
                    const Vec2d t1 = exit_.second.normalized();
                    const double turn_strength = std::abs(
                        cross2d(t0, chord.normalized()));
                    if (turn_strength > 0.25)
                        return false;

                    const Vec2d perpendicular(
                        -chord.y() / chord.norm(), chord.x() / chord.norm());
                    const std::vector<double> fractions = {
                        0.18, 0.22, 0.26, 0.30, 0.35, 0.40, 0.46,
                        0.54, 0.60, 0.65, 0.70, 0.76, 0.82};
                    const std::vector<double> offsets = {
                        -3.0, -2.0, -1.5, -1.0, -0.5,
                        0.5, 1.0, 1.5, 2.0, 2.5, 3.0};
                    const std::vector<double> alphas = {0.18, 0.24, 0.30};
                    std::vector<StrictOrdinaryCandidate> pending;
                    for (double fraction : fractions) {
                        const Vec2d base = entry.first + fraction * chord;
                        for (double offset : offsets) {
                            const Vec2d waypoint = base + offset * perpendicular;
                            if ((waypoint - entry.first).norm() < 0.30 ||
                                (exit_.first - waypoint).norm() < 0.30)
                                continue;
                            std::vector<Vec2d> mid_tangents;
                            const Vec2d to_exit = exit_.first - waypoint;
                            const Vec2d from_entry = waypoint - entry.first;
                            if (to_exit.norm() > 1e-8)
                                mid_tangents.push_back(to_exit.normalized());
                            if (from_entry.norm() > 1e-8 && to_exit.norm() > 1e-8) {
                                const Vec2d bisector =
                                    from_entry.normalized() + to_exit.normalized();
                                if (bisector.norm() > 1e-8)
                                    mid_tangents.push_back(bisector.normalized());
                            }
                            for (const Vec2d& mid_tangent : mid_tangents) {
                                for (double alpha : alphas) {
                                    BezierCurve candidate = makeCurveFromKnots(
                                        {entry.first, waypoint, exit_.first},
                                        {t0, mid_tangent, t1}, alpha);
                                    if (!isTwoSegmentOrdinaryWaypointShapeAcceptable(
                                            candidate, entry.first, t0,
                                            exit_.first, t1, chord.norm(),
                                            turn_strength) ||
                                        curveSelfIntersectsBusiness(candidate, 1.0) ||
                                        strictCandidateHasExactPhysicalRisk(candidate))
                                        continue;
                                    StrictOrdinaryCandidate item;
                                    item.curve = std::move(candidate);
                                    item.sampled = sampleStrictPairCurve(item.curve);
                                    item.score =
                                        std::abs(offset) +
                                        0.50 * std::abs(
                                            item.curve.arcLength() / chord.norm() - 1.03) +
                                        item.curve.maxCurvature(40) +
                                        0.005 * item.curve.arcLength();
                                    pending.push_back(std::move(item));
                                }
                            }
                        }
                    }
                    std::sort(pending.begin(), pending.end(),
                              [](const StrictOrdinaryCandidate& a,
                                 const StrictOrdinaryCandidate& b) {
                                  return a.score < b.score;
                              });
                    auto& pool = closure_pools[id];
                    constexpr std::size_t kMaxWaypointCandidates = 96;
                    for (auto& item : pending) {
                        if (pool.size() >= 512 + kMaxWaypointCandidates)
                            break;
                        const std::string key = strict_curve_key(item.curve);
                        bool duplicate = false;
                        for (const auto& existing : pool) {
                            if (strict_curve_key(existing.curve) == key) {
                                duplicate = true;
                                break;
                            }
                        }
                        if (!duplicate)
                            pool.push_back(std::move(item));
                    }
                    if (isgDebugPairRepair())
                        fprintf(stderr,
                                "[STRICT-WAYPOINT-POOL] id=%s added=%zu total=%zu\n",
                                id.c_str(), pool.size() -
                                    std::min<std::size_t>(pool.size(), 512),
                                pool.size());
                    return pool.size() > 512;
                };

                std::unordered_set<ConnId> waypoint_targets;
                const auto closure_contains = [&](const ConnId& id) {
                    return closure_id_set.count(id) != 0;
                };
                const auto choose_waypoint_target = [&](const ConnId& a,
                                                        const ConnId& b) {
                    const auto is_eligible = [&](const ConnId& id) {
                        if (!is_movable_closure_member(id))
                            return false;
                        const Connectivity* conn = scene.view.connectivity(id);
                        if (!conn)
                            return false;
                        const auto entry = scene.view.entryFrame(conn->entry_lane_id);
                        const auto exit_ = scene.view.exitFrame(conn->exit_lane_id);
                        const Vec2d chord = exit_.first - entry.first;
                        return chord.norm() > 15.0 && entry.second.norm() > 1e-8 &&
                            std::abs(cross2d(entry.second.normalized(),
                                              chord.normalized())) <= 0.25;
                    };
                    if (is_eligible(a))
                        return a;
                    if (is_eligible(b))
                        return b;
                    return ConnId();
                };
                for (const auto& pair : cluster_solver_.pairs()) {
                    if (pair.exempt == CrossExemption::StructuralCross)
                        continue;
                    const bool a_variable = closure_contains(pair.id_a);
                    const bool b_variable = closure_contains(pair.id_b);
                    if (!a_variable && !b_variable)
                        continue;
                    ConnId target;
                    if (a_variable && b_variable) {
                        const auto& pool_a = closure_pools[pair.id_a];
                        const auto& pool_b = closure_pools[pair.id_b];
                        bool supported = false;
                        for (const auto& candidate_a : pool_a) {
                            for (const auto& candidate_b : pool_b) {
                                if (!sampledCurvesIntersectBusiness(
                                        candidate_a.sampled, candidate_b.sampled,
                                        kSharedEndpointPairCrossTol)) {
                                    supported = true;
                                    break;
                                }
                            }
                            if (supported)
                                break;
                        }
                        if (!supported)
                            target = choose_waypoint_target(pair.id_a, pair.id_b);
                    } else {
                        const ConnId variable = a_variable ? pair.id_a : pair.id_b;
                        const ConnId fixed = a_variable ? pair.id_b : pair.id_a;
                        auto fixed_it = strict_result_samples.find(fixed);
                        if (fixed_it == strict_result_samples.end())
                            continue;
                        bool supported = false;
                        for (const auto& candidate : closure_pools[variable]) {
                            if (!sampledCurvesIntersectBusiness(
                                    candidate.sampled, fixed_it->second,
                                    kSharedEndpointPairCrossTol)) {
                                supported = true;
                                break;
                            }
                        }
                        if (!supported)
                            target = choose_waypoint_target(variable, variable);
                    }
                    if (!target.empty())
                        waypoint_targets.insert(target);
                }
                for (const ConnId& id : waypoint_targets)
                    append_strict_waypoint_candidates(id);

                if (isgDebugPairRepair() &&
                    family_item.first == "entry:1015689") {
                    fprintf(stderr, "[STRICT-CLOSURE-POOLS] ok=%d sizes=",
                            closure_pool_ok ? 1 : 0);
                    for (const ConnId& id : closure_ids) {
                        auto it = closure_pools.find(id);
                        fprintf(stderr, "%s:%zu,", id.c_str(),
                                it == closure_pools.end() ? 0 : it->second.size());
                    }
                    fprintf(stderr, "\n");
                    for (const ConnId& id : closure_ids) {
                        if (id != "5" && id != "9")
                            continue;
                        const auto it = closure_pools.find(id);
                        if (it == closure_pools.end())
                            continue;
                        fprintf(stderr, "[STRICT-CLOSURE-POOL-VIEW] id=%s handles=",
                                id.c_str());
                        for (const auto& candidate : it->second) {
                            const auto& g = candidate.curve.segs.front().ctrl;
                            fprintf(stderr, "%.3f/%.3f,",
                                    (g[1] - g[0]).norm(),
                                    (g[3] - g[2]).norm());
                        }
                        fprintf(stderr, "\n");
                    }
                }

                // 直接阻塞者自身可能在闭包外还有可变普通兄弟。只沿“某个
                // 候选实际相交”的边扩展一轮，并限制变量总数；固定形态、
                // U-turn 和不可变连接留在后续矩阵中作为硬约束。这样既能
                // 覆盖 5/9 的分裂阻塞，又不会把全路口关系图展开成全局穷举。
                // 大族往往被多个跨端点阻塞者同时挡住（43106530 的候选分别
                // 被 41 和 47 挡住），每轮只纳入一个新变量无法解锁。轮数与
                // 每轮发现数只在升级重试轮放宽，仍受 kMaxStrictClosureMembers
                // 与节点预算约束。
                const std::size_t closure_discovery_per_round =
                    (closure_escalate && family.size() >= 4) ? 2u : 1u;
                const int closure_max_rounds =
                    (closure_escalate && family.size() >= 4) ? 3 : 1;
                for (int closure_round = 0;
                     closure_pool_ok && closure_round < closure_max_rounds &&
                     closure_ids.size() < kMaxStrictClosureMembers;
                     ++closure_round) {
                    std::vector<ConnId> discovered;
                    std::unordered_set<ConnId> discovered_set;
                    for (const ConnId& variable_id : closure_ids) {
                        if (discovered.size() >= closure_discovery_per_round)
                            break;
                        const auto pool_it = closure_pools.find(variable_id);
                        if (pool_it == closure_pools.end())
                            continue;
                        const auto neighbor_it =
                            strict_pair_neighbors.find(variable_id);
                        if (neighbor_it == strict_pair_neighbors.end())
                            continue;
                        for (const auto& candidate : pool_it->second) {
                            if (discovered.size() >= closure_discovery_per_round)
                                break;
                            for (const CurvePair* pair_ptr : neighbor_it->second) {
                                if (discovered.size() >= closure_discovery_per_round)
                                    break;
                                const CurvePair& pair = *pair_ptr;
                                ConnId other_id;
                                if (pair.id_a == variable_id)
                                    other_id = pair.id_b;
                                else if (pair.id_b == variable_id)
                                    other_id = pair.id_a;
                                else
                                    continue;
                                if (closure_id_set.count(other_id) ||
                                    discovered_set.count(other_id) ||
                                    pair.exempt == CrossExemption::StructuralCross ||
                                    !is_movable_closure_member(other_id))
                                    continue;
                                auto oi = final_idx.find(other_id);
                                if (oi == final_idx.end() ||
                                    !results[oi->second].curve)
                                    continue;
                                const auto other_sample =
                                    strict_result_samples.find(other_id);
                                if (other_sample == strict_result_samples.end() ||
                                    !sampledCurvesIntersectBusiness(
                                        candidate.sampled, other_sample->second,
                                        kSharedEndpointPairCrossTol))
                                    continue;
                                // 族内候选可能分别被不同外部兄弟阻塞：没有任意
                                // 一条兄弟会独占挡住全部候选，但这些阻塞的并集
                                // 仍会清空原始族的可行组合。只要存在实际冲突候选，
                                // 就把可变普通兄弟纳入有限闭包，由矩阵统一求解。
                                if (discovered_set.insert(other_id).second)
                                    discovered.push_back(other_id);
                                break;
                            }
                        }
                    }
                    for (const ConnId& id : discovered) {
                        if (closure_ids.size() >= kMaxStrictClosureMembers)
                            break;
                        if (closure_id_set.insert(id).second)
                            closure_ids.push_back(id);
                    }
                    if (discovered.empty())
                        break;
                    build_missing_closure_pools();
                }

                struct ClosurePair {
                    std::size_t a = 0;
                    std::size_t b = 0;
                    ConnId fixed_id;
                    std::vector<unsigned char> conflict;
                };
                std::vector<ClosurePair> closure_pairs;
                std::unordered_map<ConnId, std::size_t> closure_index;
                for (std::size_t ci = 0; ci < closure_ids.size(); ++ci)
                    closure_index[closure_ids[ci]] = ci;
                if (closure_pool_ok) {
                    const auto closure_forbidden =
                        [&](const BezierCurve& curve_a, const SampledCurve& a,
                            const BezierCurve& curve_b, const SampledCurve& b) {
                            // 闭包的组合硬门禁只判断真实曲线在连接点外的
                            // 相交/接触/重叠。控制折线互穿是候选排序前兆，
                            // 不能在跨端点族的事务中当作真实几何无解；否则
                            // 合法的扇出顺序会被前兆门禁提前清空。
                            const double endpoint_tol =
                                (curve_a.numSegments() > 1 ||
                                 curve_b.numSegments() > 1)
                                ? kStrictSameClusterEndpointTol
                                : kSharedEndpointPairCrossTol;
                            return sampledCurvesIntersectBusiness(
                                       a, b, endpoint_tol);
                        };
                    for (const auto& pair : cluster_solver_.pairs()) {
                        if (pair.exempt == CrossExemption::StructuralCross)
                            continue;
                        const auto ia = closure_index.find(pair.id_a);
                        const auto ib = closure_index.find(pair.id_b);
                        const bool a_variable = ia != closure_index.end();
                        const bool b_variable = ib != closure_index.end();
                        if (!a_variable && !b_variable)
                            continue;
                        auto ra = final_idx.find(pair.id_a);
                        auto rb = final_idx.find(pair.id_b);
                        if (ra == final_idx.end() || rb == final_idx.end() ||
                            !results[ra->second].curve ||
                            !results[rb->second].curve) {
                            closure_pool_ok = false;
                            break;
                        }
                        ClosurePair item;
                        item.a = a_variable ? ia->second : ib->second;
                        if (a_variable && b_variable) {
                            item.b = ib->second;
                            // item.a/item.b 是按闭包成员索引规范化后的方向；
                            // 不能再使用 CurvePair 原始 id_a/id_b 取池，否则
                            // 记录为 38|5、闭包顺序为 5,9,38 时会把矩阵行列
                            // 交换，传播索引越界并错误清空候选域。
                            const auto& pool_a = closure_pools[
                                closure_ids[item.a]];
                            const auto& pool_b = closure_pools[
                                closure_ids[item.b]];
                            item.conflict.resize(pool_a.size() * pool_b.size(), 0);
                            for (std::size_t ai = 0; ai < pool_a.size(); ++ai)
                                for (std::size_t bi = 0; bi < pool_b.size(); ++bi)
                                    item.conflict[ai * pool_b.size() + bi] =
                                        closure_forbidden(
                                            pool_a[ai].curve, pool_a[ai].sampled,
                                            pool_b[bi].curve, pool_b[bi].sampled)
                                            ? 1 : 0;
                        } else {
                            item.b = std::numeric_limits<std::size_t>::max();
                            const ConnId& variable_id = a_variable
                                ? pair.id_a : pair.id_b;
                            const ConnId& fixed_id = a_variable
                                ? pair.id_b : pair.id_a;
                            item.fixed_id = fixed_id;
                            const auto& pool = closure_pools[variable_id];
                            const SampledCurve& fixed =
                                strict_result_samples.at(fixed_id);
                            item.conflict.resize(pool.size(), 0);
                            for (std::size_t ai = 0; ai < pool.size(); ++ai)
                                item.conflict[ai] = closure_forbidden(
                                    pool[ai].curve, pool[ai].sampled,
                                    *results[final_idx.at(fixed_id)].curve,
                                    fixed) ? 1 : 0;
                        }
                        closure_pairs.push_back(std::move(item));
                    }
                }
                if (isgDebugPairRepair() &&
                    family_item.first == "entry:1015689") {
                    for (const auto& item : closure_pairs) {
                        if (item.b == std::numeric_limits<std::size_t>::max())
                            continue;
                        const auto& pool_a = closure_pools[closure_ids[item.a]];
                        const auto& pool_b = closure_pools[closure_ids[item.b]];
                        std::size_t conflicts = 0;
                        std::size_t rows_with_support = 0;
                        std::size_t cols_with_support = 0;
                        for (std::size_t ai = 0; ai < pool_a.size(); ++ai) {
                            bool row_supported = false;
                            for (std::size_t bi = 0; bi < pool_b.size(); ++bi) {
                                if (!item.conflict[ai * pool_b.size() + bi])
                                    row_supported = true;
                                else
                                    ++conflicts;
                            }
                            rows_with_support += row_supported ? 1 : 0;
                        }
                        for (std::size_t bi = 0; bi < pool_b.size(); ++bi) {
                            for (std::size_t ai = 0; ai < pool_a.size(); ++ai)
                                if (!item.conflict[ai * pool_b.size() + bi]) {
                                    ++cols_with_support;
                                    break;
                                }
                        }
                        fprintf(stderr,
                                "[STRICT-CLOSURE-VARIABLE] %s|%s conflicts=%zu/%zu rows=%zu/%zu cols=%zu/%zu\n",
                                closure_ids[item.a].c_str(),
                                closure_ids[item.b].c_str(), conflicts,
                                pool_a.size() * pool_b.size(), rows_with_support,
                                pool_a.size(), cols_with_support, pool_b.size());
                        if ((closure_ids[item.a] == "5" &&
                             closure_ids[item.b] == "9") ||
                            (closure_ids[item.a] == "9" &&
                             closure_ids[item.b] == "5")) {
                            std::size_t exact_conflicts = 0;
                            std::size_t exact_rows = 0;
                            std::size_t exact_cols = 0;
                            for (std::size_t ai = 0; ai < pool_a.size(); ++ai) {
                                bool row_supported = false;
                                for (std::size_t bi = 0; bi < pool_b.size(); ++bi) {
                                    const bool conflict = curvesIntersectBusiness(
                                        pool_a[ai].curve, pool_b[bi].curve,
                                        kSharedEndpointPairCrossTol);
                                    exact_conflicts += conflict ? 1 : 0;
                                    row_supported = row_supported || !conflict;
                                }
                                exact_rows += row_supported ? 1 : 0;
                            }
                            for (std::size_t bi = 0; bi < pool_b.size(); ++bi) {
                                for (std::size_t ai = 0; ai < pool_a.size(); ++ai)
                                    if (!curvesIntersectBusiness(
                                            pool_a[ai].curve, pool_b[bi].curve,
                                            kSharedEndpointPairCrossTol)) {
                                        ++exact_cols;
                                        break;
                                    }
                            }
                            fprintf(stderr,
                                    "[STRICT-CLOSURE-EXACT] %s|%s conflicts=%zu/%zu rows=%zu/%zu cols=%zu/%zu\n",
                                    closure_ids[item.a].c_str(),
                                    closure_ids[item.b].c_str(), exact_conflicts,
                                    pool_a.size() * pool_b.size(), exact_rows,
                                    pool_a.size(), exact_cols, pool_b.size());
                        }
                    }
                }

                // 先做弧一致性约束传播：固定外部曲线直接删除冲突候选，
                // 变量之间删除在另一侧没有任何兼容支撑的候选。该过程只删掉
                // 不可能属于任意完整解的值，不改变闭包的可行解集合；之后的
                // 回溯只处理真正需要组合的剩余域。
                const auto closure_conflict = [&](const ClosurePair& item,
                                                  std::size_t ai,
                                                  std::size_t bi) {
                    if (item.b == std::numeric_limits<std::size_t>::max())
                        return item.conflict[ai] != 0;
                    const auto& pool_b = closure_pools[closure_ids[item.b]];
                    return item.conflict[ai * pool_b.size() + bi] != 0;
                };
                std::vector<std::vector<unsigned char>> closure_domains(
                    closure_ids.size());
                if (closure_pool_ok) {
                    for (std::size_t ci = 0; ci < closure_ids.size(); ++ci) {
                        const auto& pool = closure_pools[closure_ids[ci]];
                        closure_domains[ci].assign(pool.size(), 1);
                        if (pool.empty()) {
                            closure_pool_ok = false;
                            break;
                        }
                    }
                }
                const auto domain_has_value =
                    [&](std::size_t ci, std::size_t candidate) {
                        return ci < closure_domains.size() &&
                            candidate < closure_domains[ci].size() &&
                            closure_domains[ci][candidate] != 0;
                };
                if (closure_pool_ok) {
                    bool changed = true;
                    while (changed && closure_pool_ok) {
                        changed = false;
                        for (const auto& item : closure_pairs) {
                            if (item.b == std::numeric_limits<std::size_t>::max()) {
                                for (std::size_t ai = 0;
                                     ai < closure_domains[item.a].size(); ++ai) {
                                    if (domain_has_value(item.a, ai) &&
                                        item.conflict[ai]) {
                                        closure_domains[item.a][ai] = 0;
                                        changed = true;
                                    }
                                }
                                continue;
                            }
                            const std::size_t a_size =
                                closure_domains[item.a].size();
                            const std::size_t b_size =
                                closure_domains[item.b].size();
                            bool a_has_value = false;
                            bool b_has_value = false;
                            for (std::size_t ai = 0; ai < a_size; ++ai) {
                                if (!domain_has_value(item.a, ai))
                                    continue;
                                bool supported = false;
                                for (std::size_t bi = 0; bi < b_size; ++bi) {
                                    if (domain_has_value(item.b, bi) &&
                                        !closure_conflict(item, ai, bi)) {
                                        supported = true;
                                        break;
                                    }
                                }
                                if (supported)
                                    a_has_value = true;
                                else {
                                    closure_domains[item.a][ai] = 0;
                                    changed = true;
                                }
                            }
                            for (std::size_t bi = 0; bi < b_size; ++bi) {
                                if (!domain_has_value(item.b, bi))
                                    continue;
                                bool supported = false;
                                for (std::size_t ai = 0; ai < a_size; ++ai) {
                                    if (domain_has_value(item.a, ai) &&
                                        !closure_conflict(item, ai, bi)) {
                                        supported = true;
                                        break;
                                    }
                                }
                                if (supported)
                                    b_has_value = true;
                                else {
                                    closure_domains[item.b][bi] = 0;
                                    changed = true;
                                }
                            }
                            if (!a_has_value || !b_has_value) {
                                closure_pool_ok = false;
                                break;
                            }
                        }
                    }
                }
                if (isgDebugPairRepair() &&
                    family_item.first == "entry:1015689") {
                    fprintf(stderr, "[STRICT-CLOSURE-DOMAINS] ok=%d ",
                            closure_pool_ok ? 1 : 0);
                    for (std::size_t ci = 0; ci < closure_ids.size(); ++ci)
                        fprintf(stderr, "%s:%zu/%zu ", closure_ids[ci].c_str(),
                                std::count(closure_domains[ci].begin(),
                                           closure_domains[ci].end(),
                                           static_cast<unsigned char>(1)),
                                closure_domains[ci].size());
                    fprintf(stderr, "\n");
                    for (const auto& item : closure_pairs) {
                        if (item.b != std::numeric_limits<std::size_t>::max())
                            continue;
                        std::size_t conflicts = 0;
                        for (unsigned char value : item.conflict)
                            conflicts += value != 0 ? 1 : 0;
                        fprintf(stderr,
                                "[STRICT-CLOSURE-FIXED] %s vs %s conflicts=%zu/%zu\n",
                                closure_ids[item.a].c_str(), item.fixed_id.c_str(),
                                conflicts, item.conflict.size());
                    }
                }
                std::vector<std::size_t> closure_order;
                if (closure_pool_ok) {
                    closure_order.resize(closure_ids.size());
                    for (std::size_t ci = 0; ci < closure_ids.size(); ++ci)
                        closure_order[ci] = ci;
                    std::sort(closure_order.begin(), closure_order.end(),
                              [&](std::size_t a, std::size_t b) {
                                  const auto count_active =
                                      [&](std::size_t ci) {
                                          return std::count(
                                              closure_domains[ci].begin(),
                                              closure_domains[ci].end(),
                                              static_cast<unsigned char>(1));
                                      };
                                  const std::size_t a_count = count_active(a);
                                  const std::size_t b_count = count_active(b);
                                  if (a_count != b_count)
                                      return a_count < b_count;
                                  return a < b;
                              });
                }
                std::vector<std::size_t> closure_selected(
                    closure_ids.size(), 0);
                std::vector<bool> closure_selected_flags(
                    closure_ids.size(), false);
                std::size_t closure_nodes = 0;
                const std::size_t max_closure_nodes = 50000;
                const auto exact_closure_selection_is_clean = [&]() {
                    // 矩阵用于快速剪枝；叶节点逐对执行真实曲线相交复核，
                    // 避免采样近似或候选引用错误带入提交。
                    for (const auto& item : closure_pairs) {
                        const StrictOrdinaryCandidate& selected_a =
                            closure_pools[closure_ids[item.a]][
                                closure_selected[item.a]];
                        if (item.b == std::numeric_limits<std::size_t>::max()) {
                            auto fixed_it = final_idx.find(item.fixed_id);
                            if (fixed_it == final_idx.end() ||
                                !results[fixed_it->second].curve)
                                return false;
                            if (strictPairHasExactCross(
                                    selected_a.curve,
                                    *results[fixed_it->second].curve))
                                return false;
                            continue;
                        }
                        const StrictOrdinaryCandidate& selected_b =
                            closure_pools[closure_ids[item.b]][
                                closure_selected[item.b]];
                        if (strictPairHasExactCross(
                                selected_a.curve, selected_b.curve))
                            return false;
                    }
                    return true;
                };
                std::function<bool(std::size_t)> search_closure =
                    [&](std::size_t depth) {
                        if (!closure_pool_ok || ++closure_nodes > max_closure_nodes)
                            return false;
                        if (depth == closure_order.size()) {
                            for (std::size_t ci = 0; ci < closure_ids.size(); ++ci)
                                if (strictCandidateHasExactPhysicalRisk(
                                        closure_pools[closure_ids[ci]][
                                            closure_selected[ci]].curve))
                                    return false;
                            return exact_closure_selection_is_clean();
                        }
                        const std::size_t variable = closure_order[depth];
                        const auto& pool = closure_pools[closure_ids[variable]];
                        for (std::size_t candidate = 0;
                             candidate < pool.size(); ++candidate) {
                            if (!domain_has_value(variable, candidate))
                                continue;
                            bool valid = true;
                            for (const auto& item : closure_pairs) {
                                if (item.a != variable && item.b != variable)
                                    continue;
                                const std::size_t other =
                                    item.a == variable ? item.b : item.a;
                                if (other == std::numeric_limits<std::size_t>::max()) {
                                    continue;
                                }
                                if (!closure_selected_flags[other])
                                    continue;
                                const std::size_t ai = item.a == variable
                                    ? candidate : closure_selected[item.a];
                                const std::size_t bi = item.b == variable
                                    ? candidate : closure_selected[item.b];
                                const auto& pool_b = closure_pools[
                                    closure_ids[item.b]];
                                if (closure_conflict(item, ai, bi)) {
                                    valid = false;
                                    break;
                                }
                            }
                            if (!valid)
                                continue;
                            closure_selected[variable] = candidate;
                            closure_selected_flags[variable] = true;
                            if (search_closure(depth + 1))
                                return true;
                            closure_selected_flags[variable] = false;
                        }
                        return false;
                    };
                bool closure_found = closure_pool_ok && search_closure(0);
                if (!closure_found && closure_pool_ok &&
                    isgDebugPairRepair() &&
                    family_item.first == "entry:1015689") {
                    // 失败时求一次矩阵上的最小冲突组合，区分候选域完全
                    // 不可行与仅有少量外部冲突需要下一轮收口的情况。该诊断
                    // 不参与提交，避免为定位问题改变严格搜索结果。
                    int best_conflicts = std::numeric_limits<int>::max();
                    std::vector<std::size_t> best_selection(
                        closure_ids.size(), 0);
                    std::function<void(std::size_t, int)> find_best =
                        [&](std::size_t depth, int partial) {
                            if (partial >= best_conflicts)
                                return;
                            if (depth == closure_order.size()) {
                                best_conflicts = partial;
                                best_selection = closure_selected;
                                return;
                            }
                            const std::size_t variable = closure_order[depth];
                            const auto& pool = closure_pools[
                                closure_ids[variable]];
                            for (std::size_t candidate = 0;
                                 candidate < pool.size(); ++candidate) {
                                int added = 0;
                                for (const auto& item : closure_pairs) {
                                    if (item.a != variable && item.b != variable)
                                        continue;
                                    const std::size_t other = item.a == variable
                                        ? item.b : item.a;
                                    if (other == std::numeric_limits<std::size_t>::max()) {
                                        added += item.conflict[candidate] ? 1 : 0;
                                    } else if (closure_selected_flags[other]) {
                                        const std::size_t ai = item.a == variable
                                            ? candidate : closure_selected[item.a];
                                        const std::size_t bi = item.b == variable
                                            ? candidate : closure_selected[item.b];
                                        const auto& pool_b = closure_pools[
                                            closure_ids[item.b]];
                                        added += item.conflict[
                                            ai * pool_b.size() + bi] ? 1 : 0;
                                    }
                                }
                                closure_selected[variable] = candidate;
                                closure_selected_flags[variable] = true;
                                find_best(depth + 1, partial + added);
                                closure_selected_flags[variable] = false;
                            }
                        };
                    find_best(0, 0);
                    fprintf(stderr,
                            "[STRICT-CLOSURE-BEST] %s conflicts=%d handles=",
                            family_item.first.c_str(), best_conflicts);
                    for (std::size_t ci = 0; ci < closure_ids.size(); ++ci) {
                        const auto& g = closure_pools[closure_ids[ci]][
                            best_selection[ci]].curve.segs.front().ctrl;
                        fprintf(stderr, "%s:%.3f/%.3f,", closure_ids[ci].c_str(),
                                (g[1] - g[0]).norm(), (g[3] - g[2]).norm());
                    }
                    for (const auto& item : closure_pairs) {
                        const std::size_t ai = item.a < best_selection.size()
                            ? best_selection[item.a] : 0;
                        const std::size_t bi = item.b < best_selection.size()
                            ? best_selection[item.b] : 0;
                        const auto& pool_b = item.b < closure_ids.size()
                            ? closure_pools[closure_ids[item.b]]
                            : closure_pools[closure_ids[item.a]];
                        const bool conflict = item.b ==
                            std::numeric_limits<std::size_t>::max()
                            ? item.conflict[ai]
                            : item.conflict[ai * pool_b.size() + bi];
                        if (conflict)
                            fprintf(stderr, " [conflict=%s|%s]",
                                    closure_ids[item.a].c_str(),
                                    item.b == std::numeric_limits<std::size_t>::max()
                                        ? item.fixed_id.c_str()
                                        : closure_ids[item.b].c_str());
                    }
                    fprintf(stderr, "\n");
                }
                if (closure_found) {
                    std::vector<CurvePatchEntry> closure_patch;
                    closure_patch.reserve(closure_ids.size());
                    for (std::size_t ci = 0; ci < closure_ids.size(); ++ci)
                        closure_patch.push_back(CurvePatchEntry{
                            closure_ids[ci],
                            closure_pools[closure_ids[ci]][
                                closure_selected[ci]].curve});
                    AtomicCurvePatch().apply(
                        closure_patch, results,
                        [&](ConnectivityCurve& result,
                            const BezierCurve& candidate) {
                            setConnectivityCurveGeometry(result, candidate);
                            validate(result, input, sdf);
                        });
                    final_idx = resultIndexById(results);
                    family_pass_changed = true;
                    for (const ConnId& id : closure_ids)
                        strict_changed_this_pass.insert(id);
                    if (isgDebugPairRepair())
                        fprintf(stderr,
                                "[STRICT-CLOSURE] %s committed members=%zu nodes=%zu\n",
                                family_item.first.c_str(), closure_ids.size(),
                                closure_nodes);
                    continue;
                }
                if (isgDebugPairRepair())
                    fprintf(stderr,
                            "[STRICT-FAMILY] %s search failed nodes=%zu closure=%zu closure_nodes=%zu\n",
                            family_item.first.c_str(), search_nodes,
                            closure_ids.size(), closure_nodes);
                if (!closure_escalate)
                    strict_failed_this_pass.insert(family_item.first);
                continue;
            }

            if (isgDebugPairRepair())
                fprintf(stderr, "[STRICT-FAMILY] %s search accepted nodes=%zu\n",
                        family_item.first.c_str(), search_nodes);

            std::vector<CurvePatchEntry> patch;
            patch.reserve(family.size());
            for (std::size_t fi = 0; fi < family.size(); ++fi)
                patch.push_back(CurvePatchEntry{
                    family[fi], pools[fi][selected_index[fi]].curve});
            AtomicCurvePatch().apply(
                patch, results,
                [&](ConnectivityCurve& result, const BezierCurve& candidate) {
                    setConnectivityCurveGeometry(result, candidate);
                    validate(result, input, sdf);
                });
            final_idx = resultIndexById(results);
            family_pass_changed = true;
            for (const ConnId& id : family)
                strict_changed_this_pass.insert(id);
            }
            // 即使本轮没有任何族改写，只要还有"待升级重试"的失败族，就必须
            // 再走一轮：升级预算的唯一入口就是下一轮的 closure_escalate。
            if (!family_pass_changed && strict_failed_this_pass.empty())
                break;
            strict_changed_last_pass.swap(strict_changed_this_pass);
            strict_failed_last_pass.swap(strict_failed_this_pass);
        }
    }
    if (isgProfile()) {
        fprintf(stderr, "[ISG_PROFILE] ordinary shared family stage %.3f ms\n",
                std::chrono::duration<double, std::milli>(
                    std::chrono::steady_clock::now() - ordinary_family_stage_t0).count());
    }

    auto annotate_t0 = std::chrono::steady_clock::now();
    // 所有收口结束后重新计算单条状态；否则早期窄通道 WarnA2 或旧交叉状态会
    // 在几何已替换后继续残留，造成最终输出与当前曲线不一致。
    for (auto& result : results)
        result.violation.exempt_crosses.clear();
    for (auto& result : results)
        if (result.curve)
            validate(result, input, sdf);
    annotateClusterCrossings(results, cluster_solver_, kClusterEndpointTol);
    if (isgProfile()) {
        fprintf(stderr, "[ISG_PROFILE] annotate stage %.3f ms\n",
                std::chrono::duration<double, std::milli>(
                    std::chrono::steady_clock::now() - annotate_t0).count());
    }

    ElevationInterpolator(results, input);

    auto t1 = std::chrono::steady_clock::now();
    if (out_ms) *out_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    return results;
}

}
