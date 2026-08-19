#include "generator/connectivity_generation_session.h"
#include "curve/hermite_init.h"
#include "curve/curve_utils.h"
#include "constraints/boundary_safety.h"
#include "constraints/fence_check.h"
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
#include <string>

namespace isg {

// U型调头硬形态约束：进入/退出两侧都没有人行横道约束时，首尾直行段才使用2m保底距离进入中间单段掉头弧。
static constexpr double kUTurnNoCrosswalkMinLead = 2.0;
// 同配置家族U-turn的首/尾平齐点沿轴向法线相向微移，避免直行段完全重叠/贴合。
// 该量是家族级的分米级基线，真正施加时还会被几何上限截断：
// 1) 不能把中弧挤成高曲率尖弧；
// 2) 不能把首/尾直行段推出 G1 允许范围。
static constexpr double kUTurnSharedEndpointLeadStagger = 0.25;
static constexpr double kUTurnSharedEndpointPairLeadStagger = 0.01;
// 同簇曲线只允许在真实连接点相遇；连接点外的贴合/重叠/相交都按违规处理。
static constexpr double kClusterEndpointTol = 0.30;
// 同入口/同出口三段式U-turn首尾直行段从同一点扇出时，平齐点之间至少
// 需要达到该间距，才能在最终同簇标注里只检查后续中弧/尾段。
static constexpr double kUTurnSharedLeadFanMinSeparation = 0.10;


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


////////////////////////////////////////////////////////////
// 编译辅助函数
////////////////////////////////////////////////////////////
// 判别输入是否只包含连接、车道等几何数据，而不包含会触发物理约束的
// 障碍物、边界、人行横道或停车线。纯几何输入使用更轻量的拓扑收口路径。
static bool isPureGeometricInput(const IntersectionInput& input) {
    return input.obstacles.empty() && input.boundaries.empty() &&
           input.crosswalks.empty() && input.stop_lines.empty();
}

// ── 端点粘连 boundary 过滤 ─────────────────────────────────────────────────
// 某些路口臂的侧边边界由多段拼接而成（例如 RoadEdge + Connector + RoadEdge），
// 这些段链的首尾端点分别紧靠进入臂和退出臂的车道端点（距离约 0.3~0.5m）。
// 当曲线端点 p0/p1 位于链的端点附近且 boundary 方向与 T0/T1 对齐时，
// 曲线必然从该端点出发穿越链中间节点——这是正常转向/调头几何，不是违规。
// 过滤条件（同时满足以下全部才过滤）：
//   1. boundary 的某个端点距 p0 在 pos_tol 内，且 boundary 方向与 T0 对齐
//      （dot > dir_min_dot，即从进入点出发的侧边链）
//   OR
//   boundary 的某个端点距 p1 在 pos_tol 内，且 boundary 方向与 T1 对齐
//      （从退出点出发的侧边链）
// 端点/方向的分侧匹配（p0↔T0, p1↔T1）避免了将从退出点 p1 出发但沿进入方向
// T0 延伸的 RoadEdge 误过滤（这类 boundary 可能是退出口横向 RoadEdge 的端点
// 恰好与 p1 重合，方向却不沿 T1）。
// 返回过滤后的边界副本；输入边界顺序和对象内容保持不变，仅移除满足上述
// 端点粘连条件的边界。
static std::vector<Boundary> filterEndpointAdherentBoundaries(
    const std::vector<Boundary>& boundaries, const Vec2d& p0, const Vec2d& t0,
    const Vec2d& p1, const Vec2d& t1, double pos_tol = 0.6, double dir_min_dot = 0.70) {
    Vec2d T0 = t0.norm() > 1e-8 ? t0.normalized() : Vec2d(1, 0);
    Vec2d T1 = t1.norm() > 1e-8 ? t1.normalized() : Vec2d(1, 0);

    std::vector<Boundary> filtered;
    filtered.reserve(boundaries.size());
    for (const auto& bnd : boundaries) {
        const auto& pts = bnd.geometry.points;
        if (pts.empty()) {
            filtered.push_back(bnd);
            continue;
        }
        Vec2d first = xyOf(pts.front());
        Vec2d last  = xyOf(pts.back());

        Vec2d boundary_direction = last - first;
        if (boundary_direction.norm() > 1e-8)
            boundary_direction.normalize();

        // p0侧：boundary 端点距 p0 在 pos_tol 内，且方向与 T0 对齐
        bool near_p0 = (first - p0).norm() < pos_tol || (last - p0).norm() < pos_tol;
        bool p0_side_ok = near_p0 && boundary_direction.norm() > 1e-8 &&
                          std::abs(boundary_direction.dot(T0)) >= dir_min_dot;

        // p1侧：boundary 端点距 p1 在 pos_tol 内，且方向与 T1 对齐
        bool near_p1 = (first - p1).norm() < pos_tol || (last - p1).norm() < pos_tol;
        bool p1_side_ok = near_p1 && boundary_direction.norm() > 1e-8 &&
                          std::abs(boundary_direction.dot(T1)) >= dir_min_dot;

        if (p0_side_ok || p1_side_ok) {
            // 端点粘连且方向对齐对应侧 → 过滤（侧边链，允许穿越）
            continue;
        }
        filtered.push_back(bnd);
    }
    return filtered;
}

// 使用 BoundarySafety 的道路中心和采样规则检查曲线是否穿越边界或越出道路边缘。
// 空边界集合直接视为安全；返回值只表示物理边界风险，不包含障碍物或同簇风险。
static bool curveIntersectsBoundaries(
        const BezierCurve& curve, const std::vector<Boundary>& boundaries, const Vec2d& center) {
    if (boundaries.empty())
        return false;
    auto safety = curveBoundarySafety(curve, boundaries, center, 128, 0.15);
    return safety.intersects || safety.outside_road_edge;
}

// 对曲线进行采样后的线段级边界相交检查。
// 可选地只检查 RoadEdge，并忽略曲线端点和边界端点容差内的合法连接；
// 返回 true 表示存在一个非端点真实穿越。
static bool curveRawIntersectsBoundariesImpl(
        const BezierCurve& curve, const std::vector<Boundary>& boundaries,
        bool road_edges_only, double curve_endpoint_tol, double boundary_endpoint_tol) {
    if (curve.empty())
        return false;
    BoundingBox2d curve_box = curve.bbox();
    double bbox_pad = std::max(curve_endpoint_tol, boundary_endpoint_tol) + 0.02;
    curve_box.min_pt -= Vec2d(bbox_pad, bbox_pad);
    curve_box.max_pt += Vec2d(bbox_pad, bbox_pad);
    auto pts = curve.sampleByArcLength(std::max(
            64, std::min(240, (int)std::ceil(curve.arcLength() / 0.18) + 1)));
    if (pts.size() < 2)
        return false;
    for (const auto& bnd : boundaries) {
        if ((road_edges_only && bnd.type != Boundary::Type::RoadEdge) ||
            bnd.geometry.points.size() < 2)
            continue;
        if (!curve_box.intersects(bnd.geometry.bbox()))
            continue;
        const auto& bpts = bnd.geometry.points;
        for (int i = 0; i + 1 < (int)pts.size(); ++i) {
            for (int j = 0; j + 1 < (int)bpts.size(); ++j) {
                Vec2d isect;
                if (!segmentsIntersect(pts[i], pts[i + 1], bpts[j], bpts[j + 1], &isect))
                    continue;
                if ((isect - pts.front()).norm() <= curve_endpoint_tol ||
                    (isect - pts.back()).norm() <= curve_endpoint_tol)
                    continue;
                if ((j == 0 && (isect - bpts[j]).norm() <= boundary_endpoint_tol) ||
                    (j + 1 == (int)bpts.size() - 1 &&
                     (isect - bpts[j + 1]).norm() <= boundary_endpoint_tol))
                    continue;
                return true;
            }
        }
    }
    return false;
}

// 检查曲线是否与任意类型的边界发生非端点真实相交。
// 这是通用包装器，保留端点容差语义并委托给线段级实现。
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
        const BezierCurve& curve, const std::vector<Boundary>& boundaries, const Vec2d& center) {
    if (curve.numSegments() != 3)
        return curveIntersectsBoundaries(curve, boundaries, center);
    return curveRawIntersectsAnyBoundary(curve, boundaries, 0.15, 0.10);
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
        double ov = 0;
        // 性能优化: 上限从240降至80,降低采样数
        int n = std::max(20, std::min(80, (int)std::ceil(c.arcLength() / 0.25) + 1));
        for (auto& pt : c.sampleByArcLength(n))
            if (!polygonContains(input.area.geometry, pt))
                ov = std::max(ov, pointToPolygonDist(pt, input.area.geometry));
        audit.update_fence_overflow = true;
        audit.max_fence_overflow = ov;
    }
    audit.self_intersection = curveSelfIntersectsBusiness(c, 1.0);

    // 普通转向仍保留既有端点粘连边缘过滤；U型调头必须使用完整Boundary
    // 集合，curveBoundarySafety只裁掉首尾连续贴行段，避免中间弧重穿同边缘。
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

    std::vector<Boundary> filtered_bnds = input.boundaries;
    bool boundary_cross = false;
    if (is_conn_uturn) {
        Vec2d center_v = boundarySafetyCenter(input);
        boundary_cross = curveIntersectsBoundariesUTurn(c, input.boundaries, center_v);
    } else {
        filtered_bnds = filterEndpointAdherentBoundaries(input.boundaries, p0_v, t0_v, p1_v, t1_v);

        Vec2d center_v = boundarySafetyCenter(input);
        boundary_cross =
            curveIntersectsBoundaries(c, input.boundaries, center_v) ||
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
            c, p0_v, t0_v, p1_v, t1_v, 1e-5, true)) {
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

// 前置声明：判断同簇交叉是否属于允许的 U-turn/左右转结构性交叉。
static bool isAllowedSameClusterCrossing(
    const BezierCurve& a, const BezierCurve& b, double endpoint_tol);
// 前置声明：判断两条曲线是否存在必须修复的同簇相交或贴合。
static bool curvesHaveForbiddenSameClusterIntersection(
    const BezierCurve& a, const BezierCurve& b, double endpoint_tol);

// 查找距离给定点最近的 RoadEdge 线段，并返回距离及该线段的单位方向。
// 找不到有效 RoadEdge 时返回 false，并保留可用的默认方向。
static bool nearestRoadEdgeDirectionLocal(
        const std::vector<Boundary>& boundaries, const Vec2d& pt,
        double& out_dist, Vec2d& out_dir) {
    bool found = false;
    out_dist = std::numeric_limits<double>::infinity();
    out_dir = Vec2d(1, 0);
    for (const auto& bnd : boundaries) {
        if (bnd.type != Boundary::Type::RoadEdge ||
            bnd.geometry.points.size() < 2)
            continue;
        std::vector<Vec2d> pts = toVec2dArray(bnd.geometry.points);
        for (int i = 0; i + 1 < (int)pts.size(); ++i) {
            Vec2d d = pts[i + 1] - pts[i];
            double len2 = d.squaredNorm();
            if (len2 < 1e-12)
                continue;
            double t = std::max(0.0, std::min(1.0, (pt - pts[i]).dot(d) / len2));
            Vec2d proj = pts[i] + d * t;
            double pd = dist(pt, proj);
            if (!found || pd < out_dist) {
                found = true;
                out_dist = pd;
                out_dir = d.normalized();
            }
        }
    }
    return found;
}

// 以少量参数点检查 Bezier 段是否始终贴近某条 RoadEdge。
// 该局部判定用于识别合法的端点贴行段，不替代完整边界相交检查。
static bool segmentStaysNearRoadEdgeLocal(
        const BezierSegment& seg, const std::vector<Boundary>& boundaries,
        double tol = 0.14) {
    for (int i = 0; i <= 6; ++i) {
        double edge_dist = 0.0;
        Vec2d edge_dir{1, 0};
        if (!nearestRoadEdgeDirectionLocal(
                boundaries, seg.evaluate((double)i / 6.0), edge_dist, edge_dir) ||
            edge_dist > tol)
            return false;
    }
    return true;
}

// 判断曲线首段或尾段是否形成贴近 RoadEdge 的端点段。
// 返回 true 表示曲线可能属于道路端点拟合场景，可采用专门的形态豁免。
static bool curveHasRoadEdgeAdherentEndpointSectionLocal(
        const BezierCurve& curve, const std::vector<Boundary>& boundaries) {
    if (curve.segs.empty())
        return false;
    if (segmentStaysNearRoadEdgeLocal(curve.segs.front(), boundaries))
        return true;
    if (segmentStaysNearRoadEdgeLocal(curve.segs.back(), boundaries))
        return true;
    return false;
}

// 去除曲线首尾连续贴近 RoadEdge 的段，保留中间曲线副本用于交叉和形态审计。
// 若没有可安全去除的段，则返回原曲线；若全曲线均贴边，则保留原曲线避免空结果。
static BezierCurve curveWithoutRoadEdgeEndpointSectionsLocal(
        const BezierCurve& curve, const std::vector<Boundary>& boundaries) {
    if (curve.segs.empty())
        return curve;
    int first = 0;
    int last = (int)curve.segs.size() - 1;
    while (first <= last &&
           segmentStaysNearRoadEdgeLocal(curve.segs[first], boundaries)) {
        ++first;
    }
    while (last >= first &&
           segmentStaysNearRoadEdgeLocal(curve.segs[last], boundaries)) {
        --last;
    }
    if (first == 0 && last == (int)curve.segs.size() - 1)
        return curve;
    if (first > last)
        return curve;
    BezierCurve middle;
    for (int i = first; i <= last; ++i)
        middle.segs.push_back(curve.segs[i]);
    return middle;
}


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
static bool sampledCurvesIntersectBusiness(
    const SampledCurve& a, const SampledCurve& b, double endpoint_tol,
    Vec2d* out = nullptr, bool detect_near_overlap = false) {
    if (a.pts.size() < 2 || b.pts.size() < 2 || !a.bbox.intersects(b.bbox))
        return false;
    for (int ai = 0; ai + 1 < (int)a.pts.size(); ++ai) {
        Vec2d amid = 0.5 * (a.pts[ai] + a.pts[ai + 1]);
        for (int bi = 0; bi + 1 < (int)b.pts.size(); ++bi) {
            Vec2d bmid = 0.5 * (b.pts[bi] + b.pts[bi + 1]);
            if ((amid - bmid).squaredNorm() > 900.0)
                continue;
            Vec2d isect;
            if (!segmentsIntersect(a.pts[ai], a.pts[ai + 1], b.pts[bi], b.pts[bi + 1], &isect)) {
                if (!detect_near_overlap)
                    continue;
                Vec2d ad = a.pts[ai + 1] - a.pts[ai];
                Vec2d bd = b.pts[bi + 1] - b.pts[bi];
                if (ad.norm() < 1e-8 || bd.norm() < 1e-8)
                    continue;
                if (std::abs(ad.normalized().dot(bd.normalized())) < 0.96)
                    continue;
                double near_dist = std::min({
                    pointToSegment(a.pts[ai], b.pts[bi], b.pts[bi + 1]).first,
                    pointToSegment(a.pts[ai + 1], b.pts[bi], b.pts[bi + 1]).first,
                    pointToSegment(b.pts[bi], a.pts[ai], a.pts[ai + 1]).first,
                    pointToSegment(b.pts[bi + 1], a.pts[ai], a.pts[ai + 1]).first
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
    bool detect_near_overlap = false) {
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
                curve, sib.curve, endpoint_tol))
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

// 前置声明：供种子影响闭包使用的同簇交叉豁免判定。
static bool isAllowedSameClusterCrossing(
    const BezierCurve& a, const BezierCurve& b, double endpoint_tol);

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

// 根据弦线与入口切向的偏转判断曲线是否具有明显左右转形态。
// U-turn 和近直线曲线返回 false，供结构性交叉豁免规则区分转向类型。
static bool curveLooksLeftRightTurnForClusterExemption(const BezierCurve& curve) {
    if (curve.empty() || curveLooksUTurnForClusterExemption(curve))
        return false;
    Vec2d chord = curve.endPt() - curve.startPt();
    Vec2d st = curve.startTan();
    if (chord.norm() < 1e-8 || st.norm() < 1e-8)
        return false;
    return std::abs(cross2d(st.normalized(), chord.normalized())) > 0.35;
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
static bool isAllowedSameClusterCrossing(
    const BezierCurve& a, const BezierCurve& b, double endpoint_tol = 0.15) {
    if (curvesHaveForbiddenAdherence(a, b, endpoint_tol))
        return false;
    bool a_uturn = curveLooksUTurnForClusterExemption(a);
    bool b_uturn = curveLooksUTurnForClusterExemption(b);
    if (a_uturn == b_uturn)
        return false;
    return a_uturn
        ? curveLooksLeftRightTurnForClusterExemption(b)
        : curveLooksLeftRightTurnForClusterExemption(a);
}

// 汇总同簇相交和贴合规则，返回是否存在必须修复的非豁免冲突。
// 先排除允许的结构性交叉，再检查业务相交，因此是同簇约束的统一入口。
static bool curvesHaveForbiddenSameClusterIntersection(
    const BezierCurve& a, const BezierCurve& b, double endpoint_tol = 0.15) {
    if (curvesHaveForbiddenAdherence(a, b, endpoint_tol))
        return true;
    if (isAllowedSameClusterCrossing(a, b, endpoint_tol))
        return false;
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

// 采样曲线内部点并计算其到所有 RoadEdge 的最小距离。
// 两端指定距离内的点被跳过，以避免合法连接端点影响道路边缘净距评估。
static double minRoadEdgeDistanceAlongCurve(
    const BezierCurve& curve, const std::vector<Boundary>& boundaries, double endpoint_skip = 0.75) {
    if (curve.empty() || boundaries.empty())
        return std::numeric_limits<double>::infinity();
    double min_edge_dist = std::numeric_limits<double>::infinity();
    int n = std::max(32, std::min(160, (int)std::ceil(curve.arcLength() / 0.20) + 1));
    for (auto& pt : curve.sampleByArcLength(n)) {
        if ((pt - curve.startPt()).norm() <= endpoint_skip ||
            (pt - curve.endPt()).norm() <= endpoint_skip)
            continue;
        for (const auto& bnd : boundaries) {
            if (bnd.type != Boundary::Type::RoadEdge)
                continue;
            for (int bi = 0; bi + 1 < (int)bnd.geometry.points.size(); ++bi) {
                min_edge_dist = std::min(
                    min_edge_dist, pointToSegment(pt, bnd.geometry.points[bi],bnd.geometry.points[bi + 1]).first);
            }
        }
    }
    return min_edge_dist;
}

// 将曲线到 RoadEdge 的最小距离转换为当前 mode 下的净距违规量。
// 无需净距约束或没有有效边界时返回零。
static double roadEdgeClearanceViolation(
    const BezierCurve& curve, const std::vector<Boundary>& boundaries, int mode) {
    double clearance = roadEdgeAvoidanceClearanceForMode(mode);
    if (clearance <= 0.0)
        return 0.0;
    double min_edge_dist = minRoadEdgeDistanceAlongCurve(curve, boundaries);
    if (!std::isfinite(min_edge_dist))
        return 0.0;
    return std::max(0.0, clearance - min_edge_dist);
}

// 计算候选曲线的边界避让代价，合并道路边缘净距和真实边界穿越惩罚。
// 该值供 U-turn 搜索器排序使用，不直接修改曲线。
static double boundaryAvoidancePenalty(
    const BezierCurve& curve, const IntersectionInput& input) {
    double penalty = roadEdgeClearanceViolation(curve, input.boundaries, input.mode);
    // 所有mode都禁止与Boundary贴合、重叠或相交；mode=2在此基础上
    // 额外要求RoadEdge非端点1m净距。
    if (curveIntersectsBoundaries(curve, input.boundaries, boundarySafetyCenter(input)) ||
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

// 检查一条直线段是否真实穿越 RoadEdge，并忽略两端合法端点相交。
static bool segmentRawIntersectsRoadEdge(
    const Vec2d& a, const Vec2d& b, const std::vector<Boundary>& boundaries,
    double endpoint_tol = 0.15, double boundary_endpoint_tol = 0.10) {
    for (const auto& bnd : boundaries) {
        if (bnd.type != Boundary::Type::RoadEdge ||
            bnd.geometry.points.size() < 2)
            continue;
        std::vector<Vec2d> pts = toVec2dArray(bnd.geometry.points);
        for (int i = 0; i + 1 < (int)pts.size(); ++i) {
            Vec2d isect;
            if (!segmentsIntersect(a, b, pts[i], pts[i + 1], &isect))
                continue;
            if ((isect - a).norm() <= endpoint_tol ||
                (isect - b).norm() <= endpoint_tol)
                continue;
            if ((i == 0 && (isect - pts[i]).norm() <= boundary_endpoint_tol) ||
                (i + 1 == (int)pts.size() - 1 &&
                 (isect - pts[i + 1]).norm() <= boundary_endpoint_tol))
                continue;
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

// 采样检查曲线内部点是否离开围栏多边形，并允许小的数值边界误差。
// 空围栏不产生风险。
static bool curveLeavesFence(const BezierCurve& curve, const Polygon2d& fence) {
    if (fence.outer.empty())
        return false;
    auto pts = curve.sample(24);
    for (int i = 1; i + 1 < (int)pts.size(); ++i) {
        const auto& pt = pts[i];
        if (!polygonContains(fence, pt) && pointToPolygonDist(pt, fence) > 0.10)
            return true;
    }
    return false;
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

// 校验非 U-turn 单段转向的弧长、曲率和曲率符号约束。
// 短急弯采用受限例外；其它转向必须保持可见单拱形态且不能形成尖钩或 S 弯。
static bool isNonUTurnTurnShapeAcceptable(
    const BezierCurve& curve, double chord_len, double turn_strength) {
    if (chord_len < 1e-6 || turn_strength <= 0.35)
        return true;
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
    double min_arc = chord_len < 12.0 ? 1.02 : 1.06;
    if (arc_chord < min_arc || arc_chord > 1.35)
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

// 判断曲线是否满足基础直行的严格形态限制。
// 在宽松直行判定之上进一步限制弧弦比和最大曲率，供自然直行恢复使用。
static bool isStrictStraightBaseShapeAcceptable(
    const BezierCurve& curve, double chord_len) {
    if (!isStraightLikeShapeAcceptable(curve, chord_len))
        return false;
    return curve.arcLength() / chord_len < 1.02 &&
           curve.maxCurvature(40) < 0.10;
}

// 判断含 RoadEdge 端点贴行段的候选是否具备可接受的中段形态。
// 该检查允许端点贴边，但要求去除端点段后的中段不穿越道路边缘。
static bool roadEdgeEndpointFitShapeAcceptableLocal(
    const BezierCurve& curve, const std::vector<Boundary>& boundaries,
    double full_chord_len, double turn_strength) {
    if (full_chord_len < 1e-6 ||
        !curveHasRoadEdgeAdherentEndpointSectionLocal(curve, boundaries))
        return false;
    BezierCurve middle = curveWithoutRoadEdgeEndpointSectionsLocal(curve, boundaries);
    if (middle.empty())
        return false;
    Vec2d middle_chord = middle.endPt() - middle.startPt();
    double middle_chord_len = middle_chord.norm();
    if (middle_chord_len < 1e-6)
        return false;
    if (curve.arcLength() / full_chord_len >= 1.30 ||
        middle.maxCurvature(40) >= 1.25)
        return false;
    if (turn_strength < 0.25)
        return isStraightLikeShapeAcceptable(middle, middle_chord_len);
    return isNonUTurnTurnShapeAcceptable(
        middle, middle_chord_len, turn_strength);
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
    double min_arc = chord_len < 12.0 ? 1.03 : 1.08;
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
static CurveRisk assessCurveRisk(
    const BezierCurve& curve, const IntersectionInput& input, const SDFField& sdf,
    const std::vector<SampledSiblingCurve>& sampled_siblings, bool include_fence = true,
    bool is_uturn = false) {
    CurveRisk risk;
    double ms = minSDFAlongCurveAdaptive(curve, sdf);
    risk.obstacle = curveIntersectsObstacles(curve, input.obstacles) || (ms < 0.0);

    // 对于U-turn，使用特殊的边界检查逻辑，优先考虑弧形边界穿越
    if (is_uturn) {
        // 使用特殊的U-turn边界检查，优先穿越弧形边界
        risk.boundary = curveIntersectsBoundariesUTurn(
            curve, input.boundaries, boundarySafetyCenter(input));
    } else {
        risk.boundary = curveIntersectsBoundaries(
            curve, input.boundaries, boundarySafetyCenter(input)) ||
            curveRawIntersectsAnyBoundary(
                curve, input.boundaries, 0.15, 0.10);
    }

    if (roadEdgeClearanceViolation(curve, input.boundaries, input.mode) > 0.0)
        risk.boundary = true;
    risk.fence = include_fence && (!input.area.is_rough && curveLeavesFence(curve, input.area.geometry));
    risk.sibling_crosses =
        sampledSiblingCrossCount(curve, sampled_siblings, true, kClusterEndpointTol);
    return risk;
}

// 当当前曲线存在物理风险时，在有限把手长度集合中寻找安全的单段 G1 三次曲线。
// 候选必须消除障碍物、边界和围栏风险，并在 U-turn 拓扑场景下保持同簇无交叉；
// 成功时原地替换 curve，失败时保持输入不变并返回 false。
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
    if (curve.numSegments() <= 1)
        return false;
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
    for (double alpha : {0.42, 0.46, 0.50, 0.38, 0.34}) {
        BezierCurve candidate;
        candidate.segs.push_back(makeCubicG1(p0, T0, p1, T1, alpha));
        if (candidate.empty() || curveSelfIntersectsBusiness(candidate, 1.0))
            continue;
        if (!ordinarySingleCubicControlsValid(
                candidate, p0, T0, p1, T1, 1e-5, true))
            continue;
        if (!isNonUTurnTurnShapeAcceptable(candidate, chord_len, turn_strength))
            continue;
        double pure_min_arc = chord_len < 12.0 ? 1.03 : 1.08;
        if (candidate.arcLength() / chord_len < pure_min_arc)
            continue;
        CurveRisk risk = assessCurveRisk(
            candidate, input, sdf, sampled_siblings, include_fence);
        if (risk.physical())
            continue;
        double arc_chord = candidate.arcLength() / chord_len;
        double score = 1000.0 * risk.sibling_crosses +
                       10.0 * std::abs(arc_chord - 1.10) +
                       candidate.maxCurvature(40) + 0.02 * candidate.arcLength();
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
    if (curveHasRoadEdgeAdherentEndpointSectionLocal(curve, input.boundaries))
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
                    "[STRAIGHT-NATURAL] reject p0=(%.2f,%.2f) p1=(%.2f,%.2f) current_cross=%d total=%d shape=%d phys=%d cross=%d\n",
                    p0.x(), p0.y(), p1.x(), p1.y(),
                    current_risk.sibling_crosses,
                    dbg_total, dbg_shape, dbg_phys, dbg_cross);
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
    add_conservative_two_segment_candidates();
    if (have_best && best_cross == 0) {
        curve = best;
        return true;
    }
    add_endpoint_opening_candidates();

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
static bool naturalSingleSegmentCandidate(
    const Connectivity& conn, const IntersectionInput& input, BezierCurve& out_curve) {
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
static std::unordered_set<ConnId> collectFixedIdsBlockingNaturalShapes(
    const IntersectionInput& input, const std::vector<ConnectivityCurve>& fixed_results,
    const ClusterOrderSolver& cs, const std::unordered_set<ConnId>& preserved_fixed_ids,
    double endpoint_tol = kClusterEndpointTol) {
    std::unordered_set<ConnId> drop_ids;
    auto idx = resultIndexById(fixed_results);
    for (const auto& conn : input.connectivities) {
        BezierCurve natural;
        if (!naturalSingleSegmentCandidate(conn, input, natural))
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
static double uturnSharedEndpointStagger(
    const Connectivity& conn, const IntersectionInput& input,
    const ClusterOrderSolver& cs, UTurnAlignmentScope scope) {
    // 分档必须与平齐站位使用相同的同入/同出传递连通家族。若只在当前
    // entry 或 exit 的一跳家族中排序，连接两个家族的 U-turn 会在中弧
    // 起点/终点拿到同一分档，重新产生非连接点交叉。
    const bool require_shared_endpoint_pair =
        scope == UTurnAlignmentScope::LaneEndpoint;
    const UTurnFamilyRank rank = UTurnFamilyBuilder().radiusRank(
        conn, input, scope, &cs, require_shared_endpoint_pair);
    if (rank.family_size < 2)
        return 0.0;
    // family按短径到长径排列；错开量按逆序给N，使内层曲线更紧、
    // 外层曲线更开，同时保留同入口直行段的扇出间距。
    return (rank.family_size >= 3
                ? kUTurnSharedEndpointLeadStagger
                : kUTurnSharedEndpointPairLeadStagger) *
           static_cast<double>(rank.reverse_radius_rank);
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
        double pure_min_arc = chord_len < 12.0 ? 1.03 : 1.08;
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
                candidate, input, sdf, sampled_siblings);
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
        return boundaryAvoidancePenalty(curve, input);
    };
    backend.endpoint_side_violation = [&](const SampledCurve& curve) {
        return sharedEndpointSideViolation(
            curve, sampled_siblings, kClusterEndpointTol);
    };
    return backend;
}

// 创建三段式 U-turn 候选审计器，统一执行物理风险和同簇交叉检查。
// include_fence 由调用阶段决定是否把围栏纳入硬门禁。
static SegmentedUTurnAuditor makeSegmentedUTurnAuditor(
    const IntersectionInput& input, const SDFField& sdf,
    const std::vector<SampledSiblingCurve>& sampled_siblings) {
    return [&](const BezierCurve& curve, bool include_fence) {
        const CurveRisk risk = assessCurveRisk(
            curve, input, sdf, sampled_siblings, include_fence, true);
        SegmentedUTurnAudit audit;
        audit.physical_violation = risk.physical();
        audit.sibling_crosses = risk.sibling_crosses;
        return audit;
    };
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
    const UTurnFamilyInfo& uturn_family_info = context.uturn_family;
    std::vector<Crosswalk> uturn_clearance_crosswalks;
    if (is_uturn_geom) {
        uturn_min_lead0 = uturn_family_info.lead0;
        uturn_min_lead1 = uturn_family_info.lead1;
        uturn_alignment_station = uturn_family_info.aligned_station;
        uturn_clearance_crosswalks =
            uturn_family_info.clearance_crosswalks;
        uturn_shared_endpoint_stagger = uturnSharedEndpointStagger(
            conn, input, cluster_solver_,
            direction_cfg_.uturn_alignment_scope);
        uturn_shared_entry_stagger = uturn_shared_endpoint_stagger;
        uturn_shared_exit_stagger = uturn_shared_endpoint_stagger;
        // 同配置家族的多条U-turn按两侧最大顺序号计算统一错开距离，
        // 首/尾两个已轴向平齐的点沿轴向法线同步相向微移；
        // 移动后仍保持 q0/q1 的轴向平齐关系。
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
    CurveRisk risk = assessCurveRisk(initial, input, sdf, sampled_for_gate, enforce_fence);
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
            SegmentedUTurnCandidateSearch().search(
                p0, t0, p1, t1, input, sampled_for_gate,
                makeSegmentedUTurnAuditor(input, sdf, sampled_for_gate), enforce_fence,
                segmented, uturn_min_lead0, uturn_min_lead1,
                &uturn_clearance_crosswalks, uturn_shared_endpoint_stagger,
                uturn_lead0_extra_after_align, uturn_lead1_extra_after_align,
                uturn_alignment_station,
                uturn_shared_entry_stagger, uturn_shared_exit_stagger);
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
                assessCurveRisk(initial, input, sdf, sampled_for_gate, enforce_fence);
            CurveRisk segmented_risk =
                assessCurveRisk(segmented_arch, input, sdf, sampled_for_gate, enforce_fence);
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
            SegmentedUTurnCandidateSearch().search(
                p0, t0, p1, t1, input, sampled_for_gate,
                makeSegmentedUTurnAuditor(input, sdf, sampled_for_gate), enforce_fence,
                segmented, uturn_min_lead0, uturn_min_lead1,
                &uturn_clearance_crosswalks, uturn_shared_endpoint_stagger,
                uturn_lead0_extra_after_align, uturn_lead1_extra_after_align,
                uturn_alignment_station,
                uturn_shared_entry_stagger, uturn_shared_exit_stagger);
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
        SegmentedUTurnCandidateSearch().search(
            p0, t0, p1, t1, input, sampled_for_gate,
            makeSegmentedUTurnAuditor(input, sdf, sampled_for_gate), enforce_fence,
            final_c, uturn_min_lead0, uturn_min_lead1,
            &uturn_clearance_crosswalks, uturn_shared_endpoint_stagger,
            uturn_lead0_extra_after_align, uturn_lead1_extra_after_align,
            uturn_alignment_station,
            uturn_shared_entry_stagger, uturn_shared_exit_stagger);
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
        if (SegmentedUTurnCandidateSearch().search(
                p0, t0, p1, t1, input, sampled_for_gate,
                makeSegmentedUTurnAuditor(input, sdf, sampled_for_gate),
                enforce_fence, segmented,
                uturn_min_lead0, uturn_min_lead1,
                &uturn_clearance_crosswalks, uturn_shared_endpoint_stagger,
                uturn_lead0_extra_after_align, uturn_lead1_extra_after_align,
                uturn_alignment_station,
                uturn_shared_entry_stagger, uturn_shared_exit_stagger)) {
            final_c = segmented;
        }
    }

    if (!is_uturn_geom && !boundary_repaired)
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
            auto endpointFitAlreadySafe = [&]() {
                if (!boundary_repaired)
                    return false;
                if (!roadEdgeEndpointFitShapeAcceptableLocal(
                        final_c, input.boundaries, chord_len2, turn_strength))
                    return false;
                BezierCurve middle =
                    curveWithoutRoadEdgeEndpointSectionsLocal(
                        final_c, input.boundaries);
                return !curveRawIntersectsRoadEdge(middle, input.boundaries);
            };
            bool straight_like = turn_strength < 0.25;
            double arc_chord = final_c.arcLength() / chord_len2;
            bool shape_bad = straight_like
                ? (arc_chord > 1.08 || final_c.maxCurvature(40) > 1.0)
                : isNonUTurnTurnShapeRisk(final_c, chord_len2, turn_strength);
            if (shape_bad && !endpointFitAlreadySafe()) {
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
    CurveRisk exit_risk = assessCurveRisk(final_c, input, sdf, sampled_for_gate, enforce_fence);
    bool exit_endpoint_fit_ok = false;
    if (!is_uturn_geom && exit_risk.boundary && boundary_repaired) {
        Vec2d chord = p1 - p0;
        if (chord.norm() > 1e-8) {
            Vec2d T0 = t0.norm() > 1e-8 ? t0.normalized() : chord.normalized();
            double turn_strength = std::abs(cross2d(T0, chord.normalized()));
            BezierCurve middle =
                curveWithoutRoadEdgeEndpointSectionsLocal(
                    final_c, input.boundaries);
            exit_endpoint_fit_ok =
                roadEdgeEndpointFitShapeAcceptableLocal(
                    final_c, input.boundaries, chord.norm(), turn_strength) &&
                !curveRawIntersectsRoadEdge(middle, input.boundaries) &&
                !exit_risk.obstacle && !exit_risk.fence;
        }
    }
    if (!is_uturn_geom && exit_risk.boundary && !exit_endpoint_fit_ok) {
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
                exit_risk = assessCurveRisk(final_c, input, sdf, sampled_for_gate, enforce_fence);
            }
        }
    }
    if (!is_uturn_geom && !exit_risk.physical()) {
        tryNaturalSingleCubicIfSafe(
            p0, t0, p1, t1, input, sdf, sampled_for_gate, enforce_fence, final_c);
        tryNaturalStraightSingleCubicIfSafe(
            p0, t0, p1, t1, input, sdf, sampled_for_gate, enforce_fence, final_c);
        exit_risk = assessCurveRisk(final_c, input, sdf, sampled_for_gate, enforce_fence);
    }
    if (is_uturn_geom) {
        if (SegmentedUTurnCandidateSearch().search(
                p0, t0, p1, t1, input, sampled_for_gate,
                makeSegmentedUTurnAuditor(input, sdf, sampled_for_gate), enforce_fence,
                final_c, uturn_min_lead0, uturn_min_lead1,
                &uturn_clearance_crosswalks, uturn_shared_endpoint_stagger,
                uturn_lead0_extra_after_align, uturn_lead1_extra_after_align,
                uturn_alignment_station,
                uturn_shared_entry_stagger, uturn_shared_exit_stagger))
            exit_risk = assessCurveRisk(final_c, input, sdf, sampled_for_gate, enforce_fence);
    }
    if (out_physical_risk) {
        *out_physical_risk = exit_risk.physical() || shape_risk;
    }
    setConnectivityCurveGeometry(cc, final_c);
    validate(cc, input, sdf);
    if (pre.narrow_passage && cc.status == CurveStatus::OK)
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
        input, results, cluster_solver_, preserved_fixed_ids, kClusterEndpointTol);
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
                if (curveHasRoadEdgeAdherentEndpointSectionLocal(
                        *cc.curve, input.boundaries))
                    roadedge_boundary_repair_ids.insert(cid);
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
                double shared_endpoint_stagger = uturnSharedEndpointStagger(
                    conn, input, cluster_solver_, direction_cfg_.uturn_alignment_scope);
                double shared_entry_stagger = shared_endpoint_stagger;
                double shared_exit_stagger = shared_endpoint_stagger;
                const std::vector<double> extra_options = input.mode == 2
                    ? std::vector<double>{0.0, 2.0, 4.0, 6.0, 8.0}
                    : std::vector<double>{0.0, 1.0, 2.0, 4.0};
                const std::vector<double> alpha_options = input.mode == 2
                    ? std::vector<double>{2.0 / 3.0, 0.85, 1.0, 1.25, 1.50,
                                           1.75, 2.0, 0.50, 0.38, 0.28, 0.16}
                    : std::vector<double>{2.0 / 3.0, 0.50, 0.38, 0.28};
                const std::vector<double> stagger_options = input.mode == 2
                    ? std::vector<double>{0.0, 0.2, 0.4, 0.6}
                    : std::vector<double>{0.0, 0.2, 0.4};
                const size_t candidate_cap = input.mode == 2 ? 96 : 32;
                for (double extra0 : extra_options) {
                    for (double extra1 : extra_options) {
                        for (double arc_alpha : alpha_options) {
                            for (double extra_stagger : stagger_options) {
                            // 即使只共享一侧端点，q0/q1 也必须使用相同总错开量。
                            // 继续传旧的 entry/exit 值会覆盖
                            // buildSegmented 的通用错开量，使
                            // extra_stagger 实际失效。
                            double total_stagger = std::max(0.0, shared_endpoint_stagger + extra_stagger);
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
            double best_score = std::numeric_limits<double>::infinity();
            auto consider = [&](const BezierCurve& candidate, double offset) {
                if (candidate.empty() ||
                    curveSelfIntersectsBusiness(candidate, 1.0) ||
                    !isNonUTurnTurnShapeAcceptable(
                        candidate, chord_len, turn_strength))
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
                    candidate.maxCurvature(40) + 0.02 * candidate.arcLength();
                if (!have_best ||
                    std::max(cross_count, shared_cross) < best_cross ||
                    (std::max(cross_count, shared_cross) == best_cross &&
                     fixed_cross < best_fixed_cross) ||
                    (std::max(cross_count, shared_cross) == best_cross &&
                     fixed_cross == best_fixed_cross &&
                     score < best_score)) {
                    best = candidate;
                    best_cross = std::max(cross_count, shared_cross);
                    best_fixed_cross = fixed_cross;
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

    bool topology_repair = input.obstacles.empty();
    if (topology_repair) {
        auto topology_t0 = std::chrono::steady_clock::now();
        auto result_idx = resultIndexById(results);
        auto neighbors = constrainedNeighborMap(results, cluster_solver_);
        // 固有掉头可能同时影响多条同簇直行/左转；一轮修复后其它曲线的
        // 位置变化会改变交叉关系，必须迭代复查到稳定。
        for (int pass = 0; pass < repair_budget.pure_topology_passes; ++pass) {
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
                    double shared_endpoint_stagger = uturnSharedEndpointStagger(
                        *conn, input, cluster_solver_,
                        direction_cfg_.uturn_alignment_scope);
                    double shared_entry_stagger = shared_endpoint_stagger;
                    double shared_exit_stagger = shared_endpoint_stagger;
                    if (!SegmentedUTurnCandidateSearch().search(
                            entry.first, entry.second, exit_.first, exit_.second,
                            input, sampled,
                            makeSegmentedUTurnAuditor(input, sdf, sampled),
                            true, repaired,
                            family_info.lead0, family_info.lead1, nullptr,
                            shared_endpoint_stagger, 0.0, 0.0,
                            alignment_station, shared_entry_stagger,
                            shared_exit_stagger))
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
                double pure_turn_min_arc = chord_len < 12.0 ? 1.03 : 1.08;
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
        int base_expression_passes = complex_boundary_fast_path
            ? 0 : repair_budget.base_expression_passes;
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
                        double pure_min_arc = chord_len < 12.0 ? 1.03 : 1.08;
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
                    BezierCurve shape_curve = straight_like
                        ? curveWithoutRoadEdgeEndpointSectionsLocal(
                              accepted_candidate, input.boundaries)
                        : accepted_candidate;
                    if (shape_curve.empty())
                        shape_curve = accepted_candidate;
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
                                             bool broad_turn_candidates) {
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
            const std::vector<double> pair_alphas = broad_turn_candidates
                ? std::vector<double>{0.18, 0.24, 0.30, 0.34, 0.42, 0.50}
                : std::vector<double>{0.34, 0.42, 0.50};
            for (double alpha : pair_alphas) {
                BezierCurve c;
                c.segs.push_back(makeCubicG1(
                    entry.first, T0, exit_.first, T1, alpha));
                add(c);
            }
            const std::vector<double> pair_a0 = broad_turn_candidates
                ? std::vector<double>{0.005, 0.02, 0.04, 0.055, 0.08, 0.12, 0.16, 0.22, 0.30, 0.42, 0.55}
                : std::vector<double>{0.06, 0.10, 0.16, 0.22, 0.30};
            const std::vector<double> pair_a1 = broad_turn_candidates
                ? std::vector<double>{0.005, 0.02, 0.04, 0.08, 0.12, 0.16, 0.18, 0.25, 0.34, 0.42, 0.55, 0.70, 0.90, 1.10, 1.30}
                : std::vector<double>{0.55, 0.70, 0.90, 1.10, 1.30};
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
                     (turn_a >= 0.25 && turn_b >= 0.25));
                if (!participates)
                    continue;
                bool same_side_turn_pair =
                    turn_a >= 0.25 && turn_b >= 0.25 &&
                    signed_a * signed_b > 0.0;
                double pair_tol = same_side_turn_pair ? 0.15 : 0.15;
                if (curvesIntersectBeyondSharedEndpointOverlap(
                        curve, *results[oi->second].curve, pair_tol))
                    ++count;
            }
            return count;
        };
        std::unordered_set<std::string> failed_pair_repair_attempts;
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
                if (p.exempt == CrossExemption::StructuralCross &&
                    !same_side_turn_pair && !mixed_turn_pair)
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
                double pair_tol = 0.15;
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
                    std::string attempt_key = cid + "\n" + other_id + "\n" +
                        (same_side_turn_pair ? "same" :
                         (mixed_turn_pair ? "opposite" : "mixed")) + "\n" +
                        (straight_turn_pair ? "straight-turn" : "turn-turn");
                    if (complex_boundary_fast_path &&
                        failed_pair_repair_attempts.count(attempt_key))
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
                    int dbg_total = 0;
                    int dbg_self = 0;
                    int dbg_shape = 0;
                    int dbg_other = 0;
                    int dbg_phys = 0;
                    int dbg_guard = 0;
                    int dbg_fixed = 0;
                    for (const auto& candidate : make_pair_safe_candidates(
                             *conn, shared_turn_endpoint ||
                                    (!straight_like && straight_turn_pair))) {
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
                        CurveRisk candidate_risk =
                            assessCurveRisk(candidate, input, sdf, {}, true);
                        if (candidate_risk.physical()) {
                            ++dbg_phys;
                            continue;
                        }
                        int cand_cross = 0;
                        if (guard_strict_shared_crosses) {
                            cand_cross = strict_shared_cross_count(
                                cid, candidate, result_idx);
                            if (cand_cross > current_cross) {
                                ++dbg_guard;
                                continue;
                            }
                        }
                        int fixed_cross = constrainedFixedCrossCountForId(
                            cid, candidate, results, result_idx, pair_safe_neighbors,
                            kClusterEndpointTol);
                        if (fixed_cross > current_fixed) {
                            ++dbg_fixed;
                            continue;
                        }
                        double target = straight_like ? 1.01 : 1.10;
                        double arc_chord = candidate.arcLength() / chord_len;
                        double score = (guard_strict_shared_crosses ? 1000.0 * cand_cross : 0.0) +
                            std::abs(arc_chord - target) +
                            candidate.maxCurvature(40) + 0.02 * candidate.arcLength();
                        if (!have_best || score < best_score) {
                            best = candidate;
                            best_score = score;
                            have_best = true;
                        }
                        if (complex_boundary_fast_path &&
                            guard_strict_shared_crosses && cand_cross == 0) {
                            out = candidate;
                            return true;
                        }
                    }
                    if (isgDebugPairRepair()) {
                        fprintf(stderr,
                                "[PAIR-REPAIR] %s vs %s same_side=%d straight_turn=%d current_cross=%d current_fixed=%d total=%d self=%d shape=%d other=%d phys=%d guard=%d fixed=%d have=%d score=%.3f\n",
                                cid.c_str(), other_id.c_str(),
                                same_side_turn_pair ? 1 : 0,
                                straight_turn_pair ? 1 : 0,
                                current_cross, current_fixed, dbg_total, dbg_self,
                                dbg_shape, dbg_other, dbg_phys, dbg_guard,
                                dbg_fixed, have_best ? 1 : 0,
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
            auto rebuilt_candidates = make_pair_safe_candidates(*conn, true);
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
        if (!roadedge_boundary_repair_ids.empty()) {
            auto final_boundary_t0 = std::chrono::steady_clock::now();
            auto final_boundary_idx = resultIndexById(results);
            int final_boundary_seen = 0;
            for (const auto& group : plan.batches) {
                for (const auto& cid : group.conn_ids) {
                    ++final_boundary_seen;
                    auto ri = final_boundary_idx.find(cid);
                    if (ri == final_boundary_idx.end() || preserved_fixed_ids.count(cid))
                        continue;
                    if (!roadedge_boundary_repair_ids.count(cid))
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
                    CurveRisk current_boundary_risk =
                        assessCurveRisk(*cc.curve, input, sdf, {}, true);
                    Vec2d current_chord = exit_.first - entry.first;
                    if (current_chord.norm() < 1e-8)
                        continue;
                    Vec2d current_t0 = entry.second.norm() > 1e-8
                        ? entry.second.normalized() : current_chord.normalized();
                    double current_turn_strength = std::abs(
                        cross2d(current_t0, current_chord.normalized()));
                    BezierCurve current_middle =
                        curveWithoutRoadEdgeEndpointSectionsLocal(
                            *cc.curve, input.boundaries);
                    bool current_endpoint_fit_ok =
                        roadEdgeEndpointFitShapeAcceptableLocal(
                            *cc.curve, input.boundaries,
                            current_chord.norm(), current_turn_strength) &&
                        !curveRawIntersectsRoadEdge(
                            current_middle, input.boundaries);
                    if (current_endpoint_fit_ok &&
                        !current_boundary_risk.obstacle &&
                        !current_boundary_risk.fence) {
                        if (cc.status == CurveStatus::Degraded &&
                            cc.violation.reason == "curve intersects boundary away from endpoints") {
                            cc.status = CurveStatus::OK;
                            cc.violation.reason.clear();
                        }
                        continue;
                    }
                }
            }
            if (isgDebugBoundaryRepair()) {
                fprintf(stderr,
                        "[BOUNDARY-REPAIR] final-pass early scanned=%d results=%zu\n",
                        final_boundary_seen, results.size());
            }
            if (isgProfile()) {
                fprintf(stderr, "[ISG_PROFILE] pure final boundary stage %.3f ms\n",
                        std::chrono::duration<double, std::milli>(
                            std::chrono::steady_clock::now() - final_boundary_t0).count());
            }
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
                Vec2d ref_perp = pair.ref_perp.norm() > 1e-8
                    ? pair.ref_perp.normalized() : Vec2d(-t0.y(), t0.x());
                const double lateral_a =
                    (exit_a.first - entry.first).dot(ref_perp);
                const double lateral_b =
                    (exit_b.first - entry.first).dot(ref_perp);
                // 同侧组内离入口扇区中心更远者为外侧；该定义同时适用于
                // 左转和右转，避免把左右转混成一个全局方向。
                const bool a_is_outer = std::abs(lateral_a) > std::abs(lateral_b);
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
                Vec2d ref_perp = pair.ref_perp.norm() > 1e-8
                    ? pair.ref_perp.normalized() : Vec2d(-t1.y(), t1.x());
                auto entry_a = scene.view.entryFrame(cna->entry_lane_id);
                auto entry_b = scene.view.entryFrame(cnb->entry_lane_id);
                const double lateral_a =
                    (entry_a.first - exit_frame.first).dot(ref_perp);
                const double lateral_b =
                    (entry_b.first - exit_frame.first).dot(ref_perp);
                const bool a_is_outer = std::abs(lateral_a) > std::abs(lateral_b);
                auto& outer = a_is_outer ? ca : cb;
                auto& inner = a_is_outer ? cb : ca;
                const Connectivity* outer_conn = a_is_outer ? cna : cnb;
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

                const auto bounds = ordinarySingleCubicHandleBounds(
                    outer_entry.first, outer_entry.second,
                    exit_frame.first, exit_frame.second, true);
                const double target_handle = std::min(
                    bounds.end_max, inner_handle + kHandleOrderMargin);
                if (target_handle <= outer_handle + 1e-6 ||
                    target_handle < bounds.end_min - 1e-6)
                    continue;

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
                        candidate, *inner.curve, kClusterEndpointTol))
                    continue;

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
                bool roadedge_fit_shape_ok =
                    roadEdgeEndpointFitShapeAcceptableLocal(
                        repaired, input.boundaries,
                        natural_chord.norm(), turn_strength);
                bool repaired_shape_ok = turn_strength < 0.25
                    ? (isStraightLikeShapeAcceptable(repaired, natural_chord.norm()) ||
                       roadedge_fit_shape_ok)
                    : (isNonUTurnTurnShapeAcceptable(
                           repaired, natural_chord.norm(), turn_strength) &&
                       !curveHasCurvatureSignFlip(repaired)) ||
                      roadedge_fit_shape_ok;
                CurveRisk repaired_risk = assessCurveRisk(
                    repaired, input, sdf, sampled, true);
                if (repaired_risk.physical() &&
                    !(roadedge_fit_shape_ok &&
                      !repaired_risk.obstacle && !repaired_risk.fence))
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
                if (roadedge_fit_shape_ok &&
                    cc.status == CurveStatus::Degraded &&
                    cc.violation.reason == "curve intersects boundary away from endpoints") {
                    cc.status = CurveStatus::OK;
                    cc.violation.reason.clear();
                }
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

    // 同一进入车道的 U-turn 必须作为一个整体保持半径顺序。逐条后修复会
    // 让外层曲线绕开 Boundary 后反穿内层，因此这里原子化搜索一个全家族
    // 共用的中弧把手系数；仅当整组同时满足物理与同簇拓扑约束时才替换。
    {
        std::unordered_map<LaneId, std::vector<const Connectivity*>> families;
        for (const auto& conn : input.connectivities) {
            if (isGeometricUTurnConn(conn, input))
                families[conn.entry_lane_id].push_back(&conn);
        }
        const std::vector<double> family_arc_alphas =
            {2.0 / 3.0, 0.58, 0.50, 0.44, 0.38, 0.32, 0.28, 0.24, 0.20, 0.16};
        const std::vector<double> family_lead_extras =
            {0.0, 0.5, 1.0, 2.0};
        const std::vector<double> family_stagger_extras =
            {0.0, 0.2, 0.4, 0.6};
        for (auto& family_item : families) {
            auto& family = family_item.second;
            if (family.size() < 2)
                continue;
            std::stable_sort(
                family.begin(), family.end(),
                [&](const Connectivity* a, const Connectivity* b) {
                    double ra = uturnRadiusKey(*a, input);
                    double rb = uturnRadiusKey(*b, input);
                    if (std::abs(ra - rb) > 1e-9)
                        return ra < rb;
                    return a->id < b->id;
                });

            std::unordered_set<ConnId> family_ids;
            for (const Connectivity* conn : family)
                family_ids.insert(conn->id);

            bool family_repaired = false;
            double accepted_alpha = 0.0;
            std::vector<CurvePatchEntry> accepted_patch;
            for (double arc_alpha : family_arc_alphas) {
                for (double family_lead_extra : family_lead_extras) {
                    for (double family_stagger_extra : family_stagger_extras) {
                        std::vector<ConnectivityCurve> candidate_results = results;
                        auto candidate_idx = resultIndexById(candidate_results);
                        std::vector<CurvePatchEntry> family_patch;
                        family_patch.reserve(family.size());
                        bool valid_family = true;
                        ConnId rejected_id;
                        const char* rejected_reason = "";
                        for (const Connectivity* conn : family) {
                            auto ri = candidate_idx.find(conn->id);
                            if (ri == candidate_idx.end()) {
                                valid_family = false;
                                rejected_id = conn->id;
                                rejected_reason = "missing";
                                break;
                            }
                            auto entry = scene.view.entryFrame(conn->entry_lane_id);
                            auto exit_ = scene.view.exitFrame(conn->exit_lane_id);
                            const UTurnFamilyInfo& family_info =
                                pairedFamilySnapshot(conn->id);
                            const double min_lead0 = family_info.lead0;
                            const double min_lead1 = family_info.lead1;
                            double stagger = uturnSharedEndpointStagger(
                                *conn, input, cluster_solver_,
                                direction_cfg_.uturn_alignment_scope);
                            stagger = std::max(0.0, stagger + family_stagger_extra);
                            BezierCurve candidate = UTurnCurveInitializer().buildSegmented(
                                entry.first, entry.second, exit_.first, exit_.second,
                                min_lead0, min_lead1, arc_alpha, stagger,
                                family_lead_extra, family_lead_extra,
                                stagger, stagger);
                            const std::vector<Crosswalk>& clearance_crosswalks =
                                family_info.clearance_crosswalks;
                            double chord_len = (exit_.first - entry.first).norm();
                            CurveRisk candidate_risk = assessCurveRisk(
                                candidate, input, sdf, {}, true, true);
                            if (candidate.numSegments() != 3 || chord_len < 1e-6 ||
                                !segmentLooksStraight(candidate.segs.front()) ||
                                !segmentLooksStraight(candidate.segs.back()) ||
                                candidate.segs[1].maxCurvature(30) <= 0.03 ||
                                candidate.arcLength() / chord_len < 1.35 ||
                                curveSelfIntersectsBusiness(candidate, 1.0) ||
                                !segmentedUTurnMiddleArcClearsCrosswalks(
                                    candidate, clearance_crosswalks) ||
                                candidate_risk.physical()) {
                                valid_family = false;
                                rejected_id = conn->id;
                                rejected_reason = candidate_risk.physical()
                                    ? (candidate_risk.boundary ? "boundary" :
                                       (candidate_risk.obstacle ? "obstacle" : "fence"))
                                    : "shape-or-crosswalk";
                                break;
                            }
                            family_patch.push_back(
                                CurvePatchEntry{conn->id, candidate});
                        }
                        if (valid_family) {
                            valid_family = AtomicCurvePatch().apply(
                                family_patch, candidate_results,
                                [&](ConnectivityCurve& result,
                                    const BezierCurve& candidate) {
                                    setConnectivityCurveGeometry(result, candidate);
                                });
                            if (!valid_family)
                                rejected_reason = "atomic-prepare";
                        }
                        if (!valid_family) {
                            if (isgProfile()) {
                                fprintf(stderr,
                                        "[ISG_PROFILE] shared-entry U-turn family %s alpha=%.2f reject %s:%s\n",
                                        family_item.first.c_str(), arc_alpha,
                                        rejected_id.c_str(), rejected_reason);
                            }
                            continue;
                        }

                        // 同起点/同终点三段式 U-turn 只有在平齐点已经扇出后，
                        // 才裁掉共享连接点附近的首/尾直线继续检查中弧和剩余段。
                        // 若 q0/q1 仍贴得过近，候选直接按同簇冲突拒绝。
                        for (const auto& pair : cluster_solver_.pairs()) {
                            if (pair.exempt == CrossExemption::StructuralCross ||
                                (!family_ids.count(pair.id_a) &&
                                 !family_ids.count(pair.id_b)))
                                continue;
                            auto ia = candidate_idx.find(pair.id_a);
                            auto ib = candidate_idx.find(pair.id_b);
                            if (ia == candidate_idx.end() || ib == candidate_idx.end() ||
                                !candidate_results[ia->second].curve ||
                                !candidate_results[ib->second].curve)
                                continue;
                            BezierCurve ca = *candidate_results[ia->second].curve;
                            BezierCurve cb = *candidate_results[ib->second].curve;
                            bool same_family = family_ids.count(pair.id_a) &&
                                               family_ids.count(pair.id_b);
                            bool forbidden = false;
                            if (same_family) {
                                if (ca.numSegments() < 2 || cb.numSegments() < 2)
                                    continue;
                                bool shared_start =
                                    (ca.startPt() - cb.startPt()).norm() <=
                                    kClusterEndpointTol;
                                bool shared_end =
                                    (ca.endPt() - cb.endPt()).norm() <=
                                    kClusterEndpointTol;
                                if (shared_start) {
                                    double q0_sep =
                                        (ca.segs.front().ctrl[3] -
                                         cb.segs.front().ctrl[3]).norm();
                                    if (q0_sep <= kUTurnSharedLeadFanMinSeparation) {
                                        forbidden = true;
                                    } else {
                                        ca.segs.erase(ca.segs.begin());
                                        cb.segs.erase(cb.segs.begin());
                                    }
                                }
                                if (!forbidden && shared_end) {
                                    double q1_sep =
                                        (ca.segs.back().ctrl[0] -
                                         cb.segs.back().ctrl[0]).norm();
                                    if (q1_sep <= kUTurnSharedLeadFanMinSeparation) {
                                        forbidden = true;
                                    } else {
                                        ca.segs.pop_back();
                                        cb.segs.pop_back();
                                    }
                                }
                                if (!forbidden)
                                    forbidden = curvesIntersectBusiness(
                                        ca, cb, 1.5);
                            } else {
                                forbidden = curvesHaveForbiddenSameClusterIntersection(
                                    ca, cb, kClusterEndpointTol);
                            }
                            if (forbidden) {
                                valid_family = false;
                                rejected_id = family_ids.count(pair.id_a)
                                    ? pair.id_a : pair.id_b;
                                rejected_reason = "cluster";
                                break;
                            }
                        }
                        if (!valid_family) {
                            if (isgProfile()) {
                                fprintf(stderr,
                                        "[ISG_PROFILE] shared-entry U-turn family %s alpha=%.2f reject %s:%s\n",
                                        family_item.first.c_str(), arc_alpha,
                                        rejected_id.c_str(), rejected_reason);
                            }
                            continue;
                        }

                        accepted_alpha = arc_alpha;
                        accepted_patch = std::move(family_patch);
                        family_repaired = true;
                        break;
                    }
                }
                if (family_repaired)
                    break;
            }
            if (!family_repaired)
                continue;
            const bool committed = AtomicCurvePatch().apply(
                accepted_patch, results,
                [&](ConnectivityCurve& result,
                    const BezierCurve& candidate) {
                    setConnectivityCurveGeometry(result, candidate);
                    validate(result, input, sdf);
                });
            if (!committed)
                continue;
            if (isgProfile()) {
                fprintf(stderr,
                        "[ISG_PROFILE] shared-entry U-turn family %s size=%zu arc_alpha=%.2f\n",
                        family_item.first.c_str(), family.size(), accepted_alpha);
            }
        }
    }

    // 最终普通曲线形态审计：前面的同簇修复可能压短某条共享端点
    // 曲线的尾把手，导致短急弯曲率超过允许上限。此处只对无固定形态、
    // 非掉头的单段普通曲线做受约束重建，避免把保序修复留下非法形态。
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

    auto annotate_t0 = std::chrono::steady_clock::now();
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
