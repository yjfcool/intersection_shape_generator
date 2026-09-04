#include "initialization/avoidance_candidate_generator.h"
#include "initialization/curve_initializer_registry.h"
#include "initialization/fixed_shape_initializer.h"
#include "initialization/ordinary_curve_initializer.h"
#include "initialization/uturn_curve_initializer.h"

#include "curve/bezier.h"
#include "curve/curve_utils.h"
#include "preprocessing/uturn_family_builder.h"
#include "utils.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <limits>

namespace isg {

namespace {

// 求最小的 lead，使 |lead * dir + offset| >= target_chord。
// offset 是 U 型调头首尾平齐点的横向错开量：它垂直于 U 轴，但一般不垂直于
// 端点切向 dir，因此在 -dir 上留有分量，会让"最小直行段"的弦长读数比 lead
// 本身短。按二次方程反解出补偿后的 lead，弦长读数才能重新达到 target_chord，
// 而错开量本身仍然垂直于轴线，家族平齐站位分毫不动。
double leadForChordWithLateralOffset(const Vec2d& dir, const Vec2d& offset,
                                    double target_chord) {
    if (target_chord <= 0.0)
        return 0.0;
    const double b = dir.dot(offset);
    const double disc =
        b * b - offset.squaredNorm() + target_chord * target_chord;
    if (disc <= 0.0)
        return std::max(0.0, -b);
    return std::max(0.0, -b + std::sqrt(disc));
}

// 平齐站位处中弧走廊的下限。中弧的弦就是 q1-q0，它必须宽到能撑起
// segmentedUTurnMiddleArcLooksRound 的绝对鼓包下限 0.25m：对称拱的鼓包约为
// 0.75 * arc_alpha * gap，候选搜索的最大把手 arc_alpha = 2.0，因此几何可行的
// 最小走廊是 0.25 / (0.75 * 2.0) ≈ 0.167m。取 0.30m 留出余量，使 arc_alpha
// 从 1.25 起的四个档位都能满足鼓包门禁，而不是只剩最极端的一个。
constexpr double kMinStationCorridor = 0.30;

// 站位口径的走廊不变量：把中弧走廊在"平齐站位处"抬到 kMinStationCorridor 以上。
//
// 家族横向阶梯与内缩上限都按端点口径的 radiusKey 裁定（见 buildSegmented 内
// 的长注释与 UTurnFamilyBuilder::familyLateralLadder），这对家族排序是必须的：
// 两处口径必须同源，否则家族裁定的档位会被就地改写。但端点口径无法反映另一件
// 事——入口切向与出口反向切向一般不严格反平行，两侧各拉出 lead 后 q0/q1 会沿
// 法线漂移 (lead0+lead1) * T0·lateral。对走廊很窄的成员，这个漂移能超过整条
// 走廊本身，于是平齐站位处的有符号走廊塌成 0 甚至反号。
//
// 110000703-u 的掉头 41 就是如此：走廊 0.3099m、两侧 lead 9.32/9.24m、
// T0·lateral = 0.0173，漂移 0.3208m > 0.3099m，站位口径走廊由 +0.3099 翻成
// -0.0108（叠加家族入口档位 +0.0775 后约 -0.088）。此时 q0 与 q1 事实上互换了
// 左右，首直段 p0→q0 与尾直段 q1→p1 在离端点约 9m 处真交叉，
// curveSelfIntersectsBusiness(candidate, 1.0) 因此否决全部 832 个候选
// （实测 self_intersect 832/832），三段式只能退化成单段 ref_arch
// （nseg 3→1、maxκ 107.198），而单段掉头必然违反 shape.crosswalk.segments：
// 这正是"41 未跨越附近人行横道后再掉头"的直接成因。
//
// 修复量只施加在出口侧，沿 +side*lateral 相背拓宽（side 仍取端点口径符号，
// 与家族阶梯同源）。刻意不动入口侧，原因是入口侧承载家族的半径次序：
// 41 的入口档位 +0.0775 正好把它排在 44 的 +0.0387 里侧，一旦为了凑走廊把入口
// 档位吐回去，41 的首直段就会重新落到 44 外侧，44 的首直段从 41 的两条腿之间
// 穿过——那恰是上一轮修掉的 41|44 违约。出口侧没有这个负担：家族里更外层成员
// 的走廊按定义更宽（44 的站位走廊 3.4767m），把 41 的尾直段外扩 0.34m 仍远在
// 44 的尾直段里侧，嵌套次序 q0_lat(44) <= q0_lat(41) < q1_lat(41) < q1_lat(44)
// 保持成立。
//
// 外扩量受 G1/车道对齐余量约束：尾直段总横向偏移不超过 0.14 * lead1，即弦向
// 偏转不超过约 8°，仍满足 dot(lane tangent) >= 0.98 的门禁。出口余量不足时
// 宁可留下欠缺（退化为修复前的行为），也不把欠缺转嫁到入口侧。
//
// 拓宽会经 leadForChordWithLateralOffset 反解出更长的 lead，更长的 lead 又带来
// 更大的漂移，因此按不动点迭代若干轮；漂移增量是 lead 增量的 T0·lateral 倍
// （量级 1e-2），两三轮即收敛。健康成员的站位走廊是米级，一轮判断即退出，
// 其余数据集的形态分毫不动。
template <typename FinalizeLeads>
void enforceStationCorridorFloor(
    const Vec2d& p0, const Vec2d& p1, const Vec2d& lateral,
    const Vec2d& q0, const Vec2d& q1, Vec2d& q0_offset, Vec2d& q1_offset,
    const FinalizeLeads& finalize_leads) {
    // 诊断开关：ISG_STATION_CORRIDOR_FLOOR=0 关闭本修复，用于源码级 A/B 消融。
    static const bool enabled = [] {
        const char* v = std::getenv("ISG_STATION_CORRIDOR_FLOOR");
        return v == nullptr || v[0] != '0';
    }();
    if (!enabled)
        return;
    const double endpoint_corridor = (p1 - p0).dot(lateral);
    if (std::abs(endpoint_corridor) < 1e-9)
        return;
    const double side = endpoint_corridor >= 0.0 ? 1.0 : -1.0;
    for (int iteration = 0; iteration < 4; ++iteration) {
        const double station_corridor =
            ((q1 + q1_offset) - (q0 + q0_offset)).dot(lateral) * side;
        // 只在走廊真正塌掉（有符号走廊 <= 0，即 q0/q1 已互换左右）时介入。
        // 早期版本按 "< kMinStationCorridor" 介入，会把走廊仍为正、只是偏窄的
        // 健康成员一并外扩，实测把 intersection_cross 的 39/40 与 100000412 家族
        // 中一名成员从三段式打回单段（外扩改变了候选的曲率/圆度门禁结果）。
        // 走廊为正时首尾直段不换侧，不存在自交死结，没有介入的理由。
        if (iteration == 0 && station_corridor > 0.0)
            return;
        const double deficit = kMinStationCorridor - station_corridor;
        if (deficit <= 1e-9)
            return;
        const double lead1_len = (p1 - q1).norm();
        const double cap = 0.14 * lead1_len;
        const double current = q1_offset.dot(lateral) * side;
        const double target = std::min(cap, current + deficit);
        const double add = target - current;
        if (add <= 1e-9)
            return;  // G1 余量已用尽：保持修复前行为，不转嫁到入口侧。
        q1_offset += side * add * lateral;
        finalize_leads();
    }
}

}  // namespace

// 固有形态初始化：逐段转换输入折线，不引入额外平滑或避让。
bool FixedShapeInitializer::hasGeometry(
    const Connectivity& connectivity) const {
    return connectivity.geometry.points.size() >= 2;
}

BezierCurve FixedShapeInitializer::build(
    const Connectivity& connectivity) const {
    BezierCurve curve;
    if (!hasGeometry(connectivity))
        return curve;
    const std::vector<Vec2d> points =
        toVec2dArray(connectivity.geometry.points);
    for (std::size_t i = 0; i + 1 < points.size(); ++i) {
        const Vec2d chord = points[i + 1] - points[i];
        if (chord.norm() < 1e-8)
            continue;
        const Vec2d direction = chord.normalized();
        curve.segs.push_back(makeCubicG1(
            points[i], direction, points[i + 1], direction, 1.0 / 3.0));
    }
    return curve;
}

// 普通曲线初始化：严格沿端点切向构造单段三次 Bezier 候选。
BezierCurve OrdinaryCurveInitializer::buildPreferredSingleCubic(
    const Vec2d& entry_point, const Vec2d& entry_tangent,
    const Vec2d& exit_point, const Vec2d& exit_tangent) const {
    BezierCurve curve;
    const Vec2d chord = exit_point - entry_point;
    const double chord_len = chord.norm();
    if (chord_len < 1e-8)
        return curve;

    const OrdinarySingleCubicHandleBounds bounds =
        ordinarySingleCubicHandleBounds(
            entry_point, entry_tangent, exit_point, exit_tangent, true);
    const Vec2d chord_dir = chord / chord_len;
    const bool straight_like =
        bounds.start_dir.dot(bounds.end_dir) > 0.90 &&
        std::abs(cross2d(bounds.start_dir, bounds.end_dir)) < 0.25 &&
        std::abs(cross2d(bounds.start_dir, chord_dir)) < 0.25;

    // 直行把手采用等弦中点对应的统一长度；控制点仍分别落在
    // 首、尾切线轴上，因此切线不一致时不强制 C1=C2。
    double start_handle = chord_len / 2.0;
    double end_handle = chord_len / 2.0;
    const bool stable_turn_intersection =
        !straight_like && bounds.has_direction_intersection &&
        bounds.direction_intersection_start <= 1.5 * chord_len &&
        bounds.direction_intersection_end <= 1.5 * chord_len;
    if (stable_turn_intersection) {
        start_handle = (2.0 / 3.0) *
            bounds.direction_intersection_start;
        end_handle = (2.0 / 3.0) *
            bounds.direction_intersection_end;
    } else if (!straight_like) {
        // 平行、反向或过远交点不适合作为形态基准，回退到有界自然弧。
        start_handle = 0.4 * chord_len;
        end_handle = 0.4 * chord_len;
    }

    BezierSegment segment;
    segment.ctrl[0] = entry_point;
    segment.ctrl[1] = entry_point + bounds.start_dir * start_handle;
    segment.ctrl[2] = exit_point - bounds.end_dir * end_handle;
    segment.ctrl[3] = exit_point;
    curve.segs.push_back(segment);
    constrainOrdinarySingleCubicControls(
        curve, entry_point, bounds.start_dir,
        exit_point, bounds.end_dir, true);

    // 2/3 交点升阶是首选形态而不是无条件硬编码。短急弯若直接采用
    // 该比例可能产生过高曲率，回退到有界自然弧，后续仍可由候选搜索
    // 在不破坏同簇和物理约束的前提下微调。
    if (!straight_like) {
        const double arc_chord = curve.arcLength() / chord_len;
        const double max_curvature = curve.maxCurvature(40);
        const double max_allowed_curvature = chord_len <= 10.0 ? 6.0 : 2.5;
        if (arc_chord < 1.005 || arc_chord > 1.35 ||
            max_curvature > max_allowed_curvature) {
            curve.segs.front() = makeCubicG1(
                entry_point, bounds.start_dir,
                exit_point, bounds.end_dir, 0.4);
            constrainOrdinarySingleCubicControls(
                curve, entry_point, bounds.start_dir,
                exit_point, bounds.end_dir, true);
        }
    }
    return curve;
}

BezierCurve OrdinaryCurveInitializer::buildSingleCubic(
    const Vec2d& entry_point, const Vec2d& entry_tangent,
    const Vec2d& exit_point, const Vec2d& exit_tangent, double alpha) const {
    BezierCurve curve;
    const Vec2d fallback = exit_point - entry_point;
    const Vec2d start = entry_tangent.norm() > 1e-8
        ? entry_tangent.normalized()
        : (fallback.norm() > 1e-8 ? fallback.normalized() : Vec2d(1, 0));
    const Vec2d end = exit_tangent.norm() > 1e-8
        ? exit_tangent.normalized() : start;
    curve.segs.push_back(
        makeCubicG1(entry_point, start, exit_point, end, alpha));
    return curve;
}

std::vector<BezierCurve> OrdinaryCurveInitializer::buildAlphaCandidates(
    const Vec2d& entry_point, const Vec2d& entry_tangent,
    const Vec2d& exit_point, const Vec2d& exit_tangent,
    const std::vector<double>& alphas) const {
    std::vector<BezierCurve> candidates;
    candidates.reserve(alphas.size());
    for (double alpha : alphas) {
        candidates.push_back(buildSingleCubic(
            entry_point, entry_tangent, exit_point, exit_tangent, alpha));
    }
    return candidates;
}

// 物理避让初始化：将有序路点转换为保持首尾切向的分段曲线。
BezierCurve AvoidanceCandidateGenerator::buildWaypointCurve(
    const std::vector<Vec2d>& points, const Vec2d& start_tangent,
    const Vec2d& end_tangent) const {
    if (points.size() < 2)
        return BezierCurve();
    std::vector<Vec2d> tangents(points.size(), Vec2d(1, 0));
    tangents.front() = start_tangent.norm() > 1e-8
        ? start_tangent.normalized() : (points[1] - points[0]).normalized();
    tangents.back() = end_tangent.norm() > 1e-8
        ? end_tangent.normalized()
        : (points.back() - points[points.size() - 2]).normalized();
    for (std::size_t i = 1; i + 1 < points.size(); ++i) {
        const Vec2d direction = points[i + 1] - points[i - 1];
        tangents[i] = direction.norm() > 1e-8
            ? direction.normalized() : tangents[i - 1];
    }
    return makeCurveFromKnots(points, tangents, 0.34);
}

// U-turn 初始化：构造轴向平齐的单段或三段几何表达。

BezierCurve UTurnCurveInitializer::buildAligned(
    const Vec2d& p0, const Vec2d& t0, const Vec2d& p1, const Vec2d& t1,
    const Vec2d& offset_dir, double offset_m, double handle_scale,
    double min_lead0, double min_lead1) const {
    Vec2d T0 = t0.norm() > 1e-8 ? t0.normalized() : Vec2d(1, 0);
    Vec2d T1 = t1.norm() > 1e-8 ? t1.normalized() : -T0;
    Vec2d axis = T0 - T1;
    if (axis.norm() < 1e-8)
        axis = T0;
    axis.normalize();
    if (axis.dot(T0) < 0.0)
        axis = -axis;

    double fwd_bias = 0.0;
    double lat_bias = 0.0;
    if (offset_m > 0.0 && offset_dir.norm() > 1e-8) {
        Vec2d dir = offset_dir.normalized();
        fwd_bias = offset_m * dir.dot(axis);
        Vec2d lat_dir{-axis.y(), axis.x()};
        lat_bias = offset_m * dir.dot(lat_dir);
    }

    BezierCurve curve;
    curve.segs.push_back(makeAlignedUTurnCubic(
        p0, T0, p1, T1, handle_scale, fwd_bias, lat_bias,
        min_lead0, min_lead1));
    return curve;
}

BezierCurve UTurnCurveInitializer::buildSegmented(
    const Vec2d& p0, const Vec2d& t0, const Vec2d& p1, const Vec2d& t1,
    double min_lead0, double min_lead1, double arc_alpha,
    double aligned_point_stagger, double lead0_extra_after_align,
    double lead1_extra_after_align, double aligned_entry_stagger,
    double aligned_exit_stagger, double entry_lateral_bias,
    double exit_lateral_bias, double family_stagger_step) const {
    Vec2d T0 = t0.norm() > 1e-8 ? t0.normalized() : Vec2d(1, 0);
    Vec2d T1 = t1.norm() > 1e-8 ? t1.normalized() : -T0;
    Vec2d exit_back = -T1;
    Vec2d axis = T0 + exit_back;
    if (axis.norm() < 1e-8)
        axis = T0;
    axis.normalize();
    if (axis.dot(T0) < 0.0)
        axis = -axis;

    double s0 = p0.dot(axis);
    double s1 = p1.dot(axis);
    double common_s = std::max(s0, s1);
    double c0 = std::max(0.2, T0.dot(axis));
    double c1 = std::max(0.2, exit_back.dot(axis));
    double lead0 = std::max(0.0, (common_s - s0) / c0);
    double lead1 = std::max(0.0, (common_s - s1) / c1);
    lead0 = std::max(lead0, min_lead0);
    lead1 = std::max(lead1, min_lead1);

    double s0_new = p0.dot(axis) + lead0 * T0.dot(axis);
    double s1_new = p1.dot(axis) + lead1 * exit_back.dot(axis);
    double common_s_new = std::max(s0_new, s1_new);
    if (s0_new < common_s_new - 1e-6)
        lead0 += (common_s_new - s0_new) / c0;
    if (s1_new < common_s_new - 1e-6)
        lead1 += (common_s_new - s1_new) / c1;

    double extra_axis_advance = std::max(
        std::max(0.0, lead0_extra_after_align) * T0.dot(axis),
        std::max(0.0, lead1_extra_after_align) * exit_back.dot(axis));
    if (extra_axis_advance > 0.0) {
        lead0 += extra_axis_advance / c0;
        lead1 += extra_axis_advance / c1;
    }

    Vec2d q0 = p0 + lead0 * T0;
    Vec2d q1 = p1 + lead1 * exit_back;
    double entry_stagger = std::isfinite(aligned_entry_stagger)
        ? aligned_entry_stagger : aligned_point_stagger;
    double exit_stagger = std::isfinite(aligned_exit_stagger)
        ? aligned_exit_stagger : aligned_point_stagger;
    // 入口侧与出口侧的分档量必须各自独立生效，不能先合并成一个共享标量。
    // 掉头家族是一条交替由"共享入口端点"和"共享出口端点"连接起来的链，而内缩
    // 在入口侧沿 +side*lateral、在出口侧沿 -side*lateral 施加：把某成员压小以
    // 让开它共享入口的退化邻居，用共享标量就会同时把它在出口侧也压小，于是它
    // 相对共享出口的那个邻居反号。实测 100000598 的 58|57、60|59、64|63、
    // 15|14、17|16 与 100000699、100000643 的同型对全部由此产生。
    // 两侧的分档裁定见 UTurnFamilyBuilder::familyLateralLadder。
    // 横向错开与横向偏置都先只累积偏移向量，不立刻改写 q0/q1：偏移沿 U 轴法线，
    // 会在端点切向上留下 -dir 分量，直接平移会让最小直行段的弦长读数短于
    // min_lead，触发候选搜索的最小直段门禁，而任何正向加长又必然越过家族平齐
    // 站位，两道门禁互斥把候选集夹空，三段式只能退化成单段。因此先按偏移量反解
    // 补偿后的 lead，再统一平移。
    const Vec2d lateral{-axis.y(), axis.x()};
    Vec2d q0_offset(0.0, 0.0);
    Vec2d q1_offset(0.0, 0.0);
    // 家族错开方向：入口侧沿 +stagger_side*lateral，出口侧沿 -stagger_side*lateral。
    // stagger_side 已带上分档量的符号（分档量可为负，表示相背拓宽），因此它
    // 始终指向"分档量增大"的方向，即家族内更靠里的槽位方向。
    // 仅当分档量本身被完整实现（未被 lateral_gap / G1 上限压掉）时才认为"半径序
    // 正在由错开量承载"，此时横向偏置不得反向抵消它。
    double stagger_side = 0.0;
    bool entry_stagger_orders_family = false;
    bool exit_stagger_orders_family = false;
    double realized_entry_inset = 0.0;
    double realized_exit_inset = 0.0;
    if (std::abs(entry_stagger) > 1e-9 || std::abs(exit_stagger) > 1e-9) {
        Vec2d chord = q1 - q0;
        double gap = chord.norm();
        if (gap > 1e-6) {
            // 走廊宽度与内缩方向都必须按"车道端点之间"的横向间距取，不能按已经
            // 拉出直段后 q0/q1 之间的横向间距取：家族横向阶梯用的是端点口径的
            // radiusKey（cap = 0.25 * radiusKey），两处口径必须同源，否则家族
            // 裁定的档位会在这里被就地改写。
            //
            // 入口切向与出口反向切向一般不严格反平行（T0·lateral != 0），两侧各
            // 拉出直段后 q0/q1 会沿法线漂移 (lead0+lead1) * T0·lateral。对走廊
            // 很窄的退化成员，这个漂移能超过整条走廊本身：110000703-u 的掉头 41
            // 走廊 0.3099m、两侧 lead 9.32/9.24m、T0·lateral=0.0173，漂移
            // 0.3208m > 0.3099m，于是平齐站位处的有符号间距由 +0.3099 翻成
            // -0.0108。沿用 q0/q1 口径会同时错三件事：
            //   1) side 反号——家族分档被反着施加，小径反而被推到大径外侧；
            //   2) 走廊上限塌成 0.25*0.0108=0.0027——家族裁定的 0.0775 只落实
            //      3.5%，阶梯形同虚设；
            //   3) 因为分档没有完整落实，*_stagger_orders_family 为假，槽位钳制
            //      整体关闭，+0.2309 的出口横向偏置畅通无阻。
            // 三者叠加把 41 的入口直段推到 44 外侧（实测 -0.0027 vs 44 的
            // +0.0387），而它的出口端点仍在里侧，于是必然穿越 44：首中间控制点
            // 偏出 44 控制面 0.0421m，中段把手 seg1.P1-P2 与 44 的 seg1.P0-P1 相交。
            //
            // 端点口径下本函数的上限与家族上限完全一致，因此这里退化成幂等校验：
            // 家族裁定的档位被原样落实，只有 G1 余量（真正的局部约束）可以修剪它。
            double corridor = (p1 - p0).dot(lateral);
            if (std::abs(corridor) < 1e-6)
                corridor = chord.dot(lateral);
            double lateral_gap = std::abs(corridor);
            if (lateral_gap > 1e-6) {
                double lead0_len = (q0 - p0).norm();
                double lead1_len = (p1 - q1).norm();
                // 正向内缩会压窄中弧走廊，受走廊与 G1 双重上限；负向拓宽不会
                // 让 q0/q1 互换，只受 G1 上限约束。两侧的走廊上限各取
                // 0.25 * lateral_gap，因此两侧内缩之和最多吃掉走廊的一半，
                // 中弧始终保留可辨识的宽度（与合并成单一标量时的口径一致）。
                // G1 余量按各侧自身的直段长度取，短直段一侧不再被长的一侧放宽。
                const double g1_cap0 = 0.14 * lead0_len;
                const double g1_cap1 = 0.14 * lead1_len;
                double entry_inset = entry_stagger;
                double exit_inset = exit_stagger;
                if (entry_inset > 0.0) {
                    entry_inset = std::min(
                        entry_inset, std::min(0.25 * lateral_gap, g1_cap0));
                } else {
                    entry_inset = std::max(entry_inset, -g1_cap0);
                }
                if (exit_inset > 0.0) {
                    exit_inset = std::min(
                        exit_inset, std::min(0.25 * lateral_gap, g1_cap1));
                } else {
                    exit_inset = std::max(exit_inset, -g1_cap1);
                }
                double side = corridor >= 0.0 ? 1.0 : -1.0;
                // 相向内缩：打散同家族多条调头的首尾直行段重叠。
                q0_offset += side * entry_inset * lateral;
                q1_offset -= side * exit_inset * lateral;
                const double dominant =
                    std::abs(entry_inset) >= std::abs(exit_inset)
                        ? entry_inset : exit_inset;
                stagger_side = dominant >= 0.0 ? side : -side;
                // "档位是否完整落实"用容差判定而不是精确相等：上限与家族上限
                // 同源后两者只差浮点末位，一旦用 1e-9 精确比较，末位不同就会
                // 静默关掉整个槽位钳制（110000703-u 的 41 正是如此）。
                entry_stagger_orders_family =
                    std::abs(entry_inset - entry_stagger) < 1e-6 &&
                    std::abs(entry_inset) > 1e-9;
                exit_stagger_orders_family =
                    std::abs(exit_inset - exit_stagger) < 1e-6 &&
                    std::abs(exit_inset) > 1e-9;
                realized_entry_inset = std::abs(entry_inset);
                realized_exit_inset = std::abs(exit_inset);
            }
        }
    }
    // 独立的有符号横向偏置：把某一侧直行段整体推离共享车道中心线，用于与同簇
    // 的直行/转向曲线彻底分侧。与内缩不同，两侧符号由调用方分别指定。
    //
    // 但它不能把成员挤出自己的家族槽位。同入掉头家族的左右次序完全由
    // stagger = step * reverse_radius_rank + 家族统一平移 承载（半径小的排在
    // 里侧），相邻名次的横向站位恰好相差一个 step——平移对全体同量同向，不改变
    // 这个间距。一旦某个成员的偏置在反方向上超过半个 step，它就会越过上一名并
    // 与之相交；而是否相交还取决于两条曲线谁先生成，所以候选搜索里逐条做的
    // 同簇审计（当时后一名尚未生成或形态不同）抓不住它。
    // 110003285 的掉头 51/52 即属此类：48/49/50/51 的入口直段横向偏移是均匀的
    // 0/-0.249/-0.498/-0.746 等差列（step=0.25），52 本应到 -0.995，却被 +0.377
    // 的反向偏置（0.50 x 5% x lead=15.07）推回到 -0.620，越过 51 的 -0.746 而相交。
    //
    // 因此这里按槽位宽度钳制反向偏置：钳到 0.4 个 step 以内，成员就不可能挤进
    // 相邻名次的槽位，家族次序与生成顺序无关地成立。钳制而不是归零，是因为
    // 归零会让本来只需要厘米级反向微调的成员完全失去分侧手段。
    // 步长不到 10cm 时（纯几何输入的 1cm 分档）家族次序并不依赖横向站位，
    // 分档只是打散直段重叠的象征性微移，此时不钳制，让偏置继续承担与同簇
    // 直行/转向兄弟的分侧职责。
    //
    // 钳制曾经的代价（100000412 的 34/36 由不相交变为相交）已由家族层面裁定的
    // 横向阶梯从根源消除：走廊只有 0.34m 的退化成员不再让与它共享端点的外层邻居
    // 独享完整分档，而是由家族按共享端点拓扑分侧压掉冲突档位，阶梯因此天然单调，
    // 34 与 36 的分侧不再依赖任何一次越槽偏置。
    // 见 UTurnFamilyBuilder::familyLateralLadder。
    double entry_bias_eff = entry_lateral_bias;
    double exit_bias_eff = exit_lateral_bias;
    // 诊断时允许关闭槽位钳制，以验证“家族分档”和边界绕行偏置是否互相
    // 限制；正式生成仍保持默认钳制，只有家族级公共候选显式传入 0 步长时
    // 才允许改变整族的相对平移。
    if (family_stagger_step >= 0.10 &&
        std::getenv("ISG_NO_SLOT_CLAMP") == nullptr) {
        // 槽位宽度按"名义步长"和"本侧实际落实的档位"取小。等差列成员的档位本身
        // 就 >= 一个步长，取小后与只用步长完全一致；只有被家族压缩过或走廊退化的
        // 成员才收紧。这类成员与外层邻居的实际间距是家族级联给出的
        // sep = min(step, budget/组内规模)，可以远小于名义步长——110000703-u 的
        // 41 与 44 之间只有 0.0775-0.0387=0.0388m，按名义步长算出的 0.1m 槽位比
        // 真实间距还宽 2.6 倍，合规的反向偏置照样能把它挤过外层邻居。
        // 用本侧落实的档位（41 为 0.0775）作上界，槽位收到 0.031m < 0.0388m，
        // 越槽在几何上不再可能。
        const double entry_slot =
            0.4 * std::min(family_stagger_step, realized_entry_inset);
        const double exit_slot =
            0.4 * std::min(family_stagger_step, realized_exit_inset);
        if (entry_stagger_orders_family && entry_bias_eff * stagger_side < 0.0) {
            entry_bias_eff = std::copysign(
                std::min(std::abs(entry_bias_eff), entry_slot), entry_bias_eff);
        }
        if (exit_stagger_orders_family && exit_bias_eff * -stagger_side < 0.0) {
            exit_bias_eff = std::copysign(
                std::min(std::abs(exit_bias_eff), exit_slot), exit_bias_eff);
        }
    }
    q0_offset += entry_bias_eff * lateral;
    q1_offset += exit_bias_eff * lateral;
    // 按当前偏移量反解补偿后的 lead，并把两侧重新对齐到共同轴向站位。
    // 横向偏移垂直于轴，不改变站位本身，但补偿后的 lead 会加长，因此每次
    // 改动偏移量都必须重跑一遍——下面的走廊下限修复正是靠这一点迭代收敛。
    auto finalize_leads = [&]() {
        if (q0_offset.norm() <= 1e-12 && q1_offset.norm() <= 1e-12)
            return;
        lead0 = std::max(lead0, leadForChordWithLateralOffset(
                                    T0, q0_offset, min_lead0));
        lead1 = std::max(lead1, leadForChordWithLateralOffset(
                                    exit_back, q1_offset, min_lead1));
        // 补偿后两侧 lead 可能不再等长，需重新对齐到共同轴向站位，否则
        // q0/q1 轴向平齐门禁失效。横向偏移本身垂直于轴，不影响轴向站位。
        const double a0 = p0.dot(axis) + lead0 * T0.dot(axis);
        const double a1 = p1.dot(axis) + lead1 * exit_back.dot(axis);
        const double a_common = std::max(a0, a1);
        if (a0 < a_common - 1e-12)
            lead0 += (a_common - a0) / c0;
        if (a1 < a_common - 1e-12)
            lead1 += (a_common - a1) / c1;
        q0 = p0 + lead0 * T0;
        q1 = p1 + lead1 * exit_back;
    };
    finalize_leads();
    enforceStationCorridorFloor(p0, p1, lateral, q0, q1,
                                q0_offset, q1_offset, finalize_leads);
    q0 += q0_offset;
    q1 += q1_offset;

    BezierCurve curve;
    const double lead_eps = 0.20;
    Vec2d lead_dir0 = (q0 - p0).norm() > 1e-8 ? (q0 - p0).normalized() : T0;
    Vec2d lead_dir1 = (p1 - q1).norm() > 1e-8 ? (p1 - q1).normalized() : T1;
    if (lead0 > lead_eps)
        curve.segs.push_back(makeCubicG1(p0, lead_dir0, q0, lead_dir0, 1.0 / 3.0));
    if ((q1 - q0).norm() > 1e-6)
        curve.segs.push_back(makeCubicG1(q0, lead_dir0, q1, lead_dir1, arc_alpha));
    if (lead1 > lead_eps)
        curve.segs.push_back(makeCubicG1(q1, lead_dir1, p1, lead_dir1, 1.0 / 3.0));
    if (curve.empty())
        curve.segs.push_back(makeAlignedUTurnCubic(
            p0, T0, p1, T1, 1.0, 0.0, 0.0, min_lead0, min_lead1));
    return curve;
}

// 注册表保持 U-turn、fixed shape、普通曲线的既有优先级和候选顺序。
CurveInitializationOptions::CurveInitializationOptions()
    : ordinary_alphas(1, 0.4),
      allow_fixed_shape(true),
      uturn_min_lead0(0.0),
      uturn_min_lead1(0.0),
      uturn_arc_alpha(2.0 / 3.0),
      uturn_aligned_point_stagger(0.0),
      uturn_lead0_extra_after_align(0.0),
      uturn_lead1_extra_after_align(0.0),
      uturn_aligned_entry_stagger(std::numeric_limits<double>::quiet_NaN()),
      uturn_aligned_exit_stagger(std::numeric_limits<double>::quiet_NaN()) {}

void CurveInitializationOptions::applyUTurnFamily(
    const UTurnFamilyInfo& family) {
    if (!family.geometric_uturn)
        return;
    uturn_min_lead0 = family.lead0;
    uturn_min_lead1 = family.lead1;
}

std::vector<CurveCandidate> CurveInitializerRegistry::build(
    const CurveGenerationContext& context,
    const CurveInitializationOptions& options) const {
    std::vector<CurveCandidate> candidates;
    if (!context.connectivity)
        return candidates;

    const Vec2d& p0 = context.entry.first;
    const Vec2d& t0 = context.entry.second;
    const Vec2d& p1 = context.exit.first;
    const Vec2d& t1 = context.exit.second;
    const bool geometric_uturn =
        t0.norm() > 1e-8 && t1.norm() > 1e-8 &&
        t0.normalized().dot(t1.normalized()) < -0.5;

    if (geometric_uturn) {
        CurveCandidate candidate;
        candidate.origin = CandidateOrigin::SegmentedUTurn;
        candidate.curve = UTurnCurveInitializer().buildSegmented(
            p0, t0, p1, t1,
            options.uturn_min_lead0, options.uturn_min_lead1,
            options.uturn_arc_alpha, options.uturn_aligned_point_stagger,
            options.uturn_lead0_extra_after_align,
            options.uturn_lead1_extra_after_align,
            options.uturn_aligned_entry_stagger,
            options.uturn_aligned_exit_stagger);
        candidates.push_back(candidate);
        return candidates;
    }

    const FixedShapeInitializer fixed_initializer;
    if (options.allow_fixed_shape && context.connectivity->fixed_shape &&
        fixed_initializer.hasGeometry(*context.connectivity)) {
        CurveCandidate candidate;
        candidate.origin = CandidateOrigin::FixedShape;
        candidate.curve = fixed_initializer.build(*context.connectivity);
        candidate.preserves_fixed_shape = true;
        candidates.push_back(candidate);
        return candidates;
    }

    const OrdinaryCurveInitializer ordinary_initializer;
    std::vector<BezierCurve> ordinary;
    ordinary.push_back(ordinary_initializer.buildPreferredSingleCubic(
        p0, t0, p1, t1));
    const std::vector<BezierCurve> alpha_candidates =
        ordinary_initializer.buildAlphaCandidates(
            p0, t0, p1, t1, options.ordinary_alphas);
    ordinary.insert(ordinary.end(), alpha_candidates.begin(), alpha_candidates.end());
    candidates.reserve(ordinary.size());
    for (std::size_t i = 0; i < ordinary.size(); ++i) {
        CurveCandidate candidate;
        candidate.origin = CandidateOrigin::NaturalSingleCubic;
        candidate.curve = ordinary[i];
        candidate.quality_score = static_cast<double>(i);
        candidates.push_back(candidate);
    }
    return candidates;
}

}  // 命名空间 isg
