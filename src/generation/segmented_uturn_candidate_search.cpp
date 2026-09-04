#include "generation/segmented_uturn_candidate_search.h"

#include "curve/curve_utils.h"
#include "generation/uturn_shape.h"
#include "initialization/uturn_curve_initializer.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <utility>
#include <vector>

namespace isg {

static bool debugSegmentedUTurnSearch() {
    static const bool enabled = std::getenv("ISG_DEBUG_UTURN") != nullptr;
    return enabled;
}

// 量出兄弟曲线在某一侧直行段跨度内的"折算横向漂移"，即让该直段整体处在兄弟
// 曲线外侧所需的最小横向偏置。
//
// 直段从车道端点 p 出发、长度 lead_len，横向偏置 m 让它在轴向 d 处横向偏出
// m*d/lead_len。要在整段范围内都不被兄弟曲线越过，必须对跨度内每个采样点满足
// m*d/lead_len >= lat(d)，即 m >= max_d (lead_len/d)*lat(d)。同入/同出的兄弟
// 曲线在共享端点与车道切向一阶贴合，横向漂移随距离近似二次增长，因此最大值
// 通常落在直段末端；这里仍按逐点折算取上确界，避免对漂移形状做假设。
static double requiredLateralBiasToClearSibling(
    const BezierCurve& sibling, const Vec2d& origin, const Vec2d& axis,
    const Vec2d& lateral, double side, double lead_len) {
    if (lead_len <= 1e-6)
        return 0.0;
    double need = 0.0;
    const std::vector<Vec2d> pts = sibling.sampleByArcLength(200);
    for (const Vec2d& pt : pts) {
        const Vec2d rel = pt - origin;
        const double d = rel.dot(axis);
        // 轴向 0.5m 以内的采样点全部落在共享端点的一阶贴合区，折算系数
        // lead_len/d 会被噪声放大，跳过。
        if (d < 0.5 || d > lead_len)
            continue;
        const double lat = side * rel.dot(lateral);
        if (lat <= 0.0)
            continue;
        need = std::max(need, lat * lead_len / d);
    }
    return need;
}

bool SegmentedUTurnCandidateSearch::search(
    const Vec2d& p0, const Vec2d& t0, const Vec2d& p1, const Vec2d& t1,
    const IntersectionInput& input,
    const std::vector<SampledSiblingCurve>& sampled_siblings,
    const SegmentedUTurnAuditor& audit, bool include_fence,
    BezierCurve& curve, double min_lead0, double min_lead1,
    const std::vector<Crosswalk>* crosswalks_for_clearance,
    double aligned_point_stagger,
    double base_lead0_extra_after_align,
    double base_lead1_extra_after_align,
    double aligned_family_station,
    double aligned_entry_stagger,
    double aligned_exit_stagger,
    double family_stagger_step) const {
    double chord_len = (p1 - p0).norm();
    if (chord_len < 1e-6)
        return false;
    double segmented_max_curvature =
        segmentedUTurnMaxCurvatureLimit(p0, t0, p1, t1);
    bool debug = debugSegmentedUTurnSearch();

    // 使用完整Boundary集合评估三段式候选。curveBoundarySafety会只裁掉
    // 首尾连续贴行的局部段；不能再整条删除端点粘连边缘，否则中间弧
    // 再次穿回同一RoadEdge/Other边缘也会被错误豁免。
    SegmentedUTurnAudit current_risk = audit(curve, include_fence);

    struct LeadExtraPair {
        double extra0 = 0.0;
        double extra1 = 0.0;
    };
    const std::vector<double> lead_extras =
        {0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0};
    std::vector<LeadExtraPair> lead_pairs;
    lead_pairs.reserve(lead_extras.size() * lead_extras.size());
    for (double extra0 : lead_extras) {
        for (double extra1 : lead_extras) {
            LeadExtraPair pair;
            pair.extra0 = extra0;
            pair.extra1 = extra1;
            lead_pairs.push_back(pair);
        }
    }
    std::sort(
        lead_pairs.begin(), lead_pairs.end(),
        [](const LeadExtraPair& a, const LeadExtraPair& b) {
            double sum_a = a.extra0 + a.extra1;
            double sum_b = b.extra0 + b.extra1;
            if (std::abs(sum_a - sum_b) > 1e-9)
                return sum_a < sum_b;
            return std::abs(a.extra0 - a.extra1) <
                   std::abs(b.extra0 - b.extra1);
        });
    // 大把手用于小半径掉头形成接近半圆的单段中弧；较大的值允许在
    // 不改变首尾平齐站位的情况下把中弧向路口内侧推开 RoadEdge。
    const std::vector<double> arc_alphas =
        {2.0 / 3.0, 0.75, 0.85, 1.0, 1.25, 1.50, 1.75, 2.0,
         0.58, 0.50, 0.38, 0.28, 0.16};

    Vec2d T0 = t0.norm() > 1e-8 ? t0.normalized() : Vec2d(1, 0);
    Vec2d T1 = t1.norm() > 1e-8 ? t1.normalized() : -T0;
    Vec2d exit_back = -T1;
    Vec2d axis = T0 + exit_back;
    if (axis.norm() < 1e-8)
        axis = T0;
    axis.normalize();
    if (axis.dot(T0) < 0.0)
        axis = -axis;

    BezierCurve best;
    bool have_best = false;
    int best_cross = std::numeric_limits<int>::max();
    double best_score = std::numeric_limits<double>::infinity();
    double best_extra_sum = std::numeric_limits<double>::infinity();
    int rejected_min_lead = 0;
    int rejected_curvature = 0;
    int rejected_sibling_cross = 0;
    // 其余否决原因单独计数：候选集枯竭时必须能直接读出是哪一道门禁清空了
    // 候选，而不是只看到"no valid candidate"。仅在 ISG_DEBUG_UTURN 下打印。
    int rejected_shape = 0;      // 空候选或段数不为 3
    int rejected_self = 0;       // 自交
    int rejected_axis = 0;       // q0/q1 轴向未平齐
    int rejected_station = 0;    // 越过家族平齐站位
    int rejected_round = 0;      // 弧长比不足或中弧不圆整
    int rejected_crosswalk = 0;  // 中弧进入人行横道
    int rejected_physical = 0;   // Boundary/障碍/围栏等物理违约
    const std::vector<Crosswalk>& clearance_crosswalks =
        crosswalks_for_clearance ? *crosswalks_for_clearance : input.crosswalks;
    // 与家族分档反向的横向偏置在 buildSegmented 内被钳制在槽位内，钳制强度由
    // 家族分档步长决定。置 ISG_NO_SLOT_CLAMP 可整体关闭钳制，用于诊断时做
    // A/B：钳制是否是某条曲线形态变化的成因，一次运行即可判定。
    const double clamp_step =
        std::getenv("ISG_NO_SLOT_CLAMP") ? 0.0 : family_stagger_step;
    double effective_min_lead0 = min_lead0;
    double effective_min_lead1 = min_lead1;
    double effective_alignment_station = aligned_family_station;
    double entry_stagger = std::isfinite(aligned_entry_stagger)
        ? aligned_entry_stagger : aligned_point_stagger;
    double exit_stagger = std::isfinite(aligned_exit_stagger)
        ? aligned_exit_stagger : aligned_point_stagger;
    // 候选评估统一入口：加长量、把手比例、横向偏置全部作为入参，形态/物理/
    // 同簇门禁与评分只实现一次，供零偏置主搜索与横向偏置补充搜索共用。
    auto try_candidate = [&](double extra0, double extra1, double arc_alpha,
                             double entry_bias, double exit_bias,
                             bool require_strict_improvement) {
        const double extra_sum = extra0 + extra1;
        BezierCurve candidate = UTurnCurveInitializer().buildSegmented(
            p0, t0, p1, t1,
            effective_min_lead0, effective_min_lead1, arc_alpha,
            aligned_point_stagger,
            base_lead0_extra_after_align + extra0,
            base_lead1_extra_after_align + extra1,
            entry_stagger, exit_stagger, entry_bias, exit_bias,
            clamp_step);
        if (candidate.empty()) {
            ++rejected_shape;
            return;
        }
        if (candidate.numSegments() != 3) {
            ++rejected_shape;
            return;
        }
        if (!segmentedUTurnHasMinimumStraightLeads(
                candidate, effective_min_lead0, effective_min_lead1)) {
            ++rejected_min_lead;
            return;
        }
        // 门禁顺序按"代价递增"排列：轴向平齐与家族站位只是两次点积，却能否决
        // 绝大多数加长候选；自交检测与最大曲率需要密集采样，放在其后可以显著
        // 降低单路口生成耗时。判定结果与顺序无关，只有调试计数的归因会变化。
        Vec2d q0 = candidate.segs.front().ctrl[3];
        Vec2d q1 = candidate.segs.back().ctrl[0];
        if (std::abs((q0 - q1).dot(axis)) > 0.05) {
            ++rejected_axis;
            return;
        }
        // 家族平齐站位是同入/同出连通家族共同的精确站位。q0/q1
        // 不能作为某一条曲线单独沿轴继续前推，否则会在家族内形成
        // 台阶，破坏平齐并诱发中弧相交。横向偏置垂直于轴，不影响站位。
        if (std::isfinite(effective_alignment_station) &&
            (std::abs(q0.dot(axis) - effective_alignment_station) > 0.05 ||
             std::abs(q1.dot(axis) - effective_alignment_station) > 0.05)) {
            ++rejected_station;
            return;
        }
        if (curveSelfIntersectsBusiness(candidate, 1.0)) {
            ++rejected_self;
            return;
        }
        const double candidate_curvature = candidate.maxCurvature(40);
        if (candidate_curvature >= segmented_max_curvature) {
            ++rejected_curvature;
            return;
        }
        // 弧长与最大曲率都是逐段密集采样的重量级量，评分阶段直接复用门禁
        // 阶段已经算好的结果，不再重复计算(候选量为万级，重复计算会直接
        // 体现在单路口生成耗时上)。
        const double candidate_arc_length = candidate.arcLength();
        double arc_chord = candidate_arc_length / chord_len;
        if (arc_chord < 1.35 ||
            (chord_len >= 1.0 &&
             !segmentedUTurnMiddleArcLooksRound(candidate, axis))) {
            ++rejected_round;
            return;
        }
        if (!segmentedUTurnMiddleArcClearsCrosswalks(
                candidate, clearance_crosswalks)) {
            ++rejected_crosswalk;
            return;
        }

        SegmentedUTurnAudit candidate_risk = audit(candidate, include_fence);
        if (candidate_risk.physical_violation) {
            ++rejected_physical;
            return;
        }
        if (current_risk.sibling_crosses == 0 &&
            candidate_risk.sibling_crosses > 0) {
            // 同簇非端点不相交是硬约束：宁可保留当前形态，也不能为了换成
            // 三段式而新增同簇交叉。横向偏置搜索正是为了在保持这条硬约束的
            // 前提下把首尾直段推离相邻车道中心线，从而真正取得三段式形态。
            ++rejected_sibling_cross;
            return;
        }
        const double bias_cost = std::abs(entry_bias) + std::abs(exit_bias);
        double score = 1000.0 * candidate_risk.sibling_crosses +
                       candidate_curvature +
                       0.02 * candidate_arc_length +
                       0.05 * extra_sum +
                       0.50 * bias_cost;
        // 偏置候选只允许在严格减少同簇交叉时替换零偏置最优解。同交叉数下靠
        // 曲率/弧长微差胜出会无谓改写形态，扰动后续曲线的避让空间，实测会在
        // 其它数据集诱发新的同簇违约。
        const bool strictly_better =
            !have_best || candidate_risk.sibling_crosses < best_cross;
        const bool tie_better =
            have_best && candidate_risk.sibling_crosses == best_cross &&
            score < best_score;
        if (strictly_better ||
            (tie_better && !require_strict_improvement)) {
            best = candidate;
            best_cross = candidate_risk.sibling_crosses;
            best_score = score;
            best_extra_sum = extra_sum;
            have_best = true;
        }
    };
    // 第一阶段：零横向偏置的主搜索，保持既有形态优先级不变。
    for (const auto& lead_pair : lead_pairs) {
        double extra_sum = lead_pair.extra0 + lead_pair.extra1;
        if (best_cross == 0 && extra_sum > best_extra_sum + 1e-9)
            break;
        for (double arc_alpha : arc_alphas)
            try_candidate(lead_pair.extra0, lead_pair.extra1, arc_alpha,
                          0.0, 0.0, false);
    }
    // 第二阶段：横向偏置补充搜索。
    // 触发条件是主搜索没能取得任何合规候选，或最优候选仍带同簇交叉——这正是
    // "首尾直行段贴着共享车道中心线、与同簇直行曲线在毫米级缝隙内换侧"造成的
    // 被迫相交。首尾平齐点沿 U 轴法线偏置不改变轴向站位，家族平齐与人行横道
    // 净距(lead 按偏置反解补偿)都保持成立，只把直段整体移到直行曲线的同一侧。
    int bias_candidates = 0;
    // 只有存在"共享端点且本身不是几何掉头"的同簇兄弟时才启动偏置：这正是
    // 直行/普通转向与掉头汇入同一车道端点、直段在窄缝内换侧的几何。掉头之间
    // 的相互避让由深度/错开机制负责，不应由偏置介入。
    // 同时按兄弟贴的是哪一个端点，只在对应侧施加偏置，避免无谓的候选膨胀。
    bool bias_entry_side = false;
    bool bias_exit_side = false;
    auto sibling_qualifies_for_bias = [](const SampledSiblingCurve& sib) {
        return !sib.exempt_a1 && sib.shared_endpoint && !sib.curve.empty() &&
               !curveLooksUTurnForClusterExemption(sib.curve);
    };
    for (const auto& sib : sampled_siblings) {
        if (!sibling_qualifies_for_bias(sib))
            continue;
        const double d_entry = std::min((sib.curve.startPt() - p0).norm(),
                                       (sib.curve.endPt() - p0).norm());
        const double d_exit = std::min((sib.curve.startPt() - p1).norm(),
                                      (sib.curve.endPt() - p1).norm());
        if (d_entry < 0.5)
            bias_entry_side = true;
        if (d_exit < 0.5)
            bias_exit_side = true;
    }
    // 边界避让不能依赖兄弟曲线是否已经生成：同一批次中排在前面的 U-turn
    // 可能没有可用于推导偏置方向的共享兄弟，但其首尾直段仍可能贴近或穿过
    // RoadEdge。对存在 Boundary 的输入始终尝试有限的两侧偏置；候选最终仍
    // 必须通过严格物理审计，因此不会把“尝试偏置”变成任何接触豁免。
    const bool boundary_bias_needed = !input.boundaries.empty();
    if ((bias_entry_side || bias_exit_side || boundary_bias_needed) &&
        (!have_best || best_cross > 0)) {
        // G1 余量：首尾直段方向相对车道切向的偏差为 atan(bias/lead)，
        // shape.uturn.leads 要求方向点积 >= 0.98(约 11.5°)。常规档位以各侧
        // lead 下限的 5% 为上限(约 2.9°)，数据驱动档位的硬上限取 0.14*lead
        // (约 8°)，与 buildSegmented 里相向内缩的 g1_safe_inset 同源。
        const double cap0 = std::max(0.02, 0.05 * effective_min_lead0);
        const double cap1 = std::max(0.02, 0.05 * effective_min_lead1);
        // 折算所需偏置要用真实直段长度，而不是 lead 下限：家族平齐会把直段
        // 继续前推，用下限折算会低估。这里用零偏置探针候选量出实际首尾直段。
        BezierCurve probe = have_best ? best
            : UTurnCurveInitializer().buildSegmented(
                  p0, t0, p1, t1, effective_min_lead0, effective_min_lead1,
                  2.0 / 3.0, aligned_point_stagger,
                  base_lead0_extra_after_align, base_lead1_extra_after_align,
                  entry_stagger, exit_stagger, 0.0, 0.0,
                  clamp_step);
        double lead_len0 = std::max(0.5, effective_min_lead0);
        double lead_len1 = std::max(0.5, effective_min_lead1);
        if (probe.numSegments() == 3) {
            lead_len0 = std::max(
                0.5, (probe.segs.front().ctrl[3] - p0).norm());
            lead_len1 = std::max(
                0.5, (p1 - probe.segs.back().ctrl[0]).norm());
        }
        const Vec2d lateral{-axis.y(), axis.x()};
        // 数据驱动档位：固定比例档位对"弯曲直行"远远不够。110003449 的直行 12
        // 在掉头 1 的 12.05m 直段末端已横向漂移 0.93m(约 4.45°)，而 5% 档位只
        // 有 2.9°，直段全程仍留在直行曲线内侧，中弧唯一的出路就是横穿过去。
        // 这里按兄弟曲线的实测漂移反解所需偏置，两个符号都测，由物理与同簇
        // 审计决定哪一侧真正有效。
        double need0[2] = {0.0, 0.0};
        double need1[2] = {0.0, 0.0};
        for (const auto& sib : sampled_siblings) {
            if (!sibling_qualifies_for_bias(sib))
                continue;
            const double d_entry = std::min((sib.curve.startPt() - p0).norm(),
                                           (sib.curve.endPt() - p0).norm());
            const double d_exit = std::min((sib.curve.startPt() - p1).norm(),
                                          (sib.curve.endPt() - p1).norm());
            for (int s = 0; s < 2; ++s) {
                const double side = s == 0 ? 1.0 : -1.0;
                if (d_entry < 0.5) {
                    need0[s] = std::max(need0[s],
                        requiredLateralBiasToClearSibling(
                            sib.curve, p0, axis, lateral, side, lead_len0));
                }
                if (d_exit < 0.5) {
                    need1[s] = std::max(need1[s],
                        requiredLateralBiasToClearSibling(
                            sib.curve, p1, exit_back, lateral, side,
                            lead_len1));
                }
            }
        }
        // 留一点余量再钳到 G1 上限：偏置要严格超过兄弟漂移才算"整段在外侧"。
        auto finalize_need = [](double raw, double g1_cap) {
            if (raw <= 0.0)
                return 0.0;
            return std::min(g1_cap, raw * 1.15 + 0.03);
        };
        const double g1_cap0 = 0.14 * lead_len0;
        const double g1_cap1 = 0.14 * lead_len1;
        // 逐级放大偏置并在取得零交叉候选后立即停止：偏置阶段每个候选都要走
        // 形态门禁与物理审计，穷举整张网格会让掉头密集的路口生成时间显著变长。
        // 每一档是一组有序的 (entry_bias, exit_bias) 组合，先试形态代价小的。
        typedef std::vector<std::pair<double, double> > BiasTier;
        std::vector<BiasTier> tiers;
        for (double frac : {0.50, 1.00}) {
            BiasTier tier;
            if (bias_exit_side) {
                tier.push_back(std::make_pair(0.0, frac * cap1));
                tier.push_back(std::make_pair(0.0, -frac * cap1));
            }
            if (bias_entry_side) {
                tier.push_back(std::make_pair(frac * cap0, 0.0));
                tier.push_back(std::make_pair(-frac * cap0, 0.0));
            }
            // 双侧同时偏置只按同符号配对，不铺开混合符号：混合符号会把首尾直段
            // 扭向相反方向，等于用畸形形态换非交。
            if (bias_entry_side && bias_exit_side) {
                tier.push_back(std::make_pair(frac * cap0, frac * cap1));
                tier.push_back(std::make_pair(-frac * cap0, -frac * cap1));
            }
            tiers.push_back(tier);
        }
        {
            BiasTier tier;
            for (int s = 0; s < 2; ++s) {
                const double sign = s == 0 ? 1.0 : -1.0;
                const double m0 = finalize_need(need0[s], g1_cap0);
                const double m1 = finalize_need(need1[s], g1_cap1);
                // 双侧等量同向偏置是整条中弧的平移：弦长不变、两端切向只各转
                // atan(bias/lead)，中弧曲率基本不动，是形态代价最小的做法，
                // 优先尝试。单侧偏置会把该侧平齐点推向另一端，缩短中弧弦长并
                // 抬高曲率(实测 110003449 掉头 1 的 maxκ 由 0.59 升到 0.96)。
                // 为了做出平移，另一侧即使没有贴合兄弟也要跟着偏置，因此这里
                // 不再按侧门控，只受各侧 G1 上限约束。
                const double want = std::max(m0, m1);
                if (want > cap0 || want > cap1) {
                    const double pair0 = std::min(want, g1_cap0);
                    const double pair1 = std::min(want, g1_cap1);
                    tier.push_back(
                        std::make_pair(sign * pair0, sign * pair1));
                }
                // 平移被物理或同簇审计否决时，再退回单侧偏置。
                if (bias_entry_side && m0 > cap0)
                    tier.push_back(std::make_pair(sign * m0, 0.0));
                if (bias_exit_side && m1 > cap1)
                    tier.push_back(std::make_pair(0.0, sign * m1));
            }
            if (!tier.empty())
                tiers.push_back(tier);
        }
        if (boundary_bias_needed) {
            const double boundary_cap0 = g1_cap0;
            const double boundary_cap1 = g1_cap1;
            for (double fraction : {0.50, 1.00}) {
                BiasTier tier;
                // 双侧同向偏置优先保持中弧整体平移，四种符号组合都保留，
                // 让斜向 Boundary 的内外侧由严格审计实际筛选。
                tier.push_back(std::make_pair(
                    fraction * boundary_cap0, fraction * boundary_cap1));
                tier.push_back(std::make_pair(
                    -fraction * boundary_cap0, -fraction * boundary_cap1));
                tier.push_back(std::make_pair(
                    fraction * boundary_cap0, -fraction * boundary_cap1));
                tier.push_back(std::make_pair(
                    -fraction * boundary_cap0, fraction * boundary_cap1));
                tier.push_back(std::make_pair(fraction * boundary_cap0, 0.0));
                tier.push_back(std::make_pair(-fraction * boundary_cap0, 0.0));
                tier.push_back(std::make_pair(0.0, fraction * boundary_cap1));
                tier.push_back(std::make_pair(0.0, -fraction * boundary_cap1));
                tiers.push_back(tier);
            }
        }
        // 偏置阶段复用主搜索的完整把手序列：实测把它裁成 6 个代表值会让
        // 100000385-u / 100000443 / 110003285 找不到零交叉候选，整体违约数从
        // 98 回升到 100。候选量由"档位递增 + 零交叉早停"控制，而不是靠裁剪
        // 把手序列——零交叉一旦命中就立即跳出，命中不了的路口才会走满网格。
        const std::vector<double>& bias_arc_alphas = arc_alphas;
        for (const BiasTier& tier : tiers) {
            for (const auto& combo : tier) {
                // 逐档记录否决归因：偏置阶段的候选很少，只有把"哪一道门禁挡住了
                // 平移候选"打出来，才能区分"偏置幅度不够"与"平移撞上物理约束"。
                const int before_self = rejected_self;
                const int before_curv = rejected_curvature;
                const int before_axis = rejected_axis;
                const int before_station = rejected_station;
                const int before_round = rejected_round;
                const int before_cw = rejected_crosswalk;
                const int before_phys = rejected_physical;
                const int before_cross = rejected_sibling_cross;
                const int before_lead = rejected_min_lead;
                const int before_shape = rejected_shape;
                for (double arc_alpha : bias_arc_alphas) {
                    ++bias_candidates;
                    try_candidate(0.0, 0.0, arc_alpha, combo.first,
                                  combo.second, true);
                    if (have_best && best_cross == 0)
                        break;
                }
                if (debug)
                    fprintf(stderr,
                            "[UTURN-SEGMENTED] bias combo (%.3f,%.3f): "
                            "shape=%d lead=%d self=%d curvature=%d axis=%d "
                            "station=%d round=%d crosswalk=%d physical=%d "
                            "cross=%d best_cross=%d\n",
                            combo.first, combo.second,
                            rejected_shape - before_shape,
                            rejected_min_lead - before_lead,
                            rejected_self - before_self,
                            rejected_curvature - before_curv,
                            rejected_axis - before_axis,
                            rejected_station - before_station,
                            rejected_round - before_round,
                            rejected_crosswalk - before_cw,
                            rejected_physical - before_phys,
                            rejected_sibling_cross - before_cross,
                            have_best ? best_cross : -1);
                if (have_best && best_cross == 0)
                    break;
            }
            if (have_best && best_cross == 0)
                break;
        }
    }
    bool has_shared_uturn_sibling = false;
    for (const auto& sib : sampled_siblings) {
        if (!sib.exempt_a1 && sib.shared_endpoint && !sib.curve.empty() &&
            curveLooksUTurnForClusterExemption(sib.curve)) {
            has_shared_uturn_sibling = true;
            break;
        }
    }
    if (!have_best && !has_shared_uturn_sibling &&
        curve.numSegments() == 1 && current_risk.sibling_crosses == 0 &&
        !current_risk.physical_violation) {
        const BezierSegment& source = curve.segs.front();
        std::vector<double> t0_options =
            {0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.10,
             0.12, 0.15, 0.18, 0.20, 0.25, 0.30};
        std::vector<double> t1_options =
            {0.98, 0.97, 0.96, 0.94, 0.92, 0.90, 0.88,
             0.85, 0.82, 0.80, 0.75, 0.70};
        // The fallback must still be able to realize the requested lead
        // floors when the source curve is short; add fractions derived from
        // the actual source length instead of relying on the historical
        // 0.20 m split range.
        const double source_length = source.arcLength(64);
        if (source_length > 1e-6) {
            t0_options.push_back(std::min(0.45, min_lead0 / source_length));
            t1_options.push_back(std::max(0.55, 1.0 - min_lead1 / source_length));
        }
        for (double t0_split : t0_options) {
            for (double t1_split : t1_options) {
                if (t1_split <= t0_split + 0.20 || t0_split <= 0.0 ||
                    t1_split >= 1.0)
                    continue;
                auto first_rest = source.splitAt(t0_split);
                double local_t1 =
                    (t1_split - t0_split) / std::max(1e-6, 1.0 - t0_split);
                auto arc_last = first_rest.second.splitAt(local_t1);
                BezierCurve candidate;
                candidate.segs.push_back(first_rest.first);
                candidate.segs.push_back(arc_last.first);
                candidate.segs.push_back(arc_last.second);
                if ((candidate.segs.front().ctrl[3] -
                     candidate.segs.front().ctrl[0]).norm() + 1e-6 < min_lead0 ||
                    (candidate.segs.back().ctrl[3] -
                     candidate.segs.back().ctrl[0]).norm() + 1e-6 < min_lead1)
                    continue;
                if (!segmentLooksStraight(candidate.segs.front()) ||
                    !segmentLooksStraight(candidate.segs.back()))
                    continue;
                if (candidate.segs[1].maxCurvature(30) < 0.03)
                    continue;
                if (curveSelfIntersectsBusiness(candidate, 1.0))
                    continue;
                if (candidate.maxCurvature(40) >= segmented_max_curvature)
                    continue;
                Vec2d q0 = candidate.segs.front().ctrl[3];
                Vec2d q1 = candidate.segs.back().ctrl[0];
                if (std::abs((q0 - q1).dot(axis)) > 0.05)
                    continue;
                if (std::isfinite(effective_alignment_station) &&
                    (std::abs(q0.dot(axis) - effective_alignment_station) > 0.05 ||
                     std::abs(q1.dot(axis) - effective_alignment_station) > 0.05))
                    continue;
                if (candidate.arcLength() / chord_len < 1.35 ||
                    (chord_len >= 1.0 &&
                     !segmentedUTurnMiddleArcLooksRound(candidate, axis)))
                    continue;
                if (!segmentedUTurnMiddleArcClearsCrosswalks(
                        candidate, clearance_crosswalks))
                    continue;
                SegmentedUTurnAudit candidate_risk = audit(candidate, include_fence);
                if (candidate_risk.physical_violation ||
                    (current_risk.sibling_crosses == 0 &&
                     candidate_risk.sibling_crosses > 0))
                    continue;
                curve = candidate;
                return true;
            }
        }
    }
    if (!have_best) {
        if (debug)
            fprintf(stderr,
                    "[UTURN-SEGMENTED] no valid candidate, current_cross=%d "
                    "min_lead=%.3f/%.3f station=%.3f kmax=%.3f | shape=%d "
                    "lead=%d self=%d curvature=%d axis=%d station_gate=%d "
                    "round=%d crosswalk=%d physical=%d cross=%d bias_tried=%d\n",
                    current_risk.sibling_crosses, effective_min_lead0,
                    effective_min_lead1, effective_alignment_station,
                    segmented_max_curvature, rejected_shape, rejected_min_lead,
                    rejected_self, rejected_curvature, rejected_axis,
                    rejected_station, rejected_round, rejected_crosswalk,
                    rejected_physical, rejected_sibling_cross,
                    bias_candidates);
        return false;
    }
    if (debug && bias_candidates > 0)
        fprintf(stderr,
                "[UTURN-SEGMENTED] accepted after bias stage: cross=%d "
                "bias_tried=%d cross_rejected=%d\n",
                best_cross, bias_candidates, rejected_sibling_cross);
    curve = best;
    return true;
}

}  // 命名空间 isg
