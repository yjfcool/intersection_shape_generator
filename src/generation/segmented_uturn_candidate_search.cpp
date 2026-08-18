#include "generation/segmented_uturn_candidate_search.h"

#include "curve/curve_utils.h"
#include "generation/uturn_shape.h"
#include "initialization/uturn_curve_initializer.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>

namespace isg {

static bool debugSegmentedUTurnSearch() {
    static const bool enabled = std::getenv("ISG_DEBUG_UTURN") != nullptr;
    return enabled;
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
    double aligned_exit_stagger) const {
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
    const std::vector<Crosswalk>& clearance_crosswalks =
        crosswalks_for_clearance ? *crosswalks_for_clearance : input.crosswalks;
    double entry_stagger = std::isfinite(aligned_entry_stagger)
        ? aligned_entry_stagger : aligned_point_stagger;
    double exit_stagger = std::isfinite(aligned_exit_stagger)
        ? aligned_exit_stagger : aligned_point_stagger;
    double effective_min_lead0 = min_lead0;
    double effective_min_lead1 = min_lead1;
    double effective_alignment_station = aligned_family_station;
    for (const auto& lead_pair : lead_pairs) {
        double extra_sum = lead_pair.extra0 + lead_pair.extra1;
        if (best_cross == 0 && extra_sum > best_extra_sum + 1e-9)
            break;
        for (double arc_alpha : arc_alphas) {
            BezierCurve candidate = UTurnCurveInitializer().buildSegmented(
                p0, t0, p1, t1,
                effective_min_lead0, effective_min_lead1, arc_alpha, aligned_point_stagger,
                base_lead0_extra_after_align + lead_pair.extra0,
                base_lead1_extra_after_align + lead_pair.extra1,
                entry_stagger, exit_stagger);
            if (candidate.empty()) {
                continue;
            }
            if (candidate.numSegments() != 3) {
                continue;
            }
            if (curveSelfIntersectsBusiness(candidate, 1.0)) {
                continue;
            }
            if (candidate.maxCurvature(40) >= segmented_max_curvature) {
                continue;
            }
            Vec2d q0 = candidate.segs.front().ctrl[3];
            Vec2d q1 = candidate.segs.back().ctrl[0];
            if (std::abs((q0 - q1).dot(axis)) > 0.05) {
                continue;
            }
            // 家族平齐站位是同入/同出连通家族共同的精确站位。q0/q1
            // 不能作为某一条曲线单独沿轴继续前推，否则会在家族内形成
            // 台阶，破坏平齐并诱发中弧相交。需要避让时应进入专用路径。
            if (std::isfinite(effective_alignment_station) &&
                (std::abs(q0.dot(axis) - effective_alignment_station) > 0.05 ||
                 std::abs(q1.dot(axis) - effective_alignment_station) > 0.05)) {
                continue;
            }
            double arc_chord = candidate.arcLength() / chord_len;
            if (arc_chord < 1.35 ||
                (chord_len >= 1.0 &&
                 !segmentedUTurnMiddleArcLooksRound(candidate, axis))) {
                continue;
            }
            if (!segmentedUTurnMiddleArcClearsCrosswalks(
                    candidate, clearance_crosswalks)) {
                continue;
            }

            SegmentedUTurnAudit candidate_risk = audit(candidate, include_fence);
            if (candidate_risk.physical_violation) {
                continue;
            }
            if (current_risk.sibling_crosses == 0 &&
                candidate_risk.sibling_crosses > 0) {
                continue;
            }
            double score = 1000.0 * candidate_risk.sibling_crosses +
                           candidate.maxCurvature(40) +
                           0.02 * candidate.arcLength() +
                           0.05 * extra_sum;
            if (!have_best ||
                candidate_risk.sibling_crosses < best_cross ||
                (candidate_risk.sibling_crosses == best_cross &&
                 score < best_score)) {
                best = candidate;
                best_cross = candidate_risk.sibling_crosses;
                best_score = score;
                best_extra_sum = extra_sum;
                have_best = true;
            }
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
        const std::vector<double> t0_options =
            {0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.10,
             0.12, 0.15, 0.18, 0.20, 0.25, 0.30};
        const std::vector<double> t1_options =
            {0.98, 0.97, 0.96, 0.94, 0.92, 0.90, 0.88,
             0.85, 0.82, 0.80, 0.75, 0.70};
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
                     candidate.segs.front().ctrl[0]).norm() < 0.20 ||
                    (candidate.segs.back().ctrl[3] -
                     candidate.segs.back().ctrl[0]).norm() < 0.20)
                    continue;
                if (!segmentLooksStraight(candidate.segs.front()) ||
                    !segmentLooksStraight(candidate.segs.back()))
                    continue;
                if (candidate.segs[1].maxCurvature(30) < 0.03)
                    continue;
                if (curveSelfIntersectsBusiness(candidate, 1.0))
                    continue;
                if (!segmentedUTurnMiddleArcClearsCrosswalks(
                        candidate, clearance_crosswalks))
                    continue;
                curve = candidate;
                return true;
            }
        }
    }
    if (!have_best) {
        if (debug)
            fprintf(stderr, "[UTURN-SEGMENTED] no valid candidate, current_cross=%d\n",
                    current_risk.sibling_crosses);
        return false;
    }
    curve = best;
    return true;
}

}  // 命名空间 isg
