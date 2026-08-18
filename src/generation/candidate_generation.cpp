#include "generation/candidate_generation.h"

#include "constraints/uturn_envelope_constraint.h"
#include "curve/curve_utils.h"
#include "initialization/uturn_curve_initializer.h"

#include <cmath>
#include <limits>
#include <tuple>
#include <vector>

namespace isg {

const std::size_t CandidateSelector::npos = static_cast<std::size_t>(-1);

// 候选按硬约束、物理风险、簇交叉、形态、fixed shape、软代价和质量
// 依次作稳定词典序比较；字段和枚举顺序属于可观察生成行为。
std::size_t CandidateSelector::selectBest(
    const std::vector<CurveCandidate>& candidates) const {
    std::size_t best = npos;
    for (std::size_t i = 0; i < candidates.size(); ++i) {
        const CurveCandidate& candidate = candidates[i];
        if (!candidate.hasCurrentReport() || candidate.report.hasHardViolation())
            continue;
        if (best == npos ||
            std::tie(candidate.hard_violation_count,
                     candidate.physical_violation,
                     candidate.new_cluster_crosses,
                     candidate.shape_hard_violations,
                     candidate.fixed_shape_loss,
                     candidate.soft_penalty,
                     candidate.quality_score) <
            std::tie(candidates[best].hard_violation_count,
                     candidates[best].physical_violation,
                     candidates[best].new_cluster_crosses,
                     candidates[best].shape_hard_violations,
                     candidates[best].fixed_shape_loss,
                     candidates[best].soft_penalty,
                     candidates[best].quality_score)) {
            best = i;
        }
    }
    return best;
}

bool BoundedUTurnCandidateSearch::search(
    const Vec2d& entry, const Vec2d& entry_tangent,
    const Vec2d& exit, const Vec2d& exit_tangent,
    const BezierCurve& reference,
    const BoundedUTurnCandidateAuditor& auditor,
    BezierCurve& result, double min_lead0, double min_lead1) const {
    const Vec2d start_tangent = entry_tangent.norm() > 1e-8
        ? entry_tangent.normalized() : Vec2d(1, 0);
    const Vec2d end_tangent = exit_tangent.norm() > 1e-8
        ? exit_tangent.normalized() : -start_tangent;
    Vec2d axis = start_tangent - end_tangent;
    if (axis.norm() < 1e-8)
        axis = start_tangent;
    axis.normalize();

    std::vector<Vec2d> directions;
    const auto add_direction = [&](const Vec2d& raw) {
        if (raw.norm() < 1e-8)
            return;
        const Vec2d direction = raw.normalized();
        for (const auto& existing : directions) {
            if (existing.dot(direction) > 0.98)
                return;
        }
        directions.push_back(direction);
    };
    add_direction(axis);
    add_direction(-axis);
    add_direction(Vec2d(-start_tangent.y(), start_tangent.x()));
    add_direction(Vec2d(start_tangent.y(), -start_tangent.x()));

    bool found = false;
    BezierCurve best;
    int best_crosses = std::numeric_limits<int>::max();
    double best_score = std::numeric_limits<double>::max();
    const double reference_length = reference.arcLength();
    const Vec2d lateral(-axis.y(), axis.x());
    const double turn_gap = std::abs((exit - entry).dot(lateral));
    const double max_curvature = turn_gap < 1.0
        ? std::max(8.0, 4.0 / std::max(0.1, turn_gap)) : 3.0;
    const double chord_length = (exit - entry).norm();
    const UTurnEnvelopeConstraint envelope;
    const UTurnCurveInitializer initializer;

    // 固定网格及遍历顺序保持既有提前退出和评分结果，避免无界全局搜索。
    const double scales[] = {0.85, 0.70, 0.55, 0.40, 1.0, 0.30};
    const double offsets[] = {0.0, 1.0, 2.0, 3.0, 4.5, 6.0};
    for (std::size_t scale_index = 0; scale_index < 6; ++scale_index) {
        const double scale = scales[scale_index];
        for (const auto& direction : directions) {
            for (std::size_t offset_index = 0; offset_index < 6;
                 ++offset_index) {
                const double offset = offsets[offset_index];
                const BezierCurve candidate = initializer.buildAligned(
                    entry, entry_tangent, exit, exit_tangent,
                    direction, offset, scale, min_lead0, min_lead1);
                if (candidate.empty() ||
                    curveSelfIntersectsBusiness(candidate, 1.0) ||
                    envelope.exceeds(candidate, reference, entry, exit) ||
                    envelope.collapses(candidate, reference, entry, exit)) {
                    continue;
                }
                const double candidate_max_curvature =
                    candidate.maxCurvature(20);
                if (candidate_max_curvature > max_curvature)
                    continue;
                const double arc_chord = chord_length > 1e-6
                    ? candidate.arcLength() / chord_length : 1.0;
                if (turn_gap >= 10.0 &&
                    (arc_chord <= 1.40 || candidate_max_curvature >= 0.50)) {
                    continue;
                }

                const BoundedUTurnCandidateAudit audit = auditor(candidate);
                if (audit.physical_violation)
                    continue;
                const double score =
                    std::abs(candidate.arcLength() - reference_length) +
                    3.0 * std::abs(scale - 0.85) +
                    0.2 * offset + 0.3 * candidate_max_curvature;
                if (!found || audit.sibling_crosses < best_crosses ||
                    (audit.sibling_crosses == best_crosses &&
                     score < best_score)) {
                    found = true;
                    best = candidate;
                    best_crosses = audit.sibling_crosses;
                    best_score = score;
                }
            }
        }
    }

    if (!found)
        return false;
    result = best;
    return true;
}

}  // 命名空间 isg
