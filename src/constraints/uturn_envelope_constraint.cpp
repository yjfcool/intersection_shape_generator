#include "constraints/uturn_envelope_constraint.h"

#include "utils.h"

#include <algorithm>
#include <cmath>

namespace isg {

bool UTurnEnvelopeConstraint::exceeds(
    const BezierCurve& curve, const BezierCurve& reference,
    const Vec2d& entry, const Vec2d& exit) const {
    const double reference_length = std::max(
        reference.arcLength(), (exit - entry).norm());
    if (curve.arcLength() > std::max(80.0, reference_length * 4.0))
        return true;

    const Vec2d midpoint = 0.5 * (entry + exit);
    const double max_radius = std::max(60.0, reference_length * 3.0);
    for (const auto& point : curve.sampleByArcLength(80))
        if ((point - midpoint).norm() > max_radius)
            return true;
    return false;
}

double UTurnEnvelopeConstraint::lateralBulge(
    const BezierCurve& curve, const Vec2d& entry,
    const Vec2d& exit) const {
    const Vec2d chord = exit - entry;
    if (chord.norm() < 1e-8)
        return 0.0;
    Vec2d perpendicular(-chord.y(), chord.x());
    perpendicular.normalize();
    double best = 0.0;
    for (const auto& point : curve.sampleByArcLength(80)) {
        const double projection = (point - entry).dot(perpendicular);
        if (std::abs(projection) > std::abs(best))
            best = projection;
    }
    return best;
}

bool UTurnEnvelopeConstraint::collapses(
    const BezierCurve& curve, const BezierCurve& reference,
    const Vec2d& entry, const Vec2d& exit) const {
    const double reference_length = reference.arcLength();
    const double chord_length = (exit - entry).norm();
    const double reference_bulge = std::abs(
        lateralBulge(reference, entry, exit));
    const double candidate_bulge = std::abs(
        lateralBulge(curve, entry, exit));

    if (reference_length > 1.0 &&
        curve.arcLength() < std::max(
            chord_length * 1.08, reference_length * 0.72))
        return true;
    if (reference_bulge > 2.0 &&
        candidate_bulge < reference_bulge * 0.50)
        return true;
    return curve.maxCurvature(20) > std::max(
        1.0, reference.maxCurvature(20) * 8.0);
}

}  // 命名空间 isg
