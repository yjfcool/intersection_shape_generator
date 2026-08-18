#include "constraints/shape_constraint.h"

#include "constraints/fence_check.h"
#include "geometry/predicates.h"
#include "utils.h"

#include <algorithm>
#include <cmath>

namespace isg {
namespace {

ConstraintResult satisfied(const char* id) {
    ConstraintResult result;
    result.id = id;
    result.severity = ConstraintSeverity::Hard;
    result.state = ConstraintState::Satisfied;
    return result;
}

ConstraintResult violated(const char* id, const char* reason, double amount = 1.0) {
    ConstraintResult result;
    result.id = id;
    result.severity = ConstraintSeverity::Hard;
    result.state = ConstraintState::Violated;
    result.violation = amount;
    result.reason = reason;
    return result;
}

bool isGeometricUTurn(const CurveGenerationContext& context) {
    return context.entry.second.norm() > 1e-8 && context.exit.second.norm() > 1e-8 &&
           context.entry.second.normalized().dot(context.exit.second.normalized()) < -0.5;
}

bool signFlip(const BezierCurve& curve, double eps = 0.10) {
    int sign = 0;
    for (const auto& segment : curve.segs) {
        for (int i = 1; i < 20; ++i) {
            const double t = static_cast<double>(i) / 20.0;
            const Vec2d d1 = segment.evalDeriv1(t);
            const Vec2d d2 = segment.evalDeriv2(t);
            const double denominator = std::pow(d1.squaredNorm(), 1.5);
            if (denominator < 1e-12) continue;
            const double k = cross2d(d1, d2) / denominator;
            if (std::abs(k) <= eps) continue;
            const int current = k > 0.0 ? 1 : -1;
            if (sign != 0 && current != sign) return true;
            sign = current;
        }
    }
    return false;
}

bool segmentStraight(const BezierSegment& segment) {
    const Vec2d chord = segment.ctrl[3] - segment.ctrl[0];
    const double length = chord.norm();
    if (length < 1e-6) return false;
    return pointToSegment(segment.ctrl[1], segment.ctrl[0], segment.ctrl[3]).first <=
               std::max(0.05, length * 0.03) &&
           pointToSegment(segment.ctrl[2], segment.ctrl[0], segment.ctrl[3]).first <=
               std::max(0.05, length * 0.03);
}

bool pointInCrosswalks(const Vec2d& point, const std::vector<Crosswalk>& crosswalks) {
    for (const auto& crosswalk : crosswalks)
        if (polygonContains(crosswalk.geometry, point)) return true;
    return false;
}

double requiredCrosswalkLead(const Vec2d& origin, const Vec2d& direction,
                             const std::vector<Crosswalk>& crosswalks) {
    if (direction.norm() < 1e-8) return 0.0;
    const Vec2d forward = direction.normalized();
    const Vec2d lateral(-forward.y(), forward.x());
    double best_near = 1e18;
    double best_far = 0.0;
    for (const auto& crosswalk : crosswalks) {
        double near = 1e18;
        double far = 0.0;
        bool found = false;
        const std::vector<Vec2d> points = toVec2dArray(crosswalk.geometry.outer);
        for (size_t i = 0; i + 1 < points.size(); ++i) {
            const double a_forward = (points[i] - origin).dot(forward);
            const double b_forward = (points[i + 1] - origin).dot(forward);
            const double edge_near = std::min(a_forward, b_forward);
            const double edge_far = std::max(a_forward, b_forward);
            if (edge_far <= 0.0 || edge_near > 12.0) continue;
            const double a_lateral = (points[i] - origin).dot(lateral);
            const double b_lateral = (points[i + 1] - origin).dot(lateral);
            const double minimum_lateral = a_lateral * b_lateral <= 0.0 ? 0.0 :
                std::min(std::abs(a_lateral), std::abs(b_lateral));
            if (minimum_lateral > 4.0) continue;
            near = std::min(near, std::max(0.0, edge_near));
            far = std::max(far, edge_far);
            found = true;
        }
        if (found && (near < best_near - 1e-6 ||
                      (std::abs(near - best_near) <= 1e-6 && far < best_far))) {
            best_near = near;
            best_far = far;
        }
    }
    return best_near == 1e18 ? 0.0 : best_far + 0.30;
}

}  // namespace

ConstraintResult evaluateOrdinaryShape(const BezierCurve& curve,
                                       const CurveGenerationContext& context) {
    if (!context.connectivity || isGeometricUTurn(context))
        return satisfied("shape.ordinary.not_applicable");
    if (curve.empty()) return violated("shape.ordinary", "ordinary curve is empty");
    if (curve.numSegments() != 1)
        return violated("shape.ordinary.single_segment", "ordinary curve must be one cubic segment");
    const Vec2d chord = context.exit.first - context.entry.first;
    if (chord.norm() < 1e-6) return violated("shape.ordinary.chord", "ordinary curve has degenerate chord");
    const Vec2d t0 = context.entry.second.norm() > 1e-8 ? context.entry.second.normalized() : chord.normalized();
    const double turn_strength = std::abs(cross2d(t0, chord.normalized()));
    if (!ordinarySingleCubicControlsValid(curve, context.entry.first, t0,
                                          context.exit.first, context.exit.second))
        return violated("shape.ordinary.controls", "ordinary cubic controls leave endpoint tangent axes");
    const Vec2d start_tangent = curve.startTan();
    const Vec2d end_tangent = curve.endTan();
    if (start_tangent.norm() < 1e-8 || end_tangent.norm() < 1e-8 ||
        context.entry.second.norm() < 1e-8 || context.exit.second.norm() < 1e-8 ||
        start_tangent.normalized().dot(context.entry.second.normalized()) < 0.95 ||
        end_tangent.normalized().dot(context.exit.second.normalized()) < 0.95)
        return violated("shape.ordinary.g1", "ordinary curve endpoint tangents do not match lane frames");
    const double ratio = curve.arcLength() / chord.norm();
    const double lateral = [&]() {
        double maximum = 0.0;
        for (const auto& point : curve.sampleByArcLength(64))
            maximum = std::max(maximum, std::abs(cross2d(chord.normalized(), point - context.entry.first)));
        return maximum / chord.norm();
    }();
    if (turn_strength <= 0.25) {
        if (ratio > 1.08 || lateral > 0.08 || curve.maxCurvature(40) > 3.0)
            return violated("shape.ordinary.straight", "straight curve is not sufficiently straight");
    } else {
        const double minimum_ratio = chord.norm() < 12.0 ? 1.02 : 1.06;
        if (ratio < minimum_ratio || ratio > 1.35 || curve.maxCurvature(40) > 2.5 || signFlip(curve))
            return violated("shape.ordinary.turn", "turn curve is not a single smooth arch");
    }
    return satisfied("shape.ordinary");
}

ConstraintResult evaluateUTurnShape(const BezierCurve& curve,
                                    const CurveGenerationContext& context) {
    if (!context.connectivity || !isGeometricUTurn(context))
        return satisfied("shape.uturn.not_applicable");
    if (curve.numSegments() != 3)
        return violated("shape.uturn.segments", "U-turn must use straight-arc-straight");
    const Vec2d first_chord = curve.segs.front().ctrl[3] - curve.segs.front().ctrl[0];
    const Vec2d last_chord = curve.segs.back().ctrl[3] - curve.segs.back().ctrl[0];
    if (!segmentStraight(curve.segs.front()) || !segmentStraight(curve.segs.back()) ||
        first_chord.norm() < 0.05 || last_chord.norm() < 0.05 ||
        first_chord.normalized().dot(context.entry.second.normalized()) < 0.98 ||
        last_chord.normalized().dot(context.exit.second.normalized()) < 0.98)
        return violated("shape.uturn.leads", "U-turn lead segments are not straight and lane-aligned");
    const std::vector<Crosswalk>* crosswalks = context.scene ?
        &context.scene->view.input().crosswalks : nullptr;
    const double required0 = crosswalks ?
        requiredCrosswalkLead(context.entry.first, context.entry.second, *crosswalks) : 0.0;
    const double required1 = crosswalks ?
        requiredCrosswalkLead(context.exit.first, -context.exit.second, *crosswalks) : 0.0;
    if (required0 <= 0.0 && required1 <= 0.0 &&
        (first_chord.norm() < 2.0 || last_chord.norm() < 2.0))
        return violated("shape.uturn.leads", "U-turn without Crosswalk needs 2m straight leads");
    Vec2d axis = context.entry.second.normalized() - context.exit.second.normalized();
    if (axis.norm() < 1e-8) axis = context.entry.second.normalized();
    axis.normalize();
    if (axis.dot(context.entry.second) < 0.0) axis = -axis;
    const Vec2d q0 = curve.segs.front().ctrl[3];
    const Vec2d q1 = curve.segs.back().ctrl[0];
    if (std::abs((q0 - q1).dot(axis)) > 0.10)
        return violated("shape.uturn.alignment", "U-turn arc endpoints are not axially aligned");
    if (curve.segs[1].maxCurvature(30) <= 0.03)
        return violated("shape.uturn.arc", "U-turn middle segment is not a turn arc");
    return satisfied("shape.uturn");
}

ConstraintResult evaluateUTurnCrosswalk(const BezierCurve& curve,
                                        const CurveGenerationContext& context) {
    if (!context.connectivity || !isGeometricUTurn(context) || context.scene == nullptr ||
        context.scene->view.input().crosswalks.empty())
        return satisfied("shape.crosswalk.not_applicable");
    if (curve.numSegments() != 3)
        return violated("shape.crosswalk.segments", "Crosswalk rule requires segmented U-turn");
    const BezierCurve& arc = curve;
    for (int i = 0; i <= 48; ++i) {
        const Vec2d point = arc.segs[1].evaluate(static_cast<double>(i) / 48.0);
        if (pointInCrosswalks(point, context.scene->view.input().crosswalks))
            return violated("shape.crosswalk.arc", "U-turn middle arc enters Crosswalk");
    }
    const double required0 = requiredCrosswalkLead(
        context.entry.first, context.entry.second, context.scene->view.input().crosswalks);
    const double required1 = requiredCrosswalkLead(
        context.exit.first, -context.exit.second, context.scene->view.input().crosswalks);
    if ((curve.segs.front().ctrl[3] - curve.segs.front().ctrl[0]).norm() + 1e-6 < required0 ||
        (curve.segs.back().ctrl[3] - curve.segs.back().ctrl[0]).norm() + 1e-6 < required1)
        return violated("shape.crosswalk.lead", "U-turn straight lead does not clear Crosswalk");
    return satisfied("shape.crosswalk");
}

}  // 命名空间 isg
