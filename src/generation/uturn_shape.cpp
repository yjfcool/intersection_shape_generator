#include "generation/uturn_shape.h"

#include "constraints/fence_check.h"
#include "utils.h"

#include <algorithm>
#include <cmath>

namespace isg {

bool segmentLooksStraight(const BezierSegment& segment) {
    const Vec2d chord = segment.ctrl[3] - segment.ctrl[0];
    const double length = chord.norm();
    if (length < 1e-6)
        return false;
    return pointToSegment(segment.ctrl[1], segment.ctrl[0], segment.ctrl[3]).first <=
               std::max(0.05, length * 0.03) &&
           pointToSegment(segment.ctrl[2], segment.ctrl[0], segment.ctrl[3]).first <=
               std::max(0.05, length * 0.03);
}

bool curveLooksUTurnForClusterExemption(const BezierCurve& curve) {
    if (curve.empty())
        return false;
    const Vec2d start = curve.startTan();
    const Vec2d end = curve.endTan();
    return start.norm() > 1e-8 && end.norm() > 1e-8 &&
           start.normalized().dot(end.normalized()) < -0.5;
}

bool segmentedUTurnMiddleArcLooksRound(
    const BezierCurve& curve, const Vec2d& axis) {
    if (curve.numSegments() != 3 || axis.norm() < 1e-8)
        return false;
    const BezierSegment& arc = curve.segs[1];
    const double gap = (arc.ctrl[3] - arc.ctrl[0]).norm();
    if (gap < 1e-3 || arc.arcLength(32) / gap < 1.28)
        return false;
    const Vec2d direction = axis.normalized();
    double maximum_bulge = 0.0;
    for (int i = 0; i <= 32; ++i)
        maximum_bulge = std::max(maximum_bulge, std::abs(
            (arc.evaluate(static_cast<double>(i) / 32.0) - arc.ctrl[0]).dot(direction)));
    return maximum_bulge >= std::max(0.25, 0.30 * gap);
}

double segmentedUTurnMaxCurvatureLimit(
    const Vec2d& entry, const Vec2d& entry_tangent,
    const Vec2d& exit, const Vec2d& exit_tangent) {
    const double chord_length = (exit - entry).norm();
    const Vec2d start = entry_tangent.norm() > 1e-8
        ? entry_tangent.normalized() : Vec2d(1, 0);
    const Vec2d end = exit_tangent.norm() > 1e-8
        ? exit_tangent.normalized() : -start;
    Vec2d axis = start - end;
    if (axis.norm() < 1e-8)
        axis = start;
    axis.normalize();
    if (axis.dot(start) < 0.0)
        axis = -axis;
    const Vec2d lateral{-axis.y(), axis.x()};
    const double turn_gap = std::abs((exit - entry).dot(lateral));
    if (turn_gap < 1.0)
        return std::max(8.0, 4.0 / std::max(0.1, turn_gap));
    return chord_length < 4.0 ? 8.0 : (chord_length < 7.0 ? 6.0 : 3.0);
}

bool segmentedUTurnMiddleArcClearsCrosswalks(
    const BezierCurve& curve, const std::vector<Crosswalk>& crosswalks) {
    if (curve.numSegments() != 3)
        return false;
    for (int i = 0; i <= 64; ++i) {
        const Vec2d point = curve.segs[1].evaluate(static_cast<double>(i) / 64.0);
        for (const auto& crosswalk : crosswalks)
            if (polygonContains(crosswalk.geometry, point))
                return false;
    }
    return true;
}

}  // 命名空间 isg
