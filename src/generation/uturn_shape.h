#pragma once

#include "types.h"

namespace isg {

bool segmentLooksStraight(const BezierSegment& segment);
bool curveLooksUTurnForClusterExemption(const BezierCurve& curve);
bool segmentedUTurnMiddleArcLooksRound(
    const BezierCurve& curve, const Vec2d& axis);
bool segmentedUTurnHasMinimumStraightLeads(
    const BezierCurve& curve, double min_lead0, double min_lead1,
    double tolerance = 1e-6);
double segmentedUTurnMaxCurvatureLimit(
    const Vec2d& entry, const Vec2d& entry_tangent,
    const Vec2d& exit, const Vec2d& exit_tangent);
bool segmentedUTurnMiddleArcClearsCrosswalks(
    const BezierCurve& curve, const std::vector<Crosswalk>& crosswalks);

}  // 命名空间 isg
