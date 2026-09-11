#include "preprocessing/crosswalk_clearance_calculator.h"

#include "utils.h"

#include <algorithm>
#include <cmath>

namespace isg {
namespace {

void considerPolyline(
        const std::vector<Vec2d>& points, const Vec2d& origin,
        const Vec2d& forward, const Vec2d& lateral, CrosswalkClearanceResult& candidate) {
    for (size_t i = 1; i < points.size(); ++i) {
        const Vec2d& a = points[i - 1];
        const Vec2d& b = points[i];
        const double a_fwd = (a - origin).dot(forward);
        const double b_fwd = (b - origin).dot(forward);
        const double near = std::min(a_fwd, b_fwd);
        const double far = std::max(a_fwd, b_fwd);
        if (far <= 0.0 || near > 12.0)
            continue;
        const double a_lat = (a - origin).dot(lateral);
        const double b_lat = (b - origin).dot(lateral);
        const double min_abs_lat = a_lat * b_lat <= 0.0
            ? 0.0 : std::min(std::abs(a_lat), std::abs(b_lat));
        if (min_abs_lat > 4.0)
            continue;
        candidate.near = std::min(candidate.near, std::max(0.0, near));
        candidate.far = std::max(candidate.far, far);
        candidate.found = true;
    }
}

void considerPolygon(
        const Polygon2d& polygon, const Vec2d& origin,
        const Vec2d& forward, const Vec2d& lateral, CrosswalkClearanceResult& candidate) {
    if (polygon.outer.size() >= 2)
        considerPolyline(toVec2dArray(polygon.outer), origin, forward, lateral, candidate);
    for (const auto& hole : polygon.holes)
        considerPolyline(toVec2dArray(hole), origin, forward, lateral, candidate);
}

}  // namespace

bool crosswalkRelevantToUTurn(
        const Crosswalk& crosswalk,
        const Vec2d& turn_entry, const Vec2d& turn_exit,
        const CrosswalkClearanceResult& candidate,
        const CrosswalkClearanceResult* paired_candidate,
        double chord_distance,
        double paired_near_delta) {
    if (!candidate.found)
        return false;

    // 首选原有端点弦走廊判据：它能排除射线擦到相邻道路横道的情况。
    const std::vector<Vec2d> points = toVec2dArray(crosswalk.geometry.outer);
    for (size_t i = 1; i < points.size(); ++i) {
        const Vec2d& a = points[i - 1];
        const Vec2d& b = points[i];
        if (pointToSegment(a, turn_entry, turn_exit).first <= chord_distance ||
            pointToSegment(b, turn_entry, turn_exit).first <= chord_distance ||
            pointToSegment(turn_entry, a, b).first <= chord_distance ||
            pointToSegment(turn_exit, a, b).first <= chord_distance)
            return true;
    }

    // 长 U-turn 的端点弦可能远离横道，但若入口和出口射线在同一近端
    // 站位命中同一横道，说明该横道横跨整条掉头走廊，不能被弦距离误杀。
    // 近端投影差使用 1m 容差，覆盖车道轻微斜交，同时拒绝 100000547
    // 中“入口命中 6m、出口命中 10m”的相邻道路擦边误命中。
    // 极短掉头（端点间距约 0.3m）即使两条射线擦到同一宽横道，也不能
    // 仅凭近端投影对齐把整块横道纳入清距，否则会让中弧被横道集合清空。
    // 这类小走廊沿用 2m/平齐站位即可；110003285 的最小目标走廊约 3.44m。
    if ((turn_exit - turn_entry).norm() <= 1.0)
        return false;
    return paired_candidate != nullptr && paired_candidate->found &&
           paired_candidate->crosswalk_id == candidate.crosswalk_id &&
           std::isfinite(candidate.near) &&
           std::isfinite(paired_candidate->near) &&
           std::abs(candidate.near - paired_candidate->near) <= paired_near_delta;
}

CrosswalkClearanceResult CrosswalkClearanceCalculator::alongRay(
    const Vec2d& origin, const Vec2d& direction,
    const IntersectionInput& input) const {
    const Vec2d forward = direction.norm() > 1e-8
        ? direction.normalized() : Vec2d(1, 0);
    const Vec2d lateral{-forward.y(), forward.x()};
    CrosswalkClearanceResult best;
    for (const auto& crosswalk : input.crosswalks) {
        CrosswalkClearanceResult candidate;
        candidate.crosswalk_id = crosswalk.id;
        considerPolygon(crosswalk.geometry, origin, forward, lateral, candidate);
        if (!candidate.found)
            continue;
        if (!best.found || candidate.near < best.near - 1e-6 ||
            (std::abs(candidate.near - best.near) <= 1e-6 &&
             candidate.far < best.far))
            best = candidate;
    }
    if (best.found)
        best.clearance = best.far + 0.30;
    return best;
}

CrosswalkClearanceResult CrosswalkClearanceCalculator::ahead(
    const Vec2d& origin, const Vec2d& tangent, const IntersectionInput& input) const {
    return alongRay(origin, tangent, input);
}

CrosswalkClearanceResult CrosswalkClearanceCalculator::behind(
    const Vec2d& endpoint, const Vec2d& tangent,
    const IntersectionInput& input) const {
    const Vec2d forward = tangent.norm() > 1e-8
        ? tangent.normalized() : Vec2d(1, 0);
    return alongRay(endpoint, -forward, input);
}

}  // 命名空间 isg
