#include "preprocessing/crosswalk_clearance_calculator.h"

#include "utils.h"

#include <algorithm>

namespace isg {
namespace {

void considerPolyline(const std::vector<Vec2d>& points,
                      const Vec2d& origin, const Vec2d& forward,
                      const Vec2d& lateral,
                      CrosswalkClearanceResult& candidate) {
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

void considerPolygon(const Polygon2d& polygon, const Vec2d& origin,
                     const Vec2d& forward, const Vec2d& lateral,
                     CrosswalkClearanceResult& candidate) {
    if (polygon.outer.size() >= 2)
        considerPolyline(toVec2dArray(polygon.outer), origin, forward, lateral, candidate);
    for (const auto& hole : polygon.holes)
        considerPolyline(toVec2dArray(hole), origin, forward, lateral, candidate);
}

}  // namespace

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
    const Vec2d& origin, const Vec2d& tangent,
    const IntersectionInput& input) const {
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
