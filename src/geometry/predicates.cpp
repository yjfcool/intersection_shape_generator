#include "geometry/predicates.h"

#include "curve/curve_utils.h"
#include "utils.h"

#include <algorithm>
#include <limits>

namespace isg {
namespace {

bool pointInRing(const std::vector<Vec3d>& ring, const Vec2d& p) {
    if (ring.size() < 3) return false;
    bool inside = false;
    for (size_t i = 0, j = ring.size() - 1; i < ring.size(); j = i++) {
        const Vec2d a = ring[i].xy();
        const Vec2d b = ring[j].xy();
        if (((a.y() > p.y()) != (b.y() > p.y())) &&
            (p.x() < (b.x() - a.x()) * (p.y() - a.y()) /
                         (b.y() - a.y() + 1e-20) + a.x()))
            inside = !inside;
    }
    return inside;
}

bool pointInPolygon(const Polygon2d& polygon, const Vec2d& p) {
    if (!pointInRing(polygon.outer, p)) return false;
    for (const auto& hole : polygon.holes)
        if (pointInRing(hole, p)) return false;
    return true;
}

const Polygon2d& hardObstacleGeometry(const Obstacle& obstacle) {
    return obstacle.geometry.outer.empty() ? obstacle.buffered_geometry : obstacle.geometry;
}

}  // namespace

std::vector<Vec2d> sampleCurveForAudit(const BezierCurve& curve, int samples) {
    if (curve.empty()) return std::vector<Vec2d>();
    return curve.sample(std::max(2, samples));
}

bool curveSelfIntersectsForAudit(const BezierCurve& curve, double endpoint_tol) {
    return curveSelfIntersectsBusiness(curve, endpoint_tol);
}

bool curveIntersectsPolygonForAudit(const BezierCurve& curve, const Polygon2d& polygon,
                                    int samples, Vec2d* location) {
    if (curve.empty() || polygon.outer.size() < 3) return false;
    const std::vector<Vec2d> points = sampleCurveForAudit(curve, samples);
    if (points.size() < 2) return false;
    for (size_t i = 0; i < points.size(); ++i) {
        if (i != 0 && i + 1 != points.size() && pointInPolygon(polygon, points[i])) {
            if (location) *location = points[i];
            return true;
        }
    }
    for (size_t i = 0; i + 1 < points.size(); ++i) {
        for (size_t j = 0; j < polygon.outer.size(); ++j) {
            Vec2d hit;
            if (segmentsIntersect(points[i], points[i + 1], polygon.outer[j],
                                   polygon.outer[(j + 1) % polygon.outer.size()], &hit)) {
                if (location) *location = hit;
                return true;
            }
        }
    }
    return false;
}

bool curveIntersectsObstaclesForAudit(const BezierCurve& curve,
                                      const std::vector<Obstacle>& obstacles,
                                      Vec2d* location) {
    for (const auto& obstacle : obstacles) {
        if (curveIntersectsPolygonForAudit(curve, hardObstacleGeometry(obstacle), 64, location))
            return true;
    }
    return false;
}

bool curvesHaveForbiddenAdherenceForAudit(const BezierCurve& a,
                                          const BezierCurve& b,
                                          double endpoint_tol,
                                          double adherence_tol) {
    const std::vector<Vec2d> ap = sampleCurveForAudit(a, 64);
    const std::vector<Vec2d> bp = sampleCurveForAudit(b, 64);
    if (ap.size() < 2 || bp.size() < 2) return false;
    for (size_t i = 0; i + 1 < ap.size(); ++i) {
        const Vec2d ad = ap[i + 1] - ap[i];
        if (ad.norm() < 1e-8) continue;
        const Vec2d au = ad.normalized();
        for (size_t j = 0; j + 1 < bp.size(); ++j) {
            const Vec2d bd = bp[j + 1] - bp[j];
            if (bd.norm() < 1e-8 || std::abs(au.dot(bd.normalized())) < 0.96)
                continue;
            if (std::abs(cross2d(au, bp[j] - ap[i])) > adherence_tol ||
                std::abs(cross2d(au, bp[j + 1] - ap[i])) > adherence_tol)
                continue;
            const double s0 = (bp[j] - ap[i]).dot(au);
            const double s1 = (bp[j + 1] - ap[i]).dot(au);
            const double lo = std::max(0.0, std::min(s0, s1));
            const double hi = std::min(ad.norm(), std::max(s0, s1));
            if (lo > hi + adherence_tol) continue;
            const Vec2d witness = ap[i] + 0.5 * (lo + hi) * au;
            if ((witness - a.startPt()).norm() <= endpoint_tol ||
                (witness - a.endPt()).norm() <= endpoint_tol ||
                (witness - b.startPt()).norm() <= endpoint_tol ||
                (witness - b.endPt()).norm() <= endpoint_tol)
                continue;
            return true;
        }
    }
    return false;
}

double minimumCurveBoundaryDistanceForAudit(const BezierCurve& curve,
                                            const std::vector<Boundary>& boundaries,
                                            Boundary::Type type,
                                            int samples,
                                            double endpoint_exclusion,
                                            Vec2d* location) {
    const std::vector<Vec2d> points = sampleCurveForAudit(curve, samples);
    if (points.size() < 3) return 1e18;
    const double endpoint_tol = std::min(
        std::max(0.0, endpoint_exclusion), 1e-4);
    double minimum = 1e18;
    for (size_t i = 1; i + 1 < points.size(); ++i) {
        if ((points[i] - points.front()).norm() <= endpoint_tol ||
            (points[i] - points.back()).norm() <= endpoint_tol)
            continue;
        for (const auto& boundary : boundaries) {
            if (boundary.type != type) continue;
            const double distance = pointToPolyline(points[i], boundary.geometry.points);
            if (distance < minimum) {
                minimum = distance;
                if (location) *location = points[i];
            }
        }
    }
    return minimum;
}

RoadEdgeClearanceMeasure measureCurveRoadEdgeClearanceForAudit(
    const BezierCurve& curve, const std::vector<Boundary>& boundaries,
    Boundary::Type type, double clearance, int samples,
    double endpoint_exclusion) {
    RoadEdgeClearanceMeasure measure;
    const std::vector<Vec2d> points = sampleCurveForAudit(curve, samples);
    if (points.size() < 3) return measure;
    auto nearestEdge = [&](const Vec2d& p) -> double {
        double best = std::numeric_limits<double>::infinity();
        for (const auto& boundary : boundaries) {
            if (boundary.type != type) continue;
            best = std::min(best, pointToPolyline(p, boundary.geometry.points));
        }
        return best;
    };
    measure.start_distance = nearestEdge(points.front());
    measure.end_distance = nearestEdge(points.back());
    const Vec2d start_tan = curve.startTan().norm() > 1e-12
        ? curve.startTan().normalized() : Vec2d(1, 0);
    const Vec2d end_tan = curve.endTan().norm() > 1e-12
        ? curve.endTan().normalized() : Vec2d(1, 0);
    const double endpoint_tol = std::min(
        std::max(0.0, endpoint_exclusion), 1e-4);
    for (size_t i = 1; i + 1 < points.size(); ++i) {
        const double to_start = (points[i] - points.front()).norm();
        const double to_end = (points[i] - points.back()).norm();
        if (to_start <= endpoint_tol || to_end <= endpoint_tol)
            continue;
        const double distance = nearestEdge(points[i]);
        if (!std::isfinite(distance))
            continue;
        measure.valid = true;
        if (distance < measure.minimum) {
            measure.minimum = distance;
            measure.location = points[i];
        }
        if (distance >= clearance)
            continue;
        // 取较近的端点作为参照：它决定这一点上净距的可达下限。
        const bool near_start = to_start <= to_end;
        const double span = near_start ? to_start : to_end;
        double ray_distance = -1.0;
        if (roadEdgeClearanceNeedsFloor(span)) {
            const Vec2d origin = near_start ? points.front() : points.back();
            const Vec2d direction = near_start ? start_tan : -end_tan;
            ray_distance = nearestEdge(origin + direction * span);
        }
        const double floor_value = roadEdgeClearanceFloor(
            clearance, span, ray_distance);
        if (floor_value - distance > measure.deficit) {
            measure.deficit = floor_value - distance;
            measure.deficit_location = points[i];
        }
    }
    return measure;
}

}  // 命名空间 isg
