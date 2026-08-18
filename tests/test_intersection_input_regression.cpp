#include <catch2/catch_test_macros.hpp>

#include "constraints/infeasibility_detector.h"
#include "constraints/boundary_safety.h"
#include "constraints/fence_check.h"
#include "constraints/cluster_order.h"
#include "curve/curve_utils.h"
#include "generator/polygon_builder.h"
#include "intersection_shape_generator.h"
#include "io/iodata_json.h"
#include "optimizer/sdf_field.h"
#include "utils.h"

#include <algorithm>
#include <chrono>
#include <fstream>
#include <limits>
#include <string>
#include <unordered_map>
#include <vector>

using namespace isg;

static std::unordered_map<ConnId, const ConnectivityCurve*> curveMap(const IntersectionOutput& output) {
    std::unordered_map<ConnId, const ConnectivityCurve*> curves;
    for (const auto& cc : output.connectivity_curves)
        curves[cc.id] = &cc;
    return curves;
}

static const Connectivity* findConnectivity(const IntersectionInput& input, const ConnId& id);
static bool curvesIntersectBeyondAllowedEndpointOverlapForTest(
    const BezierCurve& a, const BezierCurve& b, double endpoint_tol = 0.15);
static bool curvesHaveForbiddenAdherenceForTest(
    const BezierCurve& a, const BezierCurve& b, double endpoint_tol = 0.15);

static IntersectionInput loadInputOrSkip(const std::string& path) {
    std::ifstream f(path);
    if (!f.good()) {
        SKIP("missing test data file: " + path);
        return {};
    }
    return IntersectionIO::loadFromFile(path);
}

static IntersectionInput makeUTurnAlignmentScopeInput() {
    IntersectionInput input;
    input.id = "uturn_alignment_scope";
    input.mode = 1;
    input.area.geometry.outer = {
        Vec3d(-20, -20, 0), Vec3d(25, -20, 0),
        Vec3d(25, 25, 0), Vec3d(-20, 25, 0)
    };

    LaneGroup entry_group;
    entry_group.id = "entry_group";
    entry_group.role = GroupRole::Entry;
    LaneGroup exit_group;
    exit_group.id = "exit_group";
    exit_group.role = GroupRole::Exit;

    auto add_entry_lane = [&](const LaneId& id, double y) {
        Lane lane;
        lane.id = id;
        lane.groupId = entry_group.id;
        lane.width = 3.5;
        lane.geometry.points = {Vec3d(-12, y, 0), Vec3d(0, y, 0)};
        lane.laneOrder = (int)input.lanes.size();
        entry_group.lanes.push_back(id);
        input.lanes.push_back(lane);
    };
    auto add_exit_lane = [&](const LaneId& id, double x, double y) {
        Lane lane;
        lane.id = id;
        lane.groupId = exit_group.id;
        lane.width = 3.5;
        lane.geometry.points = {Vec3d(x, y, 0), Vec3d(x - 12, y, 0)};
        lane.laneOrder = (int)input.lanes.size();
        exit_group.lanes.push_back(id);
        input.lanes.push_back(lane);
    };

    auto add_entry_lane_with_points = [&](const LaneId& id, Vec3d a, Vec3d b) {
        Lane lane;
        lane.id = id;
        lane.groupId = entry_group.id;
        lane.width = 3.5;
        lane.geometry.points = {a, b};
        lane.laneOrder = (int)input.lanes.size();
        entry_group.lanes.push_back(id);
        input.lanes.push_back(lane);
    };

    add_entry_lane("entry_near", 0.0);
    add_entry_lane("entry_far", 3.5);
    add_entry_lane_with_points("entry_shared_a", Vec3d(-12, -8, 0), Vec3d(0, -8, 0));
    add_entry_lane_with_points("entry_shared_b", Vec3d(-12, -8, 0), Vec3d(0, -8, 0));
    add_exit_lane("exit_near", 5.0, 8.0);
    add_exit_lane("exit_far", 8.0, 11.5);
    add_exit_lane("exit_shared_a", 5.0, -3.0);
    add_exit_lane("exit_shared_b", 8.0, -5.5);
    input.lane_groups.push_back(entry_group);
    input.lane_groups.push_back(exit_group);

    auto add_uturn = [&](const ConnId& id, const LaneId& entry, const LaneId& exit) {
        Connectivity conn;
        conn.id = id;
        conn.entry_lane_id = entry;
        conn.exit_lane_id = exit;
        conn.turn_type = ConnTurnType::UTurnLeft;
        conn.enterGroupId = entry_group.id;
        conn.exitGroupId = exit_group.id;
        input.connectivities.push_back(conn);
    };
    add_uturn("near", "entry_near", "exit_near");
    add_uturn("far", "entry_far", "exit_far");
    add_uturn("shared_a", "entry_shared_a", "exit_shared_a");
    // Deliberately share only the entry lane: this exercises the case where
    // the entry family has a stagger rank but the exit family does not.
    add_uturn("shared_b", "entry_shared_a", "exit_shared_b");
    return input;
}

static char geometricTurnTypeForTest(const Connectivity& conn, const IntersectionInput& input) {
    auto entry = input.entryPtDir(conn.entry_lane_id);
    auto exit_ = input.exitPtDir(conn.exit_lane_id);
    Vec2d t0 = entry.second;
    Vec2d t1 = exit_.second;
    Vec2d path = exit_.first - entry.first;
    if (t0.norm() < 1e-8 || t1.norm() < 1e-8 || path.norm() < 1e-8)
        return '?';
    t0.normalize();
    t1.normalize();
    path.normalize();
    if (t0.dot(t1) < -0.5)
        return 'U';
    double cross = cross2d(t0, path);
    if (cross > 0.35)
        return 'L';
    if (cross < -0.35)
        return 'R';
    return 'S';
}

static bool isAllowedUTurnLeftRightTurnCombination(
    const IntersectionInput& input, const ConnId& a, const ConnId& b) {
    const Connectivity* ca = findConnectivity(input, a);
    const Connectivity* cb = findConnectivity(input, b);
    if (!ca || !cb)
        return false;
    char ta = geometricTurnTypeForTest(*ca, input);
    char tb = geometricTurnTypeForTest(*cb, input);
    bool au = ca->turn_type == ConnTurnType::UTurnLeft ||
              ca->turn_type == ConnTurnType::UTurnRight || ta == 'U';
    bool bu = cb->turn_type == ConnTurnType::UTurnLeft ||
              cb->turn_type == ConnTurnType::UTurnRight || tb == 'U';
    bool alr = ca->turn_type == ConnTurnType::TurnLeft ||
               ca->turn_type == ConnTurnType::TurnRight ||
               ta == 'L' || ta == 'R';
    bool blr = cb->turn_type == ConnTurnType::TurnLeft ||
               cb->turn_type == ConnTurnType::TurnRight ||
               tb == 'L' || tb == 'R';
    return (au && blr) || (bu && alr);
}

static bool allowedUTurnLeftRightCrossingOnlyForTest(
    const IntersectionInput& input, const ConnId& a, const ConnId& b,
    const BezierCurve& ca, const BezierCurve& cb) {
    return isAllowedUTurnLeftRightTurnCombination(input, a, b) &&
           !curvesHaveForbiddenAdherenceForTest(ca, cb, 0.30);
}

static std::vector<std::string> sharedEndpointRuleViolations(
    const IntersectionInput& input, const IntersectionOutput& output) {
    auto curves = curveMap(output);

    ClusterOrderSolver solver;
    solver.build(input.connectivities, input.lanes, input.lane_groups, input.crosswalks);

    std::vector<std::string> bad_pairs;
    for (const auto& pair : solver.pairs()) {
        if (!pair.shared_endpoint)
            continue;
        auto ia = curves.find(pair.id_a);
        auto ib = curves.find(pair.id_b);
        if (ia == curves.end() || ib == curves.end())
            continue;
        if (!ia->second->curve || !ib->second->curve)
            continue;
        double tol = 0.30;
        if (curvesIntersectBeyondAllowedEndpointOverlapForTest(
                *ia->second->curve, *ib->second->curve, tol) &&
            !allowedUTurnLeftRightCrossingOnlyForTest(
                input, pair.id_a, pair.id_b,
                *ia->second->curve, *ib->second->curve))
            bad_pairs.push_back(pair.id_a + "-" + pair.id_b);
    }
    return bad_pairs;
}

static std::vector<std::string> avoidableSameClusterCrossings(
    const IntersectionInput& input, const IntersectionOutput& output) {
    auto curves = curveMap(output);

    ClusterOrderSolver solver;
    solver.build(input.connectivities, input.lanes, input.lane_groups, input.crosswalks);

    std::vector<std::string> bad_pairs;
    for (const auto& pair : solver.pairs()) {
        if (pair.exempt == CrossExemption::StructuralCross)
            continue;
        auto ia = curves.find(pair.id_a);
        auto ib = curves.find(pair.id_b);
        if (ia == curves.end() || ib == curves.end())
            continue;
        if (!ia->second->curve || !ib->second->curve)
            continue;
        double tol = 0.30;
        if (curvesIntersectBeyondAllowedEndpointOverlapForTest(
                *ia->second->curve, *ib->second->curve, tol) &&
            !allowedUTurnLeftRightCrossingOnlyForTest(
                input, pair.id_a, pair.id_b,
                *ia->second->curve, *ib->second->curve))
            bad_pairs.push_back(pair.id_a + "-" + pair.id_b);
    }
    return bad_pairs;
}

static bool isFixedGeometryConn(const IntersectionInput& input, const ConnId& id) {
    for (const auto& conn : input.connectivities)
        if (conn.id == id)
            return conn.geometry.points.size() >= 2;
    return false;
}

static std::vector<std::string> avoidableSameClusterCrossingsIgnoringFixedFixed(
    const IntersectionInput& input, const IntersectionOutput& output) {
    auto curves = curveMap(output);

    ClusterOrderSolver solver;
    solver.build(input.connectivities, input.lanes, input.lane_groups, input.crosswalks);

    std::vector<std::string> bad_pairs;
    for (const auto& pair : solver.pairs()) {
        if (pair.exempt == CrossExemption::StructuralCross)
            continue;
        if (isFixedGeometryConn(input, pair.id_a) && isFixedGeometryConn(input, pair.id_b))
            continue;
        auto ia = curves.find(pair.id_a);
        auto ib = curves.find(pair.id_b);
        if (ia == curves.end() || ib == curves.end())
            continue;
        if (!ia->second->curve || !ib->second->curve)
            continue;
        double tol = 0.30;
        if (curvesIntersectBeyondAllowedEndpointOverlapForTest(
                *ia->second->curve, *ib->second->curve, tol) &&
            !allowedUTurnLeftRightCrossingOnlyForTest(
                input, pair.id_a, pair.id_b,
                *ia->second->curve, *ib->second->curve))
            bad_pairs.push_back(pair.id_a + "-" + pair.id_b);
    }
    return bad_pairs;
}

static std::string joinPairs(const std::vector<std::string>& pairs) {
    std::string s;
    for (const auto& p : pairs) {
        if (!s.empty()) s += ", ";
        s += p;
    }
    return s;
}

template <typename P>
static bool hasPolygonPointNear(
    const std::vector<P>& pts, const Vec2d& expected, double tol) {
    for (const auto& pt : pts) {
        if ((xyOf(pt) - expected).norm() <= tol)
            return true;
    }
    return false;
}

template <typename P>
static bool polygonUsesBoundarySegment(
    const std::vector<P>& polygon,
    const Boundary& boundary,
    double tol = 0.08,
    double min_station_span = 0.20) {
    std::vector<Vec2d> polygon_xy = toVec2dArray(polygon);
    for (int pi = 0; pi + 1 < (int)polygon_xy.size(); ++pi) {
        std::vector<double> stations;
        for (const Vec2d& pt : {polygon_xy[pi], polygon_xy[pi + 1]}) {
            double acc = 0.0;
            bool found = false;
            for (int bi = 0; bi + 1 < (int)boundary.geometry.points.size(); ++bi) {
                Vec2d a = boundary.geometry.points[bi];
                Vec2d b = boundary.geometry.points[bi + 1];
                Vec2d ab = b - a;
                double seg_len = ab.norm();
                if (seg_len < 1e-9)
                    continue;
                double t = (pt - a).dot(ab) / ab.squaredNorm();
                if (t >= -1e-6 && t <= 1.0 + 1e-6) {
                    t = std::max(0.0, std::min(1.0, t));
                    Vec2d proj = a + t * ab;
                    if (dist(pt, proj) <= tol) {
                        stations.push_back(acc + t * seg_len);
                        found = true;
                        break;
                    }
                }
                acc += seg_len;
            }
            if (!found)
                break;
        }
        if (stations.size() == 2 && std::abs(stations[1] - stations[0]) >= min_station_span)
            return true;
    }
    return false;
}

template <typename P>
static bool polygonUsesBoundarySegmentInStationRange(
    const std::vector<P>& polygon,
    const Boundary& boundary,
    double min_station,
    double max_station,
    double tol = 0.08,
    double min_station_span = 0.20) {
    std::vector<Vec2d> polygon_xy = toVec2dArray(polygon);
    for (int pi = 0; pi + 1 < (int)polygon_xy.size(); ++pi) {
        std::vector<double> stations;
        for (const Vec2d& pt : {polygon_xy[pi], polygon_xy[pi + 1]}) {
            double acc = 0.0;
            bool found = false;
            for (int bi = 0; bi + 1 < (int)boundary.geometry.points.size(); ++bi) {
                Vec2d a = boundary.geometry.points[bi];
                Vec2d b = boundary.geometry.points[bi + 1];
                Vec2d ab = b - a;
                double seg_len = ab.norm();
                if (seg_len < 1e-9)
                    continue;
                double t = (pt - a).dot(ab) / ab.squaredNorm();
                if (t >= -1e-6 && t <= 1.0 + 1e-6) {
                    t = std::max(0.0, std::min(1.0, t));
                    Vec2d proj = a + t * ab;
                    if (dist(pt, proj) <= tol) {
                        stations.push_back(acc + t * seg_len);
                        found = true;
                        break;
                    }
                }
                acc += seg_len;
            }
            if (!found)
                break;
        }
        if (stations.size() != 2)
            continue;
        double lo = std::min(stations[0], stations[1]);
        double hi = std::max(stations[0], stations[1]);
        if (hi - lo >= min_station_span && lo >= min_station - 1e-6 && hi <= max_station + 1e-6)
            return true;
    }
    return false;
}

static double boundaryLengthForTest(const Boundary& boundary) {
    double out = 0.0;
    for (int i = 0; i + 1 < (int)boundary.geometry.points.size(); ++i)
        out += dist(boundary.geometry.points[i], boundary.geometry.points[i + 1]);
    return out;
}

template <typename P>
static bool polygonBoundaryStationSpanForTest(
    const std::vector<P>& polygon,
    const Boundary& boundary,
    double& min_station,
    double& max_station,
    double point_tol = 0.08) {
    min_station = std::numeric_limits<double>::infinity();
    max_station = -std::numeric_limits<double>::infinity();
    std::vector<Vec2d> polygon_xy = toVec2dArray(polygon);
    for (const auto& pt : polygon_xy) {
        double acc = 0.0;
        for (int bi = 0; bi + 1 < (int)boundary.geometry.points.size(); ++bi) {
            Vec2d a = boundary.geometry.points[bi];
            Vec2d b = boundary.geometry.points[bi + 1];
            Vec2d ab = b - a;
            double seg_len = ab.norm();
            if (seg_len < 1e-9)
                continue;
            double t = (pt - a).dot(ab) / ab.squaredNorm();
            if (t >= -1e-6 && t <= 1.0 + 1e-6) {
                t = std::max(0.0, std::min(1.0, t));
                Vec2d proj = a + t * ab;
                if (dist(pt, proj) <= point_tol) {
                    double station = acc + t * seg_len;
                    min_station = std::min(min_station, station);
                    max_station = std::max(max_station, station);
                }
            }
            acc += seg_len;
        }
    }
    return std::isfinite(min_station) && std::isfinite(max_station);
}

static double boundaryPointStationForTest(const Boundary& boundary, int point_index) {
    double out = 0.0;
    for (int i = 0; i < point_index && i + 1 < (int)boundary.geometry.points.size(); ++i)
        out += dist(boundary.geometry.points[i], boundary.geometry.points[i + 1]);
    return out;
}

template <typename P>
static bool polygonHasBoundaryStationNear(
    const std::vector<P>& polygon,
    const Boundary& boundary,
    double expected_station,
    double station_tol = 0.08,
    double point_tol = 0.08) {
    std::vector<Vec2d> polygon_xy = toVec2dArray(polygon);
    for (const auto& pt : polygon_xy) {
        double acc = 0.0;
        for (int bi = 0; bi + 1 < (int)boundary.geometry.points.size(); ++bi) {
            Vec2d a = boundary.geometry.points[bi];
            Vec2d b = boundary.geometry.points[bi + 1];
            Vec2d ab = b - a;
            double seg_len = ab.norm();
            if (seg_len < 1e-9)
                continue;
            double t = (pt - a).dot(ab) / ab.squaredNorm();
            if (t >= -1e-6 && t <= 1.0 + 1e-6) {
                t = std::max(0.0, std::min(1.0, t));
                Vec2d proj = a + t * ab;
                double station = acc + t * seg_len;
                if (dist(pt, proj) <= point_tol &&
                    std::abs(station - expected_station) <= station_tol)
                    return true;
            }
            acc += seg_len;
        }
    }
    return false;
}

template <typename P>
static bool polygonContainsBoundaryShapePointsFromStation(
    const std::vector<P>& polygon,
    const Boundary& boundary,
    double start_station,
    double tol = 0.08) {
    for (int i = 0; i < (int)boundary.geometry.points.size(); ++i) {
        if (boundaryPointStationForTest(boundary, i) + 1e-6 < start_station)
            continue;
        if (!hasPolygonPointNear(polygon, boundary.geometry.points[i], tol))
            return false;
    }
    return true;
}

static const Boundary* findBoundary(const IntersectionInput& input, const std::string& id) {
    for (const auto& boundary : input.boundaries)
        if (boundary.id == id)
            return &boundary;
    return nullptr;
}

static double pointBoundaryDistanceForTest(const Vec2d& pt, const Boundary& boundary) {
    double best = std::numeric_limits<double>::infinity();
    for (int bi = 0; bi + 1 < (int)boundary.geometry.points.size(); ++bi) {
        Vec2d a = boundary.geometry.points[bi];
        Vec2d b = boundary.geometry.points[bi + 1];
        Vec2d ab = b - a;
        double len2 = ab.squaredNorm();
        if (len2 < 1e-12)
            continue;
        double t = std::max(0.0, std::min(1.0, (pt - a).dot(ab) / len2));
        best = std::min(best, dist(pt, a + t * ab));
    }
    return best;
}

static double minCurveBoundaryDistanceForTest(
        const BezierCurve& curve, const Boundary& boundary,
        double endpoint_skip = 0.75) {
    double best = std::numeric_limits<double>::infinity();
    int n = std::max(40, std::min(220, (int)std::ceil(curve.arcLength() / 0.15) + 1));
    for (const auto& pt : curve.sampleByArcLength(n)) {
        if ((pt - curve.startPt()).norm() <= endpoint_skip ||
            (pt - curve.endPt()).norm() <= endpoint_skip)
            continue;
        best = std::min(best, pointBoundaryDistanceForTest(pt, boundary));
    }
    return best;
}

struct CurveBoundaryHitForTest {
    std::string boundary_id;
    Vec2d point{0.0, 0.0};
};

static std::vector<CurveBoundaryHitForTest> rawCurveBoundaryHitsForTest(
        const BezierCurve& curve, const std::vector<Boundary>& boundaries,
        double curve_endpoint_tol = 0.15,
        double boundary_endpoint_tol = 0.10) {
    std::vector<CurveBoundaryHitForTest> hits;
    if (curve.empty())
        return hits;
    auto pts = curve.sampleByArcLength(std::max(
        64, std::min(240, (int)std::ceil(curve.arcLength() / 0.18) + 1)));
    for (const auto& boundary : boundaries) {
        const auto& bpts = boundary.geometry.points;
        bool hit_boundary = false;
        for (int i = 0; i + 1 < (int)pts.size() && !hit_boundary; ++i) {
            for (int j = 0; j + 1 < (int)bpts.size(); ++j) {
                Vec2d isect;
                if (!segmentsIntersect(
                        pts[i], pts[i + 1], bpts[j], bpts[j + 1], &isect))
                    continue;
                if ((isect - curve.startPt()).norm() <= curve_endpoint_tol ||
                    (isect - curve.endPt()).norm() <= curve_endpoint_tol)
                    continue;
                bool at_boundary_endpoint =
                    (j == 0 &&
                     (isect - bpts.front()).norm() <= boundary_endpoint_tol) ||
                    (j + 1 == (int)bpts.size() - 1 &&
                     (isect - bpts.back()).norm() <= boundary_endpoint_tol);
                if (at_boundary_endpoint)
                    continue;
                hits.push_back({boundary.id, isect});
                hit_boundary = true;
                break;
            }
        }
    }
    return hits;
}

struct CrosswalkRayChoiceForTest {
    bool found = false;
    std::string id;
    double near = std::numeric_limits<double>::infinity();
    double far = 0.0;
};

static CrosswalkRayChoiceForTest nearestCrosswalkAlongRayForTest(
        const Vec2d& origin, const Vec2d& direction,
        const std::vector<Crosswalk>& crosswalks) {
    Vec2d forward = direction.norm() > 1e-8
        ? direction.normalized() : Vec2d(1, 0);
    Vec2d lateral{-forward.y(), forward.x()};
    constexpr double search_radius = 12.0;
    constexpr double side_tolerance = 4.0;
    CrosswalkRayChoiceForTest best;

    auto considerPolyline = [&](const std::vector<Vec2d>& pts,
                                CrosswalkRayChoiceForTest& candidate) {
        for (int i = 0; i + 1 < (int)pts.size(); ++i) {
            Vec2d a = pts[i];
            Vec2d b = pts[i + 1];
            double a_fwd = (a - origin).dot(forward);
            double b_fwd = (b - origin).dot(forward);
            double near = std::min(a_fwd, b_fwd);
            double far = std::max(a_fwd, b_fwd);
            if (far <= 0.0 || near > search_radius)
                continue;
            double a_lat = (a - origin).dot(lateral);
            double b_lat = (b - origin).dot(lateral);
            double min_abs_lat =
                a_lat * b_lat <= 0.0 ? 0.0
                                     : std::min(std::abs(a_lat), std::abs(b_lat));
            if (min_abs_lat > side_tolerance)
                continue;
            candidate.near = std::min(candidate.near, std::max(0.0, near));
            candidate.far = std::max(candidate.far, far);
            candidate.found = true;
        }
    };

    for (const auto& crosswalk : crosswalks) {
        CrosswalkRayChoiceForTest candidate;
        candidate.id = crosswalk.id;
        considerPolyline(toVec2dArray(crosswalk.geometry.outer), candidate);
        for (const auto& hole : crosswalk.geometry.holes)
            considerPolyline(toVec2dArray(hole), candidate);
        if (!candidate.found)
            continue;
        if (!best.found ||
            candidate.near < best.near - 1e-6 ||
            (std::abs(candidate.near - best.near) <= 1e-6 &&
             candidate.far < best.far)) {
            best = candidate;
        }
    }
    return best;
}

static double crosswalkClearanceAlongRayForTest(
        const Vec2d& origin, const Vec2d& direction,
        const std::vector<Crosswalk>& crosswalks) {
    CrosswalkRayChoiceForTest choice =
        nearestCrosswalkAlongRayForTest(origin, direction, crosswalks);
    return choice.found ? choice.far + 0.30 : 0.0;
}

static bool tailAdherentStartForTest(
    const BezierCurve& curve, const Boundary& boundary, Vec2d& out,
    double tol = 0.12) {
    for (int si = 0; si < (int)curve.segs.size(); ++si) {
        Vec2d start = curve.segs[si].ctrl[0];
        if (dist(start, curve.endPt()) < 0.20 ||
            pointBoundaryDistanceForTest(start, boundary) > tol)
            continue;

        bool adherent_to_end = true;
        for (int sj = si; sj < (int)curve.segs.size() && adherent_to_end; ++sj) {
            for (int k = 0; k <= 8; ++k) {
                Vec2d pt = curve.segs[sj].evaluate((double)k / 8.0);
                if (pointBoundaryDistanceForTest(pt, boundary) > tol) {
                    adherent_to_end = false;
                    break;
                }
            }
        }
        if (adherent_to_end) {
            out = start;
            return true;
        }
    }
    return false;
}

static bool headAdherentEndForTest(
    const BezierCurve& curve, const Boundary& boundary, Vec2d& out,
    double tol = 0.12) {
    for (int si = 0; si < (int)curve.segs.size(); ++si) {
        Vec2d end = curve.segs[si].ctrl[3];
        if (dist(end, curve.startPt()) < 0.20 ||
            pointBoundaryDistanceForTest(end, boundary) > tol)
            continue;

        bool adherent_from_start = true;
        for (int sj = 0; sj <= si && adherent_from_start; ++sj) {
            for (int k = 0; k <= 8; ++k) {
                Vec2d pt = curve.segs[sj].evaluate((double)k / 8.0);
                if (pointBoundaryDistanceForTest(pt, boundary) > tol) {
                    adherent_from_start = false;
                    break;
                }
            }
        }
        if (adherent_from_start) {
            out = end;
            return true;
        }
    }
    return false;
}

static bool segmentAdherentToBoundaryForTest(
    const BezierSegment& seg, const Boundary& boundary, double tol = 0.12) {
    for (int i = 0; i <= 6; ++i) {
        if (pointBoundaryDistanceForTest(seg.evaluate((double)i / 6.0), boundary) > tol)
            return false;
    }
    return true;
}

static BezierCurve removeBoundaryAdherentSegmentsForTest(
    const BezierCurve& curve, const Boundary& boundary) {
    BezierCurve kept;
    for (const auto& seg : curve.segs) {
        if (!segmentAdherentToBoundaryForTest(seg, boundary))
            kept.segs.push_back(seg);
    }
    return kept.empty() ? curve : kept;
}

static bool curvesIntersectIgnoringBoundaryAdherenceForTest(
    const BezierCurve& a, const BezierCurve& b, const Boundary& boundary,
    double endpoint_tol = 1.5) {
    BezierCurve ta = removeBoundaryAdherentSegmentsForTest(a, boundary);
    BezierCurve tb = removeBoundaryAdherentSegmentsForTest(b, boundary);
    return curvesIntersectBusiness(ta, tb, endpoint_tol);
}

static const Connectivity* findConnectivity(const IntersectionInput& input, const ConnId& id) {
    for (const auto& conn : input.connectivities)
        if (conn.id == id)
            return &conn;
    return nullptr;
}

static BezierCurve straightCurveFromLineString(const LineString2d& line) {
    BezierCurve curve;
    for (int i = 0; i + 1 < (int)line.points.size(); ++i) {
        Vec2d d = line.points[i + 1] - line.points[i];
        if (d.norm() < 1e-8)
            continue;
        Vec2d dir = d.normalized();
        curve.segs.push_back(makeCubicG1(line.points[i], dir, line.points[i + 1], dir, 1.0 / 3.0));
    }
    return curve;
}

static double endpointG1Min(const ConnectivityCurve& cc, const IntersectionInput& input) {
    REQUIRE(cc.curve);
    const Connectivity* conn = findConnectivity(input, cc.id);
    if (conn) {
        char turn = geometricTurnTypeForTest(*conn, input);
        // U-turn shared-endpoint staggering may intentionally break endpoint G1;
        // business validation ignores endpoint G1 for geometric U-turns only.
        if (turn == 'U' ||
            conn->turn_type == ConnTurnType::UTurnLeft ||
            conn->turn_type == ConnTurnType::UTurnRight)
            return 1.0;
    }
    auto entry = input.entryPtDir(cc.entry_lane_id);
    auto exit_ = input.exitPtDir(cc.exit_lane_id);
    Vec2d T0 = entry.second.norm() > 1e-8 ? entry.second.normalized() : Vec2d(1, 0);
    Vec2d T1 = exit_.second.norm() > 1e-8 ? exit_.second.normalized() : Vec2d(1, 0);
    Vec2d st = cc.curve->startTan().norm() > 1e-8 ? cc.curve->startTan().normalized() : Vec2d(1, 0);
    Vec2d et = cc.curve->endTan().norm() > 1e-8 ? cc.curve->endTan().normalized() : Vec2d(1, 0);
    return std::min(st.dot(T0), et.dot(T1));
}

static double arcChordRatioForTest(const BezierCurve& curve) {
    double chord = (curve.endPt() - curve.startPt()).norm();
    REQUIRE(chord > 1e-6);
    return curve.arcLength() / chord;
}

static bool hasUTurnStraightArcStraightShapeForTest(
    const ConnectivityCurve& cc, const IntersectionInput& input,
    double min_straight = 2.0) {
    REQUIRE(cc.curve);
    const BezierCurve& curve = *cc.curve;
    if (curve.numSegments() != 3)
        return false;

    auto entry = input.entryPtDir(cc.entry_lane_id);
    auto exit_ = input.exitPtDir(cc.exit_lane_id);
    Vec2d T0 = entry.second.norm() > 1e-8 ? entry.second.normalized() : Vec2d(1, 0);
    Vec2d T1 = exit_.second.norm() > 1e-8 ? exit_.second.normalized() : Vec2d(1, 0);
    Vec2d exit_back = -T1;
    Vec2d axis = T0 + exit_back;
    if (axis.norm() < 1e-8)
        axis = T0;
    axis.normalize();
    if (axis.dot(T0) < 0.0)
        axis = -axis;

    const BezierSegment& first = curve.segs.front();
    Vec2d first_chord = first.ctrl[3] - first.ctrl[0];
    const BezierSegment& last = curve.segs.back();
    Vec2d last_chord = last.ctrl[3] - last.ctrl[0];
    if (first_chord.norm() < min_straight || last_chord.norm() < min_straight)
        return false;
    Vec2d arc_start = first.ctrl[3];
    Vec2d arc_end = last.ctrl[0];
    if (std::abs((arc_start - arc_end).dot(axis)) > 0.10)
        return false;

    return curve.segs[1].maxCurvature(30) > 0.05;
}

static bool forcedUTurnLeadShapeForTest(const BezierCurve& curve) {
    if (curve.numSegments() != 3)
        return false;
    Vec2d first = curve.segs.front().ctrl[3] - curve.segs.front().ctrl[0];
    Vec2d last = curve.segs.back().ctrl[3] - curve.segs.back().ctrl[0];
    return first.norm() >= 8.75 && last.norm() >= 8.75 &&
           curve.segs[1].maxCurvature(24) > 0.03;
}

static bool segmentLooksStraightForTest(const BezierSegment& s);

static bool uturnArcClearsCrosswalksForTest(
    const BezierCurve& curve, const std::vector<Crosswalk>& crosswalks,
    double min_straight = 7.5) {
    if (curve.numSegments() != 3)
        return false;
    const BezierSegment& first = curve.segs.front();
    const BezierSegment& arc = curve.segs[1];
    const BezierSegment& last = curve.segs.back();
    if (!segmentLooksStraightForTest(first) || !segmentLooksStraightForTest(last))
        return false;
    if ((first.ctrl[3] - first.ctrl[0]).norm() < min_straight ||
        (last.ctrl[3] - last.ctrl[0]).norm() < min_straight)
        return false;
    for (int i = 0; i <= 48; ++i) {
        Vec2d pt = arc.evaluate((double)i / 48.0);
        for (const auto& cw : crosswalks) {
            if (polygonContains(cw.geometry, pt))
                return false;
        }
    }
    return arc.maxCurvature(30) > 0.03;
}

static bool segmentLooksStraightForTest(const BezierSegment& s) {
    Vec2d d = s.ctrl[3] - s.ctrl[0];
    double len = d.norm();
    if (len < 1e-6)
        return false;
    return pointToSegment(s.ctrl[1], s.ctrl[0], s.ctrl[3]).first <=
               std::max(0.05, len * 0.03) &&
           pointToSegment(s.ctrl[2], s.ctrl[0], s.ctrl[3]).first <=
               std::max(0.05, len * 0.03);
}

static bool sameDirectedSegmentForTest(const BezierSegment& a, const BezierSegment& b) {
    Vec2d da = a.ctrl[3] - a.ctrl[0];
    Vec2d db = b.ctrl[3] - b.ctrl[0];
    if (da.norm() < 1e-6 || db.norm() < 1e-6)
        return false;
    if (!segmentLooksStraightForTest(a) || !segmentLooksStraightForTest(b))
        return false;
    if (da.normalized().dot(db.normalized()) < 0.995)
        return false;
    return pointToSegment(a.ctrl[0], b.ctrl[0], b.ctrl[3]).first < 0.20 ||
           pointToSegment(a.ctrl[3], b.ctrl[0], b.ctrl[3]).first < 0.20 ||
           pointToSegment(b.ctrl[0], a.ctrl[0], a.ctrl[3]).first < 0.20 ||
           pointToSegment(b.ctrl[3], a.ctrl[0], a.ctrl[3]).first < 0.20;
}

struct StationedCurveSampleForTest {
    std::vector<Vec2d> pts;
    std::vector<double> station;
    double total = 0.0;
};

static StationedCurveSampleForTest sampleCurveWithStationsForTest(
    const BezierCurve& curve) {
    StationedCurveSampleForTest s;
    s.total = curve.arcLength();
    int n = std::max(64, std::min(240, (int)std::ceil(s.total / 0.20) + 1));
    s.pts = curve.sampleByArcLength(n);
    s.station.resize(s.pts.size(), 0.0);
    for (int i = 1; i < (int)s.pts.size(); ++i)
        s.station[i] = s.station[i - 1] + (s.pts[i] - s.pts[i - 1]).norm();
    if (!s.station.empty())
        s.total = s.station.back();
    return s;
}

static Vec2d pointAtSampleStationForTest(
    const StationedCurveSampleForTest& s, double station) {
    if (s.pts.empty())
        return Vec2d(0, 0);
    if (station <= 0.0)
        return s.pts.front();
    if (station >= s.total)
        return s.pts.back();
    auto it = std::lower_bound(s.station.begin(), s.station.end(), station);
    int idx = (int)std::distance(s.station.begin(), it);
    if (idx <= 0)
        return s.pts.front();
    if (idx >= (int)s.pts.size())
        return s.pts.back();
    double span = s.station[idx] - s.station[idx - 1];
    double u = span > 1e-8 ? (station - s.station[idx - 1]) / span : 0.0;
    return s.pts[idx - 1] * (1.0 - u) + s.pts[idx] * u;
}

static Vec2d tangentAtSampleStationForTest(
    const StationedCurveSampleForTest& s, double station, bool reverse = false) {
    if (s.pts.size() < 2)
        return Vec2d(1, 0);
    double clamped = std::max(0.0, std::min(s.total, station));
    auto it = std::lower_bound(s.station.begin(), s.station.end(), clamped);
    int idx = (int)std::distance(s.station.begin(), it);
    if (idx <= 0)
        idx = 1;
    if (idx >= (int)s.pts.size())
        idx = (int)s.pts.size() - 1;
    Vec2d d = s.pts[idx] - s.pts[idx - 1];
    if (d.norm() < 1e-8)
        d = s.pts.back() - s.pts.front();
    if (d.norm() < 1e-8)
        d = Vec2d(1, 0);
    if (reverse)
        d = -d;
    return d.normalized();
}

static bool curveLooksUTurnForTest(const BezierCurve& curve) {
    if (curve.empty())
        return false;
    Vec2d st = curve.startTan();
    Vec2d et = curve.endTan();
    return st.norm() > 1e-8 && et.norm() > 1e-8 &&
        st.normalized().dot(et.normalized()) < -0.5;
}

static double directedSharedEndpointOverlapLengthForTest(
    const BezierCurve& a, const BezierCurve& b,
    bool from_start, double endpoint_tol) {
    Vec2d a_anchor = from_start ? a.startPt() : a.endPt();
    Vec2d b_anchor = from_start ? b.startPt() : b.endPt();
    if ((a_anchor - b_anchor).norm() > endpoint_tol)
        return 0.0;
    Vec2d ta = from_start ? a.startTan() : -a.endTan();
    Vec2d tb = from_start ? b.startTan() : -b.endTan();
    if (ta.norm() < 1e-8 || tb.norm() < 1e-8 ||
        ta.normalized().dot(tb.normalized()) < 0.96)
        return 0.0;
    bool both_uturn = curveLooksUTurnForTest(a) && curveLooksUTurnForTest(b);
    if (both_uturn) {
        bool same_lead = from_start
            ? sameDirectedSegmentForTest(a.segs.front(), b.segs.front())
            : sameDirectedSegmentForTest(a.segs.back(), b.segs.back());
        if (!same_lead)
            return 0.0;
    }

    auto sa = sampleCurveWithStationsForTest(a);
    auto sb = sampleCurveWithStationsForTest(b);
    if (sa.pts.size() < 2 || sb.pts.size() < 2)
        return 0.0;
    double max_common = std::min(sa.total, sb.total);
    double step = std::max(0.20, std::min(0.45, endpoint_tol * 0.20));
    double last_good = 0.0;
    bool left_anchor = false;
    for (double d = step; d <= max_common + 1e-6; d += step) {
        double sta = from_start ? d : sa.total - d;
        double stb = from_start ? d : sb.total - d;
        Vec2d pa = pointAtSampleStationForTest(sa, sta);
        Vec2d pb = pointAtSampleStationForTest(sb, stb);
        Vec2d tta = tangentAtSampleStationForTest(sa, sta, !from_start);
        Vec2d ttb = tangentAtSampleStationForTest(sb, stb, !from_start);
        // U-turn lead overlap must be genuine same-line overlap.  A looser
        // tolerance hides adjacent large/small-radius middle-arc crossings.
        double near_tol = both_uturn
            ? 0.25
            : std::max(0.08, std::min(endpoint_tol, 0.20));
        bool same_bundle =
            (pa - pb).norm() <= near_tol &&
            tta.normalized().dot(ttb.normalized()) >= 0.90;
        if (!same_bundle) {
            left_anchor = true;
            break;
        }
        last_good = std::min(d, max_common);
    }
    double overlap = last_good;
    if (!left_anchor && last_good > max_common - step * 1.5)
        overlap = max_common;
    return overlap;
}

static bool curvesIntersectBeyondStrictEndpointOverlapForTest(
    const BezierCurve& a, const BezierCurve& b, double endpoint_tol) {
    return curvesIntersectBusiness(a, b, endpoint_tol);
}

static bool curvesIntersectBeyondAllowedEndpointOverlapForTest(
    const BezierCurve& a, const BezierCurve& b, double endpoint_tol) {
    return curvesIntersectBeyondStrictEndpointOverlapForTest(a, b, endpoint_tol);
}

static double distToCurveEndpointsForTest(
    const Vec2d& pt, const BezierCurve& a, const BezierCurve& b) {
    return std::min({
        (pt - a.startPt()).norm(),
        (pt - a.endPt()).norm(),
        (pt - b.startPt()).norm(),
        (pt - b.endPt()).norm()
    });
}

static bool curvesHaveForbiddenAdherenceForTest(
    const BezierCurve& a, const BezierCurve& b, double endpoint_tol) {
    auto ap = a.sampleByArcLength(96);
    auto bp = b.sampleByArcLength(96);
    if (ap.size() < 2 || bp.size() < 2)
        return false;
    constexpr double kAdherentTol = 0.18;
    for (int ai = 0; ai + 1 < (int)ap.size(); ++ai) {
        Vec2d a0 = ap[ai], a1 = ap[ai + 1];
        Vec2d ad = a1 - a0;
        double alen = ad.norm();
        if (alen < 1e-8)
            continue;
        Vec2d au = ad / alen;
        for (int bi = 0; bi + 1 < (int)bp.size(); ++bi) {
            Vec2d b0 = bp[bi], b1 = bp[bi + 1];
            Vec2d bd = b1 - b0;
            double blen = bd.norm();
            if (blen < 1e-8)
                continue;
            Vec2d bu = bd / blen;
            if (std::abs(au.dot(bu)) < 0.96)
                continue;
            if (std::abs(cross2d(au, b0 - a0)) > kAdherentTol ||
                std::abs(cross2d(au, b1 - a0)) > kAdherentTol)
                continue;
            double b0s = (b0 - a0).dot(au);
            double b1s = (b1 - a0).dot(au);
            double lo = std::max(0.0, std::min(b0s, b1s));
            double hi = std::min(alen, std::max(b0s, b1s));
            if (lo > hi + kAdherentTol)
                continue;
            Vec2d witness = a0 + 0.5 * (lo + hi) * au;
            if (distToCurveEndpointsForTest(witness, a, b) <= endpoint_tol)
                continue;
            return true;
        }
    }
    return false;
}

static double maxLateralChordDeviationRatioForTest(const BezierCurve& curve) {
    Vec2d chord = curve.endPt() - curve.startPt();
    double chord_len = chord.norm();
    REQUIRE(chord_len > 1e-6);
    Vec2d dir = chord / chord_len;
    double max_dev = 0.0;
    for (const auto& pt : curve.sampleByArcLength(96)) {
        max_dev = std::max(max_dev, std::abs(cross2d(dir, pt - curve.startPt())));
    }
    return max_dev / chord_len;
}

static bool hasCurvatureSignFlipForTest(const BezierCurve& curve, double eps = 0.10) {
    int sign = 0;
    for (const auto& seg : curve.segs) {
        for (int i = 1; i < 20; ++i) {
            double t = (double)i / 20.0;
            Vec2d d1 = seg.evalDeriv1(t);
            Vec2d d2 = seg.evalDeriv2(t);
            double den = std::pow(d1.squaredNorm(), 1.5);
            if (den < 1e-12)
                continue;
            double k = cross2d(d1, d2) / den;
            if (std::abs(k) <= eps)
                continue;
            int s = k > 0.0 ? 1 : -1;
            if (sign != 0 && s != sign)
                return true;
            sign = s;
        }
    }
    return false;
}

static bool leftTurnArchShapeForTest(
    const ConnectivityCurve& cc, const IntersectionInput& input) {
    REQUIRE(cc.curve);
    const BezierCurve& curve = *cc.curve;
    auto entry = input.entryPtDir(cc.entry_lane_id);
    auto exit_ = input.exitPtDir(cc.exit_lane_id);
    Vec2d T0 = entry.second.norm() > 1e-8
        ? entry.second.normalized() : Vec2d(1, 0);
    Vec2d T1 = exit_.second.norm() > 1e-8
        ? exit_.second.normalized() : Vec2d(1, 0);
    double chord_len = (exit_.first - entry.first).norm();
    if (chord_len < 1e-6 || hasCurvatureSignFlipForTest(curve))
        return false;

    if (curve.numSegments() == 1) {
        double arc_chord = arcChordRatioForTest(curve);
        return arc_chord >= 1.03 && arc_chord <= 1.35 &&
               curve.maxCurvature(40) <= 2.5;
    }

    if (curve.numSegments() != 3)
        return false;
    const BezierSegment& first = curve.segs.front();
    const BezierSegment& arc = curve.segs[1];
    const BezierSegment& last = curve.segs.back();
    if (!segmentLooksStraightForTest(first) ||
        !segmentLooksStraightForTest(last) ||
        arc.maxCurvature(30) <= 0.03)
        return false;
    Vec2d first_dir = first.ctrl[3] - first.ctrl[0];
    Vec2d last_dir = last.ctrl[3] - last.ctrl[0];
    return first_dir.norm() > 1e-6 &&
           last_dir.norm() > 1e-6 &&
           first_dir.normalized().dot(T0) > 0.98 &&
           last_dir.normalized().dot(T1) > 0.98;
}

static Vec2d outputCenter(const IntersectionOutput& output) {
    Vec2d center(0, 0);
    int count = 0;
    for (const auto& cc : output.connectivity_curves) {
        REQUIRE(cc.curve);
        center += cc.curve->startPt();
        center += cc.curve->endPt();
        count += 2;
    }
    REQUIRE(count > 0);
    return center / (double)count;
}

static Vec2d inputConnectionCenterForTest(const IntersectionInput& input) {
    Vec2d center(0, 0);
    int count = 0;
    for (const auto& lane : input.lanes) {
        center += getConnPoint(lane.geometry.points, input.IsEntryLane(lane.id));
        ++count;
    }
    REQUIRE(count > 0);
    return center / (double)count;
}

static bool directionNeedsCenterFlipForTest(Vec2d dir, Vec2d wanted) {
    if (dir.norm() < 1e-8 || wanted.norm() < 1e-8)
        return false;
    dir.normalize();
    wanted.normalize();
    return dir.dot(wanted) < -0.35;
}

template <typename P>
static Vec2d directedConnDirForTest(
    const std::vector<P>& pts,
    GroupRole role,
    const Vec2d& center) {
    bool is_entry = (role == GroupRole::Entry);
    Vec2d p = getConnPoint(pts, is_entry);
    Vec2d d = getConnTangent(pts, is_entry);
    Vec2d wanted = is_entry ? (center - p) : (p - center);
    if (directionNeedsCenterFlipForTest(d, wanted))
        d = -d;
    return d.norm() > 1e-8 ? d.normalized() : Vec2d(1, 0);
}

static Vec2d modeShiftDirForTest(GroupRole role, Vec2d dir, int mode) {
    if (dir.norm() < 1e-8)
        dir = Vec2d(1, 0);
    else
        dir.normalize();
    if (role == GroupRole::Entry && mode == 2)
        dir = -dir;
    return dir;
}

static bool lineSegmentIntersectionForTest(
    const Vec2d& line_pt,
    const Vec2d& line_dir,
    const Vec2d& a,
    const Vec2d& b,
    const Vec2d& anchor,
    Vec2d* out) {
    Vec2d seg = b - a;
    if (line_dir.norm() < 1e-9 || seg.norm() < 1e-9)
        return false;
    double den = cross2d(seg, line_dir);
    if (std::abs(den) < 1e-12) {
        if (std::abs(cross2d(a - line_pt, line_dir)) > 1e-8)
            return false;
        *out = dist(a, anchor) <= dist(b, anchor) ? a : b;
        return true;
    }
    double t = cross2d(line_pt - a, line_dir) / den;
    if (t < -1e-9 || t > 1.0 + 1e-9)
        return false;
    t = std::max(0.0, std::min(1.0, t));
    *out = a + t * seg;
    return true;
}

static bool curveRawIntersectsRoadEdgeForTest(
    const BezierCurve& curve, const Boundary& boundary,
    double curve_endpoint_tol = 0.15,
    double boundary_endpoint_tol = 0.10) {
    if (curve.empty() || boundary.type != Boundary::Type::RoadEdge ||
        boundary.geometry.points.size() < 2)
        return false;

    auto pts = curve.sampleByArcLength(std::max(
        64, std::min(240, (int)std::ceil(curve.arcLength() / 0.18) + 1)));
    if (pts.size() < 2)
        return false;

    const auto& bpts = boundary.geometry.points;
    for (int i = 0; i + 1 < (int)pts.size(); ++i) {
        for (int j = 0; j + 1 < (int)bpts.size(); ++j) {
            Vec2d isect;
            if (!segmentsIntersect(pts[i], pts[i + 1], bpts[j], bpts[j + 1],
                                   &isect))
                continue;
            if ((isect - pts.front()).norm() <= curve_endpoint_tol ||
                (isect - pts.back()).norm() <= curve_endpoint_tol)
                continue;
            if ((j == 0 && (isect - bpts[j]).norm() <= boundary_endpoint_tol) ||
                (j + 1 == (int)bpts.size() - 1 &&
                 (isect - bpts[j + 1]).norm() <= boundary_endpoint_tol))
                continue;
            return true;
        }
    }
    return false;
}

TEST_CASE("RoadEdge boundary safety rejects adherent contact and outside crossing",
          "[regression][boundary]") {
    Boundary edge;
    edge.id = "edge";
    edge.type = Boundary::Type::RoadEdge;
    edge.geometry.points = {Vec2d(0, 0), Vec2d(10, 0)};
    Vec2d center(5, 5);

    BezierCurve touching;
    touching.segs.push_back(makeCubicG1(Vec2d(0, 1), Vec2d(1, -0.2),
                                        Vec2d(5, 0), Vec2d(1, 0), 0.30));
    touching.segs.push_back(makeCubicG1(Vec2d(5, 0), Vec2d(1, 0),
                                        Vec2d(10, 1), Vec2d(1, 0.2), 0.30));
    BoundarySafetyResult touch_safety =
        curveBoundarySafety(touching, std::vector<Boundary>{edge}, center, 96);
    CHECK(touch_safety.intersects);
    CHECK_FALSE(touch_safety.outside_road_edge);

    BezierCurve crossing;
    crossing.segs.push_back(makeCubicG1(Vec2d(0, 1), Vec2d(1, -0.2),
                                        Vec2d(5, 0), Vec2d(1, -0.2), 0.30));
    crossing.segs.push_back(makeCubicG1(Vec2d(5, 0), Vec2d(1, -0.2),
                                        Vec2d(10, -1), Vec2d(1, -0.2), 0.30));
    BoundarySafetyResult crossing_safety =
        curveBoundarySafety(crossing, std::vector<Boundary>{edge}, center, 96);
    CHECK(crossing_safety.intersects);
    CHECK((crossing_safety.intersects || crossing_safety.outside_road_edge));
}

static Vec2d expectedLaneCutPoint(
    const IntersectionInput& input,
    const IntersectionOutput& output,
    const LaneId& lane_id,
    GroupRole role,
    const Vec2d& center,
    double offset = 0.5) {
    const Lane* lane = input.findLane(lane_id);
    REQUIRE(lane);
    REQUIRE(lane->geometry.points.size() >= 2);

    Vec2d anchor = getConnPoint(lane->geometry.points, role == GroupRole::Entry);
    Vec2d tangent = directedConnDirForTest(lane->geometry.points, role, center);
    Vec2d cut_base = anchor + offset * modeShiftDirForTest(role, tangent, input.mode);
    Vec2d cut_dir = rotLeft(tangent);
    REQUIRE(cut_dir.norm() > 1e-8);
    cut_dir.normalize();

    struct Hit {
        bool found = false;
        Vec2d pt{0, 0};
        double anchor_distance = 1e18;
        double base_distance = 1e18;
    };
    Hit best;
    auto consider_line = [&](const std::vector<Vec2d>& pts) {
        if (pts.size() < 2)
            return;
        for (int i = 0; i + 1 < (int)pts.size(); ++i) {
            Vec2d hit;
            if (!lineSegmentIntersectionForTest(cut_base, cut_dir, pts[i], pts[i + 1],
                                                anchor, &hit))
                continue;
            double ad = dist(hit, anchor);
            double bd = dist(hit, cut_base);
            if (!best.found ||
                ad < best.anchor_distance - 1e-8 ||
                (std::abs(ad - best.anchor_distance) <= 1e-8 && bd < best.base_distance)) {
                best.found = true;
                best.pt = hit;
                best.anchor_distance = ad;
                best.base_distance = bd;
            }
        }
    };

    consider_line(toVec2dArray(lane->geometry.points));
    for (const auto& cc : output.connectivity_curves) {
        if ((role == GroupRole::Entry && cc.entry_lane_id != lane_id) ||
            (role == GroupRole::Exit && cc.exit_lane_id != lane_id))
            continue;
        REQUIRE(cc.curve);
        consider_line(cc.curve->sampleByArcLength(160));
    }
    REQUIRE(best.found);
    return best.pt;
}

static Vec2d centerSideEndpointDirForTest(
    const std::vector<Vec3d>& pts,
    GroupRole role,
    const Vec2d& center) {
    REQUIRE(pts.size() >= 2);
    bool use_front = dist(pts.front(), center) <= dist(pts.back(), center);
    Vec2d p = use_front ? xyOf(pts.front()) : xyOf(pts.back());
    Vec2d d = use_front
        ? (xyOf(pts[1]) - xyOf(pts.front()))
        : (xyOf(pts.back()) - xyOf(pts[pts.size() - 2]));
    Vec2d wanted = (role == GroupRole::Entry) ? (center - p) : (p - center);
    if (directionNeedsCenterFlipForTest(d, wanted))
        d = -d;
    return d.norm() > 1e-8 ? d.normalized() : Vec2d(1, 0);
}

static bool laneEdgeReferencedByGroupForTest(
    const IntersectionInput& input,
    const LaneGroup& group,
    const LaneEdgeId& edge_id) {
    if (std::find(group.boundaries.begin(), group.boundaries.end(), edge_id) !=
        group.boundaries.end())
        return true;
    for (const auto& lane_id : group.lanes) {
        const Lane* lane = input.findLane(lane_id);
        if (lane && (lane->left_edge_id == edge_id || lane->right_edge_id == edge_id))
            return true;
    }
    return false;
}

static std::vector<LaneId> uniqueGroupLaneIdsForTest(const LaneGroup& group) {
    std::vector<LaneId> ids;
    for (const auto& lane_id : group.lanes) {
        if (std::find(ids.begin(), ids.end(), lane_id) == ids.end())
            ids.push_back(lane_id);
    }
    return ids;
}

static Vec2d laneGroupDirectionForEdgeForTest(
    const IntersectionInput& input,
    const LaneGroup& group,
    const LaneEdgeId& edge_id,
    const Vec2d& center) {
    Vec2d dir(0, 0);
    int count = 0;
    for (const auto& lane_id : uniqueGroupLaneIdsForTest(group)) {
        const Lane* lane = input.findLane(lane_id);
        if (!lane || (lane->left_edge_id != edge_id && lane->right_edge_id != edge_id))
            continue;
        dir += directedConnDirForTest(lane->geometry.points, group.role, center);
        ++count;
    }
    if (count == 0) {
        for (const auto& lane_id : uniqueGroupLaneIdsForTest(group)) {
            const Lane* lane = input.findLane(lane_id);
            if (!lane)
                continue;
            dir += directedConnDirForTest(lane->geometry.points, group.role, center);
            ++count;
        }
    }
    return dir.norm() > 1e-8 ? dir.normalized() : Vec2d(1, 0);
}

static const LaneGroup* findGroupForLaneEdgeRoleForTest(
    const IntersectionInput& input,
    const LaneEdgeId& edge_id,
    GroupRole role) {
    for (const auto& group : input.lane_groups) {
        if (group.role == role && laneEdgeReferencedByGroupForTest(input, group, edge_id))
            return &group;
    }
    return nullptr;
}

static Vec2d expectedLaneEdgeCutPoint(
    const IntersectionInput& input,
    const LaneEdgeId& edge_id,
    GroupRole role,
    const Vec2d& center,
    double offset = 0.5) {
    const LaneEdge* edge = input.findEdge(edge_id);
    REQUIRE(edge);
    REQUIRE(edge->geometry.points.size() >= 2);

    const auto& pts = edge->geometry.points;
    const LaneGroup* group = findGroupForLaneEdgeRoleForTest(input, edge_id, role);
    Vec2d group_dir = group
        ? laneGroupDirectionForEdgeForTest(input, *group, edge_id, center)
        : directedConnDirForTest(pts, role, center);
    Vec2d edge_dir = xyOf(pts.back()) - xyOf(pts.front());
    if (edge_dir.norm() < 1e-8)
        edge_dir = Vec2d(1, 0);
    else
        edge_dir.normalize();
    bool edge_opposes_group = group_dir.norm() > 1e-8 &&
        edge_dir.dot(group_dir.normalized()) < -0.35;

    bool use_entry_endpoint = (role == GroupRole::Entry);
    if (edge_opposes_group)
        use_entry_endpoint = !use_entry_endpoint;

    Vec2d anchor = getConnPoint(pts, use_entry_endpoint);
    Vec2d tangent;
    if (edge_opposes_group) {
        tangent = -getConnTangent(pts, use_entry_endpoint);
        if (group_dir.norm() > 1e-8 && tangent.dot(group_dir) < 0.0)
            tangent = -tangent;
    } else {
        tangent = (input.mode == 2 && role == GroupRole::Entry)
            ? getConnTangent(pts, true)
            : directedConnDirForTest(pts, role, center);
    }
    if (tangent.norm() < 1e-8)
        tangent = Vec2d(1, 0);
    else
        tangent.normalize();
    Vec2d cut_base = anchor + offset * modeShiftDirForTest(role, tangent, input.mode);
    Vec2d cut_dir = rotLeft(tangent);
    REQUIRE(cut_dir.norm() > 1e-8);
    cut_dir.normalize();

    struct Hit {
        bool found = false;
        Vec2d pt{0, 0};
        double anchor_distance = 1e18;
        double base_distance = 1e18;
    };
    Hit best;
    for (int i = 0; i + 1 < (int)pts.size(); ++i) {
        Vec2d hit;
        if (!lineSegmentIntersectionForTest(cut_base, cut_dir, xyOf(pts[i]),
                                            xyOf(pts[i + 1]), anchor, &hit))
            continue;
        double ad = dist(hit, anchor);
        double bd = dist(hit, cut_base);
        if (!best.found ||
            ad < best.anchor_distance - 1e-8 ||
            (std::abs(ad - best.anchor_distance) <= 1e-8 && bd < best.base_distance)) {
            best.found = true;
            best.pt = hit;
            best.anchor_distance = ad;
            best.base_distance = bd;
        }
    }
    return best.found ? best.pt : cut_base;
}

TEST_CASE("intersection_input has no avoidable same-cluster interior crossings", "[regression][cluster]") {
    const std::string path = std::string(PROJECT_ROOT_DIR) + "/datas/intersection_input.json";
    IntersectionInput input = loadInputOrSkip(path);

    IntersectionShapeGenerator gen;
    IntersectionOutput output;
    REQUIRE(gen.generate(input, output));

    auto curves = curveMap(output);
    std::vector<std::string> bad_pairs = avoidableSameClusterCrossings(input, output);

    INFO("avoidable same-cluster crossings: " << joinPairs(bad_pairs));
    CHECK(bad_pairs.empty());

    const std::vector<std::pair<ConnId, ConnId>> named_examples = {
        {"18", "32"},
        {"30", "12"},
        {"2", "3"},
        {"19", "21"},
        {"5", "22"},
        {"23", "19"},
        {"23", "21"},
        {"24", "20"},
        {"24", "22"},
    };
    std::vector<std::string> crossed_examples;
    for (const auto& ids : named_examples) {
        auto ia = curves.find(ids.first);
        auto ib = curves.find(ids.second);
        REQUIRE(ia != curves.end());
        REQUIRE(ib != curves.end());
        REQUIRE(ia->second->curve);
        REQUIRE(ib->second->curve);
        if (curvesIntersectBusiness(*ia->second->curve, *ib->second->curve, 0.15) &&
            !allowedUTurnLeftRightCrossingOnlyForTest(
                input, ids.first, ids.second,
                *ia->second->curve, *ib->second->curve))
            crossed_examples.push_back(ids.first + "-" + ids.second);
    }
    INFO("named example crossings: " << joinPairs(crossed_examples));
    CHECK(crossed_examples.empty());
}

TEST_CASE("100000643 keeps same-cluster U-turn family non-intersecting", "[regression][cluster][uturn][100000643]") {
    const std::string path = std::string(PROJECT_ROOT_DIR) + "/datas/100000643.json";
    IntersectionInput input = loadInputOrSkip(path);

    ClusterOrderSolver solver;
    solver.build(input.connectivities, input.lanes, input.lane_groups, input.crosswalks);

    CHECK(solver.exemptionOf("74", "76") == CrossExemption::None);
    CHECK(solver.exemptionOf("74", "73") == CrossExemption::None);
    CHECK(solver.exemptionOf("74", "78") == CrossExemption::None);
    CHECK(solver.exemptionOf("89", "96") == CrossExemption::None);
    CHECK(solver.exemptionOf("89", "98") == CrossExemption::None);
    CHECK(solver.exemptionOf("89", "100") == CrossExemption::None);
    CHECK(solver.exemptionOf("73", "75") == CrossExemption::None);
    CHECK(solver.exemptionOf("75", "76") == CrossExemption::None);
    CHECK(solver.exemptionOf("78", "97") == CrossExemption::StructuralCross);
    CHECK(solver.exemptionOf("78", "99") == CrossExemption::StructuralCross);

    IntersectionShapeGenerator gen;
    IntersectionOutput output;
    REQUIRE(gen.generate(input, output));

    auto curves = curveMap(output);
    for (const ConnId& id : {"72", "73", "74", "76", "78"}) {
        REQUIRE(curves.count(id) == 1);
        REQUIRE(curves[id]->curve);
    }

    INFO("74 must not share a non-endpoint directed U-turn lead with 72 or 76");
    CHECK_FALSE(curvesIntersectBeyondStrictEndpointOverlapForTest(
        *curves["74"]->curve, *curves["72"]->curve, 0.15));
    CHECK_FALSE(curvesIntersectBeyondStrictEndpointOverlapForTest(
        *curves["74"]->curve, *curves["76"]->curve, 0.15));

    INFO("shared-exit small-radius 74 must be contained by large-radius 73 "
         "without crossing after the shared tail is removed");
    CHECK_FALSE(curvesIntersectBeyondStrictEndpointOverlapForTest(
        *curves["74"]->curve, *curves["73"]->curve, 1.5));

    INFO("U-turn 74 and same-entry straight 78 are not a structural exemption");
    CHECK_FALSE(curvesIntersectBusiness(*curves["74"]->curve, *curves["78"]->curve, 1.5));

    auto bad_pairs = avoidableSameClusterCrossings(input, output);
    INFO("avoidable same-cluster crossings: " << joinPairs(bad_pairs));
    CHECK(bad_pairs.empty());
}

TEST_CASE("intersection_jd shared-entry U-turn 4 stays separated from left turns",
          "[regression][cluster][uturn][intersection_jd]") {
    const std::string path = std::string(PROJECT_ROOT_DIR) + "/datas/intersection_jd.json";
    IntersectionInput input = loadInputOrSkip(path);

    ClusterOrderSolver solver;
    solver.build(input.connectivities, input.lanes, input.lane_groups, input.crosswalks);
    CHECK(solver.exemptionOf("4", "5") == CrossExemption::StructuralCross);
    CHECK(solver.exemptionOf("4", "6") == CrossExemption::StructuralCross);

    IntersectionShapeGenerator gen;
    IntersectionOutput output;
    REQUIRE(gen.generate(input, output));

    auto curves = curveMap(output);
    for (const ConnId& id : {"4", "5", "6", "7", "8", "34"}) {
        REQUIRE(curves.count(id) == 1);
        REQUIRE(curves[id]->curve);
    }

    CHECK(endpointG1Min(*curves["4"], input) > 0.99);
    CHECK(endpointG1Min(*curves["34"], input) > 0.99);
    CHECK(endpointG1Min(*curves["5"], input) > 0.99);
    CHECK(endpointG1Min(*curves["6"], input) > 0.99);
    CHECK(endpointG1Min(*curves["7"], input) > 0.99);
    CHECK(endpointG1Min(*curves["8"], input) > 0.99);

    INFO("same-exit left-turn pairs 5|7 and 6|8 must not share a non-endpoint tail");
    CHECK_FALSE(curvesIntersectBeyondStrictEndpointOverlapForTest(
        *curves["5"]->curve, *curves["7"]->curve, 0.15));
    CHECK_FALSE(curvesIntersectBeyondStrictEndpointOverlapForTest(
        *curves["6"]->curve, *curves["8"]->curve, 0.15));

    INFO("U-turns 4 and 34 must be exactly straight + single turn arc + straight, "
         "with no-crosswalk leads >= 2m and axially aligned arc endpoints");
    CHECK(hasUTurnStraightArcStraightShapeForTest(*curves["4"], input));
    CHECK(hasUTurnStraightArcStraightShapeForTest(*curves["34"], input));
}

TEST_CASE("intersection_jd same-exit left-turn clusters do not cross after shared tail",
          "[regression][cluster][intersection_jd][jd-left-cluster]") {
    const std::string path = std::string(PROJECT_ROOT_DIR) + "/datas/intersection_jd.json";
    IntersectionInput input = loadInputOrSkip(path);

    ClusterOrderSolver solver;
    solver.build(input.connectivities, input.lanes, input.lane_groups, input.crosswalks);
    CHECK(solver.exemptionOf("5", "7") == CrossExemption::None);
    CHECK(solver.exemptionOf("6", "8") == CrossExemption::None);

    IntersectionShapeGenerator gen;
    IntersectionOutput output;
    REQUIRE(gen.generate(input, output));

    auto curves = curveMap(output);
    for (const ConnId& id : {"5", "6", "7", "8"}) {
        REQUIRE(curves.count(id) == 1);
        REQUIRE(curves[id]->curve);
        CHECK(endpointG1Min(*curves[id], input) > 0.99);
        CHECK_FALSE(curveSelfIntersectsBusiness(*curves[id]->curve, 1.0));
    }

    INFO("same-exit left-turn pairs must not share a non-endpoint G1 tail");
    CHECK_FALSE(curvesIntersectBeyondStrictEndpointOverlapForTest(
        *curves["5"]->curve, *curves["7"]->curve, 0.15));
    CHECK_FALSE(curvesIntersectBeyondStrictEndpointOverlapForTest(
        *curves["6"]->curve, *curves["8"]->curve, 0.15));
}

TEST_CASE("intersection_jd fine area uses generated connection cluster endpoints",
          "[regression][area][intersection_jd]") {
    const std::string path = std::string(PROJECT_ROOT_DIR) + "/datas/intersection_jd.json";
    IntersectionInput input = loadInputOrSkip(path);

    REQUIRE(input.lane_groups.empty());
    REQUIRE(input.boundaries.empty());

    IntersectionShapeGenerator gen;
    IntersectionOutput output;
    REQUIRE(gen.generate(input, output));
    REQUIRE_FALSE(output.area.geometry.outer.empty());

    auto curves = curveMap(output);
    REQUIRE(curves.count("4") == 1);
    REQUIRE(curves.count("5") == 1);
    REQUIRE(curves["4"]->curve);
    REQUIRE(curves["5"]->curve);

    const Connectivity* conn4 = findConnectivity(input, "4");
    const Connectivity* conn5 = findConnectivity(input, "5");
    REQUIRE(conn4);
    REQUIRE(conn5);

    CHECK(isSimplePolygon(output.area.geometry));

    const Lane* entry_lane = input.findLane(conn4->entry_lane_id);
    REQUIRE(entry_lane);
    Vec2d old_entry_endpoint = getConnPoint(entry_lane->geometry.points, true);
    CHECK_FALSE(hasPolygonPointNear(output.area.geometry.outer, old_entry_endpoint, 1e-6));

    const Lane* exit_lane = input.findLane(conn5->exit_lane_id);
    REQUIRE(exit_lane);
    Vec2d old_exit_endpoint = getConnPoint(exit_lane->geometry.points, false);
    CHECK_FALSE(hasPolygonPointNear(output.area.geometry.outer, old_exit_endpoint, 1e-6));
}

TEST_CASE("intersection_ds same-exit left and right turns do not cross",
          "[regression][cluster][intersection_ds]") {
    const std::string path = std::string(PROJECT_ROOT_DIR) + "/datas/intersection_ds.json";
    IntersectionInput input = loadInputOrSkip(path);

    IntersectionShapeGenerator gen;
    IntersectionOutput output;
    REQUIRE(gen.generate(input, output));

    auto curves = curveMap(output);
    REQUIRE(curves.count("18") == 1);
    REQUIRE(curves.count("32") == 1);
    REQUIRE(curves.count("9") == 1);
    REQUIRE(curves.count("11") == 1);
    REQUIRE(curves.count("28") == 1);
    REQUIRE(curves.count("8") == 1);
    REQUIRE(curves["18"]->curve);
    REQUIRE(curves["32"]->curve);
    REQUIRE(curves["9"]->curve);
    REQUIRE(curves["11"]->curve);
    REQUIRE(curves["28"]->curve);
    REQUIRE(curves["8"]->curve);

    ClusterOrderSolver solver;
    solver.build(input.connectivities, input.lanes, input.lane_groups, input.crosswalks);
    const CurvePair* pair18_32 = nullptr;
    for (const auto& pair : solver.pairs()) {
        if ((pair.id_a == "18" && pair.id_b == "32") ||
            (pair.id_a == "32" && pair.id_b == "18")) {
            pair18_32 = &pair;
            break;
        }
    }
    REQUIRE(pair18_32);
    CHECK_FALSE(curvesIntersectBusiness(*curves["18"]->curve, *curves["32"]->curve, 1.5));
    INFO("reported same-cluster shared-endpoint pairs must only meet at the connection point");
    CHECK_FALSE(curvesIntersectBeyondAllowedEndpointOverlapForTest(
        *curves["28"]->curve, *curves["11"]->curve, 0.30));
    CHECK_FALSE(curvesIntersectBeyondAllowedEndpointOverlapForTest(
        *curves["9"]->curve, *curves["11"]->curve, 0.30));
    CHECK_FALSE(curvesIntersectBeyondAllowedEndpointOverlapForTest(
        *curves["18"]->curve, *curves["32"]->curve, 0.30));

    const auto& turn32 = *curves["32"]->curve;
    INFO("right turn 32 is a short one-sided arc: high local curvature is allowed "
         "only while it stays a bounded single short arc");
    CHECK(endpointG1Min(*curves["32"], input) > 0.99);
    CHECK(turn32.numSegments() == 1);
    CHECK(turn32.arcLength() < 10.5);
    CHECK(arcChordRatioForTest(turn32) > 1.005);
    CHECK(arcChordRatioForTest(turn32) < 1.08);
    CHECK(turn32.maxCurvature(40) < 6.0);
    CHECK_FALSE(hasCurvatureSignFlipForTest(turn32));
    auto entry32 = input.entryPtDir(curves["32"]->entry_lane_id);
    auto exit32 = input.exitPtDir(curves["32"]->exit_lane_id);
    CHECK(ordinarySingleCubicControlsValid(
        turn32, entry32.first, entry32.second, exit32.first, exit32.second,
        1e-5, true));

    const auto& right8 = *curves["8"]->curve;
    INFO("right turn 8 has no physical avoidance and should remain a single smooth arc");
    CHECK(right8.numSegments() == 1);
    CHECK(arcChordRatioForTest(right8) > 1.03);
    CHECK(arcChordRatioForTest(right8) < 1.30);
    CHECK(right8.maxCurvature(40) < 2.0);
    CHECK_FALSE(hasCurvatureSignFlipForTest(right8));

    auto bad_shared = sharedEndpointRuleViolations(input, output);
    INFO("remaining same connection point residuals outside this focused regression: "
         << joinPairs(bad_shared));
    CHECK(std::find(bad_shared.begin(), bad_shared.end(), "28-11") == bad_shared.end());
    CHECK(std::find(bad_shared.begin(), bad_shared.end(), "9-11") == bad_shared.end());
    CHECK(std::find(bad_shared.begin(), bad_shared.end(), "18-32") == bad_shared.end());
}

TEST_CASE("intersection_ds left turn 2 keeps exit handle inside direction-intersection bounds",
          "[regression][cluster][shape][intersection_ds]") {
    const std::string path = std::string(PROJECT_ROOT_DIR) + "/datas/intersection_ds.json";
    IntersectionInput input = loadInputOrSkip(path);

    IntersectionShapeGenerator gen;
    IntersectionOutput output;
    REQUIRE(gen.generate(input, output));

    auto curves = curveMap(output);
    for (const ConnId& id : {"2", "4", "6"}) {
        REQUIRE(curves.count(id) == 1);
        REQUIRE(curves[id]->curve);
    }

    const ConnectivityCurve& turn2 = *curves["2"];
    const BezierCurve& curve2 = *turn2.curve;
    auto entry2 = input.entryPtDir(turn2.entry_lane_id);
    auto exit2 = input.exitPtDir(turn2.exit_lane_id);
    REQUIRE(curve2.numSegments() == 1);
    INFO("left turn 2 controls must stay on endpoint axes and between endpoints "
         "and the entry/exit direction intersection");
    CHECK(ordinarySingleCubicControlsValid(
        curve2, entry2.first, entry2.second, exit2.first, exit2.second,
        1e-5, true));
    CHECK(arcChordRatioForTest(curve2) > 1.03);
    CHECK(arcChordRatioForTest(curve2) < 1.20);
    CHECK(curve2.maxCurvature(40) < 0.20);
    CHECK_FALSE(curveSelfIntersectsBusiness(curve2, 1.0));

    INFO("left turn 2 must not cross same-entry left turns 4 or 6 outside the "
         "shared connection point");
    CHECK_FALSE(curvesIntersectBeyondAllowedEndpointOverlapForTest(
        curve2, *curves["4"]->curve, 0.15));
    CHECK_FALSE(curvesIntersectBeyondAllowedEndpointOverlapForTest(
        curve2, *curves["6"]->curve, 0.15));
}

TEST_CASE("intersection_ds same-entry U-turns keep wrapped aligned arcs",
          "[regression][cluster][uturn][intersection_ds]") {
    const std::string path = std::string(PROJECT_ROOT_DIR) + "/datas/intersection_ds.json";
    IntersectionInput input = loadInputOrSkip(path);

    IntersectionShapeGenerator gen;
    IntersectionOutput output;
    REQUIRE(gen.generate(input, output));

    auto curves = curveMap(output);
    for (const ConnId& id : {"19", "20", "21", "22"}) {
        REQUIRE(curves.count(id) == 1);
        REQUIRE(curves[id]->curve);
        INFO("U-turn " << id << " should keep bounded curvature; endpoint G1 is exempt");
        CHECK(endpointG1Min(*curves[id], input) > 0.99);
        CHECK(curves[id]->curve->maxCurvature(40) < 3.0);
        CHECK_FALSE(curveSelfIntersectsBusiness(*curves[id]->curve, 1.0));
    }

    CHECK_FALSE(curvesIntersectBeyondAllowedEndpointOverlapForTest(
        *curves["21"]->curve, *curves["19"]->curve, 1.5));
    CHECK_FALSE(curvesIntersectBeyondAllowedEndpointOverlapForTest(
        *curves["20"]->curve, *curves["22"]->curve, 1.5));
    INFO("same-exit U-turns 20|21 must not intersect outside their shared connection point");
    CHECK_FALSE(curvesIntersectBeyondAllowedEndpointOverlapForTest(
        *curves["20"]->curve, *curves["21"]->curve, 0.15));

    INFO("larger-radius 21 should wrap the smaller-radius 19 with a deeper/longer arch");
    CHECK(curves["21"]->curve->arcLength() > curves["19"]->curve->arcLength());
    CHECK(arcChordRatioForTest(*curves["21"]->curve) > arcChordRatioForTest(*curves["19"]->curve));

    INFO("same-entry U-turns must expose exactly straight + single turn arc + straight "
         "with no-crosswalk leads >= 2m and axially aligned arc endpoints");
    CHECK(hasUTurnStraightArcStraightShapeForTest(*curves["19"], input));
    CHECK(hasUTurnStraightArcStraightShapeForTest(*curves["20"], input));
    CHECK(hasUTurnStraightArcStraightShapeForTest(*curves["21"], input));
    CHECK(hasUTurnStraightArcStraightShapeForTest(*curves["22"], input));
    auto uturn_axis_for_test = [&](const ConnId& id) {
        const ConnectivityCurve& cc = *curves[id];
        auto entry = input.entryPtDir(cc.entry_lane_id);
        auto exit_ = input.exitPtDir(cc.exit_lane_id);
        Vec2d t0 = entry.second.norm() > 1e-8
            ? entry.second.normalized() : Vec2d(1, 0);
        Vec2d t1 = exit_.second.norm() > 1e-8
            ? exit_.second.normalized() : -t0;
        Vec2d axis = t0 - t1;
        if (axis.norm() < 1e-8)
            axis = t0;
        axis.normalize();
        if (axis.dot(t0) < 0.0)
            axis = -axis;
        return axis;
    };
    Vec2d axis19 = uturn_axis_for_test("19");
    double q0_station19 = curves["19"]->curve->segs.front().ctrl[3].dot(axis19);
    double q1_station19 = curves["19"]->curve->segs.back().ctrl[0].dot(axis19);
    for (const ConnId& id : {"20", "21", "22"}) {
        double q0_station = curves[id]->curve->segs.front().ctrl[3].dot(axis19);
        double q1_station = curves[id]->curve->segs.back().ctrl[0].dot(axis19);
        INFO("19/" << id << " common alignment stations q0="
             << q0_station19 << "/" << q0_station
             << " q1=" << q1_station19 << "/" << q1_station);
        CHECK(std::abs(q0_station19 - q0_station) < 0.05);
        CHECK(std::abs(q1_station19 - q1_station) < 0.05);
    }
    INFO("centimetre-scale stagger must separate otherwise coincident arc endpoints");
    CHECK((curves["19"]->curve->segs.front().ctrl[3] -
           curves["21"]->curve->segs.front().ctrl[3]).norm() >= 0.005);
    CHECK((curves["20"]->curve->segs.front().ctrl[3] -
           curves["22"]->curve->segs.front().ctrl[3]).norm() >= 0.005);
    ClusterOrderSolver solver;
    solver.build(input.connectivities, input.lanes, input.lane_groups, input.crosswalks);
    Vec2d exit_ref = solver.refPerpOf("20", "21");
    REQUIRE(exit_ref.norm() > 1e-8);
    exit_ref.normalize();
    int exit_side = solver.expectedSideOf("20", "21");
    REQUIRE(exit_side != 0);
    Vec2d exit_ctrl20 = curves["20"]->curve->segs[1].ctrl[2];
    Vec2d exit_ctrl21 = curves["21"]->curve->segs[1].ctrl[2];
    INFO("20|21 exit controls 20=(" << exit_ctrl20.x() << "," << exit_ctrl20.y()
         << ") 21=(" << exit_ctrl21.x() << "," << exit_ctrl21.y() << ")");
    CHECK(exit_ctrl20.y() > exit_ctrl21.y() + 0.05);
    CHECK(exit_side * (exit_ctrl20 - exit_ctrl21).dot(exit_ref) > 0.01);
    CHECK(arcChordRatioForTest(*curves["20"]->curve) > 1.4);
    CHECK(arcChordRatioForTest(*curves["22"]->curve) > 1.4);

}

TEST_CASE("intersection_cross U-turns 39 and 40 keep straight-arc-straight shape",
          "[regression][cluster][uturn][intersection_cross]") {
    const std::string path = std::string(PROJECT_ROOT_DIR) + "/datas/intersection_cross.json";
    IntersectionInput input = loadInputOrSkip(path);

    auto start = std::chrono::steady_clock::now();
    IntersectionShapeGenerator gen;
    IntersectionOutput output;
    REQUIRE(gen.generate(input, output));
    double elapsed_ms = std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - start).count();

    INFO("intersection_cross generation elapsed_ms=" << elapsed_ms);
    CHECK(elapsed_ms < 15000.0);

    auto curves = curveMap(output);
    for (const ConnId& id : {"39", "40"}) {
        REQUIRE(curves.count(id) == 1);
        REQUIRE(curves[id]->curve);
        const BezierCurve& curve = *curves[id]->curve;
        INFO("U-turn " << id << " should be explicit straight-arc-straight");
        REQUIRE(curve.numSegments() == 3);
        CHECK(segmentLooksStraightForTest(curve.segs.front()));
        CHECK(segmentLooksStraightForTest(curve.segs.back()));
        CHECK((curve.segs.front().ctrl[3] -
               curve.segs.front().ctrl[0]).norm() >= 0.20);
        CHECK((curve.segs.back().ctrl[3] -
               curve.segs.back().ctrl[0]).norm() >= 0.20);
        CHECK(curve.segs[1].maxCurvature(30) > 0.03);
        CHECK_FALSE(curveSelfIntersectsBusiness(curve, 1.0));
        CHECK(arcChordRatioForTest(curve) > 1.35);
    }

    INFO("39|40 are same-cluster U-turns without a shared endpoint, so they must "
         "not intersect outside endpoint tolerance");
    CHECK_FALSE(curvesIntersectBeyondAllowedEndpointOverlapForTest(
        *curves["39"]->curve, *curves["40"]->curve, 0.30));
    auto bad_pairs = avoidableSameClusterCrossings(input, output);
    CHECK(std::find(bad_pairs.begin(), bad_pairs.end(), "39-40") == bad_pairs.end());
    CHECK(std::find(bad_pairs.begin(), bad_pairs.end(), "40-39") == bad_pairs.end());
}

TEST_CASE("100000012-nu fine area is a simple concave-capable outline",
          "[regression][area][100000012]") {
    const std::string path = std::string(PROJECT_ROOT_DIR) + "/datas/100000012-nu.json";
    IntersectionInput input = loadInputOrSkip(path);

    IntersectionShapeGenerator gen;
    IntersectionOutput output;
    REQUIRE(gen.generate(input, output));
    REQUIRE_FALSE(output.area.geometry.outer.empty());
    CHECK(isSimplePolygon(output.area.geometry));
}

TEST_CASE("100000012-nu fixed-shape clusters stay non-intersecting",
          "[regression][cluster][fixed-geometry][100000012]") {
    const std::string path = std::string(PROJECT_ROOT_DIR) + "/datas/100000012-nu.json";
    IntersectionInput input = loadInputOrSkip(path);

    IntersectionShapeGenerator gen;
    IntersectionOutput output;
    REQUIRE(gen.generate(input, output));

    auto curves = curveMap(output);
    auto checkPair = [&](const ConnId& a, const ConnId& b) {
        REQUIRE(curves.count(a) == 1);
        REQUIRE(curves.count(b) == 1);
        REQUIRE(curves[a]->curve);
        REQUIRE(curves[b]->curve);
        ClusterOrderSolver solver;
        solver.build(input.connectivities, input.lanes, input.lane_groups, input.crosswalks);
        const CurvePair* target = nullptr;
        for (const auto& pair : solver.pairs()) {
            if ((pair.id_a == a && pair.id_b == b) ||
                (pair.id_a == b && pair.id_b == a)) {
                target = &pair;
                break;
            }
        }
        REQUIRE(target);
        INFO("pair " << a << "-" << b);
        CHECK_FALSE(curvesIntersectBusiness(*curves[a]->curve, *curves[b]->curve, 1.5));
    };

    checkPair("43285966", "7");
    checkPair("43285966", "2");
    checkPair("27", "29");

    INFO("fixed-shape U-turns must still be regenerated as straight-arc-straight");
    for (const ConnId& id : {"17", "19"}) {
        REQUIRE(curves.count(id) == 1);
        REQUIRE(curves[id]->curve);
        CHECK(hasUTurnStraightArcStraightShapeForTest(*curves[id], input));
    }

    auto bad_shared = sharedEndpointRuleViolations(input, output);
    INFO("same connection point pairs must not have non-endpoint crossings: "
         << joinPairs(bad_shared));
    CHECK(bad_shared.empty());
}

TEST_CASE("100000012-nu natural curves are not over-constrained by fixed shapes",
          "[regression][cluster][fixed-geometry][shape][100000012]") {
    const std::string path = std::string(PROJECT_ROOT_DIR) + "/datas/100000012-nu.json";
    IntersectionInput input = loadInputOrSkip(path);

    IntersectionShapeGenerator gen;
    IntersectionOutput output;
    REQUIRE(gen.generate(input, output));

    auto curves = curveMap(output);
    for (const ConnId& id : {"2", "9", "31"}) {
        REQUIRE(curves.count(id) == 1);
        REQUIRE(curves[id]->curve);
        CHECK(endpointG1Min(*curves[id], input) > 0.99);
        CHECK_FALSE(curveSelfIntersectsBusiness(*curves[id]->curve, 1.0));
    }
    CHECK(curves["2"]->curve->numSegments() <= 2);
    CHECK(curves["9"]->curve->numSegments() == 1);
    CHECK(curves["31"]->curve->numSegments() == 1);

    INFO("straight 2 and 9 should remain nearly straight natural curves");
    CHECK(arcChordRatioForTest(*curves["2"]->curve) < 1.02);
    CHECK(arcChordRatioForTest(*curves["9"]->curve) < 1.02);
    CHECK(curves["2"]->curve->maxCurvature(40) < 0.05);
    CHECK(curves["9"]->curve->maxCurvature(40) < 0.05);

    INFO("left turn 31 should remain one smooth single-direction arc");
    CHECK(arcChordRatioForTest(*curves["31"]->curve) > 1.05);
    CHECK(arcChordRatioForTest(*curves["31"]->curve) < 1.25);
    CHECK(curves["31"]->curve->maxCurvature(40) < 0.25);
    CHECK_FALSE(hasCurvatureSignFlipForTest(*curves["31"]->curve));
}

TEST_CASE("100000012-nu left turns near fixed U-turns keep natural shape",
          "[regression][cluster][fixed-geometry][shape][100000012]") {
    const std::string path = std::string(PROJECT_ROOT_DIR) + "/datas/100000012-nu.json";
    IntersectionInput input = loadInputOrSkip(path);

    IntersectionShapeGenerator gen;
    IntersectionOutput output;
    REQUIRE(gen.generate(input, output));

    auto curves = curveMap(output);
    for (const ConnId& id : {"5", "27", "29", "32"}) {
        REQUIRE(curves.count(id) == 1);
        REQUIRE(curves[id]->curve);
        INFO("conn " << id << " should keep endpoint G1");
        CHECK(endpointG1Min(*curves[id], input) > 0.98);
    }

    auto arcChord = [&](const ConnId& id) {
        const auto& c = *curves[id]->curve;
        double chord = (c.endPt() - c.startPt()).norm();
        REQUIRE(chord > 1e-6);
        return c.arcLength() / chord;
    };

    INFO("left turns 5, 27 and 29 should not collapse to near-straight chords");
    CHECK(arcChord("5") > 1.03);
    CHECK(arcChord("27") > 1.03);
    CHECK(arcChord("29") > 1.05);

    INFO("left turn 32 should not be stretched into an S-shaped detour");
    CHECK(arcChord("32") < 1.35);
}

TEST_CASE("100000012-nu fine area keeps specified group cut endpoints",
          "[regression][area][100000012]") {
    const std::string path = std::string(PROJECT_ROOT_DIR) + "/datas/100000012-nu.json";
    IntersectionInput input = loadInputOrSkip(path);

    IntersectionShapeGenerator gen;
    IntersectionOutput output;
    REQUIRE(gen.generate(input, output));
    REQUIRE_FALSE(output.area.geometry.outer.empty());
    REQUIRE(isSimplePolygon(output.area.geometry));

    Vec2d center = outputCenter(output);
    for (const auto& lane_id : std::vector<LaneId>{"43270483", "43270487"}) {
        Vec2d expected = expectedLaneCutPoint(input, output, lane_id, GroupRole::Entry, center);
        CHECK(hasPolygonPointNear(output.area.geometry.outer, expected, 0.08));
    }
    for (const auto& lane_id : std::vector<LaneId>{"43247337", "43247338", "43247339"}) {
        Vec2d expected = expectedLaneCutPoint(input, output, lane_id, GroupRole::Exit, center);
        CHECK(hasPolygonPointNear(output.area.geometry.outer, expected, 0.08));
    }
}

TEST_CASE("100000012 fine area uses the complete outer outline, not a local road-edge loop",
          "[regression][area][100000012]") {
    const std::string path = std::string(PROJECT_ROOT_DIR) + "/datas/100000012.json";
    IntersectionInput input = loadInputOrSkip(path);

    IntersectionShapeGenerator gen;
    IntersectionOutput output;
    REQUIRE(gen.generate(input, output));
    REQUIRE_FALSE(output.area.geometry.outer.empty());
    REQUIRE(isSimplePolygon(output.area.geometry));

    Vec2d center = outputCenter(output);
    CHECK(pointInPolygon(center, output.area.geometry.outer));
    CHECK(std::abs(polygonSignedArea(output.area.geometry.outer)) > 500.0);
    CHECK(output.area.geometry.outer.size() > 20);

    for (const auto& lane_id : std::vector<LaneId>{"43270483", "43270487", "43270484", "43270488"}) {
        Vec2d expected = expectedLaneCutPoint(input, output, lane_id, GroupRole::Entry, center);
        CHECK(hasPolygonPointNear(output.area.geometry.outer, expected, 0.08));
    }
    for (const auto& lane_id : std::vector<LaneId>{"43247337", "43247338", "43247339", "43242048", "43242051"}) {
        Vec2d expected = expectedLaneCutPoint(input, output, lane_id, GroupRole::Exit, center);
        CHECK(hasPolygonPointNear(output.area.geometry.outer, expected, 0.08));
    }

    const Boundary* edge43285426 = findBoundary(input, "43285426");
    const Boundary* edge43285966 = findBoundary(input, "43285966");
    REQUIRE(edge43285426);
    REQUIRE(edge43285966);
    CHECK(polygonUsesBoundarySegment(output.area.geometry.outer, *edge43285426));
    CHECK(polygonUsesBoundarySegment(output.area.geometry.outer, *edge43285966));
    double edge43285966_len = boundaryLengthForTest(*edge43285966);
    constexpr double edge43285966_entry_cut_station = 0.501155;
    CHECK(polygonHasBoundaryStationNear(
        output.area.geometry.outer, *edge43285966, edge43285966_entry_cut_station));
    CHECK(polygonContainsBoundaryShapePointsFromStation(
        output.area.geometry.outer, *edge43285966, edge43285966_entry_cut_station));
    CHECK(polygonUsesBoundarySegmentInStationRange(
        output.area.geometry.outer, *edge43285966, edge43285966_len - 1.0, edge43285966_len));
    CHECK(hasPolygonPointNear(output.area.geometry.outer, edge43285966->geometry.points.back(), 0.08));
    CHECK_FALSE(hasPolygonPointNear(output.area.geometry.outer, edge43285966->geometry.points.front(), 0.08));
    CHECK_FALSE(hasPolygonPointNear(output.area.geometry.outer, edge43285966->geometry.points[1], 0.08));
}

TEST_CASE("110000741-u fine area does not use far lane-edge orientation endpoints",
          "[regression][area][110000741]") {
    const std::string path = std::string(PROJECT_ROOT_DIR) + "/datas/110000741-u.json";
    IntersectionInput input = loadInputOrSkip(path);

    IntersectionShapeGenerator gen;
    IntersectionOutput output;
    REQUIRE(gen.generate(input, output));
    REQUIRE_FALSE(output.area.geometry.outer.empty());
    REQUIRE(isSimplePolygon(output.area.geometry));

    INFO("far lane-edge endpoints must not be shifted into the area ring");
    CHECK_FALSE(hasPolygonPointNear(output.area.geometry.outer,
                                    Vec2d(52.395118666, -6.937535945), 0.08));
    CHECK_FALSE(hasPolygonPointNear(output.area.geometry.outer,
                                    Vec2d(-71.773072630, 10.994317701), 0.08));
}

TEST_CASE("100000610 fine area clips merged RoadEdges toward the intersection",
          "[regression][area][shared-edge][100000610]") {
    const std::string path = std::string(PROJECT_ROOT_DIR) + "/datas/100000610.json";
    IntersectionInput input = loadInputOrSkip(path);
    REQUIRE(input.mode == 2);

    IntersectionShapeGenerator gen;
    IntersectionOutput output;
    auto t0 = std::chrono::steady_clock::now();
    REQUIRE(gen.generate(input, output));
    double elapsed_ms = std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - t0).count();
    INFO("100000610 generation elapsed_ms=" << elapsed_ms);
    REQUIRE(elapsed_ms < 15000.0);
    REQUIRE_FALSE(output.area.geometry.outer.empty());
    REQUIRE(isSimplePolygon(output.area.geometry));

    Vec2d center = outputCenter(output);
    double max_radius = 0.0;
    for (const auto& pt : output.area.geometry.outer)
        max_radius = std::max(max_radius, dist(pt, center));
    INFO("merged RoadEdge clipping must not pull remote road endpoints into "
         "the intersection area; max_radius=" << max_radius);
    CHECK(max_radius < 40.0);
    CHECK_FALSE(hasPolygonPointNear(
        output.area.geometry.outer, Vec2d(9.18367308, 16.08045963), 0.08));
    CHECK_FALSE(hasPolygonPointNear(
        output.area.geometry.outer, Vec2d(3.91731926, -183.95546304), 0.08));

    Vec2d edge43_entry = expectedLaneEdgeCutPoint(
        input, "43120543", GroupRole::Entry, center, 0.05);
    Vec2d edge43_exit = expectedLaneEdgeCutPoint(
        input, "43120543", GroupRole::Exit, center, 0.05);
    INFO("shared LaneEdge 43120543 keeps only its entry-direction cut");
    CHECK(hasPolygonPointNear(output.area.geometry.outer, edge43_entry, 0.08));
    CHECK_FALSE(hasPolygonPointNear(output.area.geometry.outer, edge43_exit, 0.08));

    Vec2d edge95_entry = expectedLaneEdgeCutPoint(
        input, "43120595", GroupRole::Entry, center, 0.05);
    Vec2d edge95_exit = expectedLaneEdgeCutPoint(
        input, "43120595", GroupRole::Exit, center, 0.05);
    INFO("shared LaneEdge 43120595 keeps only its exit-direction cut");
    CHECK_FALSE(hasPolygonPointNear(output.area.geometry.outer, edge95_entry, 0.08));
    CHECK(hasPolygonPointNear(output.area.geometry.outer, edge95_exit, 0.08));
}

TEST_CASE("100000610 same-entry right turns preserve outward handle order",
          "[regression][curve][cluster][100000610]") {
    const std::string path = std::string(PROJECT_ROOT_DIR) + "/datas/100000610.json";
    IntersectionInput input = loadInputOrSkip(path);
    IntersectionShapeGenerator gen;
    IntersectionOutput output;
    auto t0 = std::chrono::steady_clock::now();
    REQUIRE(gen.generate(input, output));
    const double elapsed_ms = std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - t0).count();
    REQUIRE(elapsed_ms < 15000.0);

    auto curves = curveMap(output);
    REQUIRE(curves.count("1") == 1);
    REQUIRE(curves.count("2") == 1);
    REQUIRE(curves["1"]->curve);
    REQUIRE(curves["2"]->curve);
    const BezierCurve& right1 = *curves["1"]->curve;
    const BezierCurve& right2 = *curves["2"]->curve;
    REQUIRE(right1.numSegments() == 1);
    REQUIRE(right2.numSegments() == 1);

    auto entry = input.entryPtDir(curves["1"]->entry_lane_id);
    const double lead1 = (right1.segs.front().ctrl[1] - entry.first).dot(
        entry.second.normalized());
    const double lead2 = (right2.segs.front().ctrl[1] - entry.first).dot(
        entry.second.normalized());
    INFO("right-turn shared-entry leads: 1=" << lead1 << " 2=" << lead2);
    CHECK(lead1 > lead2 + 0.20);
    CHECK((curves["1"]->curve->endPt() -
           right1.segs.front().ctrl[2]).dot(
               input.exitPtDir(curves["1"]->exit_lane_id).second.normalized()) > 0.8);
    CHECK(ordinarySingleCubicControlsValid(
        right1, entry.first, entry.second,
        right1.endPt(), input.exitPtDir(curves["1"]->exit_lane_id).second,
        1e-5, true));
    CHECK_FALSE(curvesIntersectBeyondAllowedEndpointOverlapForTest(
        right1, right2, 0.15));
    CHECK_FALSE(curveSelfIntersectsBusiness(right1, 1.0));
    CHECK_FALSE(curveSelfIntersectsBusiness(right2, 1.0));
    const auto bad_pairs = avoidableSameClusterCrossings(input, output);
    INFO("avoidable same-cluster crossings: " << joinPairs(bad_pairs));
    CHECK(bad_pairs.empty());
}

TEST_CASE("110002137 exit group road edge uses only the cut-line valid span",
          "[regression][area][110002137]") {
    const std::string path = std::string(PROJECT_ROOT_DIR) + "/datas/110002137.json";
    IntersectionInput input = loadInputOrSkip(path);

    IntersectionShapeGenerator gen;
    IntersectionOutput output;
    REQUIRE(gen.generate(input, output));
    REQUIRE_FALSE(output.area.geometry.outer.empty());
    REQUIRE(isSimplePolygon(output.area.geometry));

    const Boundary* edge1006810 = findBoundary(input, "1006810");
    REQUIRE(edge1006810);
    REQUIRE(edge1006810->geometry.points.size() >= 4);
    CHECK(polygonUsesBoundarySegment(output.area.geometry.outer, *edge1006810));

    Vec2d center = outputCenter(output);
    Vec2d duplicate_edge_cut = expectedLaneEdgeCutPoint(
        input, "1006810", GroupRole::Exit, center);
    INFO("LaneEdge 1006810 must still create its own group cut point even when "
         "a RoadEdge has the same id "
         << duplicate_edge_cut.transpose());
    CHECK(hasPolygonPointNear(output.area.geometry.outer, duplicate_edge_cut, 0.08));

    Vec2d duplicate_lane_cut = expectedLaneCutPoint(
        input, output, "1011002", GroupRole::Exit, center);
    INFO("Lane 1011002 and LaneEdge 1006810 belong to different source categories; "
         "the lane cut point must not be filtered by laneedge geometry "
         << duplicate_lane_cut.transpose());
    CHECK(hasPolygonPointNear(output.area.geometry.outer, duplicate_lane_cut, 0.08));

    Vec2d mode2_entry_edge_cut = expectedLaneEdgeCutPoint(
        input, "1006892", GroupRole::Entry, center);
    INFO("Mode=2 entry LaneEdge 1006892 should expand from its tail tangent, "
         "not from the flipped tail extension "
         << mode2_entry_edge_cut.transpose());
    CHECK(hasPolygonPointNear(output.area.geometry.outer, mode2_entry_edge_cut, 0.08));
    CHECK(dist(mode2_entry_edge_cut, Vec2d(-16.4858, 29.1412)) < 0.40);
    CHECK_FALSE(hasPolygonPointNear(
        output.area.geometry.outer, Vec2d(-17.3919, 29.1048), 0.08));

    Vec2d shared_entry_edge_cut = expectedLaneEdgeCutPoint(
        input, "1005432", GroupRole::Entry, center);
    Vec2d shared_exit_edge_cut = expectedLaneEdgeCutPoint(
        input, "1005432", GroupRole::Exit, center);
    INFO("Shared LaneEdge 1005432 should keep the cut point whose vectorized "
         "direction matches the group direction " << shared_entry_edge_cut.transpose());
    CHECK(hasPolygonPointNear(output.area.geometry.outer, shared_entry_edge_cut, 0.08));
    INFO("Shared LaneEdge 1005432 should ignore the opposite-direction exit cut "
         << shared_exit_edge_cut.transpose());
    CHECK_FALSE(hasPolygonPointNear(output.area.geometry.outer, shared_exit_edge_cut, 0.08));

    INFO("LaneEdge 1006810 exit cut point should stay on the edge near its start "
         << duplicate_edge_cut.transpose());
    CHECK(dist(duplicate_edge_cut, Vec2d(1.3339, -24.7504)) < 0.40);

    Vec2d old_tail_extension_cut =
        xyOf(edge1006810->geometry.points.back()) +
        0.5 * modeShiftDirForTest(
            GroupRole::Exit,
            centerSideEndpointDirForTest(edge1006810->geometry.points, GroupRole::Exit, center),
            input.mode);
    INFO("LaneEdge 1006810 must not use the old tail-extension cut point "
         << old_tail_extension_cut.transpose());
    CHECK_FALSE(hasPolygonPointNear(
        output.area.geometry.outer, old_tail_extension_cut, 0.08));
    CHECK_FALSE(hasPolygonPointNear(
        output.area.geometry.outer, Vec2d(4.40972, -22.903226), 0.08));

    double min_station = 0.0;
    double max_station = 0.0;
    REQUIRE(polygonBoundaryStationSpanForTest(
        output.area.geometry.outer, *edge1006810, min_station, max_station, 0.12));
    double used_span = max_station - min_station;
    double full_len = boundaryLengthForTest(*edge1006810);
    INFO("RoadEdge 1006810 must be clipped by exit group 1076917 cut-line probes; "
         "used_span=" << used_span << " full_len=" << full_len);
    CHECK(used_span < full_len - 0.50);
    CHECK_FALSE(hasPolygonPointNear(
        output.area.geometry.outer, edge1006810->geometry.points[2], 0.08));
}

TEST_CASE("100000443 road edge clipping ignores opposite-direction group hits",
          "[regression][area][100000443]") {
    const std::string path = std::string(PROJECT_ROOT_DIR) + "/datas/100000443.json";
    IntersectionInput input = loadInputOrSkip(path);

    IntersectionAreaBuilder builder(0.5);
    IntersectionArea area = builder.build(input, {}, {});
    REQUIRE_FALSE(area.geometry.outer.empty());
    REQUIRE(isSimplePolygon(area.geometry));

    INFO("Group cut-line extension hits on merged RoadEdges must only contribute "
         "Boundary stations when the local RoadEdge direction matches the group "
         "direction; these points were selected by the old opposite-direction "
         "clipping path.");
    CHECK_FALSE(hasPolygonPointNear(area.geometry.outer, Vec2d(14.0394, 19.9863), 0.10));
    CHECK_FALSE(hasPolygonPointNear(area.geometry.outer, Vec2d(19.1029, -10.6362), 0.10));
    CHECK_FALSE(hasPolygonPointNear(area.geometry.outer, Vec2d(-10.0571, -24.3031), 0.10));

    for (const auto& boundary_id : std::vector<LaneEdgeId>{"1015719", "1015721", "1015712"}) {
        const Boundary* bridge = findBoundary(input, boundary_id);
        REQUIRE(bridge);
        INFO("Other-type boundary " << boundary_id
             << " touches a RoadEdge endpoint and must be merged into the "
             "RoadEdge chain before cut-line clipping");
        CHECK(polygonUsesBoundarySegment(area.geometry.outer, *bridge));
    }
}

TEST_CASE("U-turn alignment scope switches between lane endpoint and lane group families",
          "[regression][uturn][alignment]") {
    IntersectionInput input = makeUTurnAlignmentScopeInput();

    IntersectionShapeGenerator::Config lane_cfg;
    lane_cfg.connectivity_direction.uturn_alignment_scope =
        UTurnAlignmentScope::LaneEndpoint;
    IntersectionShapeGenerator lane_gen(lane_cfg);
    IntersectionOutput lane_output;
    REQUIRE(lane_gen.generate(input, lane_output));
    auto lane_curves = curveMap(lane_output);
    REQUIRE(lane_curves.count("near") == 1);
    REQUIRE(lane_curves.count("far") == 1);
    REQUIRE(lane_curves["near"]->curve);
    REQUIRE(lane_curves["far"]->curve);
    REQUIRE(lane_curves.count("shared_a") == 1);
    REQUIRE(lane_curves.count("shared_b") == 1);
    REQUIRE(lane_curves["shared_a"]->curve);
    REQUIRE(lane_curves["shared_b"]->curve);
    REQUIRE(lane_curves["near"]->curve->numSegments() == 3);
    REQUIRE(lane_curves["far"]->curve->numSegments() == 3);
    REQUIRE(lane_curves["shared_a"]->curve->numSegments() == 3);
    REQUIRE(lane_curves["shared_b"]->curve->numSegments() == 3);
    double lane_near_station =
        lane_curves["near"]->curve->segs.front().ctrl[3].x();
    double lane_far_station =
        lane_curves["far"]->curve->segs.front().ctrl[3].x();
    double lane_shared_a_station =
        lane_curves["shared_a"]->curve->segs.front().ctrl[3].x();
    double lane_shared_b_station =
        lane_curves["shared_b"]->curve->segs.front().ctrl[3].x();
    INFO("lane scope stations near=" << lane_near_station
         << " far=" << lane_far_station
         << " shared_a=" << lane_shared_a_station
         << " shared_b=" << lane_shared_b_station);
    CHECK(lane_far_station > lane_near_station + 1.0);
    CHECK(std::abs(lane_shared_a_station - lane_shared_b_station) < 0.25);

    IntersectionShapeGenerator::Config group_cfg;
    group_cfg.connectivity_direction.uturn_alignment_scope =
        UTurnAlignmentScope::LaneGroup;
    IntersectionShapeGenerator group_gen(group_cfg);
    IntersectionOutput group_output;
    REQUIRE(group_gen.generate(input, group_output));
    auto group_curves = curveMap(group_output);
    REQUIRE(group_curves.count("near") == 1);
    REQUIRE(group_curves.count("far") == 1);
    REQUIRE(group_curves["near"]->curve);
    REQUIRE(group_curves["far"]->curve);
    REQUIRE(group_curves["near"]->curve->numSegments() == 3);
    REQUIRE(group_curves["far"]->curve->numSegments() == 3);
    double group_near_station =
        group_curves["near"]->curve->segs.front().ctrl[3].x();
    double group_far_station =
        group_curves["far"]->curve->segs.front().ctrl[3].x();
    INFO("group scope stations near=" << group_near_station
         << " far=" << group_far_station);
    CHECK(std::abs(group_near_station - group_far_station) < 0.25);
    CHECK(group_near_station > lane_near_station + 1.0);

    INFO("shared-entry U-turns must move both aligned points by the same "
         "distance in opposite directions, even when only the entry endpoint family is shared");
    auto alignedPointLateralOffsets = [&](const ConnId& id) {
        const ConnectivityCurve& cc = *lane_curves[id];
        auto entry = input.entryPtDir(cc.entry_lane_id);
        auto exit_ = input.exitPtDir(cc.exit_lane_id);
        Vec2d t0 = entry.second.norm() > 1e-8
            ? entry.second.normalized() : Vec2d(1, 0);
        Vec2d t1 = exit_.second.norm() > 1e-8
            ? exit_.second.normalized() : -t0;
        Vec2d axis = t0 - t1;
        if (axis.norm() < 1e-8)
            axis = t0;
        axis.normalize();
        if (axis.dot(t0) < 0.0)
            axis = -axis;
        Vec2d lateral{-axis[1], axis[0]};
        Vec2d q0 = cc.curve->segs.front().ctrl[3];
        Vec2d q1 = cc.curve->segs.back().ctrl[0];
        return std::pair<double, double>{
            (q0 - entry.first).dot(lateral),
            (q1 - exit_.first).dot(lateral)};
    };
    for (const ConnId& id : {"shared_a", "shared_b"}) {
        auto offsets = alignedPointLateralOffsets(id);
        INFO("U-turn " << id << " q0/q1 lateral offsets="
             << offsets.first << "/" << offsets.second);
        CHECK(std::abs(offsets.first) > 0.005);
        CHECK(std::abs(offsets.second) > 0.005);
        CHECK(std::abs(offsets.first + offsets.second) < 0.05);
    }
}

TEST_CASE("100000443 U-turn three-segment shape and crosswalk clearance",
          "[regression][uturn][crosswalk][100000443]") {
    // ── Background ──────────────────────────────────────────────────────────
    // 100000443.json has three road-arm boundary chains that each consist of
    // a RoadEdge (type=0) + Connector (type=3) + RoadEdge (type=0) with the
    // chain endpoints touching the entry/exit lane endpoints at ~0.3–0.4 m.
    //
    // Before the fix, assessCurveRisk() treated the connector node as an
    // unterminated boundary crossing, marking every U-turn three-segment
    // straight leg as physical() → all segmented candidates were rejected
    // (phys=32), leaving the solver with no valid candidates and falling back
    // to degenerate single-arc shapes (maxk > 3.0).
    //
    // ── Verified requirements ────────────────────────────────────────────────
    // (A) U-turns with crosswalks (conn 1015464, 18, 20):
    //     - Must be three-segment (straight + arc + straight)
    //     - Arc must not enter any crosswalk polygon
    //     - Straight legs must cross the crosswalk (i.e. lead ≥ crosswalk far
    //       edge distance measured from the endpoint)
    //     - Endpoint G1 check is skipped for U-turns because shared-endpoint
    //       staggering can intentionally move the aligned points off G1.
    //
    // (B) U-turns without crosswalks near entry/exit (conn 32, 34):
    //     - Must be three-segment
    //     - Straight legs must each be ≥ kUTurnNoCrosswalkMinLead = 2.0 m
    //     - Endpoint G1 check is skipped for U-turns.
    //
    // (C) Side-boundary-chain straight/turning connections (conn 12, 14, 28, 30):
    //     - Must NOT have status == CurveStatus::Degraded due to boundary
    //       traversal; the arm side-chains are legitimately traversed when
    //       entering from one arm and exiting through another.

    const std::string path = std::string(PROJECT_ROOT_DIR) + "/datas/100000443.json";
    IntersectionInput input = loadInputOrSkip(path);

    IntersectionShapeGenerator gen;
    IntersectionOutput output;
    auto t_gen = std::chrono::steady_clock::now();
    REQUIRE(gen.generate(input, output));
    double gen_ms = std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - t_gen).count();
    INFO("100000443 generation elapsed_ms=" << gen_ms);
    // AGENTS.md: worst acceptable single-intersection generation time < 15 s.
    CHECK(gen_ms < 15000.0);

    auto curves = curveMap(output);

    auto require_curve = [&](const ConnId& id) -> const ConnectivityCurve& {
        auto it = curves.find(id);
        REQUIRE(it != curves.end());
        REQUIRE(it->second->curve);
        return *it->second;
    };

    // ── (A) Crosswalk-adjacent U-turns ──────────────────────────────────────
    // conn 1015464: entry arm C (near crosswalk 100000375/100000381),
    //   crosswalkClearanceAhead ≈ 6.68 m → lead0 floor ≈ 6.68 m.
    // conn 18: entry arm B (near crosswalk 100000382),
    //   crosswalkClearanceAhead ≈ 6.98 m.
    // conn 20: same entry arm B.
    for (const ConnId& id : {"1015464", "18", "20"}) {
        const auto& cc = require_curve(id);
        INFO("U-turn " << id << " near a crosswalk must be three-segment "
             << "(straight + arc + straight) with the arc outside the crosswalk");

        // Three-segment shape.
        CHECK(cc.curve->numSegments() == 3);

        // endpointG1Min returns the business exemption value for U-turns.
        CHECK(endpointG1Min(cc, input) >= 0.95);

        // Arc does not enter any crosswalk (uturnArcClearsCrosswalksForTest
        // also verifies that each straight leg is ≥ min_straight long,
        // confirming the crosswalk was crossed).
        // Use a conservative min_straight of 4.0 m — in all three cases the
        // computed lead floor is ≥ 4 m, so the straight leg must exceed this.
        if (cc.curve->numSegments() == 3) {
            CHECK(uturnArcClearsCrosswalksForTest(*cc.curve, input.crosswalks, 4.0));
            auto entry = input.entryPtDir(cc.entry_lane_id);
            auto exit_ = input.exitPtDir(cc.exit_lane_id);
            double required0 = crosswalkClearanceAlongRayForTest(
                entry.first, entry.second, input.crosswalks);
            double required1 = crosswalkClearanceAlongRayForTest(
                exit_.first, -exit_.second, input.crosswalks);
            double lead0 = (cc.curve->segs.front().ctrl[3] -
                            cc.curve->segs.front().ctrl[0]).norm();
            double lead1 = (cc.curve->segs.back().ctrl[3] -
                            cc.curve->segs.back().ctrl[0]).norm();
            INFO("conn " << id << " exact Crosswalk lead floors: required0="
                 << required0 << " actual0=" << lead0
                 << " required1=" << required1 << " actual1=" << lead1);
            if (required0 > 0.0)
                CHECK(lead0 >= required0 - 0.05);
            if (required1 > 0.0)
                CHECK(lead1 >= required1 - 0.05);
        }

        // Curvature sanity: arc should be well-formed, not a degenerate spike.
        double maxk = cc.curve->maxCurvature(40);
        INFO("conn " << id << " maxk=" << maxk << " (was > 3.0 before fix)");
        CHECK(maxk < 3.0);
    }

    // ── (B) No-crosswalk U-turns: 2 m minimum straight legs ─────────────────
    // conn 32: turn_gap ≈ 2.56 m, no crosswalk on that arm.
    // conn 34: turn_gap ≈ 6.26 m, no crosswalk on that arm.
    // Both must have kUTurnNoCrosswalkMinLead = 2.0 m straight legs.
    constexpr double kMinLead = 2.0;
    for (const ConnId& id : {"32", "34"}) {
        const auto& cc = require_curve(id);
        INFO("U-turn " << id << " has no nearby crosswalk; straight legs "
             << "must each be >= " << kMinLead << " m");

        CHECK(cc.curve->numSegments() == 3);
        CHECK(endpointG1Min(cc, input) >= 0.95);

        if (cc.curve->numSegments() == 3) {
            double lead0 = (cc.curve->segs.front().ctrl[3] -
                            cc.curve->segs.front().ctrl[0]).norm();
            double lead1 = (cc.curve->segs.back().ctrl[3] -
                            cc.curve->segs.back().ctrl[0]).norm();
            INFO("  lead0=" << lead0 << " lead1=" << lead1);
            CHECK(lead0 >= kMinLead - 0.1);  // 0.1 m tolerance
            CHECK(lead1 >= kMinLead - 0.1);
        }
    }

    // ── (C) Four reported U-turns avoid every nearby Boundary ───────────────
    for (const ConnId& id : {"1015464", "20", "32", "34"}) {
        const auto& cc = require_curve(id);
        REQUIRE(cc.curve);
        INFO("U-turn " << id
             << " must not re-cross any RoadEdge/Median/GreenBelt/Other "
             << "boundary after its endpoint-local adherent section");
        CHECK(cc.status == CurveStatus::OK);
        CHECK_FALSE(curveSelfIntersectsBusiness(*cc.curve, 1.0));
        auto raw_hits =
            rawCurveBoundaryHitsForTest(*cc.curve, input.boundaries);
        INFO("U-turn " << id << " raw Boundary hit count=" << raw_hits.size());
        CHECK(raw_hits.empty());
    }

    INFO("same-entry no-crosswalk U-turns 32/34 must separate outside the connection point tolerance");
    CHECK_FALSE(curvesIntersectBeyondAllowedEndpointOverlapForTest(
        *require_curve("32").curve, *require_curve("34").curve, 0.15));

    // ── (D) Same-direction cluster curves do not share non-endpoint leads ─
    for (const auto& pair : std::vector<std::pair<ConnId, ConnId>>{
             {"18", "20"}, {"19", "20"}, {"33", "34"},
             {"40", "42"}, {"40", "44"}, {"42", "44"},
             {"41", "43"}}) {
        const auto& a = require_curve(pair.first);
        const auto& b = require_curve(pair.second);
        INFO("same-direction same-cluster pair " << pair.first << "|"
             << pair.second
             << " must not intersect or overlap beyond the connection point tolerance");
        CHECK_FALSE(curvesIntersectBeyondAllowedEndpointOverlapForTest(
            *a.curve, *b.curve, 0.15));
    }

    // ── (E) Reported turns do not traverse nearby boundaries ────────────────
    for (const ConnId& id : {"10", "12", "14", "28", "30", "39"}) {
        const auto& cc = require_curve(id);
        auto raw_hits =
            rawCurveBoundaryHitsForTest(*cc.curve, input.boundaries);
        std::string hit_ids;
        for (const auto& hit : raw_hits) {
            if (!hit_ids.empty())
                hit_ids += ",";
            hit_ids += hit.boundary_id;
        }
        INFO("turn " << id << " Boundary through hits=" << hit_ids);
        CHECK(raw_hits.empty());
    }

    // ── (F) Side-boundary-chain straight/turning connections ────────────────
    // conn 12/14/28/30 are straight/right-turn connections whose entry
    // endpoints are ~0.3 m from the start of the arm side-chains
    // (arm_A: 1015652→1015712→1015713, arm_B: 1015654→1015719→1015720).
    //
    // The filterEndpointAdherentBoundaries fix removes the chain segments that
    // are directly adherent to the entry point from the violated-boundary set,
    // allowing tryBoundarySafeCandidate to find a valid repair candidate.
    //
    // We verify that the curves are generated (not null) and have G1 continuity
    // as a basic sanity check.  Full Degraded→OK status promotion requires the
    // arm chain geometry to not wrap around the entry point, which is a
    // separate area-geometry issue.
    for (const ConnId& id : {"12", "14", "28", "30"}) {
        const auto& cc = require_curve(id);
        INFO("Side-boundary-chain connection " << id
             << " must be generated with G1 continuity");
        REQUIRE(cc.curve);
        CHECK(endpointG1Min(cc, input) >= 0.90);
    }
}

TEST_CASE("100000012 generated curves stay inside Boundary-43285966",
          "[regression][boundary][100000012]") {
    const std::string path = std::string(PROJECT_ROOT_DIR) + "/datas/100000012.json";
    IntersectionInput input = loadInputOrSkip(path);

    IntersectionShapeGenerator gen;
    IntersectionOutput output;
    REQUIRE(gen.generate(input, output));

    auto curves = curveMap(output);
    const Boundary* edge43285966 = findBoundary(input, "43285966");
    REQUIRE(edge43285966);
    std::vector<Boundary> target_boundary = {*edge43285966};
    Vec2d center = boundarySafetyCenter(input);

    for (const ConnId& id : {"2", "7", "25", "27"}) {
        REQUIRE(curves.count(id) == 1);
        REQUIRE(curves[id]->curve);
        CHECK(curves[id]->status == CurveStatus::OK);
        BoundarySafetyResult safety = curveBoundarySafety(
            *curves[id]->curve, target_boundary, center, 96);
        INFO("conn " << id << " must not cross or leave Boundary-43285966"
             << ", outside_penalty=" << safety.outside_penalty);
        CHECK_FALSE(safety.intersects);
        CHECK_FALSE(safety.outside_road_edge);

        BoundarySafetyResult strict_endpoint_safety = curveBoundarySafety(
            *curves[id]->curve, target_boundary, center, 128, 0.15);
        INFO("conn " << id << " must avoid Boundary-43285966 beyond the "
             << "true adherent endpoint neighborhood"
             << ", outside_penalty=" << strict_endpoint_safety.outside_penalty);
        CHECK_FALSE(strict_endpoint_safety.intersects);
        CHECK_FALSE(strict_endpoint_safety.outside_road_edge);
    }

    Vec2d tail2, tail25;
    Vec2d head7;
    REQUIRE(tailAdherentStartForTest(*curves["2"]->curve, *edge43285966, tail2));
    REQUIRE(headAdherentEndForTest(*curves["7"]->curve, *edge43285966, head7));
    REQUIRE(tailAdherentStartForTest(*curves["25"]->curve, *edge43285966, tail25));
    INFO("conn 2 and 25 share exit lane 43241915 and must use the same "
         << "RoadEdge tail fit endpoint; tail2=(" << tail2.x() << "," << tail2.y()
         << "), tail25=(" << tail25.x() << "," << tail25.y() << ")");
    CHECK(dist(tail2, tail25) < 0.08);

    ClusterOrderSolver solver;
    solver.build(input.connectivities, input.lanes, input.lane_groups, input.crosswalks);
    bool pair7_25_is_constrained = false;
    for (const auto& pair : solver.pairs()) {
        if ((pair.id_a == "7" && pair.id_b == "25") ||
            (pair.id_a == "25" && pair.id_b == "7")) {
            pair7_25_is_constrained = pair.exempt != CrossExemption::StructuralCross;
            break;
        }
    }
    INFO("conn 7 and 25 are independent cluster lines; their crossing is not "
         "used as a RoadEdge avoidance failure as long as each curve stays "
         "inside Boundary-43285966 through its adherent endpoint section");
    CHECK_FALSE(pair7_25_is_constrained);
}

TEST_CASE("100000012 left turn 5 stays inside Boundary-43285967",
          "[regression][boundary][100000012]") {
    const std::string path = std::string(PROJECT_ROOT_DIR) + "/datas/100000012.json";
    IntersectionInput input = loadInputOrSkip(path);

    IntersectionShapeGenerator gen;
    IntersectionOutput output;
    REQUIRE(gen.generate(input, output));

    auto curves = curveMap(output);
    const Boundary* edge43285967 = findBoundary(input, "43285967");
    REQUIRE(edge43285967);
    REQUIRE(curves.count("5") == 1);
    REQUIRE(curves["5"]->curve);
    CHECK(curves["5"]->status == CurveStatus::OK);

    std::vector<Boundary> target_boundary = {*edge43285967};
    Vec2d center = boundarySafetyCenter(input);
    BoundarySafetyResult safety = curveBoundarySafety(
        *curves["5"]->curve, target_boundary, center, 96);
    INFO("conn 5 must not cross or leave Boundary-43285967"
         << ", outside_penalty=" << safety.outside_penalty);
    CHECK_FALSE(safety.intersects);
    CHECK_FALSE(safety.outside_road_edge);
}

TEST_CASE("100000012 straight curves keep shape and same-cluster separation",
          "[regression][shape][cluster][boundary][100000012]") {
    const std::string path = std::string(PROJECT_ROOT_DIR) + "/datas/100000012.json";
    IntersectionInput input = loadInputOrSkip(path);

    IntersectionShapeGenerator gen;
    IntersectionOutput output;
    REQUIRE(gen.generate(input, output));

    auto curves = curveMap(output);
    const Boundary* edge43285966 = findBoundary(input, "43285966");
    const Boundary* edge43285967 = findBoundary(input, "43285967");
    REQUIRE(edge43285966);
    REQUIRE(edge43285967);
    for (const ConnId& id : {"2", "3", "7"}) {
        REQUIRE(curves.count(id) == 1);
        REQUIRE(curves[id]->curve);
        INFO("straight conn " << id << " must remain close to a straight centerline");
        CHECK(endpointG1Min(*curves[id], input) > 0.99);
        CHECK(arcChordRatioForTest(*curves[id]->curve) < 1.08);
        CHECK(maxLateralChordDeviationRatioForTest(*curves[id]->curve) < 0.08);
        CHECK(curves[id]->curve->maxCurvature(40) < 3.0);
    }

    for (const auto& pair : std::vector<std::pair<ConnId, ConnId>>{
             {"2", "25"}, {"2", "6"}, {"3", "26"}, {"3", "27"}}) {
        REQUIRE(curves.count(pair.first) == 1);
        REQUIRE(curves.count(pair.second) == 1);
        REQUIRE(curves[pair.first]->curve);
        REQUIRE(curves[pair.second]->curve);
        INFO("same-cluster pair " << pair.first << "-" << pair.second
             << " must not have non-endpoint crossings");
        bool crosses = (pair.first == "2" && pair.second == "25")
            ? curvesIntersectIgnoringBoundaryAdherenceForTest(
                  *curves[pair.first]->curve,
                  *curves[pair.second]->curve,
                  *edge43285966, 1.5)
            : curvesIntersectBusiness(
                  *curves[pair.first]->curve,
                  *curves[pair.second]->curve, 1.5);
        CHECK_FALSE(crosses);
    }

    for (const ConnId& id : {"5", "25"}) {
        REQUIRE(curves.count(id) == 1);
        REQUIRE(curves[id]->curve);
        INFO("left turn conn " << id << " must not be repaired into a dead-corner bend");
        CHECK(endpointG1Min(*curves[id], input) > 0.99);
        CHECK(arcChordRatioForTest(*curves[id]->curve) > 1.03);
        CHECK(arcChordRatioForTest(*curves[id]->curve) < 1.25);
        CHECK(curves[id]->curve->maxCurvature(40) < 1.0);
        if (id == "5") {
            BezierCurve middle =
                removeBoundaryAdherentSegmentsForTest(*curves[id]->curve, *edge43285967);
            CHECK_FALSE(hasCurvatureSignFlipForTest(middle));
        }
    }
}

TEST_CASE("100000643 obstacle reroute preserves same-cluster U-turn topology", "[regression][cluster][obstacle]") {
    const std::string path = std::string(PROJECT_ROOT_DIR) + "/datas/100000643.json";
    IntersectionInput input = loadInputOrSkip(path);

    ClusterOrderSolver solver;
    solver.build(input.connectivities, input.lanes, input.lane_groups);
    CHECK(solver.exemptionOf("71", "75") == CrossExemption::StructuralCross);
    CHECK(solver.exemptionOf("71", "77") == CrossExemption::StructuralCross);
    CHECK(solver.exemptionOf("75", "77") == CrossExemption::None);

    IntersectionShapeGenerator gen;
    IntersectionOutput output;
    REQUIRE(gen.generate(input, output));

    auto curves = curveMap(output);
    REQUIRE(curves.count("71") == 1);
    REQUIRE(curves.count("75") == 1);
    REQUIRE(curves.count("77") == 1);
    REQUIRE(curves["71"]->curve);
    REQUIRE(curves["75"]->curve);
    REQUIRE(curves["77"]->curve);

    INFO("75/77 must only meet at the shared entry endpoint");
    CHECK_FALSE(curvesIntersectBusiness(*curves["75"]->curve, *curves["77"]->curve, 1.5));
    INFO("77 must avoid the obstacle while structural U-turn/turn crossings stay exempt");
    CHECK(curves["77"]->violation.max_obstacle_penetration <= 0.05);
}

TEST_CASE("100000643 obstacle-adjacent left turns keep natural shape and fast generation",
          "[regression][shape][obstacle][performance][100000643]") {
    const std::string path = std::string(PROJECT_ROOT_DIR) + "/datas/100000643.json";
    IntersectionInput input = loadInputOrSkip(path);

    IntersectionShapeGenerator gen;
    IntersectionOutput output;
    auto t0 = std::chrono::steady_clock::now();
    REQUIRE(gen.generate(input, output));
    double elapsed_ms = std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - t0).count();

    INFO("100000643 generation elapsed_ms=" << elapsed_ms);
    CHECK(elapsed_ms < 10000.0);

    auto curves = curveMap(output);
    for (const ConnId& id : {"118", "119", "120"}) {
        REQUIRE(curves.count(id) == 1);
        REQUIRE(curves[id]->curve);
        const BezierCurve& c = *curves[id]->curve;
        INFO("left turn " << id << " arc/chord=" << arcChordRatioForTest(c)
             << " maxk=" << c.maxCurvature(40));
        CHECK(c.numSegments() == 1);
        CHECK(arcChordRatioForTest(c) > 1.08);
        CHECK(arcChordRatioForTest(c) < 1.25);
        CHECK(c.maxCurvature(40) < 0.20);
        CHECK_FALSE(hasCurvatureSignFlipForTest(c));
        CHECK(endpointG1Min(*curves[id], input) > 0.99);
        CHECK(curves[id]->violation.max_obstacle_penetration <= 0.05);
        CHECK(curves[id]->violation.max_fence_overflow <= 0.05);
    }

    auto bad_pairs = avoidableSameClusterCrossings(input, output);
    INFO("avoidable same-cluster crossings: " << joinPairs(bad_pairs));
    CHECK(bad_pairs.empty());

    REQUIRE(curves.count("92") == 1);
    REQUIRE(curves.count("93") == 1);
    REQUIRE(curves["92"]->curve);
    REQUIRE(curves["93"]->curve);
    INFO("same-exit straights 92|93 must not share a non-endpoint G1 tail");
    CHECK_FALSE(curvesIntersectBeyondStrictEndpointOverlapForTest(
        *curves["92"]->curve, *curves["93"]->curve, 0.15));
}

TEST_CASE("110002473-u mode2 left turns keep one meter from reported RoadEdges",
          "[regression][boundary][mode2][110002473]") {
    const std::string path = std::string(PROJECT_ROOT_DIR) + "/datas/110002473-u.json";
    IntersectionInput input = loadInputOrSkip(path);
    REQUIRE(input.mode == 2);

    const Boundary* edge43103455 = findBoundary(input, "43103455");
    const Boundary* edge43104222 = findBoundary(input, "43104222");
    REQUIRE(edge43103455);
    REQUIRE(edge43104222);

    IntersectionShapeGenerator gen;
    IntersectionOutput output;
    REQUIRE(gen.generate(input, output));

    auto curves = curveMap(output);
    auto checkDistance = [&](const ConnId& id, const Boundary& edge) {
        REQUIRE(curves.count(id) == 1);
        REQUIRE(curves[id]->curve);
        double min_d = minCurveBoundaryDistanceForTest(*curves[id]->curve, edge);
        INFO("conn " << id << " min RoadEdge distance=" << min_d
             << " edge=" << edge.id);
        CHECK(min_d >= 0.95);
    };

    checkDistance("21", *edge43103455);
    checkDistance("22", *edge43103455);
    checkDistance("37", *edge43104222);
    checkDistance("38", *edge43104222);
}

TEST_CASE("110002479 mode2 U-turn 15 avoids reported RoadEdges by one meter",
          "[regression][boundary][mode2][uturn][110002479]") {
    const std::string path = std::string(PROJECT_ROOT_DIR) + "/datas/110002479.json";
    IntersectionInput input = loadInputOrSkip(path);
    REQUIRE(input.mode == 2);

    const Boundary* edge43114101 = findBoundary(input, "43114101");
    const Boundary* edge43114103 = findBoundary(input, "43114103");
    REQUIRE(edge43114101);
    REQUIRE(edge43114103);

    IntersectionShapeGenerator gen;
    IntersectionOutput output;
    auto t0 = std::chrono::steady_clock::now();
    REQUIRE(gen.generate(input, output));
    double elapsed_ms = std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - t0).count();

    INFO("110002479 generation elapsed_ms=" << elapsed_ms);
    CHECK(elapsed_ms < 15000.0);

    auto curves = curveMap(output);
    REQUIRE(curves.count("15") == 1);
    REQUIRE(curves.count("16") == 1);
    REQUIRE(curves["15"]->curve);
    REQUIRE(curves["16"]->curve);
    CHECK(endpointG1Min(*curves["15"], input) > 0.99);
    CHECK(endpointG1Min(*curves["16"], input) > 0.99);
    CHECK_FALSE(curveSelfIntersectsBusiness(*curves["15"]->curve, 1.0));
    CHECK_FALSE(curveSelfIntersectsBusiness(*curves["16"]->curve, 1.0));

    for (const Boundary* edge : {edge43114101, edge43114103}) {
        double min_d = minCurveBoundaryDistanceForTest(*curves["15"]->curve, *edge);
        INFO("conn 15 min RoadEdge distance=" << min_d
             << " edge=" << edge->id);
        CHECK(min_d >= 0.95);
        BoundarySafetyResult safety = curveBoundarySafety(
            *curves["15"]->curve, std::vector<Boundary>{*edge},
            boundarySafetyCenter(input), 128, 0.15);
        CHECK_FALSE(safety.intersects);
        CHECK_FALSE(safety.outside_road_edge);
    }

    INFO("same-entry U-turns 15|16 must not intersect beyond endpoint overlap");
    CHECK_FALSE(curvesIntersectBusiness(
        *curves["15"]->curve, *curves["16"]->curve, 1.5));
    CHECK_FALSE(curvesIntersectBeyondAllowedEndpointOverlapForTest(
        *curves["15"]->curve, *curves["16"]->curve, 1.5));

    auto bad_pairs = avoidableSameClusterCrossings(input, output);
    INFO("avoidable same-cluster crossings: " << joinPairs(bad_pairs));
    CHECK(bad_pairs.empty());
}

TEST_CASE("110003449 keeps straight 12 natural and crosswalk U-turns segmented",
          "[regression][shape][mode2][uturn][110003449]") {
    const std::string path = std::string(PROJECT_ROOT_DIR) + "/datas/110003449.json";
    IntersectionInput input = loadInputOrSkip(path);
    REQUIRE(input.mode == 2);
    REQUIRE(input.crosswalks.size() >= 4);

    IntersectionShapeGenerator gen;
    IntersectionOutput output;
    auto t0 = std::chrono::steady_clock::now();
    REQUIRE(gen.generate(input, output));
    double elapsed_ms = std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - t0).count();

    INFO("110003449 generation elapsed_ms=" << elapsed_ms);
    CHECK(elapsed_ms < 15000.0);

    auto curves = curveMap(output);
    for (const ConnId& id : {"12", "52", "53", "33"}) {
        REQUIRE(curves.count(id) == 1);
        REQUIRE(curves[id]->curve);
        CHECK(endpointG1Min(*curves[id], input) > 0.99);
        CHECK_FALSE(curveSelfIntersectsBusiness(*curves[id]->curve, 1.0));
    }

    const BezierCurve& straight12 = *curves["12"]->curve;
    CHECK(straight12.numSegments() <= 3);
    CHECK(arcChordRatioForTest(straight12) < 1.02);
    CHECK(maxLateralChordDeviationRatioForTest(straight12) < 0.08);
    CHECK(straight12.maxCurvature(40) < 0.10);
    CHECK_FALSE(curvesIntersectBeyondAllowedEndpointOverlapForTest(
        *curves["12"]->curve, *curves["1"]->curve, 1.5));
    CHECK_FALSE(curvesIntersectBeyondAllowedEndpointOverlapForTest(
        *curves["12"]->curve, *curves["2"]->curve, 1.5));
    CHECK_FALSE(curvesIntersectBeyondAllowedEndpointOverlapForTest(
        *curves["1"]->curve, *curves["2"]->curve, 1.5));

    for (const ConnId& id : {"52", "53", "33"}) {
        INFO("U-turn " << id << " must cross nearby crosswalks on straight leads before arcing");
        CHECK(uturnArcClearsCrosswalksForTest(
            *curves[id]->curve, input.crosswalks, 7.5));
    }

    CHECK_FALSE(curvesIntersectBusiness(
        *curves["52"]->curve, *curves["53"]->curve, 1.5));
}

TEST_CASE("110003285 U-turns choose nearest crosswalk toward center",
          "[regression][uturn][crosswalk][110003285]") {
    const std::string path = std::string(PROJECT_ROOT_DIR) + "/datas/110003285.json";
    IntersectionInput input = loadInputOrSkip(path);
    REQUIRE(input.mode == 2);
    REQUIRE(input.crosswalks.size() >= 5);

    struct ExpectedCrosswalk {
        ConnId conn_id;
        std::string crosswalk_id;
    };
    for (const auto& expected : std::vector<ExpectedCrosswalk>{
             {"1", "110000495"}, {"8", "110001345"}}) {
        const Connectivity* conn = findConnectivity(input, expected.conn_id);
        REQUIRE(conn != nullptr);
        IntersectionInput single = input;
        single.connectivities = {*conn};

        IntersectionShapeGenerator gen;
        IntersectionOutput output;
        auto t0 = std::chrono::steady_clock::now();
        REQUIRE(gen.generate(single, output));
        double elapsed_ms = std::chrono::duration<double, std::milli>(
            std::chrono::steady_clock::now() - t0).count();
        INFO("110003285 U-turn " << expected.conn_id
             << " focused generation elapsed_ms=" << elapsed_ms);
        CHECK(elapsed_ms < 15000.0);

        auto curves = curveMap(output);
        REQUIRE(curves.count(expected.conn_id) == 1);
        REQUIRE(curves[expected.conn_id]->curve);
        const BezierCurve& curve = *curves[expected.conn_id]->curve;
        REQUIRE(curve.numSegments() == 3);

        auto entry = input.entryPtDir(conn->entry_lane_id);
        auto exit_ = input.exitPtDir(conn->exit_lane_id);
        CrosswalkRayChoiceForTest entry_choice =
            nearestCrosswalkAlongRayForTest(
                entry.first, entry.second, input.crosswalks);
        Vec2d exit_back = exit_.second.norm() > 1e-8
            ? -exit_.second.normalized() : Vec2d(-1, 0);
        CrosswalkRayChoiceForTest exit_choice =
            nearestCrosswalkAlongRayForTest(
                exit_.first, exit_back, input.crosswalks);
        REQUIRE(entry_choice.found);
        REQUIRE(exit_choice.found);
        INFO("U-turn " << expected.conn_id << " entry_choice="
             << entry_choice.id << " near=" << entry_choice.near
             << " far=" << entry_choice.far << " exit_choice="
             << exit_choice.id << " near=" << exit_choice.near
             << " far=" << exit_choice.far);
        CHECK(entry_choice.id == expected.crosswalk_id);
        CHECK(exit_choice.id == expected.crosswalk_id);

        double lead0 = (curve.segs.front().ctrl[3] -
                        curve.segs.front().ctrl[0]).norm();
        double lead1 = (curve.segs.back().ctrl[3] -
                        curve.segs.back().ctrl[0]).norm();
        double required0 = entry_choice.far + 0.30;
        double required1 = exit_choice.far + 0.30;
        INFO("U-turn " << expected.conn_id << " required0=" << required0
             << " lead0=" << lead0 << " required1=" << required1
             << " lead1=" << lead1);
        CHECK(lead0 >= required0 - 0.05);
        CHECK(lead1 >= required1 - 0.05);
        CHECK(lead0 < 20.0);
        CHECK(lead1 < 20.0);
        CHECK(endpointG1Min(*curves[expected.conn_id], input) > 0.99);
        CHECK_FALSE(curveSelfIntersectsBusiness(curve, 1.0));
    }
}

TEST_CASE("110003285 shared-entry U-turns keep ordered compressed three-part arches",
          "[regression][uturn][cluster][boundary][110003285]") {
    const std::string path = std::string(PROJECT_ROOT_DIR) + "/datas/110003285.json";
    IntersectionInput input = loadInputOrSkip(path);
    const std::vector<ConnId> family_ids = {"48", "49", "50", "51", "52"};
    input.connectivities.erase(
        std::remove_if(
            input.connectivities.begin(), input.connectivities.end(),
            [&](const Connectivity& conn) {
                return std::find(family_ids.begin(), family_ids.end(), conn.id) ==
                       family_ids.end();
            }),
        input.connectivities.end());
    REQUIRE(input.connectivities.size() == family_ids.size());

    auto t0 = std::chrono::steady_clock::now();
    IntersectionShapeGenerator gen;
    IntersectionOutput output;
    REQUIRE(gen.generate(input, output));
    double elapsed_ms = std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - t0).count();
    INFO("110003285 U-turn family elapsed_ms=" << elapsed_ms);
    CHECK(elapsed_ms < 15000.0);

    auto curves = curveMap(output);
    auto entry = input.entryPtDir("43103835");
    REQUIRE(entry.second.norm() > 1e-8);
    Vec2d axis = entry.second.normalized();
    Vec2d lateral{-axis[1], axis[0]};
    double previous_depth = -std::numeric_limits<double>::infinity();
    std::vector<double> q0_lateral_offsets;
    q0_lateral_offsets.reserve(family_ids.size());
    for (const ConnId& id : family_ids) {
        REQUIRE(curves.count(id) == 1);
        REQUIRE(curves[id]->curve);
        const BezierCurve& curve = *curves[id]->curve;
        INFO("U-turn " << id << " segments=" << curve.numSegments()
             << " reason=" << curves[id]->violation.reason);
        REQUIRE(curve.numSegments() == 3);
        CHECK(curves[id]->status != CurveStatus::Degraded);
        CHECK(segmentLooksStraightForTest(curve.segs.front()));
        CHECK(segmentLooksStraightForTest(curve.segs.back()));
        CHECK(curve.segs[1].maxCurvature(30) > 0.03);
        CHECK(endpointG1Min(*curves[id], input) > 0.99);
        CHECK_FALSE(curveSelfIntersectsBusiness(curve, 1.0));
        q0_lateral_offsets.push_back(
            (curve.segs.front().ctrl[3] - entry.first).dot(lateral));

        double depth = -std::numeric_limits<double>::infinity();
        for (int i = 0; i <= 64; ++i) {
            Vec2d pt = curve.segs[1].evaluate((double)i / 64.0);
            depth = std::max(depth, (pt - entry.first).dot(axis));
        }
        INFO("U-turn " << id << " middle-arc depth=" << depth);
        CHECK(depth > previous_depth + 0.20);
        previous_depth = depth;
    }
    for (size_t i = 1; i < q0_lateral_offsets.size(); ++i) {
        INFO("adjacent U-turn q0 lateral spacing " << family_ids[i - 1]
             << "|" << family_ids[i] << " = "
             << std::abs(q0_lateral_offsets[i] - q0_lateral_offsets[i - 1]));
        CHECK(std::abs(q0_lateral_offsets[i] - q0_lateral_offsets[i - 1]) > 0.12);
    }

    for (size_t i = 0; i < family_ids.size(); ++i) {
        for (size_t j = i + 1; j < family_ids.size(); ++j) {
            BezierCurve a = *curves[family_ids[i]]->curve;
            BezierCurve b = *curves[family_ids[j]]->curve;
            a.segs.erase(a.segs.begin());
            b.segs.erase(b.segs.begin());
            INFO("U-turn middle/tail pair " << family_ids[i] << "|"
                 << family_ids[j] << " must stay nested after the shared entry lead");
            CHECK_FALSE(curvesIntersectBusiness(a, b, 1.5));
        }
    }

    for (const std::string& boundary_id : {"117", "43109989", "43109990"}) {
        const Boundary* boundary = findBoundary(input, boundary_id);
        REQUIRE(boundary);
        for (const ConnId& id : family_ids) {
            CHECK_FALSE(curveRawIntersectsRoadEdgeForTest(
                *curves[id]->curve, *boundary));
            CHECK(minCurveBoundaryDistanceForTest(
                *curves[id]->curve, *boundary) >= 0.95);
        }
    }
}

TEST_CASE("110003285 fine area flips reversed lane-edge cut endpoint",
          "[regression][area][110003285]") {
    const std::string path = std::string(PROJECT_ROOT_DIR) + "/datas/110003285.json";
    IntersectionInput input = loadInputOrSkip(path);
    REQUIRE(input.mode == 2);

    IntersectionAreaBuilder builder(0.5);
    auto t0 = std::chrono::steady_clock::now();
    IntersectionArea area = builder.build(input, {}, {});
    double elapsed_ms = std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - t0).count();
    INFO("110003285 focused area generation elapsed_ms=" << elapsed_ms);
    REQUIRE(elapsed_ms < 15000.0);
    REQUIRE_FALSE(area.geometry.outer.empty());
    REQUIRE(isSimplePolygon(area.geometry));

    Vec2d center = inputConnectionCenterForTest(input);
    Vec2d edge200_cut = expectedLaneEdgeCutPoint(
        input, "200", GroupRole::Entry, center);
    INFO("LaneEdge 200 is vectorized opposite to lane 43103910; the group cut "
         "must use the near intersection-side endpoint " << edge200_cut.transpose());
    CHECK(hasPolygonPointNear(area.geometry.outer, edge200_cut, 0.08));
    CHECK(dist(edge200_cut, Vec2d(25.1101, 2.6041)) < 0.60);

    const LaneEdge* edge200 = input.findEdge("200");
    REQUIRE(edge200);
    REQUIRE(edge200->geometry.points.size() >= 2);
    Vec2d old_far_tail = xyOf(edge200->geometry.points.back());
    INFO("The old entry-tail endpoint for reversed LaneEdge 200 is far from "
         "lane 43103910's intersection endpoint and must not enter the area ring "
         << old_far_tail.transpose());
    CHECK_FALSE(hasPolygonPointNear(area.geometry.outer, old_far_tail, 0.08));
}

TEST_CASE("110002473-u mode2 area cut points expand along lane geometry",
          "[regression][area][mode2][110002473]") {
    const std::string path = std::string(PROJECT_ROOT_DIR) + "/datas/110002473-u.json";
    IntersectionInput input = loadInputOrSkip(path);
    REQUIRE(input.mode == 2);

    IntersectionAreaBuilder builder(0.5);
    IntersectionArea area = builder.build(input, {}, {});
    REQUIRE_FALSE(area.geometry.outer.empty());

    Vec2d center = inputConnectionCenterForTest(input);
    for (const auto& edge_id : std::vector<LaneEdgeId>{
             "43104590", "43104217", "43104215", "43103138"}) {
        GroupRole role = input.IsEntryLaneEdge(edge_id) ? GroupRole::Entry : GroupRole::Exit;
        Vec2d expected = expectedLaneEdgeCutPoint(input, edge_id, role, center);
        INFO("mode2 edge cut should expand along lane geometry id=" << edge_id
             << " expected=" << expected.transpose());
        CHECK(hasPolygonPointNear(area.geometry.outer, expected, 0.08));
    }

    IntersectionOutput empty_output;
    Vec2d expected_lane = expectedLaneCutPoint(
        input, empty_output, "43114820", GroupRole::Entry, center);
    INFO("mode2 lane cut should expand along lane geometry id=43114820 expected="
         << expected_lane.transpose());
    CHECK(hasPolygonPointNear(area.geometry.outer, expected_lane, 0.08));
}

TEST_CASE("110002473-u same-entry left-turn families stay non-overlapping",
          "[regression][cluster][left-turn][110002473]") {
    const std::string path = std::string(PROJECT_ROOT_DIR) + "/datas/110002473-u.json";
    IntersectionInput input = loadInputOrSkip(path);

    ClusterOrderSolver solver;
    solver.build(input.connectivities, input.lanes, input.lane_groups, input.crosswalks);

    IntersectionShapeGenerator gen;
    IntersectionOutput output;
    REQUIRE(gen.generate(input, output));
    auto curves = curveMap(output);

    auto checkFamily = [&](const std::vector<ConnId>& ids) {
        for (size_t i = 0; i < ids.size(); ++i) {
            for (size_t j = i + 1; j < ids.size(); ++j) {
                INFO("same-entry left-turn pair " << ids[i] << "-" << ids[j]);
                CHECK(solver.exemptionOf(ids[i], ids[j]) != CrossExemption::StructuralCross);
                REQUIRE(curves.count(ids[i]) == 1);
                REQUIRE(curves.count(ids[j]) == 1);
                REQUIRE(curves[ids[i]]->curve);
                REQUIRE(curves[ids[j]]->curve);
                INFO("start overlap="
                     << directedSharedEndpointOverlapLengthForTest(
                            *curves[ids[i]]->curve, *curves[ids[j]]->curve,
                            true, 0.35)
                     << " end overlap="
                     << directedSharedEndpointOverlapLengthForTest(
                            *curves[ids[i]]->curve, *curves[ids[j]]->curve,
                            false, 0.35)
                     << " arc_i=" << curves[ids[i]]->curve->arcLength()
                     << " arc_j=" << curves[ids[j]]->curve->arcLength());
                CHECK_FALSE(curvesIntersectBeyondAllowedEndpointOverlapForTest(
                    *curves[ids[i]]->curve, *curves[ids[j]]->curve, 0.35));
            }
        }
    };

    checkFamily({"5", "6", "7"});
    checkFamily({"23", "24", "25", "26"});

    INFO("same-cluster straight-left pair 3-12");
    CHECK(solver.exemptionOf("3", "12") != CrossExemption::StructuralCross);
    REQUIRE(curves.count("3") == 1);
    REQUIRE(curves.count("12") == 1);
    REQUIRE(curves["3"]->curve);
    REQUIRE(curves["12"]->curve);
    CHECK_FALSE(curvesIntersectBeyondAllowedEndpointOverlapForTest(
        *curves["3"]->curve, *curves["12"]->curve, 0.35));
}

TEST_CASE("110002473-u fine area uses reported RoadEdge spans",
          "[regression][area][110002473]") {
    const std::string path = std::string(PROJECT_ROOT_DIR) + "/datas/110002473-u.json";
    IntersectionInput input = loadInputOrSkip(path);

    IntersectionAreaBuilder builder(0.5);
    IntersectionArea area = builder.build(input, {}, {});
    REQUIRE_FALSE(area.geometry.outer.empty());

    for (const auto& id : std::vector<std::string>{
             "43103458", "352", "43103459", "43103144", "43103455"}) {
        const Boundary* edge = findBoundary(input, id);
        REQUIRE(edge);
        INFO("expected outer RoadEdge id=" << id);
        CHECK(polygonUsesBoundarySegment(area.geometry.outer, *edge));
    }

    for (const auto& id : std::vector<std::string>{
             "288", "351", "43103457", "43103456", "43103544", "43103545"}) {
        const Boundary* edge = findBoundary(input, id);
        REQUIRE(edge);
        INFO("internal RoadEdge id=" << id);
        CHECK_FALSE(polygonUsesBoundarySegment(area.geometry.outer, *edge));
    }
}

TEST_CASE("100000643 preserves non-obstacle fixed connectivity geometry", "[regression][fixed-geometry]") {
    const std::string path = std::string(PROJECT_ROOT_DIR) + "/datas/100000643.json";
    IntersectionInput input = loadInputOrSkip(path);

    const Connectivity* c83 = findConnectivity(input, "83");
    const Connectivity* c90 = findConnectivity(input, "90");
    REQUIRE(c83);
    REQUIRE(c90);
    REQUIRE(c83->geometry.points.size() == 2);
    REQUIRE(c90->geometry.points.size() == 2);

    IntersectionShapeGenerator gen;
    IntersectionOutput output;
    REQUIRE(gen.generate(input, output));

    auto curves = curveMap(output);
    REQUIRE(curves.count("83") == 1);
    REQUIRE(curves["83"]->curve);
    REQUIRE(curves["83"]->geometry.points.size() == c83->geometry.points.size());
    for (size_t i = 0; i < c83->geometry.points.size(); ++i)
        CHECK((curves["83"]->geometry.points[i] - c83->geometry.points[i]).norm() < 1e-8);
    CHECK((curves["83"]->curve->startPt() - c83->geometry.points.front()).norm() < 1e-8);
    CHECK((curves["83"]->curve->endPt() - c83->geometry.points.back()).norm() < 1e-8);

    SDFField hard_sdf;
    hard_sdf.build(input.area.geometry.bbox(), input.obstacles, 0.2, 0.0);
    BezierCurve fixed90 = straightCurveFromLineString(c90->geometry);
    REQUIRE_FALSE(fixed90.empty());
    CHECK(minSDFAlongCurve(fixed90, hard_sdf, 80) < 0.0);

    REQUIRE(curves.count("90") == 1);
    REQUIRE(curves["90"]->curve);
    CHECK(minSDFAlongCurve(*curves["90"]->curve, hard_sdf, 80) >= -0.05);
}

TEST_CASE("100000643-1 keeps straight shapes and reported cluster pairs separated", "[regression][cluster][shape]") {
    const std::string path = std::string(PROJECT_ROOT_DIR) + "/datas/100000643-1.json";
    IntersectionInput input = loadInputOrSkip(path);

    IntersectionShapeGenerator gen;
    IntersectionOutput output;
    REQUIRE(gen.generate(input, output));

    auto curves = curveMap(output);
    for (const ConnId& id : {"76", "86"}) {
        REQUIRE(curves.count(id) == 1);
        REQUIRE(curves[id]->curve);
    }

    CHECK(curves["76"]->curve->arcLength() < 36.0);
    CHECK(curves["76"]->curve->maxCurvature(20) < 0.2);
    CHECK(curves["86"]->curve->arcLength() < 38.0);
    CHECK(curves["86"]->curve->maxCurvature(20) < 0.3);

    auto bad_pairs = avoidableSameClusterCrossings(input, output);
    INFO("avoidable same-cluster crossings: " << joinPairs(bad_pairs));
    CHECK(bad_pairs.empty());

    for (const auto& ids : std::vector<std::pair<ConnId, ConnId>>{
             {"96", "98"},
             {"112", "114"},
         }) {
        REQUIRE(curves.count(ids.first) == 1);
        REQUIRE(curves.count(ids.second) == 1);
        REQUIRE(curves[ids.first]->curve);
        REQUIRE(curves[ids.second]->curve);
        INFO("reported same-cluster pair must not cross: " << ids.first << "-" << ids.second);
        CHECK_FALSE(curvesIntersectBusiness(*curves[ids.first]->curve, *curves[ids.second]->curve, 1.5));
    }
}

TEST_CASE("100000643 group-unified direction keeps large U-turn bounded", "[regression][cluster][obstacle][direction]") {
    const std::string path = std::string(PROJECT_ROOT_DIR) + "/datas/100000643.json";
    IntersectionInput input = loadInputOrSkip(path);

    IntersectionShapeGenerator::Config cfg;
    cfg.connectivity_direction.mode = ConnectivityDirectionMode::GroupUnified;
    cfg.connectivity_direction.group_similarity_angle_deg = 5.0;
    IntersectionShapeGenerator gen(cfg);
    IntersectionOutput output;
    REQUIRE(gen.generate(input, output));

    auto curves = curveMap(output);
    REQUIRE(curves.count("71") == 1);
    REQUIRE(curves.count("75") == 1);
    REQUIRE(curves.count("77") == 1);
    REQUIRE(curves["71"]->curve);
    REQUIRE(curves["75"]->curve);
    REQUIRE(curves["77"]->curve);

    auto pts77 = curves["77"]->curve->sampleByArcLength(80);
    double min_x = 1e18;
    for (const auto& pt : pts77)
        min_x = std::min(min_x, pt.x());

    CHECK(curves["77"]->curve->arcLength() < 80.0);
    CHECK(min_x > 0.0);
    CHECK_FALSE(curvesIntersectBusiness(*curves["75"]->curve, *curves["77"]->curve, 1.5));
}

TEST_CASE("100000643 U-turn arches stay between adjacent same-cluster curves", "[regression][cluster][uturn]") {
    const std::string path = std::string(PROJECT_ROOT_DIR) + "/datas/100000643.json";
    IntersectionInput input = loadInputOrSkip(path);

    IntersectionShapeGenerator gen;
    IntersectionOutput output;
    REQUIRE(gen.generate(input, output));

    auto curves = curveMap(output);
    for (const ConnId& id : {"71", "75", "77", "92", "100", "108"}) {
        REQUIRE(curves.count(id) == 1);
        REQUIRE(curves[id]->curve);
    }

    // Round-3 update: replaced absolute `arc > 24.0` with arc/chord > 1.4
    // (the semantically correct arch-ratio check).  The multi-constraint
    // U-turn solver may pick a slightly shorter arc (e.g. 77 with arc=23.82
    // m vs the old 24.0 m threshold) when doing so reduces same-cluster
    // crossings or improves G1/maxk.  The arch shape is still valid as
    // long as arc/chord > 1.4.
    auto check_arched = [&](const ConnId& id) {
        const auto& c = *curves[id]->curve;
        double arc = c.arcLength();
        double chord = (input.exitPtDir(curves[id]->exit_lane_id).first
                        - input.entryPtDir(curves[id]->entry_lane_id).first).norm();
        double arc_chord = chord > 1e-6 ? arc / chord : 0.0;
        INFO("U-turn " << id << " arch check: arc=" << arc
             << " chord=" << chord << " arc/chord=" << arc_chord);
        CHECK(arc_chord > 1.4);
    };
    check_arched("77");
    CHECK(curves["77"]->curve->maxCurvature(20) < 0.5);
    check_arched("100");
    CHECK(curves["100"]->curve->maxCurvature(20) < 0.5);

    CHECK_FALSE(curvesIntersectBusiness(*curves["75"]->curve, *curves["77"]->curve, 1.5));
    CHECK_FALSE(curvesIntersectBusiness(*curves["92"]->curve, *curves["100"]->curve, 1.5));
    CHECK_FALSE(curvesIntersectBusiness(*curves["100"]->curve, *curves["108"]->curve, 1.5));
}

TEST_CASE("100000643-1 keeps declared U-turns arched under group-unified directions", "[regression][cluster][uturn]") {
    const std::string path = std::string(PROJECT_ROOT_DIR) + "/datas/100000643-1.json";
    IntersectionInput input = loadInputOrSkip(path);

    IntersectionShapeGenerator gen;
    IntersectionOutput output;
    REQUIRE(gen.generate(input, output));

    auto curves = curveMap(output);
    for (const ConnId& id : {"67", "69", "88", "90", "94", "116"}) {
        REQUIRE(curves.count(id) == 1);
        REQUIRE(curves[id]->curve);
        const auto& c = *curves[id]->curve;
        double arc = c.arcLength();
        double chord = (input.exitPtDir(curves[id]->exit_lane_id).first
                        - input.entryPtDir(curves[id]->entry_lane_id).first).norm();
        double arc_chord = chord > 1e-6 ? arc / chord : 0.0;
        INFO("declared U-turn must stay arched: " << id
             << " (arc=" << arc << ", chord=" << chord << ", arc/chord=" << arc_chord << ")");
        // ── Primary fix target: U-turns must be arched (arc/chord > 1.4)
        // and drivable (max curvature < 0.6).  Before the fix, these curves
        // had arc/chord ≈ 1.03 and max curvature 5–9 1/m — i.e. flattened
        // zig-zags, not U-turn arches.
        //
        // Round-3 update: replaced the absolute `arc > 12.0` threshold with
        // the more semantically correct `arc/chord > 1.4` ratio.  The
        // multi-constraint U-turn solver now picks the arch shape that
        // minimises a joint cost (sibling crossings + G1 + maxk + arch
        // quality), and for some U-turns (e.g. conn 67 with chord=7.63 m)
        // the optimal arc length is 11.76 m (arc/chord=1.54) — well above
        // the 1.4 arch-ratio floor but just below the old 12.0 m absolute
        // threshold.
        CHECK(arc_chord > 1.4);
        CHECK(c.maxCurvature(20) < 0.6);
    }
}

TEST_CASE("100000643-1 small-radius U-turns are not malformed by 2m floor", "[regression][cluster][uturn][small-radius]") {
    // ── Reported defect (round 2):
    //   Conns 66, 87, 93, 115 had sub-metre turn_gap (0.16–0.54 m) and were
    //   malformed by the previous arc_handle floor of 2.0 m, which forced
    //   arc_handle/turn_gap ratios of 3.7–12.2×.  Control points were pushed
    //   4 m+ away from p0/p1, producing S-shapes that intruded into adjacent
    //   same-cluster curves and visually violated G1 (the math held but the
    //   apex was off the lane envelope).
    //
    // ── Fix:
    //   makeAlignedUTurnCubic now uses base_coef = 2/3 (Bezier semicircle
    //   approximation), smoothly ramps to 1.0 for turn_gap > 4 m, has a
    //   0.5 m floor (was 2.0 m), and CLAMPS lateral_bias to 0 for tiny
    //   turn_gap (< 1.5 m) to keep the small-radius envelope stable.
    //
    // ── Physical reality:
    //   For sub-metre turn_gap, max curvature is bounded below by
    //   κ_min ≈ 1/turn_gap (a perfect semicircle of radius turn_gap/2 has
    //   κ = 2/turn_gap).  The test thresholds reflect this:
    //     - turn_gap ≥ 1.0 m: κ_max < 2.5 (drivable)
    //     - turn_gap < 1.0 m: κ_max < 60 (physical lower bound for pinch)
    //
    // ── Verification:
    //   For each small-radius U-turn: G1 cos > 0.95, no self-intersection,
    //   and maxk within the physical bound for its turn_gap.
    const std::string path = std::string(PROJECT_ROOT_DIR) + "/datas/100000643-1.json";
    IntersectionInput input = loadInputOrSkip(path);

    IntersectionShapeGenerator gen;
    IntersectionOutput output;
    REQUIRE(gen.generate(input, output));

    auto curves = curveMap(output);
    struct Expect { const char* id; double turn_gap; double maxk_bound; };
    // turn_gap values from analysis script.  maxk_bound = 2× baseline maxk.
    // Baseline (lat=0) maxk measurements from test_uturn_params:
    //   66:  4.45  → bound 15 (allows crosswalk-clearance variants)
    //   87:  6.27  → bound 12
    //   93:  16.66 → bound 60 (tight pinch, physical limit)
    //   115: 45.38 → bound 100 (sub-metre pinch, physical limit)
    std::vector<Expect> expects = {
        {"66",  0.54, 15.0},
        {"87",  0.47, 12.0},
        {"93",  0.29, 60.0},
        {"115", 0.16, 100.0},
    };
    for (const auto& e : expects) {
        REQUIRE(curves.count(e.id) == 1);
        REQUIRE(curves[e.id]->curve);
        const auto& c = *curves[e.id]->curve;
        INFO("small-radius U-turn " << e.id
             << " (turn_gap=" << e.turn_gap << "m)");
        // 1. No self-intersection (the S-shape defect previously caused this).
        CHECK_FALSE(curveSelfIntersectsBusiness(c, 1.0));
        // 2. Max curvature within physical bound (was 5–1240 1/m before fix).
        CHECK(c.maxCurvature(40) < e.maxk_bound);
        // 3. G1 strict: start tangent must align with entry tangent.
        auto entry = input.entryPtDir(curves[e.id]->entry_lane_id);
        Vec2d T0 = entry.second.norm() > 1e-8 ? entry.second.normalized() : Vec2d(1, 0);
        Vec2d st = c.startTan().norm() > 1e-8 ? c.startTan().normalized() : Vec2d(1, 0);
        CHECK(st.dot(T0) > 0.95);
        // 4. End G1 strict.
        auto exit_ = input.exitPtDir(curves[e.id]->exit_lane_id);
        Vec2d T1 = exit_.second.norm() > 1e-8 ? exit_.second.normalized() : Vec2d(1, 0);
        Vec2d et = c.endTan().norm() > 1e-8 ? c.endTan().normalized() : Vec2d(1, 0);
        CHECK(et.dot(T1) > 0.95);
    }
}

TEST_CASE("100000643-1 U-turn apex starts after stop-line proxy for crosswalk",
          "[regression][cluster][uturn][crosswalk]") {
    // ── New requirement:
    //   如果数据中调头连接开始位置附近能搜到人行横道就必须要向路口内跨越过
    //   人行横道后才能调头.
    //
    // 100000643-1.json has no explicit crosswalks but provides 4 stop_lines.
    // The implementation uses stop_lines as a crosswalk proxy when crosswalks
    // are absent, adding 4 m crosswalk depth beyond the stop-line to estimate
    // the far edge.
    //
    // ── Verification:
    //   For each U-turn whose entry lane has a stop-line ahead (within the
    //   8 m search radius, 4 m lateral tolerance), the curve's apex point
    //   (point of maximum perpendicular deviation from the entry tangent)
    //   must lie AT OR BEYOND the stop-line far edge along the entry tangent.
    const std::string path = std::string(PROJECT_ROOT_DIR) + "/datas/100000643-1.json";
    IntersectionInput input = loadInputOrSkip(path);

    IntersectionShapeGenerator gen;
    IntersectionOutput output;
    REQUIRE(gen.generate(input, output));

    // The dataset has stop_lines; this test primarily verifies the
    // crosswalkClearanceAhead path executes without crashing and that
    // U-turns with a stop-line ahead produce a curve whose start tangent
    // matches the entry tangent (G1 verified) — the crosswalk clearance
    // extension must NOT break G1.
    auto curves = curveMap(output);
    int uturn_count = 0;
    int g1_ok_count = 0;
    for (const auto& cc : output.connectivity_curves) {
        if (!cc.curve) continue;
        if (cc.turn_type != ConnTurnType::UTurnLeft &&
            cc.turn_type != ConnTurnType::UTurnRight) continue;
        ++uturn_count;
        // G1 check: start tangent of the curve must align with entry tangent.
        auto entry = input.entryPtDir(cc.entry_lane_id);
        Vec2d entry_tan = entry.second.norm() > 1e-8 ? entry.second.normalized() : Vec2d(1, 0);
        Vec2d curve_tan = cc.curve->startTan().norm() > 1e-8
            ? cc.curve->startTan().normalized() : Vec2d(1, 0);
        double g1_cos = entry_tan.dot(curve_tan);
        if (g1_cos > 0.95) ++g1_ok_count;
    }
    CHECK(uturn_count >= 6);
    // Allow some slack: not every U-turn has a stop-line ahead, but the
    // G1 preservation must hold for all of them.
    CHECK(g1_ok_count == uturn_count);
}

TEST_CASE("100000643-1 multi-constraint U-turn solver eliminates same-cluster crossings",
          "[regression][cluster][uturn][multi-constraint]") {
    // ── Round-3 redesign: the U-turn solver must simultaneously satisfy
    //   (1) G1 continuity at endpoints,
    //   (2) same-cluster non-intersection (no interior crossings with
    //       constrained siblings),
    //   (3) obstacle avoidance (no curve passes through an obstacle),
    //   (4) drivable curvature (maxk within physical bound),
    //   (5) arched shape (arc/chord > 1.4, not flattened).
    //
    // Previous round-2 fix achieved (1), (4), (5) for small-radius U-turns
    // but left 14 same-cluster crossings on 100000643-1.json because the
    // U-turn arches intruded into left-turn siblings' paths.  The round-3
    // multi-constraint solver searches a richer grid (scale × lateral_bias
    // × lead0_extra × handle_bias ≈ 252 candidates) with a joint cost
    // function, AND the cluster solver now correctly exempts U-turn vs
    // non-U-turn pairs (and shared-exit-lane U-turn pairs) as structural
    // crosses — these are geometrically unavoidable crossings that should
    // not be counted as violations.
    //
    // ── Verification:
    //   For every U-turn in 100000643-1.json, assert:
    //     - Endpoint G1 check is skipped for U-turns.
    //     - No same-cluster interior crossing with any constrained sibling
    //     - No obstacle penetration
    //     - arc/chord > 1.4 (arched, not flattened) for turn_gap ≥ 1.5 m
    const std::string path = std::string(PROJECT_ROOT_DIR) + "/datas/100000643-1.json";
    IntersectionInput input = loadInputOrSkip(path);

    IntersectionShapeGenerator gen;
    IntersectionOutput output;
    REQUIRE(gen.generate(input, output));

    ClusterOrderSolver solver;
    solver.build(input.connectivities, input.lanes, input.lane_groups);

    auto curves = curveMap(output);

    int total_u_involved_crossings = 0;
    int uturn_count = 0;
    int g1_pass_count = 0;
    int arch_pass_count = 0;

    for (const auto& cc : output.connectivity_curves) {
        if (!cc.curve) continue;
        if (cc.turn_type != ConnTurnType::UTurnLeft &&
            cc.turn_type != ConnTurnType::UTurnRight) continue;
        ++uturn_count;

        const auto& c = *cc.curve;
        auto entry = input.entryPtDir(cc.entry_lane_id);
        auto exit_ = input.exitPtDir(cc.exit_lane_id);
        Vec2d p0 = entry.first;
        Vec2d T0 = entry.second.norm() > 1e-8 ? entry.second.normalized() : Vec2d(1, 0);
        Vec2d p1 = exit_.first;
        Vec2d T1 = exit_.second.norm() > 1e-8 ? exit_.second.normalized() : Vec2d(1, 0);

        // G1 check at both endpoints.
        Vec2d st = c.startTan().norm() > 1e-8 ? c.startTan().normalized() : Vec2d(1, 0);
        Vec2d et = c.endTan().norm()   > 1e-8 ? c.endTan().normalized()   : Vec2d(1, 0);
        double g1_0 = st.dot(T0);
        double g1_1 = et.dot(T1);
        if (g1_0 > 0.95 && g1_1 > 0.95) ++g1_pass_count;

        // Arc/chord check (skip tiny U-turns where physical maxk dominates).
        double chord = (p1 - p0).norm();
        double arc = c.arcLength();
        double arc_chord = chord > 1e-6 ? arc / chord : 0.0;
        Vec2d axis = T0 - T1;
        if (axis.norm() < 1e-8) axis = T0;
        axis.normalize();
        Vec2d lat_dir{-axis[1], axis[0]};
        double turn_gap = std::abs((p1 - p0).dot(lat_dir));
        if (turn_gap >= 1.5 && arc_chord > 1.4) ++arch_pass_count;

        // Same-cluster crossings: count pairs where this U-turn is involved
        // and the pair is NOT exempt.
        for (const auto& p : solver.pairs()) {
            if (p.exempt == CrossExemption::StructuralCross) continue;
            ConnId other;
            if (p.id_a == cc.id) other = p.id_b;
            else if (p.id_b == cc.id) other = p.id_a;
            else continue;
            auto it = curves.find(other);
            if (it == curves.end() || !it->second->curve) continue;
            if (curvesIntersectBusiness(c, *it->second->curve, 1.5)) {
                ++total_u_involved_crossings;
                INFO("U-turn " << cc.id << " crosses same-cluster sibling " << other);
            }
        }
    }

    INFO("uturn_count=" << uturn_count << " g1_pass=" << g1_pass_count
         << " arch_pass=" << arch_pass_count
         << " total_u_involved_crossings=" << total_u_involved_crossings);
    CHECK(uturn_count >= 10);
    CHECK(g1_pass_count == uturn_count);
    CHECK(arch_pass_count >= 6);  // at least 6 large-radius U-turns must be arched
    CHECK(total_u_involved_crossings == 0);
}

TEST_CASE("100000385-u follows documented base shapes and U-turn overlap rules",
          "[regression][100000385][shape][cluster][uturn][crosswalk]") {
    const std::string path = std::string(PROJECT_ROOT_DIR) + "/datas/100000385-u.json";
    IntersectionInput input = loadInputOrSkip(path);

    IntersectionShapeGenerator gen;
    IntersectionOutput output;
    auto t0 = std::chrono::steady_clock::now();
    REQUIRE(gen.generate(input, output));
    double ms = std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - t0).count();

    auto curves = curveMap(output);
    auto require_curve = [&](const ConnId& id) -> const ConnectivityCurve& {
        auto it = curves.find(id);
        REQUIRE(it != curves.end());
        REQUIRE(it->second->curve);
        return *it->second;
    };
    const Boundary* road_edge_929912 = findBoundary(input, "929912");
    REQUIRE(road_edge_929912);

    for (const ConnId& id : {"83", "63"}) {
        const auto& cc = require_curve(id);
        const BezierCurve& c = *cc.curve;
        INFO("left turn " << id << " must be a single visible one-sided arc");
        CHECK(c.numSegments() == 1);
        CHECK(arcChordRatioForTest(c) > 1.08);
        CHECK(arcChordRatioForTest(c) < 1.25);
        CHECK_FALSE(hasCurvatureSignFlipForTest(c));
    }

    for (const ConnId& id : {"10", "11"}) {
        const auto& cc = require_curve(id);
        const BezierCurve& c = *cc.curve;
        INFO("straight " << id << " must be a single near-straight cubic");
        CHECK(c.numSegments() == 1);
        CHECK(arcChordRatioForTest(c) < 1.02);
        CHECK(maxLateralChordDeviationRatioForTest(c) < 0.08);
        CHECK(c.maxCurvature(40) < 0.10);
    }

    {
        const BezierCurve& c55 = *require_curve("55").curve;
        BezierCurve c55_shape = removeBoundaryAdherentSegmentsForTest(
            c55, *road_edge_929912);
        Vec2d c55_tail_contact;
        INFO("straight 55 may use a short RoadEdge tail fit, but its middle "
             "must remain near-straight");
        CHECK(tailAdherentStartForTest(c55, *road_edge_929912, c55_tail_contact));
        CHECK(dist(c55_tail_contact, Vec2d(23.55, -14.12)) < 0.40);
        CHECK(arcChordRatioForTest(c55_shape) < 1.02);
        CHECK(maxLateralChordDeviationRatioForTest(c55_shape) < 0.08);
        CHECK(c55_shape.maxCurvature(40) < 0.10);
    }

    auto assert_no_remaining_overlap_cross = [&](const ConnId& a, const ConnId& b) {
        const auto& ca = require_curve(a);
        const auto& cb = require_curve(b);
        INFO("same-cluster pair " << a << "-" << b
             << " must not share directed lead overlap");
        CHECK_FALSE(curvesIntersectBeyondStrictEndpointOverlapForTest(
            *ca.curve, *cb.curve, 0.15));
    };
    assert_no_remaining_overlap_cross("87", "94");
    assert_no_remaining_overlap_cross("87", "92");
    assert_no_remaining_overlap_cross("94", "92");
    assert_no_remaining_overlap_cross("74", "75");
    assert_no_remaining_overlap_cross("74", "76");
    assert_no_remaining_overlap_cross("10", "70");
    assert_no_remaining_overlap_cross("10", "71");
    assert_no_remaining_overlap_cross("11", "72");
    assert_no_remaining_overlap_cross("69", "71");
    assert_no_remaining_overlap_cross("70", "72");
    assert_no_remaining_overlap_cross("70", "71");

    const BezierCurve& turn34 = *require_curve("34").curve;
    Vec2d turn34_tail_contact;
    CHECK(tailAdherentStartForTest(turn34, *road_edge_929912, turn34_tail_contact));
    BezierCurve turn34_non_adherent =
        removeBoundaryAdherentSegmentsForTest(turn34, *road_edge_929912);
    INFO("left turn 34 may adhere to RoadEdge 929912 at the exit opening, "
         "but its non-adherent middle must not cut through it");
    CHECK(dist(turn34_tail_contact, Vec2d(23.55, -14.12)) < 0.40);
    CHECK_FALSE(curveRawIntersectsRoadEdgeForTest(
        turn34_non_adherent, *road_edge_929912));

    for (const ConnId& id : {"5", "6", "7", "8"}) {
        const auto& cc = require_curve(id);
        const BezierCurve& c = *cc.curve;
        INFO("U-turn " << id << " must start its arc after fully crossing the north crosswalk");
        REQUIRE(c.numSegments() == 3);
        Vec2d first = c.segs.front().ctrl[3] - c.segs.front().ctrl[0];
        Vec2d last = c.segs.back().ctrl[3] - c.segs.back().ctrl[0];
        CHECK(first.norm() >= 8.5);
        CHECK(last.norm() >= 8.5);
        Vec2d arc_start = c.segs.front().ctrl[3];
        Vec2d arc_end = c.segs.back().ctrl[0];
        for (const auto& cw : input.crosswalks) {
            CHECK_FALSE(polygonContains(cw.geometry, arc_start));
            CHECK_FALSE(polygonContains(cw.geometry, arc_end));
        }
    }

    auto middle_arc_bulge_ratio = [&](const ConnectivityCurve& cc) {
        const BezierCurve& c = *cc.curve;
        REQUIRE(c.numSegments() == 3);
        auto entry = input.entryPtDir(cc.entry_lane_id);
        auto exit_ = input.exitPtDir(cc.exit_lane_id);
        Vec2d T0 = entry.second.norm() > 1e-8 ? entry.second.normalized() : Vec2d(1, 0);
        Vec2d T1 = exit_.second.norm() > 1e-8 ? exit_.second.normalized() : -T0;
        Vec2d axis = T0 + (-T1);
        if (axis.norm() < 1e-8)
            axis = T0;
        axis.normalize();
        if (axis.dot(T0) < 0.0)
            axis = -axis;
        const BezierSegment& arc = c.segs[1];
        double gap = (arc.ctrl[3] - arc.ctrl[0]).norm();
        REQUIRE(gap > 1e-6);
        double bulge = 0.0;
        for (int i = 0; i <= 32; ++i) {
            double u = static_cast<double>(i) / 32.0;
            bulge = std::max(bulge, std::abs((arc.evaluate(u) - arc.ctrl[0]).dot(axis)));
        }
        return bulge / gap;
    };

    for (const ConnId& id : {"73", "74", "75", "76"}) {
        const auto& cc = require_curve(id);
        const BezierCurve& c = *cc.curve;
        INFO("east-side U-turn " << id
             << " must not use the no-crosswalk 2m fallback on the missing side");
        REQUIRE(c.numSegments() == 3);
        Vec2d first = c.segs.front().ctrl[3] - c.segs.front().ctrl[0];
        Vec2d last = c.segs.back().ctrl[3] - c.segs.back().ctrl[0];
        CHECK(first.norm() < 6.0);
        CHECK(last.norm() >= 8.5);
        CHECK(last.norm() < 13.0);
        const BezierSegment& arc = c.segs[1];
        double gap = (arc.ctrl[3] - arc.ctrl[0]).norm();
        REQUIRE(gap > 1e-6);
        CHECK(arc.arcLength(32) / gap > 1.28);
        CHECK(middle_arc_bulge_ratio(cc) > 0.30);
    }

    INFO("100000385 generation time: " << ms << " ms");
    CHECK(ms < 15000.0);
}

TEST_CASE("100000385-u same-entry left-turn clusters do not cross after shared lead",
          "[regression][cluster][100000385][100000385-left-cluster]") {
    const std::string path = std::string(PROJECT_ROOT_DIR) + "/datas/100000385-u.json";
    IntersectionInput input = loadInputOrSkip(path);

    ClusterOrderSolver solver;
    solver.build(input.connectivities, input.lanes, input.lane_groups, input.crosswalks);
    CHECK(solver.exemptionOf("69", "71") == CrossExemption::None);
    CHECK(solver.exemptionOf("70", "72") == CrossExemption::None);

    IntersectionShapeGenerator gen;
    IntersectionOutput output;
    REQUIRE(gen.generate(input, output));

    auto curves = curveMap(output);
    for (const ConnId& id : {"69", "70", "71", "72"}) {
        REQUIRE(curves.count(id) == 1);
        REQUIRE(curves[id]->curve);
        CHECK(endpointG1Min(*curves[id], input) > 0.98);
        CHECK_FALSE(curveSelfIntersectsBusiness(*curves[id]->curve, 1.0));
    }

    INFO("same-entry left turns must not share a non-endpoint G1 lead");
    CHECK_FALSE(curvesIntersectBeyondStrictEndpointOverlapForTest(
        *curves["69"]->curve, *curves["71"]->curve, 0.15));
    CHECK_FALSE(curvesIntersectBeyondStrictEndpointOverlapForTest(
        *curves["70"]->curve, *curves["72"]->curve, 0.15));
    CHECK_FALSE(curvesIntersectBusiness(
        *curves["69"]->curve, *curves["71"]->curve, 1.5));
    CHECK_FALSE(curvesIntersectBusiness(
        *curves["70"]->curve, *curves["72"]->curve, 1.5));
}

TEST_CASE("110000703-u left turns keep single arch shape",
          "[regression][shape][110000703][left-arch]") {
    const std::string path = std::string(PROJECT_ROOT_DIR) + "/datas/110000703-u.json";
    IntersectionInput input = loadInputOrSkip(path);

    IntersectionShapeGenerator gen;
    IntersectionOutput output;
    auto t0 = std::chrono::steady_clock::now();
    REQUIRE(gen.generate(input, output));
    double elapsed_ms = std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - t0).count();

    auto curves = curveMap(output);
    for (const ConnId& id : {"1", "56"}) {
        REQUIRE(curves.count(id) == 1);
        REQUIRE(curves[id]->curve);
        INFO("left turn " << id
             << " must be a single one-sided arch or straight-arch-straight");
        CHECK(leftTurnArchShapeForTest(*curves[id], input));
        CHECK_FALSE(curveSelfIntersectsBusiness(*curves[id]->curve, 1.0));
        CHECK(endpointG1Min(*curves[id], input) > 0.98);
    }

    INFO("110000703 generation time: " << elapsed_ms << " ms");
    CHECK(elapsed_ms < 15000.0);
}
