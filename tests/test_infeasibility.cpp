#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_all.hpp>
#include "constraints/infeasibility_detector.h"
#include "constraints/fence_check.h"
#include "optimizer/sdf_field.h"

using namespace isg;

static Polygon2d makeSquarePoly(Vec2d c, double half) {
    Polygon2d p;
    p.outer = {{c[0]-half,c[1]-half},{c[0]+half,c[1]-half},
               {c[0]+half,c[1]+half},{c[0]-half,c[1]+half}};
    return p;
}

// ──────────────────────────────────────────────────────────────
TEST_CASE("BFS: open corridor is reachable", "[infeasibility]") {
    BoundingBox2d roi; roi.min_pt={-1,-5}; roi.max_pt={11,5};
    SDFField sdf; sdf.build(roi, {}, 0.5, 0.0);

    Polygon2d fence = makeSquarePoly(Vec2d(5,0), 8.0);
    bool ok = bfsReachable(sdf, fence, Vec2d(0,0), Vec2d(10,0), 0.2);
    REQUIRE(ok);
}

TEST_CASE("BFS: blocked by obstacle spanning full width", "[infeasibility]") {
    // Obstacle fills the entire corridor
    Polygon2d obs_poly;
    obs_poly.outer = {{3,-5},{7,-5},{7,5},{3,5}};
    Obstacle obs; obs.geometry = obs_poly;

    BoundingBox2d roi; roi.min_pt={-1,-6}; roi.max_pt={11,6};
    SDFField sdf; sdf.build(roi, {obs}, 0.5, 0.0);

    Polygon2d fence = makeSquarePoly(Vec2d(5,0), 6.0);
    bool ok = bfsReachable(sdf, fence, Vec2d(0,0), Vec2d(10,0), 0.3);
    REQUIRE_FALSE(ok);
}

TEST_CASE("detectSandwich: obstacle pressed against fence", "[infeasibility]") {
    // Obstacle touching right fence edge
    Polygon2d obs_poly;
    obs_poly.outer = {{8,-1},{10,-1},{10,1},{8,1}};
    Obstacle obs; obs.geometry = obs_poly;

    BoundingBox2d roi; roi.min_pt={-1,-6}; roi.max_pt={11,6};
    SDFField sdf; sdf.build(roi, {obs}, 0.5, 0.0);

    Polygon2d fence = makeSquarePoly(Vec2d(5,0), 5.0);  // right edge at x=10
    bool sandwich = detectSandwich(sdf, fence, 0.5);
    // Obstacle butts against fence → sandwich
    REQUIRE(sandwich);
}

TEST_CASE("preCheck: C-type topological block detected", "[infeasibility]") {
    Polygon2d obs_poly;
    obs_poly.outer = {{3,-5},{7,-5},{7,5},{3,5}};
    Obstacle obs; obs.geometry = obs_poly;

    BoundingBox2d roi; roi.min_pt={-1,-6}; roi.max_pt={11,6};
    SDFField sdf; sdf.build(roi, {obs}, 0.5, 0.0);

    Polygon2d fence = makeSquarePoly(Vec2d(5,0), 5.0);
    auto res = preCheck(sdf, fence, Vec2d(0,0), Vec2d(10,0), 3.5);
    REQUIRE(res.topological_block);
    REQUIRE(res.type == ViolationInfo::InfeasibilityType::TopologicalBlock);
}

TEST_CASE("preCheck: A-type narrow passage detected", "[infeasibility]") {
    // Small obstacle on one side leaving narrow gap
    Polygon2d obs_poly;
    obs_poly.outer = {{3,-4},{7,-4},{7,-0.5},{3,-0.5}};
    Obstacle obs; obs.geometry = obs_poly;

    BoundingBox2d roi; roi.min_pt={-1,-5}; roi.max_pt={11,5};
    SDFField sdf; sdf.build(roi, {obs}, 0.5, 0.0);

    Polygon2d fence = makeSquarePoly(Vec2d(5,0), 4.5);
    // Required width 6m (3 lanes × 2m): corridor only ~0.5m
    auto res = preCheck(sdf, fence, Vec2d(0,0), Vec2d(10,0), 6.0);
    REQUIRE(res.narrow_passage);
}

TEST_CASE("preCheck: no issue for clear corridor", "[infeasibility]") {
    BoundingBox2d roi; roi.min_pt={-1,-6}; roi.max_pt={11,6};
    SDFField sdf; sdf.build(roi, {}, 0.5, 0.0);

    Polygon2d fence = makeSquarePoly(Vec2d(5,0), 5.0);
    auto res = preCheck(sdf, fence, Vec2d(0,0), Vec2d(10,0), 3.0);
    REQUIRE_FALSE(res.topological_block);
    REQUIRE_FALSE(res.fence_sandwich);
    REQUIRE(res.type == ViolationInfo::InfeasibilityType::None);
}

TEST_CASE("makeFallbackCurve: topological block → Infeasible", "[infeasibility]") {
    PreCheckResult pre;
    pre.topological_block = true;
    pre.type = ViolationInfo::InfeasibilityType::TopologicalBlock;

    Connectivity conn;
    conn.id="X"; conn.entry_lane_id="L1"; conn.exit_lane_id="R1";
    conn.turn_type = ConnTurnType::Straight;

    auto cc = makeFallbackCurve(pre, conn, Vec2d(0,0), Vec2d(5,0));
    REQUIRE(cc.status == CurveStatus::Infeasible);
    REQUIRE(!cc.curve);
}

TEST_CASE("makeFallbackCurve: sandwich → Degraded straight line", "[infeasibility]") {
    PreCheckResult pre;
    pre.fence_sandwich = true;
    pre.type = ViolationInfo::InfeasibilityType::Sandwich;

    Connectivity conn;
    conn.id="Y"; conn.entry_lane_id="L1"; conn.exit_lane_id="R1";
    conn.turn_type = ConnTurnType::Straight;

    Vec2d p0(0,0), p1(5,0);
    auto cc = makeFallbackCurve(pre, conn, p0, p1);
    REQUIRE(cc.status == CurveStatus::Degraded);
    REQUIRE(cc.curve);

    // Straight line: endpoints match
    REQUIRE_THAT((cc.curve->startPt()-p0).norm(),
                 Catch::Matchers::WithinAbs(0.0, 0.1));
    REQUIRE_THAT((cc.curve->endPt()  -p1).norm(),
                 Catch::Matchers::WithinAbs(0.0, 0.1));
}

TEST_CASE("tryRelaxFence: expansion blocked by road edge", "[infeasibility]") {
    Polygon2d fence = makeSquarePoly(Vec2d(5,0), 5.0);  // ±5 in both dims

    // Road edge exactly at fence boundary x=10
    LineString2d edge; edge.points = {{10,-10},{10,10}};
    Boundary bnd; bnd.type = Boundary::Type::RoadEdge; bnd.geometry = edge;

    Vec2d conflict(9, 0);
    auto res = tryRelaxFence(fence, {bnd}, conflict, 1.5);
    // Expansion would cross road edge → must fail
    REQUIRE_FALSE(res.success);
}

TEST_CASE("polygonContains basic", "[fence_check]") {
    Polygon2d square = makeSquarePoly(Vec2d(0,0), 2.0);
    REQUIRE( polygonContains(square, Vec2d(0,0)));
    REQUIRE( polygonContains(square, Vec2d(1,1)));
    REQUIRE_FALSE(polygonContains(square, Vec2d(3,0)));
    REQUIRE_FALSE(polygonContains(square, Vec2d(0,5)));
}
