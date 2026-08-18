#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "optimizer/sdf_field.h"

using Catch::Matchers::WithinAbs;

using namespace isg;

static Polygon2d makeSquare(Vec2d centre, double half) {
    Polygon2d p;
    p.outer = {
        {centre[0]-half, centre[1]-half},
        {centre[0]+half, centre[1]-half},
        {centre[0]+half, centre[1]+half},
        {centre[0]-half, centre[1]+half}
    };
    return p;
}

TEST_CASE("SDFField: point far from obstacle is safe", "[sdf]") {
    Obstacle obs;
    obs.geometry = makeSquare(Vec2d(0,0), 1.0);

    BoundingBox2d roi;
    roi.min_pt = Vec2d(-10,-10);
    roi.max_pt = Vec2d( 10, 10);

    SDFField sdf;
    sdf.build(roi, {obs}, 0.2, 0.0);

    std::pair<double, Vec2d> _q = sdf.queryWithGrad(Vec2d(5, 0));
    REQUIRE(_q.first > 3.0);   // ~4m from obstacle
    REQUIRE(sdf.isSafe(Vec2d(5,0), 0.1));
}

TEST_CASE("SDFField: point inside obstacle has negative SDF", "[sdf]") {
    Obstacle obs;
    obs.geometry = makeSquare(Vec2d(0,0), 2.0);

    BoundingBox2d roi;
    roi.min_pt = Vec2d(-10,-10);
    roi.max_pt = Vec2d( 10, 10);

    SDFField sdf;
    sdf.build(roi, {obs}, 0.2, 0.0);

    std::pair<double, Vec2d> _q = sdf.queryWithGrad(Vec2d(0,0));
    REQUIRE(_q.first < 0.0);
    REQUIRE_FALSE(sdf.isSafe(Vec2d(0,0), 0.1));
}

TEST_CASE("SDFField: gradient points away from obstacle", "[sdf]") {
    Obstacle obs;
    obs.geometry = makeSquare(Vec2d(0,0), 1.0);

    BoundingBox2d roi;
    roi.min_pt = Vec2d(-10,-10);
    roi.max_pt = Vec2d( 10, 10);

    SDFField sdf;
    sdf.build(roi, {obs}, 0.2, 0.0);

    // Point to the right of obstacle
    Vec2d pt(2.5, 0);
    std::pair<double, Vec2d> _q = sdf.queryWithGrad(pt);
    REQUIRE(_q.first > 0.0);
    // Gradient should point in +x direction (away from obstacle at x=1)
    REQUIRE(_q.second[0] > 0);
}

TEST_CASE("SDFField: obstacle penalty is zero far away", "[sdf]") {
    Obstacle obs;
    obs.geometry = makeSquare(Vec2d(0,0), 1.0);

    BoundingBox2d roi;
    roi.min_pt = Vec2d(-10,-10);
    roi.max_pt = Vec2d( 10, 10);

    SDFField sdf;
    sdf.build(roi, {obs}, 0.2, 0.0);

    double pen = sdf.obstaclePenalty(Vec2d(5,0), 0.5);
    REQUIRE_THAT(pen, WithinAbs(0.0, 1e-6));
}

TEST_CASE("SDFField: obstacle penalty positive when too close", "[sdf]") {
    Obstacle obs;
    obs.geometry = makeSquare(Vec2d(0,0), 1.0);

    BoundingBox2d roi;
    roi.min_pt = Vec2d(-10,-10);
    roi.max_pt = Vec2d( 10, 10);

    SDFField sdf;
    sdf.build(roi, {obs}, 0.2, 0.0);

    // Point at (1.2, 0): distance from obstacle surface ≈ 0.2m < clearance 0.5m
    double pen = sdf.obstaclePenalty(Vec2d(1.2, 0), 0.5);
    REQUIRE(pen > 0.0);
}

TEST_CASE("SDFField: buffer makes obstacle appear larger", "[sdf]") {
    Obstacle obs;
    obs.geometry = makeSquare(Vec2d(0,0), 1.0);  // 2×2m square

    BoundingBox2d roi;
    roi.min_pt = Vec2d(-10,-10);
    roi.max_pt = Vec2d( 10, 10);

    SDFField sdf_no_buf, sdf_buf;
    sdf_no_buf.build(roi, {obs}, 0.2, 0.0);
    sdf_buf.build   (roi, {obs}, 0.2, 0.5);  // 0.5m buffer

    Vec2d pt(1.3, 0);  // just outside raw obstacle
    std::pair<double, Vec2d> _qn = sdf_no_buf.queryWithGrad(pt);
    std::pair<double, Vec2d> _qb = sdf_buf.queryWithGrad(pt);

    REQUIRE(_qn.first > 0.0);   // outside raw obstacle
    REQUIRE(_qb.first < _qn.first); // buffered obstacle brings boundary closer
    REQUIRE(_qb.first < 0.0);   // inside the 0.5m avoidance buffer
}

TEST_CASE("SDFField: worldToCell / cellToWorld round-trip", "[sdf]") {
    Obstacle obs;
    obs.geometry = makeSquare(Vec2d(0,0), 1.0);

    BoundingBox2d roi;
    roi.min_pt = Vec2d(-5,-5);
    roi.max_pt = Vec2d( 5, 5);

    SDFField sdf;
    sdf.build(roi, {obs}, 0.5, 0.0);

    Vec2d world(2.3, -1.7);
    std::pair<int,int> _wc = sdf.worldToCell(world);
    Vec2d back = sdf.cellToWorld(_wc.first, _wc.second);

    // Should round-trip to within one cell
    REQUIRE((world-back).norm() < sdf.cellSize()*1.5);
}
