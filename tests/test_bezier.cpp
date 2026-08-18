#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "curve/bezier.h"
#include "curve/curve_utils.h"
#include "utils.h"

#include "io/iodata_shapefile.h"
#include "utils/shapefile.hpp"

using Catch::Matchers::WithinAbs;

using namespace isg;

// ──────────────────────────────────────────────────────────────
TEST_CASE("BezierSegment evaluate endpoints", "[bezier]") {
    BezierSegment s;
    s.ctrl[0] = Vec2d(0,0);
    s.ctrl[1] = Vec2d(1,0);
    s.ctrl[2] = Vec2d(2,1);
    s.ctrl[3] = Vec2d(3,1);

    auto p0 = s.evaluate(0.0);
    auto p1 = s.evaluate(1.0);
    REQUIRE_THAT(p0[0], WithinAbs(0.0, 1e-9));
    REQUIRE_THAT(p0[1], WithinAbs(0.0, 1e-9));
    REQUIRE_THAT(p1[0], WithinAbs(3.0, 1e-9));
    REQUIRE_THAT(p1[1], WithinAbs(1.0, 1e-9));
}

TEST_CASE("BezierSegment tangent at endpoints", "[bezier]") {
    BezierSegment s;
    s.ctrl[0] = Vec2d(0,0);
    s.ctrl[1] = Vec2d(1,0);  // tangent at P0 ∝ (1,0)
    s.ctrl[2] = Vec2d(2,1);
    s.ctrl[3] = Vec2d(3,1);

    auto t0 = s.tangent(0.0);
    REQUIRE_THAT(t0.normalized()[0], WithinAbs(1.0, 1e-6));
    REQUIRE_THAT(t0.normalized()[1], WithinAbs(0.0, 1e-6));
}

TEST_CASE("BezierSegment de Casteljau split preserves endpoints", "[bezier]") {
    BezierSegment s;
    s.ctrl[0] = Vec2d(0,0);
    s.ctrl[1] = Vec2d(1,2);
    s.ctrl[2] = Vec2d(3,2);
    s.ctrl[3] = Vec2d(4,0);

    std::pair<BezierSegment,BezierSegment> _sp = s.splitAt(0.5);
    BezierSegment L = _sp.first, R = _sp.second;
    // Split point is shared
    REQUIRE_THAT((L.ctrl[3] - R.ctrl[0]).norm(), WithinAbs(0.0, 1e-9));
    // Original endpoints preserved
    REQUIRE_THAT((L.ctrl[0] - s.ctrl[0]).norm(), WithinAbs(0.0, 1e-9));
    REQUIRE_THAT((R.ctrl[3] - s.ctrl[3]).norm(), WithinAbs(0.0, 1e-9));
}

TEST_CASE("makeCubicG1 endpoint positions and tangents", "[bezier]") {
    Vec2d p0(0,0), t0(1,0), p1(5,0), t1(1,0);
    auto seg = makeCubicG1(p0, t0, p1, t1);

    REQUIRE_THAT((seg.ctrl[0]-p0).norm(), WithinAbs(0.0, 1e-9));
    REQUIRE_THAT((seg.ctrl[3]-p1).norm(), WithinAbs(0.0, 1e-9));

    // ctrl[1] must lie along t0 from p0
    Vec2d dir1 = (seg.ctrl[1]-p0).normalized();
    REQUIRE_THAT(dir1[0], WithinAbs(1.0, 1e-6));
    REQUIRE_THAT(dir1[1], WithinAbs(0.0, 1e-6));

    // ctrl[2] must lie along -t1 from p1
    Vec2d dir2 = (p1 - seg.ctrl[2]).normalized();
    REQUIRE_THAT(dir2[0], WithinAbs(1.0, 1e-6));
    REQUIRE_THAT(dir2[1], WithinAbs(0.0, 1e-6));
}

TEST_CASE("ordinary single cubic controls stay on directional axes and inside turn intersection",
          "[bezier][constraints]") {
    Vec2d p0(0, 0), p1(10, 10);
    Vec2d t0(1, 0), t1(0, 1);
    BezierCurve c;
    c.segs.push_back(makeCubicG1(p0, t0, p1, t1, 0.9));
    c.segs.front().ctrl[1] += Vec2d(2.0, 1.0);
    c.segs.front().ctrl[2] += Vec2d(-1.0, 2.0);

    REQUIRE_FALSE(ordinarySingleCubicControlsValid(c, p0, t0, p1, t1));
    constrainOrdinarySingleCubicControls(c, p0, t0, p1, t1);

    CHECK(ordinarySingleCubicControlsValid(c, p0, t0, p1, t1));
    CHECK_THAT(c.segs.front().ctrl[1].x(), WithinAbs(10.0, 1e-6));
    CHECK_THAT(c.segs.front().ctrl[1].y(), WithinAbs(0.0, 1e-6));
    CHECK_THAT(c.segs.front().ctrl[2].x(), WithinAbs(10.0, 1e-6));
    CHECK_THAT(c.segs.front().ctrl[2].y(), WithinAbs(0.0, 1e-6));
}

TEST_CASE("ordinary straight single cubic clamps lateral and oversized handles",
          "[bezier][constraints]") {
    Vec2d p0(0, 0), p1(8, 0), t0(1, 0), t1(1, 0);
    BezierCurve c;
    c.segs.push_back(makeCubicG1(p0, t0, p1, t1, 0.5));
    c.segs.front().ctrl[1] = Vec2d(20, 4);
    c.segs.front().ctrl[2] = Vec2d(-5, -3);

    constrainOrdinarySingleCubicControls(c, p0, t0, p1, t1);
    CHECK(ordinarySingleCubicControlsValid(c, p0, t0, p1, t1));
    CHECK_THAT(c.segs.front().ctrl[1].x(), WithinAbs(8.0, 1e-6));
    CHECK_THAT(c.segs.front().ctrl[1].y(), WithinAbs(0.0, 1e-6));
    CHECK_THAT(c.segs.front().ctrl[2].x(), WithinAbs(0.0, 1e-6));
    CHECK_THAT(c.segs.front().ctrl[2].y(), WithinAbs(0.0, 1e-6));
}

TEST_CASE("BezierCurve G1 at join points", "[bezier]") {
    std::vector<Vec2d> pts = {Vec2d(0,0), Vec2d(3,2), Vec2d(6,0)};
    std::vector<Vec2d> tans = {Vec2d(1,0), Vec2d(1,0.5), Vec2d(1,0)};
    for (auto& t : tans) t.normalize();

    BezierCurve c = makeCurveFromKnots(pts, tans);
    REQUIRE(c.segs.size() == 2);

    // Join point: P3 of seg[0] == P0 of seg[1]
    REQUIRE_THAT((c.segs[0].ctrl[3]-c.segs[1].ctrl[0]).norm(),
                 WithinAbs(0.0, 1e-9));

    // G1 at join: tangent direction is continuous
    Vec2d t_end_s0 = (c.segs[0].ctrl[3]-c.segs[0].ctrl[2]).normalized();
    Vec2d t_beg_s1 = (c.segs[1].ctrl[1]-c.segs[1].ctrl[0]).normalized();
    REQUIRE_THAT(cross2d(t_end_s0, t_beg_s1), WithinAbs(0.0, 1e-4));
}

TEST_CASE("BezierCurve arcLength is positive", "[bezier]") {
    std::vector<Vec2d> pts = {Vec2d(0,0), Vec2d(5,0)};
    std::vector<Vec2d> tans = {Vec2d(1,0), Vec2d(1,0)};
    auto c = makeCurveFromKnots(pts, tans);

    std::cout << toWKT(genVectorline({pts[0][0],pts[0][1],0}, {tans[0][0],tans[0][1],0}, 1.0), "in") << std::endl;
    std::cout << toWKT(genVectorline({pts[1][0],pts[1][1],0}, {tans[1][0],tans[1][1],0}, 1.0), "out") << std::endl;
    std::cout << toWKT(toArray(c.sample(50)), "c") << std::endl;

    double len = c.arcLength();
    REQUIRE(len > 4.5);
    REQUIRE(len < 5.5);
}

TEST_CASE("curveToParams / curveFromParams round-trip", "[bezier]") {
    std::vector<Vec2d> pts = {Vec2d(0,0), Vec2d(4,3), Vec2d(8,0)};
    std::vector<Vec2d> tans = {Vec2d(1,0), Vec2d(1,0), Vec2d(1,0)};
    auto c = makeCurveFromKnots(pts, tans);

    VecXd params = curveToParams(c);
    BezierCurve c2 = curveFromParams(params, c);

    // Endpoints must be preserved
    REQUIRE_THAT((c.startPt()-c2.startPt()).norm(), WithinAbs(0.0, 1e-9));
    REQUIRE_THAT((c.endPt()  -c2.endPt()  ).norm(), WithinAbs(0.0, 1e-9));
}

TEST_CASE("localCurvature of straight line is zero", "[curve_utils]") {
    double k = localCurvature(Vec2d(0,0), Vec2d(1,0), Vec2d(2,0));
    REQUIRE_THAT(k, WithinAbs(0.0, 1e-6));
}

TEST_CASE("localCurvature of unit circle arc", "[curve_utils]") {
    // Three points on unit circle, angle step 30°
    double a = 0, b = M_PI/6, c_ = M_PI/3;
    Vec2d pa(std::cos(a),std::sin(a));
    Vec2d pb(std::cos(b),std::sin(b));
    Vec2d pc(std::cos(c_),std::sin(c_));
    double k = localCurvature(pa, pb, pc);
    REQUIRE_THAT(k, WithinAbs(1.0, 0.05));  // curvature of unit circle = 1
}

TEST_CASE("segmentsIntersect basic", "[curve_utils]") {
    Vec2d a(0,0), b(2,0), c(1,-1), d(1,1);
    Vec2d out;
    REQUIRE(segmentsIntersect(a,b,c,d,&out));
    REQUIRE_THAT(out[0], WithinAbs(1.0, 1e-6));
    REQUIRE_THAT(out[1], WithinAbs(0.0, 1e-6));
}

TEST_CASE("segmentsIntersect parallel lines do not intersect", "[curve_utils]") {
    Vec2d a(0,0), b(2,0), c(0,1), d(2,1);
    REQUIRE_FALSE(segmentsIntersect(a,b,c,d));
}
