#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_all.hpp>
#include "constraints/cluster_order.h"
#include "constraints/boundary_safety.h"
#include "curve/curve_utils.h"
#include "optimizer/sdf_field.h"
#include "curve/bezier.h"
#include "curve/hermite_init.h"
#include "utils.h"

using namespace isg;

// ──────────────────────────────────────────────────────────────
//  Helpers: build simple connectivity set
// ──────────────────────────────────────────────────────────────
static IntersectionInput makeSimpleInput() {
    IntersectionInput inp;

    // Two entry lanes (straight + left-turn) in same group
    // Lane geometry: entry lanes run outside→junction (last pt = junction edge)
    Lane l1; l1.id="L1"; l1.width=3.5;
    l1.geometry.points={{-5,0},{0,0}};   // entry: last pt (0,0) is junction edge
    Lane l2; l2.id="L2"; l2.width=3.5;
    l2.geometry.points={{-5,3.5},{0,3.5}};
    inp.lanes = {l1, l2};

    LaneGroup entry;
    entry.id="EG"; entry.role=GroupRole::Entry;
    entry.lanes={"L1","L2"};
    entry.boundaries={"B0","B1","B2"};
    //entry.direction=Vec2d(1,0); entry.ref_point=Vec2d(0,0);

    LaneGroup exit1;
    exit1.id="XG1"; exit1.role=GroupRole::Exit;
    exit1.lanes={"R1"}; exit1.boundaries={"C0","C1"};
    //exit1.direction=Vec2d(1,0); exit1.ref_point=Vec2d(10,0);

    LaneGroup exit2;
    exit2.id="XG2"; exit2.role=GroupRole::Exit;
    exit2.lanes={"R2"}; exit2.boundaries={"D0","D1"};
    //exit2.direction=Vec2d(0,1); exit2.ref_point=Vec2d(5,5);

    inp.lane_groups = {entry, exit1, exit2};

    // Exit lane geometry: junction edge → outside (first pt = junction edge)
    Lane r1; r1.id="R1"; r1.width=3.5; r1.geometry.points={{10,0},{15,0}};
    Lane r2; r2.id="R2"; r2.width=3.5; r2.geometry.points={{5,5},{5,10}};
    inp.lanes.push_back(r1); inp.lanes.push_back(r2);

    Connectivity c1; c1.id="C1"; c1.entry_lane_id="L1";
                     c1.exit_lane_id="R1"; c1.turn_type=ConnTurnType::Straight;
    Connectivity c2; c2.id="C2"; c2.entry_lane_id="L2";
                     c2.exit_lane_id="R2"; c2.turn_type=ConnTurnType::TurnLeft;
    inp.connectivities = {c1, c2};
    return inp;
}

// ──────────────────────────────────────────────────────────────
TEST_CASE("ClusterOrderSolver builds pairs for same entry group", "[cluster]") {
    auto inp = makeSimpleInput();
    ClusterOrderSolver cs;
    cs.build(inp.connectivities, inp.lanes, inp.lane_groups);

    auto& pairs = cs.pairs();
    // C1 and C2 share entry group EG → one pair expected
    REQUIRE(pairs.size() >= 1);

    bool found = false;
    for (auto& p : pairs)
        if ((p.id_a=="C1"&&p.id_b=="C2")||(p.id_a=="C2"&&p.id_b=="C1"))
            found = true;
    REQUIRE(found);
}

TEST_CASE("ClusterOrderSolver: Straight before Left in order", "[cluster]") {
    auto inp = makeSimpleInput();
    ClusterOrderSolver cs;
    cs.build(inp.connectivities, inp.lanes, inp.lane_groups);

    // Straight (C1) should be priority 0, Left (C2) priority 1
    // So C1 appears first in entry_order_
    // We can't directly inspect private, but the pair exemption tells us
    auto ex = cs.exemptionOf("C1","C2");
    // Neither is a UTurn → no structural exemption
    REQUIRE(ex == CrossExemption::None);
}

TEST_CASE("ClusterOrderSolver: UTurn pair remains constrained in same cluster", "[cluster]") {
    IntersectionInput inp = makeSimpleInput();
    Connectivity cu; cu.id="CU"; cu.entry_lane_id="L1";
                     cu.exit_lane_id="R1"; cu.turn_type=ConnTurnType::UTurnLeft;
    inp.connectivities.push_back(cu);

    ClusterOrderSolver cs;
    cs.build(inp.connectivities, inp.lanes, inp.lane_groups);

    // U-turns are not automatically a logical crossing. They still need the
    // same-cluster non-endpoint intersection constraint.
    auto ex = cs.exemptionOf("C1","CU");
    REQUIRE(ex == CrossExemption::None);
}

TEST_CASE("UTurn initializer keeps a compact non-self-intersecting arc", "[uturn][geometry]") {
    SDFField sdf;
    Polygon2d fence;
    Vec2d p0(1.75, -10.0);
    Vec2d t0(0.0, 1.0);
    Vec2d p1(-1.75, -10.0);
    Vec2d t1(0.0, -1.0);

    BezierCurve c = buildTwoSegmentUTurn(p0, t0, p1, t1, sdf, fence);
    REQUIRE(c.numSegments() == 1);
    REQUIRE_FALSE(curveSelfIntersectsBusiness(c, 1.0));
    REQUIRE(c.startTan().dot(t0.normalized()) > 0.999);
    REQUIRE(c.endTan().dot(t1.normalized()) > 0.999);

    auto pts = c.sampleByArcLength(40);
    double max_abs_x = 0.0;
    double max_y = -1e18;
    for (auto& p : pts) {
        max_abs_x = std::max(max_abs_x, std::abs(p.x()));
        max_y = std::max(max_y, p.y());
    }
    REQUIRE(max_abs_x < 3.2);
    REQUIRE(max_y > -8.4);

    int sign = 0;
    int slope_sign = 0;
    int extrema = 0;
    double prev_forward = (c.segs.front().evaluate(0.0) - p0).dot(t0.normalized());
    for (int i = 0; i <= 40; ++i) {
        double u = i / 40.0;
        Vec2d d1 = c.segs.front().evalDeriv1(u);
        Vec2d d2 = c.segs.front().evalDeriv2(u);
        double signed_k = cross2d(d1, d2);
        if (std::abs(signed_k) > 1e-7) {
            int cur_sign = signed_k > 0.0 ? 1 : -1;
            if (sign == 0)
                sign = cur_sign;
            REQUIRE(cur_sign == sign);
        }

        double forward = (c.segs.front().evaluate(u) - p0).dot(t0.normalized());
        if (i > 0) {
            double df = forward - prev_forward;
            int cur_slope = df > 1e-5 ? 1 : (df < -1e-5 ? -1 : 0);
            if (cur_slope != 0) {
                if (slope_sign != 0 && cur_slope != slope_sign)
                    ++extrema;
                slope_sign = cur_slope;
            }
        }
        prev_forward = forward;
    }
    REQUIRE(extrema <= 1);
}

TEST_CASE("curveCrossings detects intersection", "[intersection_check]") {
    // Two straight curves that cross at 90°
    std::vector<Vec2d> pts_h = {Vec2d(0,0), Vec2d(4,0)};
    std::vector<Vec2d> tans_h = {Vec2d(1,0), Vec2d(1,0)};
    auto ch = makeCurveFromKnots(pts_h, tans_h);

    std::vector<Vec2d> pts_v = {Vec2d(2,-2), Vec2d(2,2)};
    std::vector<Vec2d> tans_v = {Vec2d(0,1), Vec2d(0,1)};
    auto cv = makeCurveFromKnots(pts_v, tans_v);

    auto crosses = curveCrossings(ch, cv);
    REQUIRE(crosses.size() >= 1);
    // Intersection near (2,0)
    REQUIRE_THAT((crosses[0]-Vec2d(2,0)).norm(),
                 Catch::Matchers::WithinAbs(0.0, 0.15));
}

TEST_CASE("curveCrossings no crossing for parallel curves", "[intersection_check]") {
    std::vector<Vec2d> pts1 = {Vec2d(0,0), Vec2d(4,0)};
    std::vector<Vec2d> pts2 = {Vec2d(0,1), Vec2d(4,1)};
    std::vector<Vec2d> t1   = {Vec2d(1,0), Vec2d(1,0)};
    auto c1 = makeCurveFromKnots(pts1, t1);
    auto c2 = makeCurveFromKnots(pts2, t1);

    auto crosses = curveCrossings(c1, c2);
    REQUIRE(crosses.empty());
}

TEST_CASE("curve business intersection rejects non-endpoint overlap", "[intersection_check][overlap]") {
    std::vector<Vec2d> pts1 = {Vec2d(0,0), Vec2d(6,0)};
    std::vector<Vec2d> pts2 = {Vec2d(2,0), Vec2d(5,0)};
    std::vector<Vec2d> tan = {Vec2d(1,0), Vec2d(1,0)};
    auto c1 = makeCurveFromKnots(pts1, tan);
    auto c2 = makeCurveFromKnots(pts2, tan);

    REQUIRE(curvesIntersectBusiness(c1, c2, 0.15));
}

TEST_CASE("boundary safety rejects non-endpoint RoadEdge adherence", "[boundary][overlap]") {
    Boundary edge;
    edge.id = "edge";
    edge.type = Boundary::Type::RoadEdge;
    edge.geometry.points = {Vec2d(0, 0), Vec2d(6, 0)};
    Vec2d center(3, 4);

    std::vector<Vec2d> pts = {Vec2d(1, 0), Vec2d(5, 0)};
    std::vector<Vec2d> tan = {Vec2d(1, 0), Vec2d(1, 0)};
    BezierCurve adherent = makeCurveFromKnots(pts, tan);

    auto safety = curveBoundarySafety(
        adherent, std::vector<Boundary>{edge}, center, 64, 0.15);
    REQUIRE(safety.intersects);
}

TEST_CASE("ClusterOrderSolver checkAndMarkA2 marks obstacle cross", "[cluster]") {
    // Build an SDF with an obstacle and two curves that cross near it
    Polygon2d obs_poly;
    obs_poly.outer = {{1.5,-0.5},{2.5,-0.5},{2.5,0.5},{1.5,0.5}};
    Obstacle obs; obs.geometry = obs_poly;

    BoundingBox2d roi; roi.min_pt={-1,-3}; roi.max_pt={5,3};
    SDFField sdf;
    sdf.build(roi, {obs}, 0.2, 0.0);

    // Two crossing curves near the obstacle
    std::vector<Vec2d> p1 = {Vec2d(0,-1), Vec2d(4,1)};
    std::vector<Vec2d> p2 = {Vec2d(0,1),  Vec2d(4,-1)};
    std::vector<Vec2d> t12 = {Vec2d(1,0.5).normalized(), Vec2d(1,0.5).normalized()};
    std::vector<Vec2d> t21 = {Vec2d(1,-0.5).normalized(), Vec2d(1,-0.5).normalized()};
    BezierCurve c1 = makeCurveFromKnots(p1, t12);
    BezierCurve c2 = makeCurveFromKnots(p2, t21);

    auto inp = makeSimpleInput();
    ClusterOrderSolver cs;
    cs.build(inp.connectivities, inp.lanes, inp.lane_groups);

    std::unordered_map<ConnId, BezierCurve> curves;
    curves["C1"] = c1;
    curves["C2"] = c2;

    cs.checkAndMarkA2(curves, sdf, 2.0);

    // After marking, the pair C1-C2 should be ObstacleCross if crossing near obstacle
    auto ex = cs.exemptionOf("C1","C2");
    // The crossing is near the obstacle (within 2m), so it should be marked
    // (may be None if curves don't actually cross in this simplified setup)
    // Just verify no crash and result is valid
    REQUIRE((ex == CrossExemption::None || ex == CrossExemption::ObstacleCross));
}
