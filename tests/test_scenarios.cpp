// ═════════════════════════════════════════════════════════════════════════════
//  test_scenarios.cpp
//  Comprehensive intersection scenario tests:
//    - 井字(4-way), T字(3-way), Y字(3-way angled), 环岛(roundabout)
//    - Obstacle variants: none / central / narrow passage / sandwich / topo-block
//    - Connectivity: 1:1, 1:N, N:1; all turn types incl. U-turn
// ═════════════════════════════════════════════════════════════════════════════
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "intersection_shape_generator.h"
#include "constraints/fence_check.h"
#include "scenario_builder.h"
#include "optimizer/sdf_field.h"
#include <cmath>
#include <vector>
#include <iostream>

#include "io/iodata_shapefile.h"
#include "io/iodata_json.h"

using Catch::Matchers::WithinAbs;

using namespace isg;

//#define TEST_4WAY
//#define TEST_T
//#define TEST_Y
//#define TEST_ROUNDABOUT
#define TEST_OBSTACLE
//#define TEST_FILE

// ─────────────────────────────────────────────────────────────────────────────
//  Common validation helpers
// ─────────────────────────────────────────────────────────────────────────────
static void checkCurves(const IntersectionOutput& out,
                        int expected_count,
                        const std::string& tag = "")
{
    INFO("Checking curves for: " << tag);
    REQUIRE((int)out.connectivity_curves.size() == expected_count);
    int feasible = 0;
    for (auto& cc : out.connectivity_curves) {
        if (cc.status == CurveStatus::Infeasible) continue;
        feasible++;
        REQUIRE(cc.curve);
        // G1: start/end points correspond to lane geometry endpoints
        // (we only check curve is non-degenerate)
        double len = cc.curve->arcLength();
        REQUIRE(len > 0.1);
    }
    // At least one feasible curve expected (unless all declared infeasible)
    (void)feasible;
}

static bool curveAvoidsPenetration(const ConnectivityCurve& cc,
                                    const IntersectionInput& inp)
{
    if (!cc.curve) return true;
    SDFField sdf;
    auto roi = inp.area.geometry.bbox();
    if (roi.empty()) return true;
    sdf.build(roi, inp.obstacles, 0.2, 0.0);
    for (auto& pt : cc.curve->sampleByArcLength(50)) {
        std::pair<double,Vec2d> _q = sdf.queryWithGrad(pt);
        if (_q.first < -0.25) return false;   // >25cm penetration = failure
    }
    return true;
}

static bool allCurvesAvoidObstacles(const IntersectionOutput& out,
                                     const IntersectionInput& inp)
{
    for (auto& cc : out.connectivity_curves) {
        if (!curveAvoidsPenetration(cc, inp)) return false;
    }
    return true;
}

static const ConnectivityCurve* findCurveById(
    const IntersectionOutput& out, const ConnId& id)
{
    for (auto& cc : out.connectivity_curves) {
        if (cc.id == id)
            return &cc;
    }
    return nullptr;
}

static double minRawObstacleDistance(
    const ConnectivityCurve& cc, const IntersectionInput& inp)
{
    if (!cc.curve || inp.obstacles.empty())
        return 1e18;
    SDFField sdf;
    auto roi = inp.area.geometry.bbox();
    if (roi.empty())
        return 1e18;
    sdf.build(roi, inp.obstacles, 0.2, 0.0);
    double min_d = 1e18;
    for (auto& pt : cc.curve->sampleByArcLength(80)) {
        std::pair<double, Vec2d> _q = sdf.queryWithGrad(pt);
        min_d = std::min(min_d, _q.first);
    }
    return min_d;
}

// ─────────────────────────────────────────────────────────────────────────────
//  Reset global connectivity counter between tests
// ─────────────────────────────────────────────────────────────────────────────
namespace { struct ResetConn { ResetConn(){gConnIdx=0;} } _reset_; }


// ═════════════════════════════════════════════════════════════════════════════
//  SECTION 1 — 井字路口 (4-way cross)
// ═════════════════════════════════════════════════════════════════════════════

// Build a standard 4-arm input.
// Arms: N(north entry from south), S, E, W
// Returns arm lane IDs: arms[arm][0]=entry_ids, arms[arm][1]=exit_ids
static IntersectionInput build4Way(int lanes_per_arm = 1,
                                    double lane_w = 3.5)
{
    gConnIdx = 0;
    IntersectionInput inp;
    addStandardArea(inp, 14.0);

    // N: inbound from south → north, junction pt at (0, -10)
    ArmDef arms[4];
    arms[0] = {"N", {0,-10}, {0, 1}, lanes_per_arm, lanes_per_arm, lane_w};
    arms[1] = {"S", {0, 10}, {0,-1}, lanes_per_arm, lanes_per_arm, lane_w};
    arms[2] = {"E", {-10,0}, {1, 0}, lanes_per_arm, lanes_per_arm, lane_w};
    arms[3] = {"W", { 10,0}, {-1,0}, lanes_per_arm, lanes_per_arm, lane_w};

    std::vector<LaneId> entry[4], exitl[4];
    for (int a = 0; a < 4; ++a) {
        std::pair<std::vector<LaneId>,std::vector<LaneId>> _arm = buildArm(arms[a], inp);
        entry[a] = _arm.first; exitl[a] = _arm.second;
    }

    // Road-edge boundaries (top/bottom/left/right of fence)
    addRoadEdgeBoundary(inp,"TOP",  {-14,14},{14,14});
    addRoadEdgeBoundary(inp,"BOT",  {-14,-14},{14,-14});
    addRoadEdgeBoundary(inp,"LEFT", {-14,-14},{-14,14});
    addRoadEdgeBoundary(inp,"RIGHT",{14,-14},{14,14});

    // N→S straight, N→E left, N→W right (for single-lane arm)
    // Repeat for each arm; full 12 connectivity 1:1
    struct ConnSpec { int from,to; ConnTurnType tt; };
    std::vector<ConnSpec> specs = {
        // N inbound
        {0,1,ConnTurnType::Straight},{0,2,ConnTurnType::TurnLeft},{0,3,ConnTurnType::TurnRight},
        // S inbound
        {1,0,ConnTurnType::Straight},{1,3,ConnTurnType::TurnLeft},{1,2,ConnTurnType::TurnRight},
        // E inbound
        {2,3,ConnTurnType::Straight},{2,1,ConnTurnType::TurnLeft},{2,0,ConnTurnType::TurnRight},
        // W inbound
        {3,2,ConnTurnType::Straight},{3,0,ConnTurnType::TurnLeft},{3,1,ConnTurnType::TurnRight},
    };
    for (auto& s : specs)
        addConn(inp, entry[s.from][0], exitl[s.to][0], s.tt);

    return inp;
}
#ifdef TEST_4WAY
TEST_CASE("4-way: clear no obstacles — all 12 connections feasible", "[4way][clear]")
{
    std::cout << "TEST_CASE: " << "4-way: clear no obstacles — all 12 connections feasible" << std::endl;

    auto inp = build4Way();
    IntersectionShapeGenerator gen;
    IntersectionOutput out;
    REQUIRE(gen.generate(inp, out));

    save(inp, "./4way-clear/input/", "");
    save(out, "./4way-clear/output/", "");

    checkCurves(out, 12, "4-way clear");
    int ok = 0;
    for (auto& cc : out.connectivity_curves)
        if (cc.status != CurveStatus::Infeasible) ok++;
    REQUIRE(ok >= 10); // at least 10 of 12 feasible
    REQUIRE(allCurvesAvoidObstacles(out, inp));
}

TEST_CASE("4-way: central square obstacle — curves reroute around it", "[4way][obstacle]")
{
    std::cout << "TEST_CASE: " << "4-way: central square obstacle — curves reroute around it" << std::endl;

    auto inp = build4Way();
    addObstacle(inp, "CENTRE", makeRectObstacle({0,0}, 2.0, 2.0));
    IntersectionShapeGenerator gen;
    IntersectionOutput out;
    REQUIRE(gen.generate(inp, out));

    save(inp, "./4way-obstacle/input/", "");
    save(out, "./4way-obstacle/output/", "");

    REQUIRE((int)out.connectivity_curves.size() == 12);
    REQUIRE(allCurvesAvoidObstacles(out, inp));
}

TEST_CASE("4-way: multi-lane (2 lanes per arm) — 1:1 connectivity", "[4way][multilane]")
{
    std::cout << "TEST_CASE: " << "4-way: multi-lane (2 lanes per arm) — 1:1 connectivity" << std::endl;

    gConnIdx = 0;
    IntersectionInput inp;
    addStandardArea(inp, 18.0);
    ArmDef arms[4] = {
        {"N",{0,-12},{0, 1},2,2,3.5},{"S",{0, 12},{0,-1},2,2,3.5},
        {"E",{-12,0},{1, 0},2,2,3.5},{"W",{12, 0},{-1,0},2,2,3.5}
    };
    std::vector<LaneId> entry[4], exitl[4];
    for (int a=0;a<4;++a){ std::pair<std::vector<LaneId>,std::vector<LaneId>> _arm=buildArm(arms[a],inp); entry[a]=_arm.first; exitl[a]=_arm.second; }

    // Per-arm: lane[0]=inner (left-turn dedicated), lane[1]=outer (straight/right)
    // N→S straight (outer-outer), N→E left (inner-inner), N→W right (outer-outer)
    // 8 connections per quadrant × 4 arms = 16 total (simplified)
    for (int a=0;a<4;++a) {
        int b=(a+1)%4, c=(a+2)%4, d=(a+3)%4;
        addConn(inp, entry[a][0], exitl[b][0], ConnTurnType::TurnLeft);
        addConn(inp, entry[a][1], exitl[c][1], ConnTurnType::Straight);
        addConn(inp, entry[a][1], exitl[d][1], ConnTurnType::TurnRight);
    }
    IntersectionShapeGenerator gen;
    IntersectionOutput out;
    REQUIRE(gen.generate(inp, out));

    save(inp, "./4way-11/input/", "");
    save(out, "./4way-11/output/", "");

    REQUIRE((int)out.connectivity_curves.size() == 12);
    REQUIRE(allCurvesAvoidObstacles(out, inp));
}

TEST_CASE("4-way: 1:N connectivity — one entry to two exits", "[4way][1toN]")
{
    std::cout << "TEST_CASE: " << "4-way: 1:N connectivity — one entry to two exits" << std::endl;

    gConnIdx = 0;
    IntersectionInput inp;
    addStandardArea(inp, 18.0);
    ArmDef arms[4] = {
        {"N",{0,-12},{0,1},1,1,3.5},{"S",{0,12},{0,-1},1,1,3.5},
        {"E",{-12,0},{1,0},1,1,3.5},{"W",{12,0},{-1,0},1,1,3.5}
    };
    std::vector<LaneId> entry[4], exitl[4];
    for (int a=0;a<4;++a){ std::pair<std::vector<LaneId>,std::vector<LaneId>> _arm=buildArm(arms[a],inp); entry[a]=_arm.first; exitl[a]=_arm.second; }

    // N lane splits to both S (straight) and E (left)
    addConn(inp, entry[0][0], exitl[1][0], ConnTurnType::Straight);
    addConn(inp, entry[0][0], exitl[2][0], ConnTurnType::TurnLeft);
    // S lane splits to N and W
    addConn(inp, entry[1][0], exitl[0][0], ConnTurnType::Straight);
    addConn(inp, entry[1][0], exitl[3][0], ConnTurnType::TurnLeft);
    // E,W: straight only
    addConn(inp, entry[2][0], exitl[3][0], ConnTurnType::Straight);
    addConn(inp, entry[3][0], exitl[2][0], ConnTurnType::Straight);

    IntersectionShapeGenerator gen;
    IntersectionOutput out;
    REQUIRE(gen.generate(inp, out));

    save(inp, "./4way-1N/input/", "");
    save(out, "./4way-1N/output/", "");

    REQUIRE((int)out.connectivity_curves.size() == 6);
    REQUIRE(allCurvesAvoidObstacles(out, inp));
}

TEST_CASE("4-way: N:1 connectivity — two entries merge to one exit", "[4way][Nto1]")
{
    std::cout << "TEST_CASE: " << "4-way: N:1 connectivity — two entries merge to one exit" << std::endl;

    gConnIdx = 0;
    IntersectionInput inp;
    addStandardArea(inp, 18.0);
    ArmDef arms[4] = {
        {"N",{0,-12},{0,1},2,1,3.5},{"S",{0,12},{0,-1},2,1,3.5},
        {"E",{-12,0},{1,0},2,1,3.5},{"W",{12,0},{-1,0},2,1,3.5}
    };
    std::vector<LaneId> entry[4], exitl[4];
    for (int a=0;a<4;++a){ std::pair<std::vector<LaneId>,std::vector<LaneId>> _arm=buildArm(arms[a],inp); entry[a]=_arm.first; exitl[a]=_arm.second; }

    // N lane[0] + N lane[1] → S exit[0] (merge: straight + straight-inner)
    addConn(inp, entry[0][0], exitl[1][0], ConnTurnType::Straight);
    addConn(inp, entry[0][1], exitl[1][0], ConnTurnType::Straight);
    // S → N
    addConn(inp, entry[1][0], exitl[0][0], ConnTurnType::Straight);
    addConn(inp, entry[1][1], exitl[0][0], ConnTurnType::Straight);
    // E → W, W → E
    addConn(inp, entry[2][0], exitl[3][0], ConnTurnType::Straight);
    addConn(inp, entry[3][0], exitl[2][0], ConnTurnType::Straight);

    IntersectionShapeGenerator gen;
    IntersectionOutput out;
    REQUIRE(gen.generate(inp, out));

    save(inp, "./4way-N1/input/", "");
    save(out, "./4way-N1/output/", "");

    REQUIRE((int)out.connectivity_curves.size() == 6);
}

TEST_CASE("4-way: U-turn connections", "[4way][uturn]")
{
    std::cout << "TEST_CASE: " << "4-way: U-turn connections" << std::endl;

    gConnIdx = 0;
    IntersectionInput inp;
    addStandardArea(inp, 16.0);
    ArmDef arms[4] = {
        {"N",{0,-10},{0,1},2,2,3.5},{"S",{0,10},{0,-1},2,2,3.5},
        {"E",{-10,0},{1,0},2,2,3.5},{"W",{10,0},{-1,0},2,2,3.5}
    };
    std::vector<LaneId> entry[4], exitl[4];
    for (int a=0;a<4;++a){ std::pair<std::vector<LaneId>,std::vector<LaneId>> _arm=buildArm(arms[a],inp); entry[a]=_arm.first; exitl[a]=_arm.second; }

    for (int a=0;a<4;++a) {
        int opp=(a+2)%4;
        // Inner lane → U-turn back to same arm exit inner lane
        addConn(inp, entry[a][0], exitl[a][0], ConnTurnType::UTurnLeft);
        // Outer lane → straight to opposite
        addConn(inp, entry[a][1], exitl[opp][1], ConnTurnType::Straight);
    }
    IntersectionShapeGenerator gen;
    IntersectionOutput out;
    REQUIRE(gen.generate(inp, out));

    save(inp, "./4way-uturn/input/", "");
    save(out, "./4way-uturn/output/", "");

    REQUIRE((int)out.connectivity_curves.size() == 8);
    REQUIRE(allCurvesAvoidObstacles(out, inp));
}

TEST_CASE("4-way: topo-block — obstacle spans full width → Infeasible", "[4way][infeasible][topo_block]")
{
    std::cout << "TEST_CASE: " << "4-way: topo-block — obstacle spans full width → Infeasible" << std::endl;

    auto inp = build4Way();
    // Obstacle blocks N→S corridor completely
    addObstacle(inp, "BLOCKER", makeRectObstacle({0,0}, 4.0, 2.5));
    IntersectionShapeGenerator gen;
    IntersectionOutput out;
    REQUIRE(gen.generate(inp, out));

    save(inp, "./4way-topo-block/input/", "");
    save(out, "./4way-topo-block/output/", "");

    // N→S and S→N should be Infeasible; others may succeed
    int infeasible = 0;
    for (auto& cc : out.connectivity_curves)
        if (cc.status == CurveStatus::Infeasible) infeasible++;
    REQUIRE(infeasible >= 1);
}

TEST_CASE("4-way: narrow passage — only one lane fits through gap", "[4way][narrow]")
{
    std::cout << "TEST_CASE: " << "4-way: narrow passage — only one lane fits through gap" << std::endl;

    auto inp = build4Way();
    // Two obstacles leaving ~3m gap (just enough for one lane)
    addObstacle(inp, "OBS_L", makeRectObstacle({-6,0}, 2.0, 8.0));
    addObstacle(inp, "OBS_R", makeRectObstacle({ 6,0}, 2.0, 8.0));
    IntersectionShapeGenerator gen;
    IntersectionOutput out;
    REQUIRE(gen.generate(inp, out));

    save(inp, "./4way-narrow/input/", "");
    save(out, "./4way-narrow/output/", "");

    // Pipeline must not crash; some curves may be degraded/warn
    REQUIRE((int)out.connectivity_curves.size() == 12);
}

TEST_CASE("4-way: fence-obstacle sandwich on east side", "[4way][sandwich]")
{
    std::cout << "TEST_CASE: " << "4-way: fence-obstacle sandwich on east side" << std::endl;

    auto inp = build4Way();
    // Obstacle pressed against right fence edge
    addObstacle(inp, "SANDWICH", makeRectObstacle({6,0}, 2, 3.0));
    IntersectionShapeGenerator gen;
    IntersectionOutput out;
    REQUIRE(gen.generate(inp, out));

    save(inp, "./4way-sandwich/input/", "");
    save(out, "./4way-sandwich/output/", "");

    REQUIRE((int)out.connectivity_curves.size() == 12);
}
#endif

// ═════════════════════════════════════════════════════════════════════════════
//  SECTION 2 — T字路口 (3-way / T-junction)
// ═════════════════════════════════════════════════════════════════════════════

static IntersectionInput buildTJunction(int lanes_per_arm = 1, double lane_w = 3.5)
{
    gConnIdx = 0;
    IntersectionInput inp;
    addStandardArea(inp, 14.0);

    // Main road: E ↔ W; branch: N (only enters from south)
    ArmDef arms[3] = {
        {"N", {0,-10}, {0, 1}, lanes_per_arm, lanes_per_arm, lane_w},  // branch
        {"E", {-10,0}, {1, 0}, lanes_per_arm, lanes_per_arm, lane_w},  // main east
        {"W", { 10,0}, {-1,0}, lanes_per_arm, lanes_per_arm, lane_w},  // main west
    };
    std::vector<LaneId> entry[3], exitl[3];
    for (int a=0;a<3;++a){ std::pair<std::vector<LaneId>,std::vector<LaneId>> _arm=buildArm(arms[a],inp); entry[a]=_arm.first; exitl[a]=_arm.second; }

    addRoadEdgeBoundary(inp,"TOP", {-14,14},{14,14});
    addRoadEdgeBoundary(inp,"BOT", {-14,-14},{14,-14});

    // 6 connections: each arm → each other arm
    struct CS{int f,t;ConnTurnType tt;};
    std::vector<CS> specs = {
        {0,1,ConnTurnType::TurnLeft},{0,2,ConnTurnType::TurnRight},     // N→E, N→W
        {1,2,ConnTurnType::Straight},{1,0,ConnTurnType::TurnRight},     // E→W, E→N
        {2,1,ConnTurnType::Straight},{2,0,ConnTurnType::TurnLeft},      // W→E, W→N
    };
    for (auto& s:specs) addConn(inp, entry[s.f][0], exitl[s.t][0], s.tt);
    return inp;
}
#ifdef TEST_T
TEST_CASE("T-junction: clear, all 6 connections", "[T][clear]")
{
    std::cout << "TEST_CASE: " << "T-junction: clear, all 6 connections" << std::endl;

    auto inp = buildTJunction();
    IntersectionShapeGenerator gen;
    IntersectionOutput out;
    REQUIRE(gen.generate(inp, out));

    save(inp, "./T-clear/input/", "");
    save(out, "./T-clear/output/", "");

    checkCurves(out, 6, "T clear");
    REQUIRE(allCurvesAvoidObstacles(out, inp));
}

TEST_CASE("T-junction: central obstacle forcing reroute", "[T][obstacle]")
{
    std::cout << "TEST_CASE: " << "T-junction: central obstacle forcing reroute" << std::endl;

    auto inp = buildTJunction();
    addObstacle(inp, "C", makeRectObstacle({-2,0}, 1.5, 1.5));
    IntersectionShapeGenerator gen;
    IntersectionOutput out;
    REQUIRE(gen.generate(inp, out));

    save(inp, "./T-obstacle/input/", "");
    save(out, "./T-obstacle/output/", "");

    REQUIRE((int)out.connectivity_curves.size() == 6);
    REQUIRE(allCurvesAvoidObstacles(out, inp));
}

TEST_CASE("T-junction: 2-lane arms with 1:N split at branch", "[T][1toN][multilane]")
{
    std::cout << "TEST_CASE: " << "T-junction: 2-lane arms with 1:N split at branch" << std::endl;

    gConnIdx = 0;
    IntersectionInput inp;
    addStandardArea(inp, 16.0);
    ArmDef arms[3] = {
        {"N",{0,-10},{0,1},1,2,3.5},  // branch: 1 entry, 2 exit
        {"E",{-10,0},{1,0},2,2,3.5},  // main east
        {"W",{10, 0},{-1,0},2,2,3.5}, // main west
    };
    std::vector<LaneId> entry[3], exitl[3];
    for (int a=0;a<3;++a){ std::pair<std::vector<LaneId>,std::vector<LaneId>> _arm=buildArm(arms[a],inp); entry[a]=_arm.first; exitl[a]=_arm.second; }

    // Branch N → left+right (1:2)
    addConn(inp, entry[0][0], exitl[1][0], ConnTurnType::TurnLeft);
    addConn(inp, entry[0][0], exitl[2][0], ConnTurnType::TurnRight);
    // E → W (both lanes straight), E → N (inner lane only)
    addConn(inp, entry[1][0], exitl[2][0], ConnTurnType::Straight);
    addConn(inp, entry[1][1], exitl[2][1], ConnTurnType::Straight);
    addConn(inp, entry[1][0], exitl[0][0], ConnTurnType::TurnRight);
    // W → E, W → N
    addConn(inp, entry[2][0], exitl[1][0], ConnTurnType::Straight);
    addConn(inp, entry[2][1], exitl[1][1], ConnTurnType::Straight);
    addConn(inp, entry[2][0], exitl[0][1], ConnTurnType::TurnLeft);

    IntersectionShapeGenerator gen;
    IntersectionOutput out;
    REQUIRE(gen.generate(inp, out));

    save(inp, "./T-1N/input/", "");
    save(out, "./T-1N/output/", "");

    REQUIRE((int)out.connectivity_curves.size() == 8);
    REQUIRE(allCurvesAvoidObstacles(out, inp));
}

TEST_CASE("T-junction: topo-block on branch arm", "[T][topo_block]")
{
    std::cout << "TEST_CASE: " << "T-junction: topo-block on branch arm" << std::endl;

    auto inp = buildTJunction();
    // Block the branch corridor completely
    addObstacle(inp, "BLOCK", makeRectObstacle({0,-4}, 5.0, 2.5));
    IntersectionShapeGenerator gen;
    IntersectionOutput out;
    REQUIRE(gen.generate(inp, out));

    save(inp, "./T-topoblock/input/", "");
    save(out, "./T-topoblock/output/", "");

    int infeasible=0;
    for (auto& cc : out.connectivity_curves)
        if (cc.status==CurveStatus::Infeasible) infeasible++;
    REQUIRE(infeasible >= 2); // N→E and N→W both blocked
}
#endif

// ═════════════════════════════════════════════════════════════════════════════
//  SECTION 3 — Y字路口 (3-way Y-junction with angled arms)
// ═════════════════════════════════════════════════════════════════════════════

static IntersectionInput buildYJunction(double lane_w = 3.5)
{
    gConnIdx = 0;
    IntersectionInput inp;
    addStandardArea(inp, 14.0);

    // Three arms at ~120° apart
    double a0 = M_PI/2;               // arm A: points north (270° = up)
    double a1 = a0 + 2*M_PI/3;        // arm B: points lower-left
    double a2 = a0 - 2*M_PI/3;        // arm C: points lower-right
    double R  = 10.0;                  // junction radius

    auto mkArm = [&](const std::string& nm, double ang) -> ArmDef {
        Vec2d dir(-std::sin(ang), std::cos(ang));   // pointing INTO junction
        Vec2d jpt = -dir * R;                        // junction pt at edge
        return {nm, jpt, dir, 1, 1, lane_w};
    };
    ArmDef arms[3] = {mkArm("A",a0), mkArm("B",a1), mkArm("C",a2)};
    std::vector<LaneId> entry[3], exitl[3];
    for (int a=0;a<3;++a){ std::pair<std::vector<LaneId>,std::vector<LaneId>> _arm=buildArm(arms[a],inp); entry[a]=_arm.first; exitl[a]=_arm.second; }

    // 6 connections
    struct CS{int f,t;ConnTurnType tt;};
    std::vector<CS> specs={
        {0,1,ConnTurnType::TurnLeft},{0,2,ConnTurnType::TurnRight},
        {1,0,ConnTurnType::TurnRight},{1,2,ConnTurnType::TurnLeft},
        {2,0,ConnTurnType::TurnLeft},{2,1,ConnTurnType::TurnRight},
    };
    for (auto& s:specs) addConn(inp, entry[s.f][0], exitl[s.t][0], s.tt);
    return inp;
}
#ifdef TEST_Y
TEST_CASE("Y-junction: clear angled arms — all 6 connections", "[Y][clear]")
{
    std::cout << "TEST_CASE: " << "Y-junction: clear angled arms — all 6 connections" << std::endl;

    auto inp = buildYJunction();
    IntersectionShapeGenerator gen;
    IntersectionOutput out;
    REQUIRE(gen.generate(inp, out));

    save(inp, "./Y-clear/input/", "");
    save(out, "./Y-clear/output/", "");

    checkCurves(out, 6, "Y clear");
    REQUIRE(allCurvesAvoidObstacles(out, inp));
}

TEST_CASE("Y-junction: central circular obstacle", "[Y][obstacle]")
{
    std::cout << "TEST_CASE: " << "Y-junction: central circular obstacle" << std::endl;

    auto inp = buildYJunction();
    addObstacle(inp, "CIRC", makeCircleObstacle({0,0}, 2.0));
    IntersectionShapeGenerator gen;
    IntersectionOutput out;
    REQUIRE(gen.generate(inp, out));

    save(inp, "./Y-obstacle/input/", "");
    save(out, "./Y-obstacle/output/", "");

    REQUIRE((int)out.connectivity_curves.size() == 6);
    REQUIRE(allCurvesAvoidObstacles(out, inp));
}

TEST_CASE("Y-junction: 1:N — one arm splits to both others", "[Y][1toN]")
{
    std::cout << "TEST_CASE: " << "Y-junction: 1:N — one arm splits to both others" << std::endl;

    gConnIdx = 0;
    IntersectionInput inp;
    addStandardArea(inp, 14.0);
    double R=10.0, lw=3.5;
    double a0=M_PI/2, a1=a0+2*M_PI/3, a2=a0-2*M_PI/3;
    auto mkArm=[&](const std::string& nm,double ang)->ArmDef{
        Vec2d dir(-std::sin(ang),std::cos(ang));
        return{nm,-dir*R,dir,1,1,lw};
    };
    ArmDef arms[3]={mkArm("A",a0),mkArm("B",a1),mkArm("C",a2)};
    std::vector<LaneId> entry[3],exitl[3];
    for(int a=0;a<3;++a){std::pair<std::vector<LaneId>,std::vector<LaneId>> _arm=buildArm(arms[a],inp);entry[a]=_arm.first;exitl[a]=_arm.second;}

    // A→B, A→C (1:2 split)
    addConn(inp,entry[0][0],exitl[1][0],ConnTurnType::TurnLeft);
    addConn(inp,entry[0][0],exitl[2][0],ConnTurnType::TurnRight);
    // B→C, C→B
    addConn(inp,entry[1][0],exitl[2][0],ConnTurnType::TurnLeft);
    addConn(inp,entry[2][0],exitl[1][0],ConnTurnType::TurnRight);

    IntersectionShapeGenerator gen;
    IntersectionOutput out;
    REQUIRE(gen.generate(inp, out));

    save(inp, "./Y-1N/input/", "");
    save(out, "./Y-1N/output/", "");

    REQUIRE((int)out.connectivity_curves.size() == 4);
    REQUIRE(allCurvesAvoidObstacles(out, inp));
}

TEST_CASE("Y-junction: narrow passage between two obstacles", "[Y][narrow]")
{
    std::cout << "TEST_CASE: " << "Y-junction: narrow passage between two obstacles" << std::endl;

    auto inp = buildYJunction();
    // Two obstacles creating a ~4m gap in the centre
    addObstacle(inp, "OL", makeRectObstacle({-5,2}, 2.5, 2.0));
    addObstacle(inp, "OR", makeRectObstacle({ 5,2}, 2.5, 2.0));
    IntersectionShapeGenerator gen;
    IntersectionOutput out;
    REQUIRE(gen.generate(inp, out));

    save(inp, "./Y-narrow/input/", "");
    save(out, "./Y-narrow/output/", "");

    REQUIRE((int)out.connectivity_curves.size() == 6);
}
#endif

// ═════════════════════════════════════════════════════════════════════════════
//  SECTION 4 — 环岛 (roundabout)
//  Modelled as 4 arms connecting to a one-way circular road.
//  Inside the roundabout we add a circular obstacle (central island).
//  Connectivity: each arm has one entry lane and one exit lane;
//  they connect through the ring (simplified as direct curves).
// ═════════════════════════════════════════════════════════════════════════════

static IntersectionInput buildRoundabout(double ring_r = 8.0,
                                          double lane_w = 3.5,
                                          bool   add_island = true)
{
    gConnIdx = 0;
    IntersectionInput inp;
    double fence_half = ring_r + lane_w * 3 + 5.0;
    addStandardArea(inp, fence_half);

    // 4 arms at N/E/S/W attaching to the ring
    struct Arm4 { std::string nm; Vec2d jpt; Vec2d dir; };
    double off = ring_r + lane_w;
    std::vector<Arm4> raw = {
        {"N",{0,-off},{0, 1}},
        {"S",{0, off},{0,-1}},
        {"E",{-off,0},{1, 0}},
        {"W",{ off,0},{-1,0}},
    };
    std::vector<LaneId> entry(4), exitl(4);
    for (int a=0;a<4;++a) {
        ArmDef def{raw[a].nm, raw[a].jpt, raw[a].dir, 1, 1, lane_w, 12.0};
        std::pair<std::vector<LaneId>,std::vector<LaneId>> _arm = buildArm(def, inp);
        entry[a] = _arm.first[0]; exitl[a] = _arm.second[0];
    }

    // Central island obstacle
    if (add_island)
        addObstacle(inp, "ISLAND", makeCircleObstacle({0,0}, ring_r - lane_w*0.5, 24));

    // Clockwise connections (right-hand traffic on roundabout):
    //  N→E (right turn), N→S (straight through ring), N→W (left through ring)
    //  + mirror for E, S, W
    struct CS{int f,t;ConnTurnType tt;};
    std::vector<CS> specs = {
        {0,2,ConnTurnType::TurnRight},{0,1,ConnTurnType::Straight},{0,3,ConnTurnType::TurnLeft},
        {2,1,ConnTurnType::TurnRight},{2,0,ConnTurnType::Straight},{2,3,ConnTurnType::TurnLeft},
        {1,3,ConnTurnType::TurnRight},{1,0,ConnTurnType::Straight},{1,2,ConnTurnType::TurnLeft},
        {3,0,ConnTurnType::TurnRight},{3,1,ConnTurnType::Straight},{3,2,ConnTurnType::TurnLeft},
    };
    for (auto& s:specs) addConn(inp, entry[s.f], exitl[s.t], s.tt);

    return inp;
}
#ifdef TEST_ROUNDABOUT
TEST_CASE("Roundabout: with central island — all 12 connections", "[roundabout][island]")
{
    std::cout << "TEST_CASE: " << "Roundabout: with central island — all 12 connections" << std::endl;

    auto inp = buildRoundabout(8.0, 3.5, true);
    IntersectionShapeGenerator gen;
    IntersectionOutput out;
    REQUIRE(gen.generate(inp, out));

    save(inp, "./roundabout-island/input/", "");
    save(out, "./roundabout-island/output/", "");

    REQUIRE((int)out.connectivity_curves.size() == 12);
    REQUIRE(allCurvesAvoidObstacles(out, inp));
}

TEST_CASE("Roundabout: no island — basic connectivity", "[roundabout][no_island]")
{
    std::cout << "TEST_CASE: " << "Roundabout: no island — basic connectivity" << std::endl;

    auto inp = buildRoundabout(8.0, 3.5, false);
    IntersectionShapeGenerator gen;
    IntersectionOutput out;
    REQUIRE(gen.generate(inp, out));

    save(inp, "./roundabout-noisland/input/", "");
    save(out, "./roundabout-noisland/output/", "");

    REQUIRE((int)out.connectivity_curves.size() == 12);
}

TEST_CASE("Roundabout: oversized island causes narrow passage", "[roundabout][narrow_island]")
{
    std::cout << "TEST_CASE: " << "Roundabout: oversized island causes narrow passage" << std::endl;

    // Island radius = ring_r = very large, leaving almost no room
    auto inp = buildRoundabout(8.0, 3.5, false);
    addObstacle(inp, "BIG_ISLAND", makeCircleObstacle({0,0}, 9.5, 24));
    IntersectionShapeGenerator gen;
    IntersectionOutput out;
    REQUIRE(gen.generate(inp, out));

    save(inp, "./roundabout-narrowisland/input/", "");
    save(out, "./roundabout-narrowisland/output/", "");

    // Some connections may degrade or warn; pipeline must not crash
    REQUIRE((int)out.connectivity_curves.size() == 12);
}

TEST_CASE("Roundabout: satellite obstacle on ring path", "[roundabout][obstacle]")
{
    std::cout << "TEST_CASE: " << "Roundabout: satellite obstacle on ring path" << std::endl;

    auto inp = buildRoundabout(8.0, 3.5, true);
    // Extra obstacle on the ring path between N and E
    addObstacle(inp, "DEBRIS", makeRectObstacle({-5,-5}, 1.2, 1.2));
    IntersectionShapeGenerator gen;
    IntersectionOutput out;
    REQUIRE(gen.generate(inp, out));

    save(inp, "./roundabout-obstacle/input/", "");
    save(out, "./roundabout-obstacle/output/", "");

    REQUIRE(allCurvesAvoidObstacles(out, inp));
}

TEST_CASE("Roundabout: topo-block — central island completely fills fence", "[roundabout][topo_block]")
{
    std::cout << "TEST_CASE: " << "Roundabout: topo-block — central island completely fills fence" << std::endl;

    auto inp = buildRoundabout(8.0, 3.5, false);
    // Monster obstacle fills the whole junction area
    addObstacle(inp, "FULLBLOCK", makeCircleObstacle({0,0}, 6.0, 24));
    IntersectionShapeGenerator gen;
    IntersectionOutput out;
    REQUIRE(gen.generate(inp, out));

    save(inp, "./roundabout-topoblock/input/", "");
    save(out, "./roundabout-topoblock/output/", "");

    int infeasible=0;
    for (auto& cc : out.connectivity_curves)
        if (cc.status==CurveStatus::Infeasible) infeasible++;
    REQUIRE(infeasible >= 4); // at least cross-arm routes blocked
}
#endif

// ═════════════════════════════════════════════════════════════════════════════
//  SECTION 5 — Cross-type multi-obstacle & edge scenarios
// ═════════════════════════════════════════════════════════════════════════════
#ifdef TEST_OBSTACLE
TEST_CASE("4-way: multiple obstacles — curves avoid all", "[4way][multi_obs]")
{
    std::cout << "TEST_CASE: " << "4-way: multiple obstacles — curves avoid all" << std::endl;

    auto inp = build4Way();
    addObstacle(inp, "O1", makeRectObstacle({-4, 4}, 1.0, 1.0));
    addObstacle(inp, "O2", makeRectObstacle({ 4,-4}, 1.0, 1.0));
    addObstacle(inp, "O3", makeCircleObstacle({-4,-4}, 0.8));
    IntersectionShapeGenerator gen;
    IntersectionOutput out;
    REQUIRE(gen.generate(inp, out));

    save(inp, "./4way-multiobs/input/", "");
    save(out, "./4way-multiobs/output/", "");

    REQUIRE((int)out.connectivity_curves.size() == 12);
    REQUIRE(allCurvesAvoidObstacles(out, inp));

    for (const auto& id : {ConnId("C1"), ConnId("C7"), ConnId("C10")}) {
        auto* cc = findCurveById(out, id);
        REQUIRE(cc != nullptr);
        REQUIRE(minRawObstacleDistance(*cc, inp) > 1.0);
    }
}

// TEST_CASE("4-way: sandwich on both sides — forced narrow centre", "[4way][double_sandwich]")
// {
//     auto inp = build4Way();
//     // Obstacles pressed near top and bottom fence edges
//     addObstacle(inp, "TOP_OBS", makeRectObstacle({0, 12}, 10.0, 1.5));
//     addObstacle(inp, "BOT_OBS", makeRectObstacle({0,-12}, 10.0, 1.5));
//     IntersectionShapeGenerator gen;
//     IntersectionOutput out;
//     REQUIRE(gen.generate(inp, out));
//
//     save(inp, "./4way-doublesandwich/input/", "");
//     save(out, "./4way-doublesandwich/output/", "");
//
//     REQUIRE((int)out.connectivity_curves.size() == 12);
// }

TEST_CASE("T-junction: U-turn on branch arm", "[T][uturn]")
{
    std::cout << "TEST_CASE: " << "T-junction: U-turn on branch arm" << std::endl;

    gConnIdx = 0;
    IntersectionInput inp;
    addStandardArea(inp, 14.0);
    ArmDef arms[3] = {
        {"N",{0,-10},{0,1},2,2,3.5},
        {"E",{-10,0},{1,0},1,1,3.5},
        {"W",{ 10,0},{-1,0},1,1,3.5},
    };
    std::vector<LaneId> entry[3], exitl[3];
    for (int a=0;a<3;++a){std::pair<std::vector<LaneId>,std::vector<LaneId>> _arm=buildArm(arms[a],inp);entry[a]=_arm.first;exitl[a]=_arm.second;}

    // Branch inner lane U-turns back
    addConn(inp, entry[0][0], exitl[0][0], ConnTurnType::UTurnLeft);
    // Branch outer lane → main road
    addConn(inp, entry[0][1], exitl[1][0], ConnTurnType::TurnLeft);
    addConn(inp, entry[0][1], exitl[2][0], ConnTurnType::TurnRight);
    // Main road through
    addConn(inp, entry[1][0], exitl[2][0], ConnTurnType::Straight);
    addConn(inp, entry[2][0], exitl[1][0], ConnTurnType::Straight);
    // Main road turns onto branch
    addConn(inp, entry[1][0], exitl[0][1], ConnTurnType::TurnRight);
    addConn(inp, entry[2][0], exitl[0][1], ConnTurnType::TurnLeft);

    IntersectionShapeGenerator gen;
    IntersectionOutput out;
    REQUIRE(gen.generate(inp, out));

    save(inp, "./T-uturn/input/", "");
    save(out, "./T-uturn/output/", "");

    REQUIRE((int)out.connectivity_curves.size() == 7);
    REQUIRE(allCurvesAvoidObstacles(out, inp));
}

TEST_CASE("T-junction: fence-obstacle sandwich on branch side", "[T][sandwich]")
{
    std::cout << "TEST_CASE: " << "T-junction: fence-obstacle sandwich on branch side" << std::endl;

    auto inp = buildTJunction();
    // Obstacle near top fence edge (above branch)
    addObstacle(inp, "SAND", makeRectObstacle({0,-5}, 1, 1));
    IntersectionShapeGenerator gen;
    IntersectionOutput out;
    REQUIRE(gen.generate(inp, out));

    save(inp, "./T-sandwich/input/", "");
    save(out, "./T-sandwich/output/", "");

    REQUIRE((int)out.connectivity_curves.size() == 6);
}

TEST_CASE("Y-junction: all U-turns (unusual but valid)", "[Y][uturn]")
{
    std::cout << "TEST_CASE: " << "Y-junction: all U-turns (unusual but valid)" << std::endl;

    gConnIdx = 0;
    IntersectionInput inp;
    addStandardArea(inp, 16.0);
    double R=10.0,lw=3.5;
    double a0=M_PI/2,a1=a0+2*M_PI/3,a2=a0-2*M_PI/3;
    auto mkArm=[&](const std::string& nm,double ang)->ArmDef{
        Vec2d dir(-std::sin(ang),std::cos(ang));
        return{nm,-dir*R,dir,2,2,lw};
    };
    ArmDef arms[3]={mkArm("A",a0),mkArm("B",a1),mkArm("C",a2)};
    std::vector<LaneId> entry[3],exitl[3];
    for(int a=0;a<3;++a){std::pair<std::vector<LaneId>,std::vector<LaneId>> _arm=buildArm(arms[a],inp);entry[a]=_arm.first;exitl[a]=_arm.second;}

    // Inner lanes U-turn, outer lanes do cross-arm turns
    for(int a=0;a<3;++a){
        int b=(a+1)%3,c=(a+2)%3;
        addConn(inp,entry[a][0],exitl[a][0],ConnTurnType::UTurnLeft);
        addConn(inp,entry[a][1],exitl[b][1],ConnTurnType::TurnLeft);
        addConn(inp,entry[a][1],exitl[c][1],ConnTurnType::TurnRight);
    }
    IntersectionShapeGenerator gen;
    IntersectionOutput out;
    REQUIRE(gen.generate(inp, out));

    save(inp, "./Y-uturn/input/", "");
    save(out, "./Y-uturn/output/", "");

    REQUIRE((int)out.connectivity_curves.size() == 9);
    REQUIRE(allCurvesAvoidObstacles(out, inp));
}

TEST_CASE("4-way: all turn types present simultaneously", "[4way][all_turns]")
{
    std::cout << "TEST_CASE: " << "4-way: all turn types present simultaneously" << std::endl;

    gConnIdx = 0;
    IntersectionInput inp;
    addStandardArea(inp, 18.0);
    ArmDef arms[4]={
        {"N",{0,-12},{0,1},3,3,3.5},{"S",{0,12},{0,-1},3,3,3.5},
        {"E",{-12,0},{1,0},3,3,3.5},{"W",{12,0},{-1,0},3,3,3.5},
    };
    std::vector<LaneId> entry[4],exitl[4];
    for(int a=0;a<4;++a){std::pair<std::vector<LaneId>,std::vector<LaneId>> _arm=buildArm(arms[a],inp);entry[a]=_arm.first;exitl[a]=_arm.second;}

    for(int a=0;a<4;++a){
        int opp=(a+2)%4,left=(a+1)%4,right=(a+3)%4;
        addConn(inp,entry[a][0],exitl[a][0],   ConnTurnType::UTurnLeft);   // inner: U-turn
        addConn(inp,entry[a][1],exitl[left][1], ConnTurnType::TurnLeft);   // middle: left
        addConn(inp,entry[a][2],exitl[opp][2],  ConnTurnType::Straight);   // outer: straight
        addConn(inp,entry[a][1],exitl[right][1],ConnTurnType::TurnRight);  // middle: right
    }
    IntersectionShapeGenerator gen;
    IntersectionOutput out;
    REQUIRE(gen.generate(inp, out));

    save(inp, "./4way-alluturn/input/", "");
    save(out, "./4way-alluturn/output/", "");

    REQUIRE((int)out.connectivity_curves.size() == 16);
    REQUIRE(allCurvesAvoidObstacles(out, inp));
}

TEST_CASE("4-way: topology validates — correct boundary counts", "[4way][validation]")
{
    std::cout << "TEST_CASE: " << "4-way: topology validates — correct boundary counts" << std::endl;

    auto inp = build4Way();
    IntersectionShapeGenerator gen;
    IntersectionOutput out;
    REQUIRE(gen.generate(inp, out));

    save(inp, "./4way-validation/input/", "");
    save(out, "./4way-validation/output/", "");

    REQUIRE(gen.lastReport().is_valid());
    REQUIRE(gen.lastReport().errors.empty());
}

// TEST_CASE("4-way: topology rejects wrong boundary count", "[4way][validation_fail]")
// {
//     auto inp = build4Way();
//     // Corrupt one group's boundaries
//     inp.lane_groups[0].boundaries.pop_back();
//     IntersectionShapeGenerator gen;
//     IntersectionOutput out;
//     REQUIRE_FALSE(gen.generate(inp, out));
//
//     save(inp, "./4way-validationfail/input/", "");
//     save(out, "./4way-validationfail/output/", "");
//
//     REQUIRE(!gen.lastReport().is_valid());
// }

TEST_CASE("All scenario types: perf stats populated and sane", "[perf]")
{
    std::cout << "TEST_CASE: " << "All scenario types: perf stats populated and sane" << std::endl;

    std::vector<IntersectionInput> inputs = {
        build4Way(), buildTJunction(), buildYJunction(),
        buildRoundabout()
    };
    for (auto& inp : inputs) {
        IntersectionShapeGenerator gen;
        IntersectionOutput out;
        if (!gen.generate(inp, out)) continue;

        save(inp, "./perf/input/", "");
        save(out, "./perf/output/", "");

        REQUIRE(out.perf.sdf_build_ms >= 0.0);
        REQUIRE(out.perf.optimize_ms  >= 0.0);
        REQUIRE(out.perf.edge_gen_ms  >= 0.0);
        REQUIRE(out.perf.area_gen_ms  >= 0.0);
    }
}
#endif

#ifdef TEST_FILE
TEST_CASE("Actual intersection of osm", "[actual]")
{
    std::cout << "TEST_CASE: " << "Actual intersection of osm" << std::endl;

    std::string fpth = std::string(PROJECT_ROOT_DIR) + "/intersection_input.json";
    std::vector<IntersectionInput> inputs = {IntersectionIO::loadFromFile(fpth)};
    for (auto& inp : inputs) {
        IntersectionShapeGenerator gen;
        IntersectionOutput out;
        if (!gen.generate(inp, out)) continue;

        save(inp, "./actual/input/", "");
        save(out, "./actual/output/", "");

        REQUIRE(out.perf.sdf_build_ms >= 0.0);
        REQUIRE(out.perf.optimize_ms  >= 0.0);
        REQUIRE(out.perf.edge_gen_ms  >= 0.0);
        REQUIRE(out.perf.area_gen_ms  >= 0.0);
    }
}
#endif
