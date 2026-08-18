#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "intersection_shape_generator.h"
#include "constraints/fence_check.h"
#include "constraints/infeasibility_detector.h"
#include "generator/polygon_builder.h"
#include "optimizer/sdf_field.h"
#include "io/iodata_json.h"

#include "io/iodata_shapefile.h"

#include <cmath>
#include <sstream>

using Catch::Matchers::WithinAbs;

using namespace isg;

template <typename P>
static bool hasPolygonPointNear(
    const std::vector<P>& pts, const Vec2d& expected, double tol = 1e-6) {
    for (const auto& pt : pts) {
        if ((xyOf(pt) - expected).norm() <= tol)
            return true;
    }
    return false;
}

template <typename P>
static bool hasPolygonPointNearZ(
    const std::vector<P>& pts, const Vec2d& expected, double expected_z,
    double xy_tol = 1e-6, double z_tol = 1e-6) {
    for (const auto& pt : pts) {
        if ((xyOf(pt) - expected).norm() <= xy_tol &&
            std::abs(xyzOf(pt).z() - expected_z) <= z_tol)
            return true;
    }
    return false;
}

template <typename P>
static std::string pointListString(const std::vector<P>& pts) {
    std::ostringstream os;
    for (const auto& pt : pts)
        os << "(" << pt[0] << "," << pt[1] << ") ";
    return os.str();
}

TEST_CASE("Integration: road edge boundaries merge by identical endpoints", "[integration][polygon]") {
    std::vector<std::vector<Vec2d>> fragments = {
        {Vec2d(0, 0), Vec2d(1, 0)},
        {Vec2d(1, 0), Vec2d(2, 0)},
        {Vec2d(4, 0), Vec2d(3, 0)},
        {Vec2d(3, 0), Vec2d(2, 0)},
        {Vec2d(2, 0), Vec2d(2, 1)}
    };

    auto merged = IntersectionAreaBuilder::MergeBoundaries(fragments);
    REQUIRE(merged.size() == 2);
    REQUIRE(merged.front().size() == 5);
    CHECK((merged.front()[0] - Vec2d(0, 0)).norm() < 1e-9);
    CHECK((merged.front()[1] - Vec2d(1, 0)).norm() < 1e-9);
    CHECK((merged.front()[2] - Vec2d(2, 0)).norm() < 1e-9);
    CHECK((merged.front()[3] - Vec2d(3, 0)).norm() < 1e-9);
    CHECK((merged.front()[4] - Vec2d(4, 0)).norm() < 1e-9);
    REQUIRE(merged.back().size() == 2);
    CHECK((merged.back()[0] - Vec2d(2, 0)).norm() < 1e-9);
    CHECK((merged.back()[1] - Vec2d(2, 1)).norm() < 1e-9);

    std::vector<std::vector<Vec2d>> corner = {
        {Vec2d(0, 0), Vec2d(1, 0)},
        {Vec2d(1, 0), Vec2d(1, 1)}
    };
    auto merged_corner = IntersectionAreaBuilder::MergeBoundaries(corner);
    REQUIRE(merged_corner.size() == 1);
    REQUIRE(merged_corner.front().size() == 3);
    CHECK((merged_corner.front()[0] - Vec2d(0, 0)).norm() < 1e-9);
    CHECK((merged_corner.front()[1] - Vec2d(1, 0)).norm() < 1e-9);
    CHECK((merged_corner.front()[2] - Vec2d(1, 1)).norm() < 1e-9);
}

// ──────────────────────────────────────────────────────────────
//  Minimal 2-direction intersection builder
//  East-entry → West-exit (straight)
// ──────────────────────────────────────────────────────────────
static IntersectionInput makeTwoDirectionInput(bool add_obstacle = false) {
    IntersectionInput inp;

    // Lanes with geometry:
    //   Entry lane (LE): outside→junction, last point = junction entry
    //   Exit  lane (LR): junction→outside, first point = junction exit
    Lane le; le.id="LE"; le.width=3.5;
    le.geometry.points = {Vec2d(-3,0), Vec2d(0,0)};   // entry at (0,0), tangent=(1,0)
    Lane lr; lr.id="LR"; lr.width=3.5;
    lr.geometry.points = {Vec2d(10,0), Vec2d(13,0)};  // exit at (10,0), tangent=(1,0)
    inp.lanes = {le, lr};

    // Lane edges (stubs — geometry not used by generator directly)
    LaneEdge be0; be0.id="BE0"; be0.geometry.points={{0,-1.75},{0,0}};
    LaneEdge be1; be1.id="BE1"; be1.geometry.points={{0,0},{0,1.75}};
    LaneEdge be2; be2.id="BE2"; be2.geometry.points={{0,1.75},{0,3.5}};
    LaneEdge br0; br0.id="BR0";
    LaneEdge br1; br1.id="BR1";
    LaneEdge br2; br2.id="BR2";
    inp.lane_edges = {be0,be1,be2,br0,br1,br2};

    // Entry group (east side, going west → x direction = (1,0))
    LaneGroup entry;
    entry.id="EG"; entry.role=GroupRole::Entry;
    entry.lanes={"LE"};
    entry.boundaries={"BE0","BE1"};
    //entry.direction=Vec2d(1,0);
    //entry.ref_point=Vec2d(0,0);

    // Exit group (west side, continuing west)
    LaneGroup exitg;
    exitg.id="XG"; exitg.role=GroupRole::Exit;
    exitg.lanes={"LR"};
    exitg.boundaries={"BR0","BR1"};
    //exitg.direction=Vec2d(1,0);
    //exitg.ref_point=Vec2d(10,0);

    inp.lane_groups = {entry, exitg};

    // Connectivity: straight through
    Connectivity c;
    c.id="CONN1"; c.entry_lane_id="LE"; c.exit_lane_id="LR";
    c.turn_type=ConnTurnType::Straight;
    inp.connectivities = {c};

    // Coarse area
    Polygon2d fence;
    fence.outer = {{-1,-5},{11,-5},{11,5},{-1,5}};
    inp.area.geometry = fence;

    // Boundary (road edges)
    Boundary top_edge;
    top_edge.id="TOP"; top_edge.type=Boundary::Type::RoadEdge;
    top_edge.geometry.points={{-1,5},{11,5}};

    Boundary bot_edge;
    bot_edge.id="BOT"; bot_edge.type=Boundary::Type::RoadEdge;
    bot_edge.geometry.points={{-1,-5},{11,-5}};

    inp.boundaries = {top_edge, bot_edge};

    if (add_obstacle) {
        Obstacle obs;
        obs.geometry.outer = {{4,-1},{6,-1},{6,1},{4,1}};
        inp.obstacles = {obs};
    }

    return inp;
}

static IntersectionInput makeTwoDirectionNearObstacleInput(int mode) {
    IntersectionInput inp = makeTwoDirectionInput(false);
    inp.mode = mode;
    Obstacle obs;
    obs.geometry.outer = {{3, 0.2}, {7, 0.2}, {7, 1.2}, {3, 1.2}};
    inp.obstacles = {obs};
    return inp;
}

static double minRawObstacleSdfAlongCurve(const ConnectivityCurve& cc,
                                          const IntersectionInput& inp) {
    if (!cc.curve || inp.obstacles.empty())
        return 1e18;
    SDFField raw_sdf;
    raw_sdf.build(inp.area.geometry.bbox(), inp.obstacles, 0.05, 0.0);
    double min_d = 1e18;
    for (auto& pt : cc.curve->sampleByArcLength(120)) {
        std::pair<double, Vec2d> _q = raw_sdf.queryWithGrad(pt);
        min_d = std::min(min_d, _q.first);
    }
    return min_d;
}

// ──────────────────────────────────────────────────────────────
TEST_CASE("Integration: simple straight, no obstacles, OK status", "[integration]") {
    auto inp = makeTwoDirectionInput(false);

    IntersectionShapeGenerator gen;
    IntersectionOutput out;
    bool ok = gen.generate(inp, out);

    save(inp, "./straight-noobstacles/input/", "");
    save(out, "./straight-noobstacles/output/", "");

    REQUIRE(ok);
    REQUIRE(out.connectivity_curves.size() == 1);

    auto& cc = out.connectivity_curves[0];
    REQUIRE(cc.curve);
    REQUIRE((cc.status == CurveStatus::OK ||
             cc.status == CurveStatus::WarnA2));

    // Endpoints near entry/exit reference points
    REQUIRE_THAT((cc.curve->startPt()-Vec2d(0,0)).norm(), WithinAbs(0.0, 0.5));
    REQUIRE_THAT((cc.curve->endPt()  -Vec2d(10,0)).norm(),WithinAbs(0.0, 0.5));
}

TEST_CASE("Integration: start/end tangent G1 continuity", "[integration]") {
    auto inp = makeTwoDirectionInput(false);

    IntersectionShapeGenerator gen;
    IntersectionOutput out;
    gen.generate(inp, out);

    save(inp, "./G1/input/", "");
    save(out, "./G1/output/", "");

    REQUIRE(!out.connectivity_curves.empty());

    auto& cc = out.connectivity_curves[0];
    REQUIRE(cc.curve);

    // Tangent at start should point roughly east (+x)
    Vec2d start_tan = cc.curve->startTan();
    REQUIRE(start_tan[0] > 0.7);

    // Tangent at end should also point roughly east
    Vec2d end_tan = cc.curve->endTan();
    REQUIRE(end_tan[0] > 0.7);
}

static IntersectionInput makeGroupedDirectionInput() {
    IntersectionInput inp;

    Lane e0; e0.id = "E0"; e0.width = 3.5; e0.groupId = "EG"; e0.laneOrder = 0;
    e0.geometry.points = {Vec2d(-5, 0), Vec2d(0, 0)};
    Lane e1; e1.id = "E1"; e1.width = 3.5; e1.groupId = "EG"; e1.laneOrder = 1;
    e1.geometry.points = {Vec2d(-5, 2), Vec2d(0, 2)};
    Lane e2; e2.id = "E2"; e2.width = 3.5; e2.groupId = "EG"; e2.laneOrder = 2;
    e2.geometry.points = {Vec2d(-2, -6), Vec2d(0, 4)};

    Lane x0; x0.id = "X0"; x0.width = 3.5; x0.groupId = "XG"; x0.laneOrder = 0;
    x0.geometry.points = {Vec2d(10, 0), Vec2d(15, 0)};
    Lane x1; x1.id = "X1"; x1.width = 3.5; x1.groupId = "XG"; x1.laneOrder = 1;
    x1.geometry.points = {Vec2d(10, 2), Vec2d(15, 2)};
    Lane x2; x2.id = "X2"; x2.width = 3.5; x2.groupId = "XG"; x2.laneOrder = 2;
    x2.geometry.points = {Vec2d(10, 4), Vec2d(15, 4)};
    inp.lanes = {e0, e1, e2, x0, x1, x2};

    LaneGroup eg;
    eg.id = "EG"; eg.role = GroupRole::Entry; eg.lanes = {"E0", "E1", "E2"};
    LaneGroup xg;
    xg.id = "XG"; xg.role = GroupRole::Exit; xg.lanes = {"X0", "X1", "X2"};
    inp.lane_groups = {eg, xg};

    for (int i = 0; i < 3; ++i) {
        Connectivity c;
        c.id = "C" + std::to_string(i);
        c.entry_lane_id = "E" + std::to_string(i);
        c.exit_lane_id = "X" + std::to_string(i);
        c.enterGroupId = "EG";
        c.exitGroupId = "XG";
        c.turn_type = ConnTurnType::Straight;
        inp.connectivities.push_back(c);
    }

    inp.area.geometry.outer = {{-1, -8}, {12, -8}, {12, 8}, {-1, 8}};
    return inp;
}

static const ConnectivityCurve* findCurveByEntryLane(
    const IntersectionOutput& out, const LaneId& entry_lane_id) {
    for (const auto& cc : out.connectivity_curves)
        if (cc.entry_lane_id == entry_lane_id)
            return &cc;
    return nullptr;
}

TEST_CASE("Integration: connectivity direction can be unified by lane group", "[integration]") {
    auto inp = makeGroupedDirectionInput();

    IntersectionShapeGenerator default_gen;
    IntersectionOutput default_out;
    REQUIRE(default_gen.generate(inp, default_out));
    auto* default_curve = findCurveByEntryLane(default_out, "E2");
    REQUIRE(default_curve);
    REQUIRE(default_curve->curve);
    REQUIRE(default_curve->curve->startTan().normalized().dot(Vec2d(1, 0)) < 0.5);

    IntersectionShapeGenerator::Config cfg;
    cfg.connectivity_direction.mode = ConnectivityDirectionMode::GroupUnified;
    cfg.connectivity_direction.group_similarity_angle_deg = 5.0;
    IntersectionShapeGenerator group_gen(cfg);
    IntersectionOutput group_out;
    REQUIRE(group_gen.generate(inp, group_out));
    auto* group_curve = findCurveByEntryLane(group_out, "E2");
    REQUIRE(group_curve);
    REQUIRE(group_curve->curve);
    REQUIRE(group_curve->curve->startTan().normalized().dot(Vec2d(1, 0)) > 0.95);
}

TEST_CASE("Integration: CON_NUM=1 lane keeps its own direction under group unification",
          "[integration][direction][cluster]") {
    auto inp = makeGroupedDirectionInput();
    for (auto& lane : inp.lanes) {
        if (lane.id == "E2")
            lane.attrs["CON_NUM"] = "1";
    }

    IntersectionShapeGenerator::Config cfg;
    cfg.connectivity_direction.mode = ConnectivityDirectionMode::GroupUnified;
    IntersectionShapeGenerator gen(cfg);
    IntersectionOutput out;
    REQUIRE(gen.generate(inp, out));

    auto* curve = findCurveByEntryLane(out, "E2");
    REQUIRE(curve);
    REQUIRE(curve->curve);
    Vec2d own_dir = Vec2d(2, 10).normalized();
    Vec2d start_tan = curve->curve->startTan().normalized();
    CHECK(start_tan.dot(own_dir) > 0.95);
    CHECK(start_tan.dot(Vec2d(1, 0)) < 0.5);
}

TEST_CASE("Integration: group-unified direction falls back to leftmost lane without straight",
          "[integration][direction][cluster]") {
    auto inp = makeGroupedDirectionInput();
    for (auto& conn : inp.connectivities)
        conn.turn_type = ConnTurnType::TurnLeft;

    IntersectionShapeGenerator::Config cfg;
    cfg.connectivity_direction.mode = ConnectivityDirectionMode::GroupUnified;
    IntersectionShapeGenerator gen(cfg);
    IntersectionOutput out;
    REQUIRE(gen.generate(inp, out));

    auto* curve = findCurveByEntryLane(out, "E0");
    REQUIRE(curve);
    REQUIRE(curve->curve);
    Vec2d leftmost_dir = Vec2d(2, 10).normalized();
    CHECK(curve->curve->startTan().normalized().dot(leftmost_dir) > 0.95);
}

TEST_CASE("Integration: curve stays inside fence", "[integration]") {
    auto inp = makeTwoDirectionInput(false);

    IntersectionShapeGenerator gen;
    IntersectionOutput out;
    gen.generate(inp, out);

    save(inp, "./fence/input/", "");
    save(out, "./fence/output/", "");

    REQUIRE(!out.connectivity_curves.empty());
    auto& cc = out.connectivity_curves[0];
    REQUIRE(cc.curve);

    auto& fence = inp.area.geometry;
    for (auto& pt : cc.curve->sampleByArcLength(40)) {
        // Allow small numerical overshoot
        double dist = pointToPolygonDist(pt, fence);
        bool inside = polygonContains(fence, pt);
        if (!inside) {
            REQUIRE(dist < 0.3);  // no more than 30cm overshoot
        }
    }
}

TEST_CASE("Integration: obstacle avoidance, no penetration", "[integration]") {
    auto inp = makeTwoDirectionInput(true);   // add central obstacle

    IntersectionShapeGenerator::Config cfg;
    cfg.lbfgs.max_iter = 200;
    IntersectionShapeGenerator gen(cfg);
    IntersectionOutput out;
    bool ok = gen.generate(inp, out);

    save(inp, "./obstacle-nopenetration/input/", "");
    save(out, "./obstacle-nopenetration/output/", "");

    REQUIRE(ok);
    REQUIRE(!out.connectivity_curves.empty());

    auto& cc = out.connectivity_curves[0];
    if (cc.status == CurveStatus::Infeasible) return;  // acceptable

    REQUIRE(cc.curve);

    // Build SDF to check penetration
    SDFField check_sdf;
    check_sdf.build(inp.area.geometry.bbox(), inp.obstacles, 0.2, 0.0);

    for (auto& pt : cc.curve->sampleByArcLength(50)) {
        std::pair<double,Vec2d> _q = check_sdf.queryWithGrad(pt);
        // Allow tiny numerical penetration
        REQUIRE(_q.first > -0.2);
    }
}

TEST_CASE("Integration: non-mode2 obstacle avoidance only forbids penetration",
          "[integration][obstacle][mode]") {
    auto inp = makeTwoDirectionNearObstacleInput(1);

    IntersectionShapeGenerator gen;
    IntersectionOutput out;
    REQUIRE(gen.generate(inp, out));
    REQUIRE(out.connectivity_curves.size() == 1);

    const auto& cc = out.connectivity_curves[0];
    REQUIRE(cc.curve);

    double min_raw_sdf = minRawObstacleSdfAlongCurve(cc, inp);
    INFO("mode=1 min raw obstacle sdf=" << min_raw_sdf);
    CHECK(min_raw_sdf > -0.05);
    CHECK(min_raw_sdf < 0.8);
}

TEST_CASE("Integration: mode2 obstacle avoidance keeps one meter clearance",
          "[integration][obstacle][mode]") {
    auto inp = makeTwoDirectionNearObstacleInput(2);

    IntersectionShapeGenerator gen;
    IntersectionOutput out;
    REQUIRE(gen.generate(inp, out));
    REQUIRE(out.connectivity_curves.size() == 1);

    const auto& cc = out.connectivity_curves[0];
    REQUIRE(cc.curve);

    double min_raw_sdf = minRawObstacleSdfAlongCurve(cc, inp);
    INFO("mode=2 min raw obstacle sdf=" << min_raw_sdf);
    CHECK(min_raw_sdf > 0.8);
    CHECK(cc.violation.max_obstacle_penetration <= 0.05);
}

TEST_CASE("Integration: topology validation catches mismatched boundaries", "[integration]") {
    auto inp = makeTwoDirectionInput(false);
    // Break boundary count: remove one boundary from entry group
    inp.lane_groups[0].boundaries.pop_back();

    IntersectionShapeGenerator gen;
    IntersectionOutput out;
    bool ok = gen.generate(inp, out);

    save(inp, "./topology-validation/input/", "");
    save(out, "./topology-validation/output/", "");

    REQUIRE_FALSE(ok);
    REQUIRE(!gen.lastReport().is_valid());
    REQUIRE(!gen.lastReport().errors.empty());
}

TEST_CASE("Integration: fine area is non-empty polygon", "[integration]") {
    auto inp = makeTwoDirectionInput(false);

    IntersectionShapeGenerator gen;
    IntersectionOutput out;
    gen.generate(inp, out);

    save(inp, "./fine-area/input/", "");
    save(out, "./fine-area/output/", "");

    REQUIRE_FALSE(out.area.geometry.outer.empty());
    REQUIRE(out.area.geometry.outer.size() >= 3);
}

TEST_CASE("Integration: lane group cut area is used without road boundaries", "[integration]") {
    IntersectionInput inp;
    inp.id = "no_road_boundary_group_cuts";
    inp.area.geometry.outer = {{-2, -2}, {2, -2}, {2, 2}, {-2, 2}};

    auto add_group = [&](const std::string& name,
                         GroupRole role,
                         const Vec2d& lane_outer,
                         const Vec2d& lane_conn,
                         const Vec2d& edge0_outer,
                         const Vec2d& edge0_conn,
                         const Vec2d& edge1_outer,
                         const Vec2d& edge1_conn) {
        Lane lane;
        lane.id = name + "_L";
        lane.groupId = name;
        lane.laneOrder = 0;
        lane.geometry.points = toVec3dArray((role == GroupRole::Entry)
            ? std::vector<Vec2d>{lane_outer, lane_conn}
            : std::vector<Vec2d>{lane_conn, lane_outer});

        LaneEdge e0;
        e0.id = name + "_E0";
        e0.groupId = name;
        e0.lineOrder = 0;
        e0.geometry.points = toVec3dArray((role == GroupRole::Entry)
            ? std::vector<Vec2d>{edge0_outer, edge0_conn}
            : std::vector<Vec2d>{edge0_conn, edge0_outer});

        LaneEdge e1;
        e1.id = name + "_E1";
        e1.groupId = name;
        e1.lineOrder = 2;
        e1.geometry.points = toVec3dArray((role == GroupRole::Entry)
            ? std::vector<Vec2d>{edge1_outer, edge1_conn}
            : std::vector<Vec2d>{edge1_conn, edge1_outer});

        LaneGroup group;
        group.id = name;
        group.role = role;
        group.lanes = {lane.id};
        group.boundaries = {e0.id, e1.id};

        inp.lanes.push_back(lane);
        inp.lane_edges.push_back(e0);
        inp.lane_edges.push_back(e1);
        inp.lane_groups.push_back(group);
    };

    add_group("WEST_ENTRY", GroupRole::Entry,
              {-10, 0}, {-1, 0}, {-10, -1}, {-1, -1}, {-10, 1}, {-1, 1});
    add_group("SOUTH_ENTRY", GroupRole::Entry,
              {0, -10}, {0, -1}, {-1, -10}, {-1, -1}, {1, -10}, {1, -1});
    add_group("EAST_EXIT", GroupRole::Exit,
              {10, 0}, {1, 0}, {10, -1}, {1, -1}, {10, 1}, {1, 1});
    add_group("NORTH_EXIT", GroupRole::Exit,
              {0, 10}, {0, 1}, {-1, 10}, {-1, 1}, {1, 10}, {1, 1});

    IntersectionAreaBuilder builder(0.5);
    IntersectionArea area = builder.build(inp, {}, {});

    REQUIRE(inp.boundaries.empty());
    REQUIRE_FALSE(area.geometry.outer.empty());
    REQUIRE(area.geometry.outer.size() >= 4);
    REQUIRE(isSimplePolygon(area.geometry));

    inp.mode = 2;
    IntersectionArea jd_area = builder.build(inp, {}, {});
    REQUIRE_FALSE(jd_area.geometry.outer.empty());
    REQUIRE(jd_area.geometry.outer.size() >= 4);
    REQUIRE(isSimplePolygon(jd_area.geometry));
}

TEST_CASE("Integration: lane group cut keeps every lane endpoint as an ordered polyline", "[integration]") {
    IntersectionInput inp;
    inp.id = "multi_lane_group_cut_polyline";
    inp.area.geometry.outer = {{-2, -2}, {2, -2}, {2, 2}, {-2, 2}};

    auto add_lane = [&](const std::string& id,
                        const std::string& group_id,
                        GroupRole role,
                        int order,
                        Vec2d conn,
                        Vec2d outer) {
        Lane lane;
        lane.id = id;
        lane.groupId = group_id;
        lane.laneOrder = order;
        lane.geometry.points = toVec3dArray((role == GroupRole::Entry)
            ? std::vector<Vec2d>{outer, conn}
            : std::vector<Vec2d>{conn, outer});
        inp.lanes.push_back(lane);
    };

    LaneGroup west;
    west.id = "WEST_ENTRY";
    west.role = GroupRole::Entry;
    west.lanes = {"W0", "W1", "W2"};
    add_lane("W0", west.id, west.role, 0, {-1.0, -1.0}, {-6.0, -1.0});
    add_lane("W1", west.id, west.role, 1, {-0.6, 0.0}, {-6.0, 0.0});
    add_lane("W2", west.id, west.role, 2, {-1.0, 1.0}, {-6.0, 1.0});

    LaneGroup east;
    east.id = "EAST_EXIT";
    east.role = GroupRole::Exit;
    east.lanes = {"E0", "E1", "E2"};
    add_lane("E0", east.id, east.role, 0, {1.0, -1.0}, {6.0, -1.0});
    add_lane("E1", east.id, east.role, 1, {1.0, 0.0}, {6.0, 0.0});
    add_lane("E2", east.id, east.role, 2, {1.0, 1.0}, {6.0, 1.0});

    inp.lane_groups = {west, east};

    Boundary bottom;
    bottom.id = "bottom";
    bottom.type = Boundary::Type::RoadEdge;
    bottom.geometry.points = {{-2.0, -1.5}, {2.0, -1.5}};
    Boundary top;
    top.id = "top";
    top.type = Boundary::Type::RoadEdge;
    top.geometry.points = {{2.0, 1.5}, {-2.0, 1.5}};
    inp.boundaries = {bottom, top};

    IntersectionAreaBuilder builder(0.5, 0.01, 0.25, "clockwise");
    IntersectionArea area = builder.build(inp, {}, {});

    REQUIRE_FALSE(area.geometry.outer.empty());
    INFO(pointListString(area.geometry.outer));
    REQUIRE(isSimplePolygon(area.geometry));
    CHECK(hasPolygonPointNear(area.geometry.outer, Vec2d(-0.5, -1.0)));
    CHECK(hasPolygonPointNear(area.geometry.outer, Vec2d(-0.1, 0.0)));
    CHECK(hasPolygonPointNear(area.geometry.outer, Vec2d(-0.5, 1.0)));

    inp.mode = 2;
    IntersectionArea jd_area = builder.build(inp, {}, {});
    REQUIRE_FALSE(jd_area.geometry.outer.empty());
    INFO(pointListString(jd_area.geometry.outer));
    REQUIRE(isSimplePolygon(jd_area.geometry));
    CHECK(hasPolygonPointNear(jd_area.geometry.outer, Vec2d(-1.5, -1.0)));
    CHECK(hasPolygonPointNear(jd_area.geometry.outer, Vec2d(-1.1, 0.0)));
    CHECK(hasPolygonPointNear(jd_area.geometry.outer, Vec2d(-1.5, 1.0)));
    CHECK(hasPolygonPointNear(jd_area.geometry.outer, Vec2d(1.5, 0.0)));
    CHECK(hasPolygonPointNear(jd_area.geometry.outer, Vec2d(1.5, -1.0)));
}

TEST_CASE("Integration: road edges are trimmed by extended group cut lines", "[integration]") {
    IntersectionInput inp;
    inp.id = "road_edge_trim_by_group_cuts";
    inp.area.geometry.outer = {{-3, -3}, {3, -3}, {3, 3}, {-3, 3}};

    auto add_group = [&](const std::string& name,
                         GroupRole role,
                         const Vec2d& lane_outer,
                         const Vec2d& lane_conn,
                         const Vec2d& edge0_outer,
                         const Vec2d& edge0_conn,
                         const Vec2d& edge1_outer,
                         const Vec2d& edge1_conn) {
        Lane lane;
        lane.id = name + "_L";
        lane.groupId = name;
        lane.laneOrder = 0;
        lane.geometry.points = toVec3dArray((role == GroupRole::Entry)
            ? std::vector<Vec2d>{lane_outer, lane_conn}
            : std::vector<Vec2d>{lane_conn, lane_outer});

        LaneEdge e0;
        e0.id = name + "_E0";
        e0.groupId = name;
        e0.lineOrder = 0;
        e0.geometry.points = toVec3dArray((role == GroupRole::Entry)
            ? std::vector<Vec2d>{edge0_outer, edge0_conn}
            : std::vector<Vec2d>{edge0_conn, edge0_outer});

        LaneEdge e1;
        e1.id = name + "_E1";
        e1.groupId = name;
        e1.lineOrder = 2;
        e1.geometry.points = toVec3dArray((role == GroupRole::Entry)
            ? std::vector<Vec2d>{edge1_outer, edge1_conn}
            : std::vector<Vec2d>{edge1_conn, edge1_outer});

        LaneGroup group;
        group.id = name;
        group.role = role;
        group.lanes = {lane.id};
        group.boundaries = {e0.id, e1.id};

        inp.lanes.push_back(lane);
        inp.lane_edges.push_back(e0);
        inp.lane_edges.push_back(e1);
        inp.lane_groups.push_back(group);
    };

    add_group("WEST_ENTRY", GroupRole::Entry,
              {-10, 0}, {-1, 0}, {-10, -1}, {-1, -1}, {-10, 1}, {-1, 1});
    add_group("SOUTH_ENTRY", GroupRole::Entry,
              {0, -10}, {0, -1}, {-1, -10}, {-1, -1}, {1, -10}, {1, -1});
    add_group("EAST_EXIT", GroupRole::Exit,
              {10, 0}, {1, 0}, {10, -1}, {1, -1}, {10, 1}, {1, 1});
    add_group("NORTH_EXIT", GroupRole::Exit,
              {0, 10}, {0, 1}, {-1, 10}, {-1, 1}, {1, 10}, {1, 1});

    auto add_boundary = [&](const std::string& id, Boundary::Type type,
                            Vec2d a, Vec2d b, double za = 0.0, double zb = 0.0) {
        Boundary boundary;
        boundary.id = id;
        boundary.type = type;
        boundary.geometry.points = {Vec3d(a, za), Vec3d(b, zb)};
        inp.boundaries.push_back(boundary);
    };

    add_boundary("bottom", Boundary::Type::RoadEdge, {-2, -2}, {2, -2}, 10.0, 14.0);
    add_boundary("right", Boundary::Type::RoadEdge, {2, -2}, {2, 2}, 20.0, 24.0);
    add_boundary("top", Boundary::Type::RoadEdge, {2, 2}, {-2, 2}, 30.0, 34.0);
    add_boundary("left", Boundary::Type::RoadEdge, {-2, 2}, {-2, -2}, 40.0, 44.0);
    add_boundary("median", Boundary::Type::MedianStrip, {10, -2}, {10, 2});

    IntersectionAreaBuilder builder(0.5, 2.0, 0.5, "clockwise");
    IntersectionArea area = builder.build(inp, {}, {});

    REQUIRE_FALSE(area.geometry.outer.empty());
    INFO(pointListString(area.geometry.outer));
    REQUIRE(hasPolygonPointNear(area.geometry.outer, Vec2d(1.5, -2.0)));
    REQUIRE(isSimplePolygon(area.geometry));
    REQUIRE_FALSE(hasPolygonPointNear(area.geometry.outer, Vec2d(-2.0, -2.0), 1e-4));
    REQUIRE_FALSE(hasPolygonPointNear(area.geometry.outer, Vec2d(10.0, 0.0), 1e-4));
}

TEST_CASE("Integration: merged road edge does not keep internal bottom cut", "[integration]") {
    IntersectionInput inp;
    inp.id = "merged_road_edge_internal_bottom_cut";
    inp.area.geometry.outer = {{-3, -3}, {3, -3}, {3, 3}, {-3, 3}};

    auto add_group = [&](const std::string& name,
                         GroupRole role,
                         const Vec2d& lane_outer,
                         const Vec2d& lane_conn,
                         const Vec2d& edge0_outer,
                         const Vec2d& edge0_conn,
                         const Vec2d& edge1_outer,
                         const Vec2d& edge1_conn) {
        Lane lane;
        lane.id = name + "_L";
        lane.groupId = name;
        lane.laneOrder = 0;
        lane.geometry.points = toVec3dArray((role == GroupRole::Entry)
            ? std::vector<Vec2d>{lane_outer, lane_conn}
            : std::vector<Vec2d>{lane_conn, lane_outer});

        LaneEdge e0;
        e0.id = name + "_E0";
        e0.groupId = name;
        e0.lineOrder = 0;
        e0.geometry.points = toVec3dArray((role == GroupRole::Entry)
            ? std::vector<Vec2d>{edge0_outer, edge0_conn}
            : std::vector<Vec2d>{edge0_conn, edge0_outer});

        LaneEdge e1;
        e1.id = name + "_E1";
        e1.groupId = name;
        e1.lineOrder = 2;
        e1.geometry.points = toVec3dArray((role == GroupRole::Entry)
            ? std::vector<Vec2d>{edge1_outer, edge1_conn}
            : std::vector<Vec2d>{edge1_conn, edge1_outer});

        LaneGroup group;
        group.id = name;
        group.role = role;
        group.lanes = {lane.id};
        group.boundaries = {e0.id, e1.id};

        inp.lanes.push_back(lane);
        inp.lane_edges.push_back(e0);
        inp.lane_edges.push_back(e1);
        inp.lane_groups.push_back(group);
    };

    add_group("SOUTH_ENTRY_L", GroupRole::Entry,
              {-2, -10}, {-2, -1}, {-3, -10}, {-3, -1}, {-1, -10}, {-1, -1});
    add_group("SOUTH_ENTRY_M", GroupRole::Entry,
              {0, -10}, {0, -1}, {-1, -10}, {-1, -1}, {1, -10}, {1, -1});
    add_group("SOUTH_ENTRY_R", GroupRole::Entry,
              {2, -10}, {2, -1}, {1, -10}, {1, -1}, {3, -10}, {3, -1});
    add_group("EAST_EXIT", GroupRole::Exit,
              {10, 0}, {1, 0}, {10, -1}, {1, -1}, {10, 1}, {1, 1});
    add_group("WEST_EXIT", GroupRole::Exit,
              {-10, 0}, {-1, 0}, {-10, -1}, {-1, -1}, {-10, 1}, {-1, 1});

    auto add_boundary = [&](const std::string& id, Boundary::Type type,
                            Vec2d a, Vec2d b, double za = 0.0, double zb = 0.0) {
        Boundary boundary;
        boundary.id = id;
        boundary.type = type;
        boundary.geometry.points = {Vec3d(a, za), Vec3d(b, zb)};
        inp.boundaries.push_back(boundary);
    };

    add_boundary("bottom", Boundary::Type::RoadEdge, {-3, -2}, {3, -2}, 10.0, 16.0);
    add_boundary("right", Boundary::Type::RoadEdge, {3, -2}, {3, 2}, 20.0, 24.0);
    add_boundary("top", Boundary::Type::RoadEdge, {3, 2}, {-3, 2}, 30.0, 36.0);
    add_boundary("left", Boundary::Type::RoadEdge, {-3, 2}, {-3, -2}, 40.0, 46.0);

    IntersectionAreaBuilder builder(0.5, 2.0, 0.5, "clockwise");
    IntersectionArea area = builder.build(inp, {}, {});

    REQUIRE_FALSE(area.geometry.outer.empty());
    INFO(pointListString(area.geometry.outer));
    CHECK(hasPolygonPointNear(area.geometry.outer, Vec2d(1.5, -2.0), 0.4));
    CHECK_FALSE(hasPolygonPointNear(area.geometry.outer, Vec2d(0.0, -2.0), 0.35));
}

TEST_CASE("Integration: perf stats are populated", "[integration]") {
    auto inp = makeTwoDirectionInput(false);

    IntersectionShapeGenerator gen;
    IntersectionOutput out;
    gen.generate(inp, out);

    save(inp, "./perf-stats/input/", "");
    save(out, "./perf-stats/output/", "");

    REQUIRE(out.perf.sdf_build_ms >= 0.0);
    REQUIRE(out.perf.optimize_ms  >= 0.0);
    REQUIRE(out.perf.edge_gen_ms  >= 0.0);
    REQUIRE(out.perf.area_gen_ms  >= 0.0);
}

// ──────────────────────────────────────────────────────────────
//  4-way intersection: 4 entry + 4 exit groups, 8 connectivities
// ──────────────────────────────────────────────────────────────
static IntersectionInput makeFourWayInput() {
    IntersectionInput inp;

    // Coarse area: 20×20m square centred at origin
    Polygon2d fence;
    fence.outer = {{-10,-10},{10,-10},{10,10},{-10,10}};
    inp.area.geometry = fence;

    // 4 directions: North(N), South(S), East(E), West(W)
    // Each direction: 1 entry lane + 1 exit lane
    struct Dir { std::string name; Vec2d ref_in; Vec2d dir_in;
                                   Vec2d ref_out; Vec2d dir_out; };
    std::vector<Dir> dirs = {
        {"N", {0,-10}, {0,1},  {0,10},  {0,1}},
        {"S", {0, 10}, {0,-1}, {0,-10}, {0,-1}},
        {"E", {-10,0}, {1,0},  {10,0},  {1,0}},
        {"W", {10, 0}, {-1,0}, {-10,0}, {-1,0}},
    };

    for (auto& d : dirs) {
        // Entry lane: outside → ref_in (last pt = junction edge)
        // Exit  lane: ref_out → outside (first pt = junction edge)
        Lane le; le.id="LE_"+d.name; le.width=3.5;
        le.geometry.points=toVec3dArray(std::vector<Vec2d>{d.ref_in - d.dir_in*5.0, d.ref_in});
        Lane lr; lr.id="LR_"+d.name; lr.width=3.5;
        lr.geometry.points=toVec3dArray(std::vector<Vec2d>{d.ref_out, d.ref_out + d.dir_out*5.0});
        inp.lanes.push_back(le);
        inp.lanes.push_back(lr);

        LaneGroup entry;
        entry.id="EG_"+d.name; entry.role=GroupRole::Entry;
        entry.lanes={"LE_"+d.name};
        entry.boundaries={"B0_"+d.name,"B1_"+d.name};
        //entry.direction=d.dir_in;
        //entry.ref_point=d.ref_in;

        LaneGroup exitg;
        exitg.id="XG_"+d.name; exitg.role=GroupRole::Exit;
        exitg.lanes={"LR_"+d.name};
        exitg.boundaries={"C0_"+d.name,"C1_"+d.name};
        //exitg.direction=d.dir_out;
        //exitg.ref_point=d.ref_out;

        inp.lane_groups.push_back(entry);
        inp.lane_groups.push_back(exitg);
    }

    // Straight-through connectivities for each direction
    std::vector<std::pair<std::string,std::string>> straights = {
        {"N","S"},{"S","N"},{"E","W"},{"W","E"}
    };
    int idx = 0;
    for (auto& kv : straights) {
        auto& from = kv.first; auto& to = kv.second;
        Connectivity c;
        c.id = "C"+std::to_string(idx++);
        c.entry_lane_id = "LE_"+from;
        c.exit_lane_id  = "LR_"+to;
        c.turn_type     = ConnTurnType::Straight;
        inp.connectivities.push_back(c);
    }

    return inp;
}

TEST_CASE("Integration: 4-way intersection generates all curves", "[integration]") {
    auto inp = makeFourWayInput();

    IntersectionShapeGenerator gen;
    IntersectionOutput out;
    bool ok = gen.generate(inp, out);

    save(inp, "./4-way/input/", "");
    save(out, "./4-way/output/", "");

    REQUIRE(ok);

    // All 4 connectivities should have curves
    REQUIRE(out.connectivity_curves.size() == 4);
    for (auto& cc : out.connectivity_curves) {
        bool feasible = (cc.status != CurveStatus::Infeasible);
        if (feasible) REQUIRE(cc.curve);
    }
}

TEST_CASE("Generated connectivity geometry interpolates endpoint elevations", "[integration][elevation]") {
    IntersectionInput inp;
    inp.id = "elevation";

    Lane entry;
    entry.id = "E";
    entry.geometry.points = {Vec3d(-3, 0, 7), Vec3d(0, 0, 10)};
    entry.groupId = "EG";

    Lane exit;
    exit.id = "X";
    exit.geometry.points = {Vec3d(10, 0, 20), Vec3d(13, 0, 23)};
    exit.groupId = "XG";

    LaneGroup eg;
    eg.id = "EG";
    eg.role = GroupRole::Entry;
    eg.lanes = {"E"};

    LaneGroup xg;
    xg.id = "XG";
    xg.role = GroupRole::Exit;
    xg.lanes = {"X"};

    Connectivity conn;
    conn.id = "C";
    conn.entry_lane_id = "E";
    conn.exit_lane_id = "X";
    conn.enterGroupId = "EG";
    conn.exitGroupId = "XG";
    conn.turn_type = ConnTurnType::Straight;

    inp.lanes = {entry, exit};
    inp.lane_groups = {eg, xg};
    inp.connectivities = {conn};
    inp.area.geometry.outer = {Vec3d(-5, -5, 0), Vec3d(15, -5, 0), Vec3d(15, 5, 0), Vec3d(-5, 5, 0)};

    IntersectionShapeGenerator gen;
    IntersectionOutput out;
    REQUIRE(gen.generate(inp, out));
    REQUIRE(out.connectivity_curves.size() == 1);
    const auto& pts = out.connectivity_curves.front().geometry.points;
    REQUIRE(pts.size() >= 3);
    CHECK_THAT(pts.front().z(), WithinAbs(10.0, 1e-9));
    CHECK_THAT(pts.back().z(), WithinAbs(20.0, 1e-9));
    CHECK(pts[pts.size() / 2].z() > 10.0);
    CHECK(pts[pts.size() / 2].z() < 20.0);

    nlohmann::json j = lineString2dToJson(out.connectivity_curves.front().geometry);
    REQUIRE(j["points"][0].contains("z"));
    CHECK_THAT(j["points"][0]["z"].get<double>(), WithinAbs(10.0, 1e-9));
}
