#include <catch2/catch_test_macros.hpp>
#include "generator/edge_line_generator.h"
#include "utils.h"
#include "curve/bezier.h"

using namespace isg;

static LaneEdge makeEdge(const std::string& id, const std::vector<Vec2d>& pts) {
    LaneEdge e;
    e.id = id;
    e.geometry.points = toVec3dArray(pts);
    return e;
}

TEST_CASE("EdgeLineGenerator preserves left/right at exit endpoint", "[edge_line]") {
    IntersectionInput inp;

    Lane entry;
    entry.id = "E";
    entry.geometry.points = {Vec2d(0, 0), Vec2d(10, 0)};
    entry.left_edge_id = "E_L";
    entry.right_edge_id = "E_R";
    entry.groupId = "EG";

    Lane exit;
    exit.id = "X";
    exit.geometry.points = {Vec2d(20, 0), Vec2d(30, 0)};
    exit.left_edge_id = "X_L";
    exit.right_edge_id = "X_R";
    exit.groupId = "XG";

    LaneGroup eg;
    eg.id = "EG";
    eg.role = GroupRole::Entry;
    eg.lanes = {"E"};
    eg.boundaries = {"E_L", "E_R"};

    LaneGroup xg;
    xg.id = "XG";
    xg.role = GroupRole::Exit;
    xg.lanes = {"X"};
    xg.boundaries = {"X_L", "X_R"};

    inp.lanes = {entry, exit};
    inp.lane_groups = {eg, xg};
    inp.lane_edges = {
        makeEdge("E_L", {Vec2d(0, 1), Vec2d(10, 1)}),
        makeEdge("E_R", {Vec2d(0, -1), Vec2d(10, -1)}),
        makeEdge("X_L", {Vec2d(20, 1), Vec2d(30, 1)}),
        makeEdge("X_R", {Vec2d(20, -1), Vec2d(30, -1)})
    };

    Connectivity conn;
    conn.id = "C1";
    conn.entry_lane_id = "E";
    conn.exit_lane_id = "X";
    conn.enterGroupId = "EG";
    conn.exitGroupId = "XG";
    inp.connectivities = {conn};

    ConnectivityCurve cc;
    cc.id = "C1";
    cc.entry_lane_id = "E";
    cc.exit_lane_id = "X";
    BezierCurve center;
    center.segs.push_back(makeCubicG1(Vec2d(10, 0), Vec2d(1, 0),
                                      Vec2d(20, 0), Vec2d(1, 0), 0.33));
    cc.curve = std::make_shared<BezierCurve>(center);
    std::vector<ConnectivityCurve> centerlines = {cc};

    EdgeLineGenerator gen;
    auto edges = gen.generate(inp, centerlines);

    const ConnectivityLaneEdge* left = nullptr;
    const ConnectivityLaneEdge* right = nullptr;
    for (auto& e : edges) {
        if (e.id == "gen_el_left_C1") left = &e;
        if (e.id == "gen_el_right_C1") right = &e;
    }

    REQUIRE(left);
    REQUIRE(right);
    REQUIRE(left->geometry.points.front().y() > 0.9);
    REQUIRE(left->geometry.points.back().y() > 0.9);
    REQUIRE(right->geometry.points.front().y() < -0.9);
    REQUIRE(right->geometry.points.back().y() < -0.9);

    Vec2d left_tail = xyOf(left->geometry.points.back()) - xyOf(left->geometry.points[left->geometry.points.size() - 2]);
    Vec2d right_tail = xyOf(right->geometry.points.back()) - xyOf(right->geometry.points[right->geometry.points.size() - 2]);
    REQUIRE(left_tail.normalized().dot(Vec2d(1, 0)) > 0.8);
    REQUIRE(right_tail.normalized().dot(Vec2d(1, 0)) > 0.8);
}
