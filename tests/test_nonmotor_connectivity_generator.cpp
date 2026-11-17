#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <algorithm>
#include <array>
#include <string>
#include <chrono>
#include <iostream>
#include <map>
#include <limits>
#include <set>

#include "io/iodata_json.h"
#include "nonmotor/nonmotor_backbone_extractor.h"
#include "nonmotor/nonmotor_connectivity_generator.h"
#include "nonmotor/nonmotor_corner_rounding.h"
#include "nonmotor/nonmotor_surface_constraints.h"
#include "toolkits/toolkits.h"
#include "curve/curve_utils.h"

namespace {

isg::Lane makeLane(const std::string& id, const isg::Vec2d& a, const isg::Vec2d& b) {
    isg::Lane lane;
    lane.id = id;
    lane.geometry.points.push_back(isg::Vec3d(a, 1.0));
    lane.geometry.points.push_back(isg::Vec3d(b, 2.0));
    return lane;
}

bool hasPoint(const isg::LineString2d& line, const isg::Vec2d& point) {
    for (const auto& candidate : line.points)
        if ((candidate.xy() - point).norm() < 1e-7) return true;
    return false;
}

}  // namespace

TEST_CASE("nonmotor connectivity follows the counter-clockwise backbone") {
    isg::IntersectionInput input;
    input.lanes.push_back(makeLane("entry", isg::Vec2d(-10, 5), isg::Vec2d(-1, 5)));
    input.lanes.push_back(makeLane("exit", isg::Vec2d(11, 5), isg::Vec2d(20, 5)));
    isg::Connectivity connection;
    connection.id = "nm-1";
    connection.entry_lane_id = "entry";
    connection.exit_lane_id = "exit";
    connection.lane_type = isg::MotorLaneType::NonMotorLane;
    input.connectivities.push_back(connection);

    isg::NonMotorBackbone backbone;
    backbone.closed = true;
    backbone.ring_geometry.points = {
        isg::Vec3d(0, 0, 0), isg::Vec3d(10, 0, 0), isg::Vec3d(10, 10, 0),
        isg::Vec3d(0, 10, 0), isg::Vec3d(0, 0, 0)};

    isg::NonMotorConnectivityGenerator generator;
    const std::vector<isg::ConnectivityCurve> curves = generator.generate(input, {}, backbone);
    REQUIRE(curves.size() == 1);
    REQUIRE(curves.front().geometry.points.size() >= 5);
    CHECK(curves.front().geometry.points.front().xy().isApprox(isg::Vec2d(-1, 5)));
    CHECK(curves.front().geometry.points.back().xy().isApprox(isg::Vec2d(11, 5)));
    CHECK(curves.front().geometry.points.front().z() == Catch::Approx(2.0));
    CHECK(curves.front().geometry.points.back().z() == Catch::Approx(1.0));
}

TEST_CASE("nonmotor U-turn selects the shorter backbone direction") {
    isg::IntersectionInput input;
    input.lanes.push_back(makeLane("entry-u", isg::Vec2d(-10, 2), isg::Vec2d(-1, 2)));
    input.lanes.push_back(makeLane("exit-u", isg::Vec2d(11, 9), isg::Vec2d(20, 9)));
    isg::Connectivity connection;
    connection.id = "nm-u-short";
    connection.entry_lane_id = "entry-u";
    connection.exit_lane_id = "exit-u";
    connection.turn_type = isg::ConnTurnType::UTurnLeft;
    connection.lane_type = isg::MotorLaneType::NonMotorLane;
    input.connectivities.push_back(connection);

    isg::NonmotorConifg config;
    config.overlap_offset = 0.0;
    config.corner_radius = 0.0;
    isg::NonMotorBackbone backbone;
    backbone.closed = true;
    backbone.ring_geometry.points = {
        isg::Vec3d(0, 0), isg::Vec3d(10, 0), isg::Vec3d(10, 10),
        isg::Vec3d(0, 10), isg::Vec3d(0, 0)};

    const auto curves = isg::NonMotorConnectivityGenerator(config).generate(input, {}, backbone);
    REQUIRE(curves.size() == 1);
    const auto& path = curves.front().geometry;
    CHECK(hasPoint(path, isg::Vec2d(0, 10)));
    CHECK(hasPoint(path, isg::Vec2d(10, 10)));
    CHECK_FALSE(hasPoint(path, isg::Vec2d(10, 0)));
}

TEST_CASE("nonmotor attachment uses the nearest branch intersection when offsets are disabled") {
    isg::IntersectionInput input;
    input.lanes.push_back(makeLane("entry-branch", isg::Vec2d(20, 4), isg::Vec2d(16, 4)));
    input.lanes.push_back(makeLane("exit-ring", isg::Vec2d(5, 11), isg::Vec2d(5, 20)));
    isg::Connectivity connection;
    connection.id = "nm-branch-attachment";
    connection.entry_lane_id = "entry-branch";
    connection.exit_lane_id = "exit-ring";
    connection.lane_type = isg::MotorLaneType::NonMotorLane;
    input.connectivities.push_back(connection);

    isg::NonmotorConifg config;
    config.overlap_offset = 0.0;
    config.corner_radius = 0.0;
    isg::NonMotorBackbone backbone;
    backbone.closed = true;
    backbone.ring_nodes = {"n0", "n1", "n2", "n3"};
    backbone.ring_geometry.points = {
        isg::Vec3d(0, 0), isg::Vec3d(10, 0), isg::Vec3d(10, 10),
        isg::Vec3d(0, 10), isg::Vec3d(0, 0)};
    isg::NonMotorBackboneEdge branch;
    branch.id = "branch";
    branch.from = "n1";
    branch.to = "n4";
    branch.branch = true;
    branch.geometry.points = {isg::Vec3d(10, 0), isg::Vec3d(15, 5)};
    backbone.edges.push_back(branch);

    const auto curves = isg::NonMotorConnectivityGenerator(config).generate(input, {}, backbone);
    REQUIRE(curves.size() == 1);
    const auto& path = curves.front().geometry.points;
    REQUIRE(path.size() >= 4);
    CHECK(path.front().xy().isApprox(isg::Vec2d(16, 4)));
    CHECK(path[1].xy().isApprox(isg::Vec2d(14, 4)));
    CHECK(path[path.size() - 2].xy().isApprox(isg::Vec2d(5, 10)));
    CHECK(path.back().xy().isApprox(isg::Vec2d(5, 11)));
}

TEST_CASE("nonmotor attachment ignores backbone intersections behind the lane endpoint") {
    isg::IntersectionInput input;
    // The entry ray points +X and the exit ray points +X out of the junction.
    input.lanes.push_back(makeLane("entry-ray", isg::Vec2d(-20, 5), isg::Vec2d(-10, 5)));
    input.lanes.push_back(makeLane("exit-ray", isg::Vec2d(20, 5), isg::Vec2d(30, 5)));
    isg::Connectivity connection;
    connection.id = "nm-ray-direction";
    connection.entry_lane_id = "entry-ray";
    connection.exit_lane_id = "exit-ray";
    connection.lane_type = isg::MotorLaneType::NonMotorLane;
    input.connectivities.push_back(connection);

    isg::NonmotorConifg config;
    config.overlap_offset = 0.0;
    config.corner_radius = 0.0;
    isg::NonMotorBackbone backbone;
    backbone.closed = true;
    backbone.ring_nodes = {"n0", "n1", "n2", "n3"};
    backbone.ring_geometry.points = {
        isg::Vec3d(0, 0), isg::Vec3d(10, 0), isg::Vec3d(10, 10),
        isg::Vec3d(0, 10), isg::Vec3d(0, 0)};
    isg::NonMotorBackboneEdge behind_entry;
    behind_entry.from = "n0";
    behind_entry.to = "be";
    behind_entry.branch = true;
    behind_entry.geometry.points = {isg::Vec3d(0, 0), isg::Vec3d(-20, 5)};
    backbone.edges.push_back(behind_entry);
    isg::NonMotorBackboneEdge ahead_exit;
    ahead_exit.from = "n1";
    ahead_exit.to = "ae";
    ahead_exit.branch = true;
    ahead_exit.geometry.points = {isg::Vec3d(10, 0), isg::Vec3d(30, 5)};
    backbone.edges.push_back(ahead_exit);

    const auto curves = isg::NonMotorConnectivityGenerator(config).generate(input, {}, backbone);
    REQUIRE(curves.size() == 1);
    const auto& points = curves.front().geometry.points;
    REQUIRE(points.size() >= 4);
    CHECK(points[1].xy().isApprox(isg::Vec2d(0, 5)));
    CHECK(points[points.size() - 2].xy().isApprox(isg::Vec2d(10, 5)));
}

TEST_CASE("216863072 connectivity 43 uses both nearest branch intersections") {
    const std::string input_path = std::string(PROJECT_ROOT_DIR) +
                                    "/datas/216863072_100002096.json";
    const isg::IntersectionInput input = isg::IntersectionIO::loadFromFile(input_path);
    const isg::IntersectionInput nonmotor = isg::ConnectivityMotorLogicFilter(input, false);
    isg::NonMotorBackboneExtractor extractor(std::string(PROJECT_ROOT_DIR)+"/cmake-build-debug");
    isg::NonMotorBackbone backbone = extractor.run(input);
    // 关闭偏移和圆角验证原始挂接来源；圆角会截去原始公共点，不再用点下标验证挂接来源。
    isg::NonmotorConifg config;
    config.overlap_offset = 0.0;
    config.corner_radius = 0.0;
    const auto curves = isg::NonMotorConnectivityGenerator(config).generate(nonmotor, {}, backbone);
    const auto connection = std::find_if(
        input.connectivities.begin(), input.connectivities.end(),
        [](const isg::Connectivity& item) { return item.id == "43"; });
    const auto curve = std::find_if(
        curves.begin(), curves.end(),
        [](const isg::ConnectivityCurve& item) { return item.id == "43"; });
    REQUIRE(connection != input.connectivities.end());
    REQUIRE(curve != curves.end());
    REQUIRE(curve->geometry.points.size() >= 4);
    const isg::Vec2d entry = input.entryPtDir(connection->entry_lane_id).first;
    const isg::Vec2d entry_direction = input.entryPtDir(connection->entry_lane_id).second;
    const isg::Vec2d exit = input.exitPtDir(connection->exit_lane_id).first;
    const isg::Vec2d exit_direction = input.exitPtDir(connection->exit_lane_id).second;
    CHECK((curve->geometry.points.front().xy() - entry).norm() < 1e-6);
    CHECK((curve->geometry.points[1].xy() - isg::Vec2d(31.3557180103, -17.1196837422)).norm() < 1e-5);
    CHECK((curve->geometry.points[curve->geometry.points.size() - 2].xy() -
           isg::Vec2d(-18.8567008966, -32.6911107410)).norm() < 1e-4);
    CHECK((curve->geometry.points.back().xy() - exit).norm() < 1e-6);
    CHECK(std::abs(isg::cross2d(curve->geometry.points[1].xy() - curve->geometry.points.front().xy(),
                                entry_direction)) < 1e-5);
    CHECK(std::abs(isg::cross2d(curve->geometry.points.back().xy() -
                                curve->geometry.points[curve->geometry.points.size() - 2].xy(),
                                exit_direction)) < 1e-5);
    CHECK((curve->geometry.points[1].xy() - curve->geometry.points.front().xy()).dot(entry_direction) > 0.0);
    CHECK((curve->geometry.points.back().xy() -
           curve->geometry.points[curve->geometry.points.size() - 2].xy()).dot(exit_direction) > 0.0);
}

TEST_CASE("nonmotor connectivity reports obstacle and fence violations") {
    isg::IntersectionInput input;
    input.lanes.push_back(makeLane("entry", isg::Vec2d(-10, 5), isg::Vec2d(-1, 5)));
    input.lanes.push_back(makeLane("exit", isg::Vec2d(11, 5), isg::Vec2d(20, 5)));
    isg::Connectivity connection;
    connection.id = "nm-2";
    connection.entry_lane_id = "entry";
    connection.exit_lane_id = "exit";
    connection.lane_type = isg::MotorLaneType::NonMotorLane;
    input.connectivities.push_back(connection);
    isg::Obstacle obstacle;
    obstacle.geometry.outer = {isg::Vec3d(4, -1), isg::Vec3d(6, -1), isg::Vec3d(6, 1), isg::Vec3d(4, 1)};
    input.obstacles.push_back(obstacle);
    input.area.geometry.outer = {isg::Vec3d(-2, -2), isg::Vec3d(12, -2), isg::Vec3d(12, 12), isg::Vec3d(-2, 12)};

    isg::NonMotorBackbone backbone;
    backbone.ring_geometry.points = {
        isg::Vec3d(0, 0), isg::Vec3d(10, 0), isg::Vec3d(10, 10),
        isg::Vec3d(0, 10), isg::Vec3d(0, 0)};
    const auto curves = isg::NonMotorConnectivityGenerator().generate(input, {}, backbone);
    REQUIRE(curves.size() == 1);
    CHECK(curves.front().status == isg::CurveStatus::WarnA2);
    CHECK(curves.front().violation.type == isg::ViolationInfo::InfeasibilityType::ForcedCross);
}

TEST_CASE("nonmotor connectivity preserves an explicitly fixed shape") {
    isg::IntersectionInput input;
    input.lanes.push_back(makeLane("entry", isg::Vec2d(-10, 5), isg::Vec2d(-1, 5)));
    input.lanes.push_back(makeLane("exit", isg::Vec2d(11, 5), isg::Vec2d(20, 5)));
    isg::Connectivity connection;
    connection.id = "nm-fixed";
    connection.entry_lane_id = "entry";
    connection.exit_lane_id = "exit";
    connection.lane_type = isg::MotorLaneType::NonMotorLane;
    connection.fixed_shape = true;
    connection.geometry.points = {isg::Vec3d(-1, 5, 7), isg::Vec3d(11, 5, 8)};
    input.connectivities.push_back(connection);

    isg::NonMotorBackbone backbone;
    backbone.ring_geometry.points = {isg::Vec3d(0, 0), isg::Vec3d(10, 0),
                                     isg::Vec3d(10, 10), isg::Vec3d(0, 10), isg::Vec3d(0, 0)};
    const auto curves = isg::NonMotorConnectivityGenerator().generate(input, {}, backbone);
    REQUIRE(curves.size() == 1);
    REQUIRE(curves.front().geometry.points.size() == 2);
    CHECK(curves.front().geometry.points[0] == connection.geometry.points[0]);
    CHECK(curves.front().geometry.points[1] == connection.geometry.points[1]);
}

namespace {
// 独立按共线投影交集验收实际输出，覆盖短段落在长段任意位置的部分重叠。
double sharedLength(const isg::LineString2d& a, const isg::LineString2d& b) {
    double sum = 0;
    for (size_t i=1;i<a.points.size();++i) for (size_t j=1;j<b.points.size();++j) {
        const auto p=a.points[i-1].xy(), q=a.points[i].xy(), r=b.points[j-1].xy(), s=b.points[j].xy();
        const isg::Vec2d d=q-p; const double length=d.norm(); if (length<1e-7 || (s-r).norm()<1e-7) continue;
        const isg::Vec2d t=d/length;
        if (std::abs(isg::cross2d(t,(s-r).normalized()))>1e-7 || std::abs(isg::cross2d(t,r-p))>1e-7 || std::abs(isg::cross2d(t,s-p))>1e-7) continue;
        const double x=(r-p).dot(t), y=(s-p).dot(t);
        sum+=std::max(0.0,std::min(length,std::max(x,y))-std::max(0.0,std::min(x,y)));
    } return sum;
}
isg::LineString2d effectiveBody(const isg::ConnectivityCurve& curve) {
    isg::LineString2d body=curve.geometry;
    if (!curve.fixed_shape && body.points.size()>=4) { body.points.erase(body.points.begin()); body.points.pop_back(); }
    return body;
}
isg::NonMotorBackbone layoutRing() {
    isg::NonMotorBackbone b; b.closed=true; b.ring_nodes={"n0","n1","n2","n3"};
    b.ring_geometry.points={isg::Vec3d(0,0),isg::Vec3d(20,0),isg::Vec3d(20,20),isg::Vec3d(0,20),isg::Vec3d(0,0)}; return b;
}
void addBottomConnection(isg::IntersectionInput& input,const std::string& id,double start,double end) {
    const auto en="entry-"+id, ex="exit-"+id;
    input.lanes.push_back(makeLane(en,{start,-10},{start,-1}));
    input.lanes.push_back(makeLane(ex,{end,-1},{end,-10}));
    isg::Connectivity c; c.id=id; c.entry_lane_id=en; c.exit_lane_id=ex; c.lane_type=isg::MotorLaneType::NonMotorLane;
    if (start>end) c.turn_type=isg::ConnTurnType::UTurnLeft; input.connectivities.push_back(c);
}

// 在固定X剖面插值实际输出，不依赖轨位分配实现或采样点下标。
double ordinateAtX(const isg::LineString2d& line, double x) {
    for (size_t i = 1; i < line.points.size(); ++i) {
        const auto a = line.points[i-1].xy(), b = line.points[i].xy();
        if (std::abs(b.x()-a.x()) < 1e-9 || x < std::min(a.x(), b.x())-1e-9 || x > std::max(a.x(), b.x())+1e-9) continue;
        return a.y()+(b.y()-a.y())*((x-a.x())/(b.x()-a.x()));
    }
    return std::numeric_limits<double>::quiet_NaN();
}
}

// 长直共享骨干中部必须保持左右居中，不能被首尾切向整段拉向路口内侧。
TEST_CASE("nonmotor shared straight backbone stays balanced after smoothing") {
    isg::IntersectionInput input;
    for (const auto& id : {"a", "b", "c"}) addBottomConnection(input, id, 4, 16);
    auto backbone = layoutRing();
    const auto curves = isg::NonMotorConnectivityGenerator().generate(input, {}, backbone);
    REQUIRE(curves.size() == 3);
    std::vector<double> offsets;
    for (const auto& curve : curves) {
        const double y = ordinateAtX(curve.geometry, 10);
        REQUIRE(std::isfinite(y)); offsets.push_back(y);
        CHECK(curve.status == isg::CurveStatus::OK);
    }
    std::sort(offsets.begin(), offsets.end());
    for (size_t i = 0; i < offsets.size(); ++i)
        CHECK(offsets[i] == Catch::Approx((static_cast<double>(i)-1)*0.3).margin(1e-6));
}

// 连通分量有3条路径，但两个局部共享区间各只有2条；两处都应使用±0.15米，而非全局序号。
TEST_CASE("nonmotor partial occupancy is centered independently on each shared interval") {
    isg::IntersectionInput input;
    addBottomConnection(input, "a", 2, 10);
    addBottomConnection(input, "b", 6, 14);
    addBottomConnection(input, "c", 12, 18);
    auto backbone = layoutRing();
    isg::NonmotorConifg config; config.corner_radius = 0;
    const auto curves = isg::NonMotorConnectivityGenerator(config).generate(input, {}, backbone);
    REQUIRE(curves.size() == 3);
    CHECK(ordinateAtX(curves[0].geometry, 8) == Catch::Approx(-0.15));
    CHECK(ordinateAtX(curves[1].geometry, 8) == Catch::Approx(0.15));
    CHECK(ordinateAtX(curves[1].geometry, 13) == Catch::Approx(-0.15));
    CHECK(ordinateAtX(curves[2].geometry, 13) == Catch::Approx(0.15));
    for (size_t i = 0; i < curves.size(); ++i) for (size_t j = i+1; j < curves.size(); ++j)
        CHECK(sharedLength(effectiveBody(curves[i]), effectiveBody(curves[j])) < 1e-6);
}

// 下侧/上侧路口边及相反行驶方向均按进出车道组所在侧安排掉头；单独掉头同样可向外偏移。
TEST_CASE("nonmotor Uturn tracks prefer the lane group side independently of direction") {
    for (bool rotate : {false, true}) for (double radius : {0.0, 2.0}) {
        isg::IntersectionInput input;
        addBottomConnection(input, "a", 4, 16);
        addBottomConnection(input, "z", 16, 4);
        if (rotate) for (auto& lane : input.lanes) for (auto& point : lane.geometry.points)
            point = isg::Vec3d(20-point.x(), 20-point.y(), point.z());
        auto backbone = layoutRing();
        isg::NonmotorConifg config; config.corner_radius = radius;
        const auto curves = isg::NonMotorConnectivityGenerator(config).generate(input, {}, backbone);
        REQUIRE(curves.size() == 2);
        const double ordinary = ordinateAtX(curves[0].geometry, 10);
        const double uturn = ordinateAtX(curves[1].geometry, 10);
        CHECK((rotate ? uturn-20 : -uturn) > 0);
        CHECK((rotate ? ordinary-20 : -ordinary) < 0);
        for (const auto& curve : curves) CHECK(curve.status == isg::CurveStatus::OK);
        std::reverse(input.connectivities.begin(), input.connectivities.end());
        const auto reversed = isg::NonMotorConnectivityGenerator(config).generate(input, {}, backbone);
        for (size_t i = 0; i < curves.size(); ++i)
            CHECK(curves[i].geometry.points == reversed[1-i].geometry.points);
    }
    isg::IntersectionInput input; addBottomConnection(input, "u", 16, 4);
    auto backbone = layoutRing();
    isg::NonmotorConifg config; config.corner_radius = 0;
    const auto curves = isg::NonMotorConnectivityGenerator(config).generate(input, {}, backbone);
    REQUIRE(curves.size() == 1);
    CHECK(ordinateAtX(curves[0].geometry, 10) < 0);
}

// 外侧被障碍完全阻断时可选择内侧可行轨位，不能为满足均衡或掉头偏好穿越障碍。
TEST_CASE("nonmotor distribution preferences yield to obstacle avoidance") {
    isg::IntersectionInput input;
    addBottomConnection(input, "a", 4, 16); addBottomConnection(input, "z", 16, 4);
    isg::Obstacle obstacle; obstacle.id = "outside-block";
    obstacle.geometry.outer = {isg::Vec3d(8,-5), isg::Vec3d(12,-5), isg::Vec3d(12,0), isg::Vec3d(8,0)};
    input.obstacles.push_back(obstacle);
    auto backbone = layoutRing();
    const auto curves = isg::NonMotorConnectivityGenerator().generate(input, {}, backbone);
    REQUIRE(curves.size() == 2);
    for (const auto& curve : curves) {
        CHECK(curve.status == isg::CurveStatus::OK);
        CHECK(ordinateAtX(curve.geometry, 10) > 0);
    }
    CHECK(sharedLength(effectiveBody(curves[0]), effectiveBody(curves[1])) < 1e-6);
}

TEST_CASE("nonmotor layout separates three identical two-anchor backbone paths") {
    isg::IntersectionInput input; for (const auto& id:{"a","b","c"}) addBottomConnection(input,id,4,16);
    auto backbone=layoutRing(); const auto curves=isg::NonMotorConnectivityGenerator().generate(input,{},backbone);
    REQUIRE(curves.size()==3);
    for (const auto& c:curves) { const auto& p=c.geometry.points; REQUIRE(p.size()>6);
        CHECK(p.front().xy().isApprox(isg::Vec2d(4,-1))); CHECK(p.back().xy().isApprox(isg::Vec2d(16,-1)));
        const auto first=(p[1].xy()-p[0].xy()).normalized(),last=(p.back().xy()-p[p.size()-2].xy()).normalized();
        // 0.2米采样弦只近似解析切向；解析G1由圆角控制点专项验证。
        CHECK(first.dot(isg::Vec2d(0,1))>0.0); CHECK(last.dot(isg::Vec2d(0,-1))>0.0);
    }
    for (size_t i=0;i<curves.size();++i) for (size_t j=i+1;j<curves.size();++j) CHECK(sharedLength(effectiveBody(curves[i]),effectiveBody(curves[j]))<1e-6);
}

// 两条共享路径的轨位为±0.15米：骨干首尾交点与中间点采用相同偏移，仅车道端点固定。
TEST_CASE("nonmotor layout offsets attachment anchors together with the entire backbone") {
    isg::IntersectionInput input;
    addBottomConnection(input, "a", 4, 16);
    addBottomConnection(input, "b", 4, 16);
    auto backbone = layoutRing();
    for (double radius : {0.0, 2.0}) {
        isg::NonmotorConifg config; config.corner_radius = radius;
        const auto curves = isg::NonMotorConnectivityGenerator(config).generate(input, {}, backbone);
        REQUIRE(curves.size() == 2);
        for (size_t i = 0; i < curves.size(); ++i) {
            INFO("corner_radius=" << radius << " path=" << curves[i].id);
            const double offset = i == 0 ? -0.15 : 0.15;
            const auto& geometry = curves[i].geometry;
            CHECK(geometry.points.front() == isg::Vec3d(4, -1, 2));
            CHECK(geometry.points.back() == isg::Vec3d(16, -1, 1));
            if (radius == 0) {
                CHECK(hasPoint(geometry, isg::Vec2d(4, offset)));
                CHECK(hasPoint(geometry, isg::Vec2d(16, offset)));
            } else {
                // 进出公共点按统一圆角截去；骨干中部仍保持原定偏移量。
                CHECK(ordinateAtX(geometry, 10) == Catch::Approx(offset));
            }
            CHECK_FALSE(hasPoint(geometry, isg::Vec2d(4, 0)));
            CHECK_FALSE(hasPoint(geometry, isg::Vec2d(16, 0)));
            CHECK(curves[i].status == isg::CurveStatus::OK);
            if (radius == 0)
                for (const auto& point : effectiveBody(curves[i]).points)
                    CHECK(point.y() == Catch::Approx(offset));
        }
        CHECK(sharedLength(effectiveBody(curves[0]), effectiveBody(curves[1])) < 1e-6);
    }
}

// 斜向车道连接偏移后的骨干交点，首尾坐标与高程固定，0.2米采样弦接近车道切向。
TEST_CASE("nonmotor layout connects shifted anchors to fixed oblique lane endpoints") {
    isg::IntersectionInput input;
    for (const auto& id : {"a", "b"}) addBottomConnection(input, id, 4, 16);
    for (size_t i = 0; i < 2; ++i) {
        input.lanes[2*i].geometry.points = {isg::Vec3d(1,-14,1), isg::Vec3d(4,-5,2)};
        input.lanes[2*i+1].geometry.points = {isg::Vec3d(16,-5,1), isg::Vec3d(19,-14,2)};
    }
    auto backbone = layoutRing();
    const auto curves = isg::NonMotorConnectivityGenerator().generate(input, {}, backbone);
    REQUIRE(curves.size() == 2);
    for (size_t i = 0; i < curves.size(); ++i) {
        INFO(curves[i].id);
        const auto& points = curves[i].geometry.points;
        REQUIRE(points.size() > 20);
        const double offset = i == 0 ? -0.15 : 0.15;
        CHECK(points.front() == isg::Vec3d(4,-5,2));
        CHECK(points.back() == isg::Vec3d(16,-5,1));
        CHECK(ordinateAtX(curves[i].geometry, 10) == Catch::Approx(offset));
        const auto incoming = input.entryPtDir(curves[i].entry_lane_id).second.normalized();
        const auto outgoing = input.exitPtDir(curves[i].exit_lane_id).second.normalized();
        CHECK((points[1].xy()-points.front().xy()).normalized().dot(incoming) > std::cos(2*isg::DEG2RAD));
        CHECK((points.back().xy()-points[points.size()-2].xy()).normalized().dot(outgoing) > std::cos(2*isg::DEG2RAD));
        // 长直段可以只保留两端；非共线圆角采样必须不超过0.2米，不能以整段曲线替代直线。
        for (size_t j = 1; j+1 < points.size(); ++j) {
            const auto before = points[j].xy()-points[j-1].xy();
            const auto after = points[j+1].xy()-points[j].xy();
            if (std::abs(isg::cross2d(before.normalized(), after.normalized())) > 1e-7 && before.norm() < 1)
                CHECK(before.norm() <= 0.20001);
        }
        CHECK(curves[i].status == isg::CurveStatus::OK);
    }
    // 车道共用端点的直线连接允许重用，有效骨干与转角仍不能正长度重叠。
    CHECK(sharedLength(effectiveBody(curves[0]), effectiveBody(curves[1])) < 1e-6);
}

// 单条路径自交仍禁止；与非机动不同曲线之间允许同簇交叉的规则区分。
TEST_CASE("nonmotor layout reports self intersection without altering fixed geometry") {
    isg::IntersectionInput input;
    addBottomConnection(input, "self", 4, 16);
    auto& connection = input.connectivities.front();
    connection.fixed_shape = true;
    connection.geometry.points = {isg::Vec3d(-2,-2), isg::Vec3d(2,2),
                                  isg::Vec3d(-2,2), isg::Vec3d(2,-2)};
    auto backbone = layoutRing();
    const auto curves = isg::NonMotorConnectivityGenerator().generate(input, {}, backbone);
    REQUIRE(curves.size() == 1);
    CHECK(curves[0].geometry.points == connection.geometry.points);
    CHECK(curves[0].status == isg::CurveStatus::Infeasible);
    CHECK(curves[0].violation.reason.find("non-motor path self-intersection") != std::string::npos);
}

TEST_CASE("nonmotor layout separates partial and reverse occupancy deterministically") {
    isg::IntersectionInput input; addBottomConnection(input,"a",2,14); addBottomConnection(input,"b",4,18); addBottomConnection(input,"u",17,3);
    auto backbone=layoutRing(); const auto curves=isg::NonMotorConnectivityGenerator().generate(input,{},backbone);
    REQUIRE(curves.size()==3);
    for (size_t i=0;i<curves.size();++i) for (size_t j=i+1;j<curves.size();++j) CHECK(sharedLength(effectiveBody(curves[i]),effectiveBody(curves[j]))<1e-6);
    std::reverse(input.connectivities.begin(),input.connectivities.end());
    const auto reverse=isg::NonMotorConnectivityGenerator().generate(input,{},backbone);
    for (size_t i=0;i<curves.size();++i) { const auto& a=curves[i].geometry.points; const auto& b=reverse[curves.size()-1-i].geometry.points;
        REQUIRE(a.size()==b.size()); for (size_t j=0;j<a.size();++j) REQUIRE((a[j]-b[j]).norm()<1e-10); }
}

TEST_CASE("nonmotor layout reserves immutable geometry and reports immutable overlap") {
    isg::IntersectionInput input; addBottomConnection(input,"fixed",2,18); addBottomConnection(input,"movable",4,16);
    input.connectivities[0].fixed_shape=true; input.connectivities[0].geometry.points={isg::Vec3d(6,0,7),isg::Vec3d(14,0,9)};
    auto backbone=layoutRing(); auto curves=isg::NonMotorConnectivityGenerator().generate(input,{},backbone);
    REQUIRE(curves.size()==2); CHECK(curves[0].geometry.points==input.connectivities[0].geometry.points);
    CHECK(sharedLength(effectiveBody(curves[0]),effectiveBody(curves[1]))<1e-6);
    input.connectivities[1].fixed_shape=true; input.connectivities[1].geometry.points={isg::Vec3d(10,0),isg::Vec3d(18,0)};
    curves=isg::NonMotorConnectivityGenerator().generate(input,{},backbone);
    for (const auto& c:curves) { CHECK(c.status==isg::CurveStatus::Infeasible); CHECK(c.violation.reason.find("effective backbone overlap")!=std::string::npos); }
}

TEST_CASE("nonmotor layout audits residual overlap when offsets are disabled") {
    isg::IntersectionInput input; addBottomConnection(input,"a",2,18); addBottomConnection(input,"b",4,16);
    isg::NonmotorConifg config; config.overlap_offset=0; config.corner_radius=0; auto backbone=layoutRing();
    const auto curves=isg::NonMotorConnectivityGenerator(config).generate(input,{},backbone);
    for (const auto& c:curves) { CHECK(c.status==isg::CurveStatus::Infeasible); CHECK(c.violation.reason.find("effective backbone overlap")!=std::string::npos); }
}

TEST_CASE("216863072 all effective nonmotor backbone segments are separated") {
    const auto input=isg::IntersectionIO::loadFromFile(std::string(PROJECT_ROOT_DIR)+"/datas/216863072_100002096.json");
    isg::NonMotorBackboneExtractor extractor(std::string(PROJECT_ROOT_DIR)+"/cmake-build-debug"); auto backbone=extractor.run(input);
    const auto nonmotor=isg::ConnectivityMotorLogicFilter(input,false);
    const auto started=std::chrono::steady_clock::now();
    const auto curves=isg::NonMotorConnectivityGenerator().generate(nonmotor,{},backbone);
    const double elapsed=std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-started).count();
    std::cout << "Layout generation: " << elapsed << " ms\n";
    CHECK(elapsed<15000);
    for (const auto& c:curves) std::cout << "Layout status " << c.id << " " << static_cast<int>(c.status) << " " << c.violation.reason << "\n";
    REQUIRE(curves.size()==10);
    for (const auto& curve : curves) {
        INFO("fixed lane endpoints: " << curve.id);
        const auto* entry_lane = nonmotor.findLane(curve.entry_lane_id);
        const auto* exit_lane = nonmotor.findLane(curve.exit_lane_id);
        REQUIRE(entry_lane != nullptr);
        REQUIRE(exit_lane != nullptr);
        CHECK(curve.geometry.points.front() == entry_lane->geometry.points.back());
        CHECK(curve.geometry.points.back() == exit_lane->geometry.points.front());
        CHECK(curve.violation.reason.find("self-intersection") == std::string::npos);
        CHECK(curve.violation.reason.find("intersects") == std::string::npos);
    }
    for (size_t i=0;i<curves.size();++i) for (size_t j=i+1;j<curves.size();++j) {
        INFO(curves[i].id << "/" << curves[j].id); CHECK(sharedLength(effectiveBody(curves[i]),effectiveBody(curves[j]))<1e-6);
    }

    // 退出43191714的两条曲线在cw_100000309_spine末端必须单调离开偏移骨干，
    // 不能在轨位事件短段处先向骨干回退再转向车道。
    const isg::Vec2d branch_origin(-16.9197850411, -25.9174886856);
    const isg::Vec2d branch_tangent =
        (isg::Vec2d(-19.8204175878, -36.0612574968)-branch_origin).normalized();
    for (const auto* id : {"13", "43"}) {
        const auto curve = std::find_if(curves.begin(), curves.end(), [&](const auto& item) { return item.id == id; });
        REQUIRE(curve != curves.end());
        size_t last_parallel = 0;
        for (size_t i=1; i<curve->geometry.points.size(); ++i) {
            const auto d = curve->geometry.points[i].xy()-curve->geometry.points[i-1].xy();
            if (d.norm()>0.05 && d.normalized().dot(branch_tangent)>1.0-1e-3) last_parallel=i;
        }
        REQUIRE(last_parallel>0);
        const auto signed_distance = [&](const isg::Vec3d& point) {
            return isg::cross2d(branch_tangent, point.xy()-branch_origin);
        };
        double previous=signed_distance(curve->geometry.points[last_parallel]);
        for (size_t i=last_parallel+1; i<curve->geometry.points.size(); ++i) {
            const double distance=signed_distance(curve->geometry.points[i]);
            CHECK(distance>=previous-1e-6);
            previous=distance;
        }

        // 反鼓肚修正不能以删除有效骨干为代价：两条43191714挂接都必须
        // 经过cw_100000309_spine及其上游si_100000506_spine_2的截取范围。
        const auto projection_range = [&](const isg::Vec2d& a, const isg::Vec2d& b) {
            double min_station=std::numeric_limits<double>::infinity(), max_station=-std::numeric_limits<double>::infinity();
            double nearest=std::numeric_limits<double>::infinity();
            const auto direction=b-a; const double length2=direction.squaredNorm();
            for (const auto& point : curve->geometry.points) {
                const double station=(point.xy()-a).dot(direction)/length2;
                const auto projected=a+direction*station;
                nearest=std::min(nearest,(point.xy()-projected).norm());
                if ((point.xy()-projected).norm()<0.5) {
                    min_station=std::min(min_station,station); max_station=std::max(max_station,station);
                }
            }
            return std::array<double,3>{{min_station,max_station,nearest}};
        };
        const auto cw309=projection_range(isg::Vec2d(-16.9197850411,-25.9174886856),
            isg::Vec2d(-19.8204175878,-36.0612574968));
        const auto si506=projection_range(isg::Vec2d(-16.8592182048,-17.3944587992),
            isg::Vec2d(-16.9197850411,-25.9174886856));
        CHECK(cw309[2] < 0.5);
        CHECK(cw309[0] < 0.0);
        CHECK(cw309[1] > 0.5);
        CHECK(si506[2] < 0.5);
        CHECK(si506[0] < 0.3);
        CHECK(si506[1] > 0.9);
    }

    // si_100000458_spine_1/cw_100000304_spine附近的真实折点必须由圆角连续
    // 通过；短站距不能留下单个大于30度的采样折角。
    for (const auto* id : {"14", "19", "28"}) {
        const auto curve = std::find_if(curves.begin(), curves.end(), [&](const auto& item) { return item.id == id; });
        REQUIRE(curve != curves.end());
        double maximum_turn=0;
        for (size_t i=1; i+1<curve->geometry.points.size(); ++i) {
            const auto a=curve->geometry.points[i].xy()-curve->geometry.points[i-1].xy();
            const auto b=curve->geometry.points[i+1].xy()-curve->geometry.points[i].xy();
            if (a.norm()<0.05 || b.norm()<0.05) continue;
            maximum_turn=std::max(maximum_turn, std::acos(std::max(-1.0,std::min(1.0,a.normalized().dot(b.normalized()))))*180.0/isg::PI);
        }
        CHECK(maximum_turn<30.0);
    }
}

TEST_CASE("216863736 left turns recover the missing backbone ring") {
    const auto input = isg::IntersectionIO::loadFromFile(
        std::string(PROJECT_ROOT_DIR) + "/datas/216863736_100000844.json");
    const auto nonmotor = isg::ConnectivityMotorLogicFilter(input, false);
    isg::NonMotorBackboneExtractor extractor(std::string(PROJECT_ROOT_DIR) + "/cmake-build-debug");
    auto backbone = extractor.run(input);
    const auto curves = isg::NonMotorConnectivityGenerator().generate(nonmotor, {}, backbone);
    CHECK(curves.size() == 8);
    CHECK(backbone.closed);
    REQUIRE(backbone.ring_geometry.points.size() >= 4);
    REQUIRE(backbone.ring_nodes.size() + 1 == backbone.ring_geometry.points.size());
    double ring_area_twice = 0.0;
    for (size_t i = 1; i < backbone.ring_geometry.points.size(); ++i)
        ring_area_twice += isg::cross2d(backbone.ring_geometry.points[i - 1].xy(),
                                        backbone.ring_geometry.points[i].xy());
    CHECK(ring_area_twice > 0.0);
    CHECK(backbone.ring_edges.size() == backbone.edges.size());
    CHECK(backbone.branch_edges.empty());
    for (const auto& edge : backbone.edges) {
        CHECK(edge.ring);
        CHECK_FALSE(edge.branch);
    }
    const auto distanceToSegment = [](const isg::Vec2d& point,
                                      const isg::Vec2d& a, const isg::Vec2d& b) {
        return isg::pointToSegment(point, a, b).first;
    };
    for (const auto* id : {"14", "18", "35", "52"}) {
        const auto curve = std::find_if(curves.begin(), curves.end(),
            [&](const isg::ConnectivityCurve& item) { return item.id == id; });
        REQUIRE(curve != curves.end());
        CHECK(curve->status == isg::CurveStatus::OK);
        REQUIRE(curve->geometry.points.size() > 100);
        size_t covered_edges = 0;
        for (size_t i = 1; i < backbone.ring_geometry.points.size(); ++i) {
            const auto a = backbone.ring_geometry.points[i - 1].xy();
            const auto b = backbone.ring_geometry.points[i].xy();
            const bool covered = std::any_of(curve->geometry.points.begin(),
                                             curve->geometry.points.end(),
                [&](const isg::Vec3d& point) {
                    return distanceToSegment(point, a, b) < 0.8;
                });
            if (covered) ++covered_edges;
        }
        CHECK(covered_edges >= 4);
    }
}

TEST_CASE("216784982 preserves explicit nonmotor declarations when lane metadata is absent") {
    const auto input = isg::IntersectionIO::loadFromFile(
        std::string(PROJECT_ROOT_DIR) + "/datas/216784982_110004764.json");
    const auto nonmotor = isg::ConnectivityMotorLogicFilter(input, false);
    REQUIRE(nonmotor.connectivities.size() == 13);
    CHECK(std::all_of(nonmotor.connectivities.begin(), nonmotor.connectivities.end(),
        [](const isg::Connectivity& connection) {
            return connection.lane_type == isg::MotorLaneType::NonMotorLane;
        }));
    isg::NonMotorBackboneExtractor extractor(std::string(PROJECT_ROOT_DIR) + "/cmake-build-debug");
    auto backbone = extractor.run(input);
    REQUIRE(backbone.closed);
    REQUIRE(backbone.ring_geometry.points.size() >= 4);
    const auto curves = isg::NonMotorConnectivityGenerator().generate(nonmotor, {}, backbone);
    REQUIRE(curves.size() == nonmotor.connectivities.size());
    CHECK(std::all_of(curves.begin(), curves.end(), [](const isg::ConnectivityCurve& curve) {
        return curve.logic_type == isg::MotorLogicType::NonMotorLogic &&
               curve.geometry.points.size() >= 2;
    }));
}

TEST_CASE("nonmotor layout preserves multi-segment branches without detouring to the ring") {
    auto backbone=layoutRing();
    isg::NonMotorBackboneEdge branch;
    branch.id="multi-branch"; branch.from="n0"; branch.to="tip"; branch.branch=true;
    branch.geometry.points={isg::Vec3d(0,0),isg::Vec3d(-5,0),isg::Vec3d(-10,0)};
    backbone.edges.push_back(branch);
    isg::IntersectionInput input;
    for (const auto& id:{"a","b","c"}) addBottomConnection(input,id,-8,-2);
    const auto curves=isg::NonMotorConnectivityGenerator().generate(input,{},backbone);
    REQUIRE(curves.size()==3);
    for (const auto& c:curves) {
        REQUIRE(c.geometry.points.size()>4);
        CHECK(c.geometry.points.front().xy().isApprox(isg::Vec2d(-8,-1)));
        CHECK(c.geometry.points.back().xy().isApprox(isg::Vec2d(-2,-1)));
        for (const auto& p:effectiveBody(c).points) CHECK(p.x() < -1.9);
    }
    for (size_t i=0;i<curves.size();++i) for (size_t j=i+1;j<curves.size();++j)
        CHECK(sharedLength(effectiveBody(curves[i]),effectiveBody(curves[j]))<1e-6);
}

TEST_CASE("nonmotor layout validates fence obstacle and RoadEdge after separation") {
    isg::IntersectionInput input;
    addBottomConnection(input,"a",2,18); addBottomConnection(input,"b",4,16);
    auto backbone=layoutRing();
    input.area.geometry.outer={isg::Vec3d(-2,-2),isg::Vec3d(22,-2),isg::Vec3d(22,4),isg::Vec3d(-2,4)};
    isg::Boundary edge; edge.id="outside"; edge.type=isg::Boundary::Type::RoadEdge;
    edge.geometry.points={isg::Vec3d(-1,-2),isg::Vec3d(21,-2)};
    input.boundaries.push_back(edge);
    isg::Obstacle obstacle; obstacle.id="outside";
    obstacle.geometry.outer={isg::Vec3d(6,2),isg::Vec3d(12,2),isg::Vec3d(12,3),isg::Vec3d(6,3)};
    input.obstacles.push_back(obstacle);
    auto curves=isg::NonMotorConnectivityGenerator().generate(input,{},backbone);
    for (const auto& c:curves) CHECK(c.status==isg::CurveStatus::OK);
    CHECK(sharedLength(effectiveBody(curves[0]),effectiveBody(curves[1]))<1e-6);

    // 固定端点被横穿的RoadEdge与不可避让障碍包围时必须报告，不能以偏移完成冒充成功。
    input.boundaries[0].geometry.points={isg::Vec3d(-1,-0.5),isg::Vec3d(21,-0.5)};
    curves=isg::NonMotorConnectivityGenerator().generate(input,{},backbone);
    for (const auto& c:curves) CHECK(c.violation.reason.find("RoadEdge")!=std::string::npos);
    input.boundaries.clear();
    input.area.geometry.outer={isg::Vec3d(-2,0.5),isg::Vec3d(22,0.5),isg::Vec3d(22,4),isg::Vec3d(-2,4)};
    curves=isg::NonMotorConnectivityGenerator().generate(input,{},backbone);
    for (const auto& c:curves) CHECK(c.violation.reason.find("fence")!=std::string::npos);
    input.area.geometry.outer.clear();
    // 覆盖允许的全部横向轨位和圆角范围，保证该场景确实不可避让。
    input.obstacles[0].geometry.outer={isg::Vec3d(8,-10),isg::Vec3d(12,-10),isg::Vec3d(12,10),isg::Vec3d(8,10)};
    curves=isg::NonMotorConnectivityGenerator().generate(input,{},backbone);
    for (const auto& c:curves) CHECK(c.violation.reason.find("obstacle")!=std::string::npos);
}

TEST_CASE("nonmotor layout allows shared lane endpoints without cluster intersection validation") {
    isg::IntersectionInput input; addBottomConnection(input,"a",4,16); addBottomConnection(input,"b",4,16);
    for (auto& c:input.connectivities) c.enterGroupId="same-entry";
    auto backbone=layoutRing();
    const auto curves=isg::NonMotorConnectivityGenerator().generate(input,{},backbone);
    CHECK(sharedLength(effectiveBody(curves[0]),effectiveBody(curves[1]))<1e-6);
    for (const auto& c:curves) {
        CHECK(c.status==isg::CurveStatus::OK);
        CHECK(c.violation.reason.find("same-cluster")==std::string::npos);
    }
}

// 非机动同簇曲线可在内部交叉；零长度交集不能误报为骨干重叠。
TEST_CASE("nonmotor layout permits interior crossings of same-cluster fixed curves") {
    isg::IntersectionInput input;
    addBottomConnection(input, "horizontal", 2, 18);
    addBottomConnection(input, "vertical", 4, 16);
    for (auto& c : input.connectivities) {
        c.fixed_shape = true;
        c.enterGroupId = "shared-entry";
        c.exitGroupId = "shared-exit";
    }
    input.connectivities[0].geometry.points = {isg::Vec3d(-3, 0, 7), isg::Vec3d(3, 0, 8)};
    input.connectivities[1].geometry.points = {isg::Vec3d(0, -3, 9), isg::Vec3d(0, 3, 10)};
    auto backbone = layoutRing();
    const auto curves = isg::NonMotorConnectivityGenerator().generate(input, {}, backbone);
    REQUIRE(curves.size() == 2);
    CHECK(sharedLength(curves[0].geometry, curves[1].geometry) == Catch::Approx(0));
    for (size_t i = 0; i < curves.size(); ++i) {
        CHECK(curves[i].status == isg::CurveStatus::OK);
        CHECK(curves[i].violation.reason.empty());
        CHECK(curves[i].geometry.points == input.connectivities[i].geometry.points);
    }
}

// 固定形态参与机动车重叠审计，但不能为了消重叠移动任何固定输入点。
TEST_CASE("nonmotor layout preserves fixed geometry when auditing overlapping motor curves") {
    isg::IntersectionInput input;
    addBottomConnection(input, "fixed", 4, 16);
    input.connectivities[0].fixed_shape = true;
    input.connectivities[0].geometry.points = {isg::Vec3d(4, 0, 7), isg::Vec3d(16, 0, 8)};
    isg::ConnectivityCurve motor;
    motor.id = "motor-reserved";
    motor.geometry.points = {isg::Vec3d(8, 0), isg::Vec3d(18, 0)};
    auto backbone = layoutRing();
    const auto curves = isg::NonMotorConnectivityGenerator().generate(input, {motor}, backbone);
    REQUIRE(curves.size() == 1);
    CHECK(curves[0].geometry.points == input.connectivities[0].geometry.points);
    CHECK(curves[0].status == isg::CurveStatus::Infeasible);
    CHECK(curves[0].violation.reason.find("motor overlap conflict: motor-reserved") != std::string::npos);
}

// 独立比较圆角前后线段净距：允许在合规范围内变化，但不能新增净距违规或穿越告警。
TEST_CASE("nonmotor corner rounding preserves RoadEdge clearance on 216863072") {
    const auto input = isg::IntersectionIO::loadFromFile(std::string(PROJECT_ROOT_DIR) + "/datas/216863072_100002096.json");
    isg::NonMotorBackboneExtractor extractor(std::string(PROJECT_ROOT_DIR) + "/cmake-build-debug");
    auto backbone = extractor.run(input);
    const auto nonmotor = isg::ConnectivityMotorLogicFilter(input, false);
    std::map<std::string, double> raw_clearance;
    std::map<std::string, isg::LineString2d> raw_geometry;
    for (double radius : {0.0, 2.0}) {
        isg::NonmotorConifg config; config.corner_radius = radius;
        const auto curves = isg::NonMotorConnectivityGenerator(config).generate(nonmotor, {}, backbone);
        for (const auto& curve : curves) {
            INFO("corner_radius=" << radius << " connectivity=" << curve.id);
            double nearest = 1e100;
            for (size_t i = 1; i < curve.geometry.points.size(); ++i) {
                const auto a = curve.geometry.points[i-1], b = curve.geometry.points[i];
                for (const auto& boundary : input.boundaries) {
                    if (boundary.type != isg::Boundary::Type::RoadEdge) continue;
                    for (size_t j = 1; j < boundary.geometry.points.size(); ++j) {
                        const auto c = boundary.geometry.points[j-1], d = boundary.geometry.points[j];
                        double distance = std::min({isg::pointToSegment(a.xy(), c, d).first,
                            isg::pointToSegment(b.xy(), c, d).first, isg::pointToSegment(c.xy(), a, b).first,
                            isg::pointToSegment(d.xy(), a, b).first});
                        // 非平行段按有向参数检查内部相交；共线重叠已由端点投影给出零距离。
                        const auto ab = b.xy()-a.xy(), cd = d.xy()-c.xy(), ca = c.xy()-a.xy();
                        const double determinant = isg::cross2d(ab, cd);
                        if (std::abs(determinant) > 1e-12) {
                            const double t = isg::cross2d(ca, cd)/determinant;
                            const double u = isg::cross2d(ca, ab)/determinant;
                            if (t >= 0 && t <= 1 && u >= 0 && u <= 1) distance = 0;
                        }
                        nearest = std::min(nearest, distance);
                    }
                }
            }
            if (radius == 0) {
                raw_clearance[curve.id] = nearest;
                raw_geometry[curve.id] = curve.geometry;
                continue;
            }
            if (raw_clearance.at(curve.id) + 1e-9 >= config.road_edge_clearance)
                CHECK(nearest + 1e-9 >= config.road_edge_clearance);
            const auto& raw = raw_geometry.at(curve.id).points;
            const auto& rounded = curve.geometry.points;
            REQUIRE(rounded.size() >= 4);
            CHECK(rounded.front() == raw.front());
            CHECK(rounded.back() == raw.back());
            if (nearest + 1e-9 < config.road_edge_clearance) {
                CHECK(curve.status == isg::CurveStatus::Degraded);
                CHECK(curve.violation.reason == "non-motor path is too close to RoadEdge");
            } else CHECK(curve.status == isg::CurveStatus::OK);
        }
    }
}

// 独立验证90度转角的截距、三次圆弧形态、解析G1和0.2米弧长采样。
TEST_CASE("nonmotor corner radius trims both incident segments and samples at 0.2 metres") {
    const std::vector<isg::Vec3d> source = {
        isg::Vec3d(0, 0, 2), isg::Vec3d(10, 0, 12), isg::Vec3d(10, 10, 22)};
    const auto rounded = isg::nonmotor_detail::roundPolylineCorners(source, 2.0);
    REQUIRE(rounded.points.size() > 10);
    CHECK(rounded.points.front() == source.front());
    CHECK(rounded.points.back() == source.back());
    CHECK((rounded.points[1] - isg::Vec3d(8, 0, 10)).norm() < 1e-9);
    CHECK((rounded.points[rounded.points.size() - 2] - isg::Vec3d(10, 2, 14)).norm() < 1e-9);
    CHECK_FALSE(hasPoint(rounded, isg::Vec2d(10, 0)));
    const auto curve = isg::nonmotor_detail::fitCornerBezier(
        isg::Vec3d(8, 0, 10), source[1], isg::Vec3d(10, 2, 14));
    const double handle = 8.0 / 3.0 * std::tan(isg::PI / 8.0);
    CHECK(curve.ctrl[1].x() == Catch::Approx(8 + handle));
    CHECK(curve.ctrl[1].y() == Catch::Approx(0));
    CHECK(curve.ctrl[2].x() == Catch::Approx(10));
    CHECK(curve.ctrl[2].y() == Catch::Approx(2 - handle));
    CHECK(curve.tangent(0).xy().normalized().isApprox(isg::Vec2d(1, 0)));
    CHECK(curve.tangent(1).xy().normalized().isApprox(isg::Vec2d(0, 1)));
    double arc_chords = 0;
    for (size_t i = 2; i + 1 < rounded.points.size(); ++i) {
        const double chord = (rounded.points[i].xy() - rounded.points[i - 1].xy()).norm();
        CHECK(chord > 0);
        CHECK(chord <= 0.20001);
        if (i + 2 < rounded.points.size()) CHECK(chord > 0.1995);
        arc_chords += chord;
        // 这是标准圆弧三次近似，不是以等参数间隔采样的折线。
        CHECK((rounded.points[i].xy() - isg::Vec2d(8, 2)).norm() == Catch::Approx(2).margin(0.0007));
        CHECK(rounded.points[i].z() >= 10);
        CHECK(rounded.points[i].z() <= 14);
    }
    CHECK(arc_chords == Catch::Approx(isg::PI).margin(0.003));
}

// 左右转、浅转与锐转都保留两端方向和坐标，高程随圆角插值。
TEST_CASE("nonmotor corner fitting respects signed turn directions and elevation") {
    for (double angle : {30.0, -30.0, 90.0, -90.0, 150.0, -150.0}) {
        const isg::Vec3d vertex(0, 0, 10);
        const isg::Vec3d before(-2, 0, 8);
        const isg::Vec2d outgoing(std::cos(angle * isg::DEG2RAD), std::sin(angle * isg::DEG2RAD));
        const isg::Vec3d after(outgoing * 2, 12);
        const auto curve = isg::nonmotor_detail::fitCornerBezier(before, vertex, after);
        CHECK(curve.tangent(0).xy().normalized().isApprox(isg::Vec2d(1, 0)));
        CHECK(curve.tangent(1).xy().normalized().isApprox(outgoing));
        const auto samples = isg::nonmotor_detail::sampleCornerBezier(curve);
        REQUIRE(samples.size() > 2);
        CHECK(samples.front() == before);
        CHECK(samples.back() == after);
        for (size_t i = 1; i < samples.size(); ++i) {
            CHECK(std::isfinite(samples[i].x()));
            CHECK(std::isfinite(samples[i].y()));
            CHECK(std::isfinite(samples[i].z()));
            CHECK((samples[i].xy() - samples[i-1].xy()).norm() <= 0.20001);
            CHECK(samples[i].z() >= samples[i-1].z());
        }
    }
}

// 相邻圆角共享短段时同比缩短截距，避免重叠；禁用、重复点和直线确定性退化。
TEST_CASE("nonmotor corner rounding bounds adjacent trims and handles degenerate input") {
    const std::vector<isg::Vec3d> short_segments = {
        isg::Vec3d(0, 0), isg::Vec3d(2, 0), isg::Vec3d(2, 1), isg::Vec3d(4, 1)};
    const auto rounded = isg::nonmotor_detail::roundPolylineCorners(short_segments, 2.0);
    CHECK(rounded.points.front() == short_segments.front());
    CHECK(rounded.points.back() == short_segments.back());
    CHECK(hasPoint(rounded, isg::Vec2d(2, 0.5)));
    CHECK_FALSE(hasPoint(rounded, isg::Vec2d(2, 0)));
    CHECK_FALSE(hasPoint(rounded, isg::Vec2d(2, 1)));
    for (size_t i = 1; i < rounded.points.size(); ++i) {
        CHECK(std::isfinite(rounded.points[i].x()));
        CHECK(std::isfinite(rounded.points[i].y()));
        CHECK(std::isfinite(rounded.points[i].z()));
        CHECK(rounded.points[i].x() >= rounded.points[i-1].x() - 1e-9);
        CHECK(rounded.points[i].y() >= rounded.points[i-1].y() - 1e-9);
        CHECK((rounded.points[i] - rounded.points[i-1]).norm() > 1e-9);
    }
    CHECK(isg::nonmotor_detail::roundPolylineCorners(short_segments, 0).points == short_segments);
    CHECK(isg::nonmotor_detail::roundPolylineCorners({}, 2).points.empty());
    const std::vector<isg::Vec3d> straight = {
        isg::Vec3d(0, 0), isg::Vec3d(0, 0), isg::Vec3d(1, 0), isg::Vec3d(2, 0)};
    const auto line = isg::nonmotor_detail::roundPolylineCorners(straight, 2);
    REQUIRE(line.points.size() == 2);
    CHECK(line.points.front() == straight.front());
    CHECK(line.points.back() == straight.back());
    const std::vector<isg::Vec3d> reversal = {
        isg::Vec3d(0, 0), isg::Vec3d(1, 0), isg::Vec3d(0, 0)};
    CHECK(isg::nonmotor_detail::roundPolylineCorners(reversal, 2).points == reversal);
}

// 切向过渡保持调用者给定的端点坐标，可用于固定车道端点或已偏移的骨干交点。
TEST_CASE("nonmotor corner transitions preserve supplied positions and tangents") {
    const isg::Vec3d first(0, 0, 3), last(10, 0, 8);
    const isg::Vec2d incoming(0, 1), outgoing(0, -1);
    const auto curve = isg::nonmotor_detail::fitAttachmentTransition(first, incoming, last, outgoing, 2);
    CHECK(curve.tangent(0).xy().normalized().isApprox(incoming));
    CHECK(curve.tangent(1).xy().normalized().isApprox(outgoing));
    const auto rounded = isg::nonmotor_detail::roundPolylineCorners({first, last}, 2, &incoming, &outgoing);
    REQUIRE(rounded.points.size() > 3);
    CHECK(rounded.points.front() == first);
    CHECK(rounded.points.back() == last);
    for (size_t i = 1; i < rounded.points.size(); ++i) {
        CHECK((rounded.points[i].xy()-rounded.points[i-1].xy()).norm() <= 0.20001);
        CHECK(rounded.points[i].y() >= 0);
    }
}

// 切向过渡局部化：20米直骨干仅前后2米发生过渡，中间保持原直线轨位。
TEST_CASE("nonmotor corner endpoint transitions preserve long straight interiors") {
    const isg::Vec2d first(0,1), last(0,-1);
    const auto line = isg::nonmotor_detail::roundPolylineCorners(
        {isg::Vec3d(0,0), isg::Vec3d(20,0)}, 2, &first, &last);
    REQUIRE(hasPoint(line, isg::Vec2d(10,0)));
    for (const auto& point : line.points)
        if (point.x() > 2+1e-7 && point.x() < 18-1e-7) CHECK(point.y() == Catch::Approx(0).margin(1e-9));
    CHECK(line.points.front() == isg::Vec3d(0,0));
    CHECK(line.points.back() == isg::Vec3d(20,0));
}

// 出骨干转角在偏移骨干侧相切：从最后一段平行直线开始，法向距离只能朝退出车道单调变化。
// 旧实现强制骨干末点取车道切向，会出现先减小法向距离再返回偏移轨位的局部鼓肚。
TEST_CASE("216863072 conn19 exit corner leaves the offset backbone without a bulge") {
    const auto input = isg::IntersectionIO::loadFromFile(std::string(PROJECT_ROOT_DIR)+"/datas/216863072_100002096.json");
    const auto nonmotor = isg::ConnectivityMotorLogicFilter(input, false);
    isg::NonMotorBackboneExtractor extractor(std::string(PROJECT_ROOT_DIR)+"/cmake-build-debug");
    auto backbone = extractor.run(input);
    const auto edge = std::find_if(backbone.edges.begin(), backbone.edges.end(), [](const isg::NonMotorBackboneEdge& item) {
        return item.id == "smooth_cw_100000307_spine_si_100000711_spine_1";
    });
    REQUIRE(edge != backbone.edges.end());
    REQUIRE(edge->geometry.points.size() == 2);
    const auto origin = edge->geometry.points.front().xy();
    const isg::Vec2d tangent = (edge->geometry.points.back().xy()-origin).normalized();
    const auto curves = isg::NonMotorConnectivityGenerator().generate(nonmotor, {}, backbone);
    const auto curve = std::find_if(curves.begin(), curves.end(), [](const isg::ConnectivityCurve& item) { return item.id == "19"; });
    REQUIRE(curve != curves.end());
    CHECK(curve->exit_lane_id == "43191717");
    CHECK(curve->status == isg::CurveStatus::OK);
    const auto& points = curve->geometry.points;
    size_t last_parallel = 0;
    for (size_t i=1; i<points.size(); ++i) {
        const isg::Vec2d direction = points[i].xy()-points[i-1].xy();
        if (direction.norm() > 0.05 && direction.normalized().dot(tangent) > 1-1e-9)
            last_parallel = i;
    }
    REQUIRE(last_parallel > 0);
    REQUIRE(last_parallel+4 < points.size());
    const auto signed_distance = [&](const isg::Vec3d& point) { return isg::cross2d(tangent, point.xy()-origin); };
    double previous = signed_distance(points[last_parallel]);
    const double sign = signed_distance(points.back()) > previous ? 1 : -1;
    const double offset = previous;
    CHECK(std::abs(offset) > 0.05);
    for (size_t i=last_parallel+1; i<points.size(); ++i) {
        INFO("conn19 exit corner sample=" << i << " offset=" << offset);
        const double distance = signed_distance(points[i]);
        CHECK(sign*(distance-previous) >= -1e-6);
        CHECK(sign*(distance-offset) >= -1e-6);
        previous = distance;
    }
    const auto* exit_lane = nonmotor.findLane(curve->exit_lane_id);
    REQUIRE(exit_lane != nullptr);
    CHECK(points.back() == exit_lane->geometry.points.front());
}

// 左右/旋转几何和正负偏移都不允许末端过渡向骨干鼓出再回到轨位。
TEST_CASE("nonmotor attached corners remain tangent and monotone on both offset sides") {
    for (double offset : {-0.3, 0.3}) for (double sign : {-1.0, 1.0}) for (bool rotate : {false, true}) {
        const auto transform = [&](const isg::Vec3d& point) {
            const isg::Vec2d xy = rotate ? isg::Vec2d(-point.y()*sign, point.x()) : isg::Vec2d(point.x(), point.y()*sign);
            return isg::Vec3d(xy, point.z());
        };
        const auto first = transform(isg::Vec3d(0,offset,2)), last = transform(isg::Vec3d(12,offset,2));
        const auto entry = transform(isg::Vec3d(0,-6,1)), exit = transform(isg::Vec3d(16,4,3));
        const isg::Vec2d incoming = (first.xy()-entry.xy()).normalized();
        const isg::Vec2d outgoing = (exit.xy()-transform(isg::Vec3d(12,0,2)).xy()).normalized();
        const auto line = isg::nonmotor_detail::roundAttachedPolyline({first,last},entry,incoming,exit,outgoing,2);
        REQUIRE(line.points.size() > 6);
        CHECK(line.points.front() == entry); CHECK(line.points.back() == exit);
        const isg::Vec2d tangent = (last.xy()-first.xy()).normalized();
        size_t straight_end = 0;
        for (size_t i=1; i<line.points.size(); ++i) {
            const isg::Vec2d d = line.points[i].xy()-line.points[i-1].xy();
            if (d.norm()>0.05 && d.normalized().dot(tangent)>1-1e-9) straight_end=i;
        }
        REQUIRE(straight_end>0);
        REQUIRE(straight_end+2<line.points.size());
        double previous=0;
        const double side = isg::cross2d(tangent,exit.xy()-first.xy())>0 ? 1 : -1;
        for (size_t i=straight_end; i<line.points.size(); ++i) {
            const double distance=side*isg::cross2d(tangent,line.points[i].xy()-first.xy());
            CHECK(distance>=previous-1e-6);
            previous=distance;
        }
    }
}

// 内部骨干折角以偏移支撑线求交，不能用两端点均值破坏两侧等距平行关系。
TEST_CASE("nonmotor offset segment joins preserve both parallel supports with bounded miters") {
    for (double degrees : {-150.0,-90.0,-30.0,30.0,90.0,150.0}) for (double offset : {-0.3,0.3}) {
        const double angle=degrees*isg::DEG2RAD;
        const isg::Vec2d incoming(1,0), outgoing(std::cos(angle),std::sin(angle));
        const isg::Vec3d vertex(0,0,7), a(0,offset,7), b(isg::Vec2d(-outgoing.y(),outgoing.x())*offset,7);
        const auto joined=isg::nonmotor_detail::joinOffsetSegments(isg::Vec3d(a.xy()-incoming*10,6),a,b,
            isg::Vec3d(b.xy()+outgoing*10,8),vertex);
        CHECK(isg::cross2d(incoming,joined.xy()-a.xy())==Catch::Approx(0).margin(1e-9));
        CHECK(isg::cross2d(outgoing,joined.xy()-b.xy())==Catch::Approx(0).margin(1e-9));
        CHECK((joined.xy()-vertex.xy()).norm()<=4*std::abs(offset)+1e-9);
        CHECK(joined.z()==Catch::Approx(7));
    }
    const auto bounded=isg::nonmotor_detail::joinOffsetSegments(isg::Vec3d(-10,0.3),isg::Vec3d(0,0.3),
        isg::Vec3d(0,-0.3),isg::Vec3d(-10,-0.29999),isg::Vec3d(0,0));
    CHECK(bounded.xy().norm()<1);
}

// 用独立的Simpson速度积分验收实际站距，避免仅用相邻采样弦长冒充0.2米弧长。
TEST_CASE("nonmotor corner samples match integrated Bezier arc length") {
    for (double scale : {0.5, 2.0, 10.0}) {
        for (double degrees : {30.0, -30.0, 90.0, -90.0, 150.0, -150.0}) {
            INFO("trim=" << scale << " angle=" << degrees);
            const double angle = degrees * isg::DEG2RAD;
            const auto curve = isg::nonmotor_detail::fitCornerBezier(
                isg::Vec3d(-scale, 0, 1), isg::Vec3d(0, 0, 2),
                isg::Vec3d(scale*std::cos(angle), scale*std::sin(angle), 3));
            const auto samples = isg::nonmotor_detail::sampleCornerBezier(curve);
            const auto arc_length = [&](double t) {
                constexpr int steps = 512;
                double sum = curve.tangent(0).xy().norm() + curve.tangent(t).xy().norm();
                for (int j = 1; j < steps; ++j)
                    sum += (j % 2 ? 4 : 2) * curve.tangent(t*j/steps).xy().norm();
                return sum*t/(3*steps);
            };
            REQUIRE(samples.size() >= 2);
            double previous_arc = 0;
            for (size_t i = 1; i+1 < samples.size(); ++i) {
                // 本组有符号转角的Y单调，二分反解参数，与生产采样查找表独立。
                double low = 0, high = 1;
                const double sign = degrees > 0 ? 1 : -1;
                for (int j = 0; j < 50; ++j) {
                    const double middle = (low+high)*0.5;
                    if (sign*curve.evaluate(middle).y() < sign*samples[i].y()) low = middle;
                    else high = middle;
                }
                const double t = (low+high)*0.5, arc = arc_length(t);
                CHECK((curve.evaluate(t)-samples[i]).norm() < 0.0002);
                CHECK(arc == Catch::Approx(i*0.2).margin(0.0002));
                previous_arc = arc;
            }
            CHECK(arc_length(1)-previous_arc > 0);
            CHECK(arc_length(1)-previous_arc <= 0.2002);
        }
    }
}

// 图节点在多点Gap之后的几何段序号不同；分支挂接必须按节点坐标映射到主环。
TEST_CASE("nonmotor branch anchors preserve a multi-point road gap ring") {
    auto backbone = layoutRing();
    backbone.ring_geometry.points = {isg::Vec3d(0, 0), isg::Vec3d(5, -2), isg::Vec3d(15, -2),
        isg::Vec3d(20, 0), isg::Vec3d(20, 20), isg::Vec3d(0, 20), isg::Vec3d(0, 0)};
    for (size_t i = 0; i < 4; ++i) {
        isg::NonMotorBackboneNode node; node.id = "n" + std::to_string(i);
        node.point = i == 0 ? isg::Vec3d(0, 0) : i == 1 ? isg::Vec3d(20, 0) :
                     i == 2 ? isg::Vec3d(20, 20) : isg::Vec3d(0, 20);
        backbone.nodes.push_back(node);
    }
    isg::NonMotorBackboneEdge branch;
    branch.id = "branch"; branch.from = "n1"; branch.to = "tip"; branch.branch = true;
    branch.geometry.points = {isg::Vec3d(20, 0), isg::Vec3d(25, 5)};
    backbone.edges.push_back(branch);
    isg::IntersectionInput input;
    input.lanes.push_back(makeLane("entry", {30, 4}, {26, 4}));
    input.lanes.push_back(makeLane("exit", {10, 21}, {10, 30}));
    isg::Connectivity c; c.id = "multi-gap"; c.entry_lane_id = "entry"; c.exit_lane_id = "exit";
    c.lane_type = isg::MotorLaneType::NonMotorLane; input.connectivities.push_back(c);
    for (auto turn : {isg::ConnTurnType::TurnLeft, isg::ConnTurnType::UTurnLeft}) {
        input.connectivities.front().turn_type = turn;
        isg::NonmotorConifg config; config.corner_radius = 0; config.overlap_offset = 0;
        const auto curves = isg::NonMotorConnectivityGenerator(config).generate(input, {}, backbone);
        REQUIRE(curves.size() == 1);
        CHECK(hasPoint(curves.front().geometry, {24, 4}));
        CHECK(hasPoint(curves.front().geometry, {20, 0}));
        CHECK(hasPoint(curves.front().geometry, {20, 20}));
        CHECK(hasPoint(curves.front().geometry, {10, 20}));
        CHECK_FALSE(hasPoint(curves.front().geometry, {5, -2}));
        CHECK(curves.front().geometry.points.front() == input.lanes.front().geometry.points.back());
        CHECK(curves.front().geometry.points.back() == input.lanes.back().geometry.points.front());
    }
}

// 回归：216784982中用户标记的九条挂接曾由偏移miter/圆角二次处理形成局部尖折或自交。
TEST_CASE("216784982 selected nonmotor attachments have bounded local turns") {
    const auto input = isg::IntersectionIO::loadFromFile(
        std::string(PROJECT_ROOT_DIR) + "/datas/216784982_110004764.json");
    const auto nonmotor = isg::ConnectivityMotorLogicFilter(input, false);
    isg::NonMotorBackboneExtractor extractor(std::string(PROJECT_ROOT_DIR) + "/cmake-build-debug");
    auto backbone = extractor.run(input);
    const auto curves = isg::NonMotorConnectivityGenerator().generate(nonmotor, {}, backbone);
    for (const auto* id : {"36", "22", "30", "59", "9", "74", "67", "19", "106"}) {
        const auto curve = std::find_if(curves.begin(), curves.end(),
            [&](const isg::ConnectivityCurve& value) { return value.id == id; });
        REQUIRE(curve != curves.end());
        CHECK(curve->geometry.points.size() >= 2);
        CHECK(curve->violation.reason.find("self-intersection") == std::string::npos);
        double maximum_turn = 0.0;
        for (size_t i = 1; i + 1 < curve->geometry.points.size(); ++i) {
            const auto first = curve->geometry.points[i].xy() - curve->geometry.points[i - 1].xy();
            const auto second = curve->geometry.points[i + 1].xy() - curve->geometry.points[i].xy();
            if (first.norm() < 0.05 || second.norm() < 0.05) continue;
            maximum_turn = std::max(maximum_turn,
                std::acos(std::max(-1.0, std::min(1.0,
                    first.normalized().dot(second.normalized())))) * 180.0 / isg::PI);
        }
        INFO("attachment=" << id << " max_turn=" << maximum_turn);
        CHECK(maximum_turn < 50.0);
    }
}

// 真实回归：挂接窗口须合并临近骨干转角，不能用加密采样掩盖短弦上的近反向尖头。
TEST_CASE("216798546 110003449 attachments merge adjacent turns without a needle U-turn") {
    const auto input = isg::IntersectionIO::loadFromFile(
        std::string(PROJECT_ROOT_DIR) + "/datas/216798546_110003449.json");
    const auto nonmotor = isg::ConnectivityMotorLogicFilter(input, false);
    auto backbone = isg::NonMotorBackboneExtractor(std::string(PROJECT_ROOT_DIR) + "/cmake-build-debug").run(input);
    const auto start = std::chrono::steady_clock::now();
    const auto curves = isg::NonMotorConnectivityGenerator().generate(nonmotor, {}, backbone);
    CHECK(std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count() < 15.0);
    for (const auto* id : {"25", "3", "7", "28", "51"}) {
        const auto curve = std::find_if(curves.begin(),curves.end(),
            [&](const isg::ConnectivityCurve& c) { return c.id == id; });
        REQUIRE(curve != curves.end());
        const auto& points = curve->geometry.points;
        REQUIRE(points.size() > 2);
        const auto* entry = nonmotor.findLane(curve->entry_lane_id);
        const auto* exit = nonmotor.findLane(curve->exit_lane_id);
        REQUIRE(entry); REQUIRE(exit);
        CHECK(points.front() == entry->geometry.points.back());
        CHECK(points.back() == exit->geometry.points.front());
        double length = 0, maximum_turn = 0, maximum_curvature = 0;
        std::vector<double> stations(points.size(),0);
        for (size_t i=1; i<points.size(); ++i)
            stations[i]=stations[i-1]+(points[i].xy()-points[i-1].xy()).norm();
        length=stations.back();
        for (size_t i=1; i+1<points.size(); ++i) {
            const bool entering = std::string(id)=="25" || std::string(id)=="3";
            if ((entering ? stations[i] : length-stations[i]) > 18.0) continue;
            const isg::Vec2d a=points[i].xy()-points[i-1].xy(), b=points[i+1].xy()-points[i].xy();
            if (a.norm()<1e-7 || b.norm()<1e-7 || (a+b).norm()<1e-7) continue;
            maximum_turn=std::max(maximum_turn,std::acos(std::max(-1.0,std::min(1.0,a.normalized().dot(b.normalized()))))*180/isg::PI);
            maximum_curvature=std::max(maximum_curvature,2*std::abs(a.x()*b.y()-a.y()*b.x())/(a.norm()*b.norm()*(a+b).norm()));
        }
        INFO("attachment=" << id << " turn=" << maximum_turn << " curvature=" << maximum_curvature);
        CHECK(maximum_turn < 15.0);
        CHECK(maximum_curvature < 1.05);
        CHECK(curve->violation.reason.find("self-intersection") == std::string::npos);
        CHECK(curve->violation.reason.find("source surface") == std::string::npos);
        if (std::string(id)=="3" || std::string(id)=="7") {
            CHECK(curve->status == isg::CurveStatus::Infeasible);
            CHECK(curve->violation.reason.find("Boundary avoidance clearance") != std::string::npos);
        }
    }
}

// 同一骨干同一分段上的短非掉头挂接必须收敛为一个公共交点，避免两端短段分别偏移后形成异常折角。
TEST_CASE("216784982 short same-segment nonmotor attachments share one backbone point") {
    const auto input = isg::IntersectionIO::loadFromFile(
        std::string(PROJECT_ROOT_DIR) + "/datas/216784982_110004764.json");
    const auto nonmotor = isg::ConnectivityMotorLogicFilter(input, false);
    isg::NonMotorBackboneExtractor extractor(std::string(PROJECT_ROOT_DIR) + "/cmake-build-debug");
    auto backbone = extractor.run(input);
    isg::NonmotorConifg config;
    config.corner_radius = 0.0;
    config.overlap_offset = 0.0;
    const auto curves = isg::NonMotorConnectivityGenerator(config).generate(nonmotor, {}, backbone);
    config.short_same_segment_max_length = 0.0;
    const auto unmerged = isg::NonMotorConnectivityGenerator(config).generate(nonmotor, {}, backbone);
    for (const auto* id : {"30", "67"}) {
        const auto curve = std::find_if(curves.begin(), curves.end(),
            [&](const isg::ConnectivityCurve& value) { return value.id == id; });
        const auto original = std::find_if(unmerged.begin(), unmerged.end(),
            [&](const isg::ConnectivityCurve& value) { return value.id == id; });
        REQUIRE(curve != curves.end());
        REQUIRE(original != unmerged.end());
        REQUIRE(curve->geometry.points.size() == 3);
        CHECK(original->geometry.points.size() > curve->geometry.points.size());
    }
}

// 局部尖折修复只能替换受限窗口，不能遗漏此前的弯曲前缀、两个窗口之间的骨干或后缀。
TEST_CASE("nonmotor sharp corner repair preserves all geometry outside local windows") {
    isg::LineString2d source;
    source.points = {isg::Vec3d(0,0,1), isg::Vec3d(5,0.5,2), isg::Vec3d(10,0,3),
        isg::Vec3d(20,0,4), isg::Vec3d(20,10,5), isg::Vec3d(21,15,6),
        isg::Vec3d(20,20,7), isg::Vec3d(30,20,8), isg::Vec3d(35,21,9)};
    const auto rounded = isg::nonmotor_detail::smoothSharpCorners(source, 0.8);
    REQUIRE(rounded.points.size() > source.points.size());
    CHECK(rounded.points.front() == source.points.front());
    CHECK(rounded.points.back() == source.points.back());
    for (size_t i : {size_t(1),size_t(2),size_t(4),size_t(5),size_t(7)}) {
        INFO("unmodified source index=" << i);
        CHECK(std::find(rounded.points.begin(), rounded.points.end(), source.points[i]) != rounded.points.end());
    }
    // 长直线内精确截取0.8米，而不是吞掉整条10米邻段；切线段和插值高程均保留。
    CHECK(hasPoint(rounded, {19.2,0}));
    CHECK(hasPoint(rounded, {20,0.8}));
    CHECK(hasPoint(rounded, {20.8,20}));
    CHECK_FALSE(hasPoint(rounded, {20,0}));
    CHECK_FALSE(hasPoint(rounded, {20,20}));
    CHECK(isg::nonmotor_detail::smoothSharpCorners(source, 0).points == source.points);
    for (size_t i = 1; i+1 < rounded.points.size(); ++i) {
        const auto& point = rounded.points[i];
        if ((point.xy()-isg::Vec2d(19.2,0)).norm() < 1e-9) {
            CHECK(point.z() == Catch::Approx(3.92));
            const isg::Vec2d tangent = (rounded.points[i+1].xy()-point.xy()).normalized();
            CHECK(tangent.dot(isg::Vec2d(1,0)) > std::cos(15*isg::DEG2RAD));
        }
        if ((point.xy()-isg::Vec2d(20,0.8)).norm() < 1e-9) {
            CHECK(point.z() == Catch::Approx(4.08));
            const isg::Vec2d tangent = (point.xy()-rounded.points[i-1].xy()).normalized();
            CHECK(tangent.dot(isg::Vec2d(0,1)) > std::cos(15*isg::DEG2RAD));
        }
    }
}

// conn59需经过南侧、东侧和东北侧长骨干，不能以一条入口到末端转角的直连线替代这些段。
TEST_CASE("216784982 conn59 follows the offset backbone through its middle sections") {
    const auto input = isg::IntersectionIO::loadFromFile(
        std::string(PROJECT_ROOT_DIR)+"/datas/216784982_110004764.json");
    const auto nonmotor = isg::ConnectivityMotorLogicFilter(input, false);
    isg::NonMotorBackboneExtractor extractor(std::string(PROJECT_ROOT_DIR)+"/cmake-build-debug");
    auto backbone = extractor.run(input);
    const auto started = std::chrono::steady_clock::now();
    const auto curves = isg::NonMotorConnectivityGenerator().generate(nonmotor, {}, backbone);
    const double elapsed = std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-started).count();
    std::cout << "conn59 route regression generation: " << elapsed << " ms\n";
    CHECK(elapsed < 15000);
    const auto curve = std::find_if(curves.begin(), curves.end(), [](const auto& c) { return c.id == "59"; });
    REQUIRE(curve != curves.end());
    REQUIRE(curve->geometry.points.size() > 3);
    INFO(curve->violation.reason);
    CHECK(curve->status == isg::CurveStatus::OK);
    CHECK(curve->geometry.points.front() == nonmotor.findLane(curve->entry_lane_id)->geometry.points.back());
    CHECK(curve->geometry.points.back() == nonmotor.findLane(curve->exit_lane_id)->geometry.points.front());
    // 从原始有效骨干长段取独立检验点，远离圆角窗口；2米容差覆盖轨位分配变化。
    for (const auto& point : {isg::Vec2d(25.2446663212,-27.0936451814),
                             isg::Vec2d(32.0190442820,-2.4913494480),
                             isg::Vec2d(31.2237951241,19.2529395952)}) {
        double nearest = std::numeric_limits<double>::infinity();
        for (size_t i = 1; i < curve->geometry.points.size(); ++i)
            nearest = std::min(nearest, isg::pointToSegment(point,
                curve->geometry.points[i-1],curve->geometry.points[i]).first);
        INFO("backbone checkpoint=" << point.transpose() << " distance=" << nearest);
        CHECK(nearest < 2.0);
    }
    CHECK(curve->violation.reason.find("self-intersection") == std::string::npos);
    for (size_t i = 0; i < curves.size(); ++i) for (size_t j = i+1; j < curves.size(); ++j)
        CHECK(sharedLength(effectiveBody(curves[i]),effectiveBody(curves[j])) < 1e-6);

}

// 来源面凹口/细洞必须按精确交点切开有效范围，不能只检查采样点或两端点。
TEST_CASE("nonmotor source surface intervals reject holes and retain shifted boundary intersections") {
    isg::Polygon2d polygon;
    polygon.outer = {{0,-1},{10,-1},{10,1},{0,1}};
    isg::nonmotor_detail::SourceSurfaceSpan span;
    span.polygon = &polygon; span.origin = {0,0}; span.tangent = {1,0};
    span.first = 2; span.last = 8; span.clearance = 0.3;
    CHECK(isg::nonmotor_detail::sourceOffsetFits(span,0.6));
    CHECK_FALSE(isg::nonmotor_detail::sourceOffsetFits(span,0.8));
    auto intervals = isg::nonmotor_detail::sourceSurfaceIntervals({0,0.6},{1,0},polygon);
    REQUIRE(intervals.size() == 1);
    CHECK(intervals[0].first == Catch::Approx(0));
    CHECK(intervals[0].second == Catch::Approx(10));
    polygon.holes.push_back({{4.501,0.4},{4.509,0.4},{4.509,0.7},{4.501,0.7}});
    CHECK_FALSE(isg::nonmotor_detail::sourceOffsetFits(span,0.6));
    CHECK(isg::nonmotor_detail::sourceOffsetFits(span,0));
    polygon.holes.clear();
    polygon.outer = {{0,-1},{10,-1},{10,1},{5.01,1},{5.01,0.2},{4.99,0.2},{4.99,1},{0,1}};
    CHECK_FALSE(isg::nonmotor_detail::sourceOffsetFits(span,0.6));
}

// 窄来源面内整组缩小轨距，左右分布且不重叠；面外Gap可以存在，不能以面范围裁掉它们。
TEST_CASE("nonmotor tracks fit a narrow owned surface without clamping onto the same boundary") {
    isg::IntersectionInput input;
    for (const auto& id : {"a","b","c"}) addBottomConnection(input,id,2,18);
    isg::Crosswalk crosswalk; crosswalk.id = "narrow";
    crosswalk.geometry.outer = {{0,-0.2},{20,-0.2},{20,0.2},{0,0.2}};
    input.crosswalks.push_back(crosswalk);
    auto backbone = layoutRing();
    isg::NonMotorBackboneEdge edge;
    edge.id = "owned"; edge.from = "n0"; edge.to = "n1"; edge.ring = true;
    edge.source_surface_ids = {"narrow"}; edge.geometry.points = {{0,0},{20,0}};
    backbone.edges.push_back(edge);
    isg::NonmotorConifg config; config.source_surface_clearance = 0.05;
    const auto curves = isg::NonMotorConnectivityGenerator(config).generate(input,{},backbone);
    REQUIRE(curves.size() == 3);
    std::vector<double> offsets;
    for (const auto& c : curves) {
        const double offset = ordinateAtX(c.geometry,10); offsets.push_back(offset);
        CHECK(std::isfinite(offset)); CHECK(std::abs(offset) <= 0.15001);
        INFO(c.violation.reason);
        CHECK(c.status == isg::CurveStatus::OK);
    }
    std::sort(offsets.begin(),offsets.end());
    CHECK(offsets[1]-offsets[0] > 0.1);
    CHECK(offsets[2]-offsets[1] > 0.1);
    for (size_t i = 0; i < curves.size(); ++i) for (size_t j = i+1; j < curves.size(); ++j)
        CHECK(sharedLength(effectiveBody(curves[i]),effectiveBody(curves[j])) < 1e-6);
    std::reverse(input.connectivities.begin(),input.connectivities.end());
    const auto reordered = isg::NonMotorConnectivityGenerator(config).generate(input,{},backbone);
    for (size_t i = 0; i < curves.size(); ++i)
        CHECK(curves[i].geometry.points == reordered[curves.size()-1-i].geometry.points);
    // 固定形态继续保留原始XYZ并占零轨位，其余轨位在剩余面宽内同比缩小。
    input.connectivities[1].fixed_shape = true;
    input.connectivities[1].geometry.points = {{4,0,7},{16,0,8}};
    const auto fixed = isg::NonMotorConnectivityGenerator(config).generate(input,{},backbone);
    CHECK(fixed[1].geometry.points == input.connectivities[1].geometry.points);
    for (const auto& c : fixed) {
        INFO(c.violation.reason);
        CHECK(c.status == isg::CurveStatus::OK);
        CHECK(std::abs(ordinateAtX(c.geometry,10)) <= 0.15001);
    }
}

// 22/36/39在弯曲横道110003312中部必须在面内、保持边界净距且靠近原骨干。
TEST_CASE("216784982 nonmotor tracks remain inside source crosswalk 110003312") {
    const auto input = isg::IntersectionIO::loadFromFile(
        std::string(PROJECT_ROOT_DIR)+"/datas/216784982_110004764.json");
    const auto nm = isg::ConnectivityMotorLogicFilter(input,false);
    isg::NonMotorBackboneExtractor extractor(std::string(PROJECT_ROOT_DIR)+"/cmake-build-debug");
    auto backbone = extractor.run(input);
    CHECK(std::any_of(backbone.edges.begin(),backbone.edges.end(),[](const auto& edge) {
        return std::find(edge.source_surface_ids.begin(),edge.source_surface_ids.end(),"110003312") != edge.source_surface_ids.end();
    }));
    const auto surface = std::find_if(input.crosswalks.begin(),input.crosswalks.end(),[](const auto& cw) {return cw.id == "110003312";});
    REQUIRE(surface != input.crosswalks.end());
    const auto started = std::chrono::steady_clock::now();
    const auto curves = isg::NonMotorConnectivityGenerator().generate(nm,{},backbone);
    const double elapsed = std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-started).count();
    std::cout << "source surface regression generation: " << elapsed << " ms\n";
    CHECK(elapsed < 15000);
    const isg::Vec2d origin(-35.268492190613,12.575000561881);
    const isg::Vec2d tangent = (isg::Vec2d(-32.297546653549,-34.678553996475)-origin).normalized();
    for (const auto* id : {"22","36","39"}) {
        const auto c = std::find_if(curves.begin(),curves.end(),[&](const auto& curve) {return curve.id == id;});
        REQUIRE(c != curves.end()); INFO("conn=" << id);
        INFO(c->violation.reason);
        CHECK(c->status == isg::CurveStatus::OK);
        for (double station : {4.,12.,20.,28.,36.,42.}) {
            double nearest = 1e100; isg::Vec2d point(0,0);
            for (size_t i = 1; i < c->geometry.points.size(); ++i) {
                const auto a = c->geometry.points[i-1].xy(), b = c->geometry.points[i].xy();
                const double sa = (a-origin).dot(tangent), sb = (b-origin).dot(tangent);
                if (std::abs(sa-sb) < 1e-9 || station < std::min(sa,sb) || station > std::max(sa,sb)) continue;
                const isg::Vec2d p = a+(b-a)*((station-sa)/(sb-sa));
                const double distance = std::abs(isg::cross2d(tangent,p-origin));
                if (distance < nearest) {nearest = distance; point = p;}
            }
            INFO("station=" << station << " offset=" << nearest);
            REQUIRE(nearest < 1.3);
            CHECK(isg::nonmotor_detail::insideSourceSurface(point,surface->geometry));
            CHECK(isg::nonmotor_detail::surfaceBoundaryDistance(point,surface->geometry) >= 0.299);
        }
    }
}

// 面内靠近原骨干的偏好不能阻止合法避障：首轮受阻时才使用仍在原面内的横向候选。
TEST_CASE("nonmotor source surface tracks can avoid an obstacle within their safe width") {
    isg::IntersectionInput input; addBottomConnection(input,"owned",2,18);
    isg::Crosswalk crosswalk; crosswalk.id = "wide";
    crosswalk.geometry.outer = {{0,-2},{20,-2},{20,2},{0,2}};
    input.crosswalks.push_back(crosswalk);
    isg::Obstacle obstacle; obstacle.geometry.outer = {{8,-0.1},{12,-0.1},{12,0.1},{8,0.1}};
    input.obstacles.push_back(obstacle);
    auto backbone = layoutRing();
    isg::NonMotorBackboneEdge edge; edge.id = "owned"; edge.ring = true;
    edge.source_surface_ids = {"wide"}; edge.geometry.points = {{0,0},{20,0}};
    backbone.edges.push_back(edge);
    const auto curves = isg::NonMotorConnectivityGenerator().generate(input,{},backbone);
    REQUIRE(curves.size() == 1); INFO(curves[0].violation.reason);
    CHECK(curves[0].status == isg::CurveStatus::OK);
    CHECK(std::abs(ordinateAtX(curves[0].geometry,10)) > 0.1);
    CHECK(std::abs(ordinateAtX(curves[0].geometry,10)) < 1.7);
}

namespace {
// 独立验收所有Boundary类型，既检查横穿，也检查共线重叠；只允许精确车道端点接触。
std::set<std::string> crossedNonmotorBoundaries(const isg::LineString2d& line,
    const std::vector<isg::Boundary>& boundaries) {
    std::set<std::string> result;
    if (line.points.size() < 2) return result;
    for (const auto& boundary : boundaries)
        for (size_t i = 1; i < line.points.size(); ++i)
            for (size_t j = 1; j < boundary.geometry.points.size(); ++j) {
                isg::Vec2d hit;
                if (isg::segmentsIntersect(line.points[i-1].xy(),line.points[i].xy(),
                    boundary.geometry.points[j-1].xy(),boundary.geometry.points[j].xy(),&hit) &&
                    (hit-line.points.front().xy()).norm() > 1e-5 &&
                    (hit-line.points.back().xy()).norm() > 1e-5) result.insert(boundary.id);
            }
    return result;
}

// 独立线段间最短距离：检查内部交点及双方端点，避免仅测曲线采样点漏掉贴边。
double nonmotorBoundaryMinimum(const isg::LineString2d& line, const std::vector<isg::Boundary>& boundaries) {
    double nearest=std::numeric_limits<double>::infinity();
    for (const auto& edge:boundaries) for (size_t i=1;i<line.points.size();++i)
        for (size_t j=1;j<edge.geometry.points.size();++j) {
            const auto& a=line.points[i-1]; const auto& b=line.points[i];
            const auto& c=edge.geometry.points[j-1]; const auto& d=edge.geometry.points[j];
            if (isg::segmentsIntersect(a.xy(),b.xy(),c.xy(),d.xy(),nullptr)) return 0;
            nearest=std::min(nearest,std::min({isg::pointToSegment(a.xy(),c,d).first,
                isg::pointToSegment(b.xy(),c,d).first,isg::pointToSegment(c.xy(),a,b).first,
                isg::pointToSegment(d.xy(),a,b).first}));
        }
    return nearest;
}
}

// 原数据目标边缘全部为Other，不能只以RoadEdge枚举过滤；所有13条曲线均独立验收全部边缘。
TEST_CASE("216784982 nonmotor attachments never cross any boundary type") {
    const auto input = isg::IntersectionIO::loadFromFile(
        std::string(PROJECT_ROOT_DIR)+"/datas/216784982_110004764.json");
    const auto nm = isg::ConnectivityMotorLogicFilter(input,false);
    auto backbone = isg::NonMotorBackboneExtractor(std::string(PROJECT_ROOT_DIR)+"/cmake-build-debug").run(input);
    const auto start = std::chrono::steady_clock::now();
    const auto curves = isg::NonMotorConnectivityGenerator().generate(nm,{},backbone);
    const double elapsed = std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-start).count();
    std::cout << "boundary avoidance regression generation: " << elapsed << " ms\n";
    CHECK(elapsed < 15000); REQUIRE(curves.size() == 13);
    for (const auto& c : curves) {
        INFO("attachment=" << c.id << " violation=" << c.violation.reason);
        REQUIRE(c.geometry.points.size() >= 2);
        CHECK(crossedNonmotorBoundaries(c.geometry,input.boundaries).empty());
        const double clearance=nonmotorBoundaryMinimum(c.geometry,input.boundaries);
        INFO("minimum Boundary clearance="<<clearance);
        std::cout<<"attachment "<<c.id<<" minimum Boundary clearance: "<<clearance<<" m\n";
        CHECK(clearance > 1.0);
        CHECK(c.status == isg::CurveStatus::OK);
        CHECK(c.geometry.points.front() == nm.findLane(c.entry_lane_id)->geometry.points.back());
        CHECK(c.geometry.points.back() == nm.findLane(c.exit_lane_id)->geometry.points.front());
        const auto entry = nm.entryPtDir(c.entry_lane_id), exit = nm.exitPtDir(c.exit_lane_id);
        CHECK((c.geometry.points[1].xy()-c.geometry.points.front().xy()).normalized().dot(entry.second.normalized()) > std::cos(5*isg::DEG2RAD));
        CHECK((c.geometry.points.back().xy()-c.geometry.points[c.geometry.points.size()-2].xy()).normalized().dot(exit.second.normalized()) > std::cos(5*isg::DEG2RAD));
    }
    for (const auto* id : {"43101498","43101499","43100890","43100892","43100992","43100993","43100472"}) {
        const auto edge = std::find_if(input.boundaries.begin(),input.boundaries.end(),[&](const auto& value) {return value.id == id;});
        REQUIRE(edge != input.boundaries.end()); CHECK(edge->type == isg::Boundary::Type::Other);
    }
}

// 进入、退出和骨干中部同时受阻；绕行保留来源面、XYZ端点及进出切向，不允许跨洞抄近路。
TEST_CASE("nonmotor local boundary avoidance covers entry exit and owned backbone") {
    isg::IntersectionInput input; addBottomConnection(input,"local",4,16);
    input.lanes[0].geometry.points.back() = {4,-6,3};
    input.lanes[1].geometry.points.front() = {16,-6,5};
    input.area.geometry.outer = {{-5,-10},{25,-10},{25,25},{-5,25}};
    isg::Crosswalk cw; cw.id = "owned"; cw.geometry.outer = {{0,-2},{20,-2},{20,2},{0,2}};
    input.crosswalks.push_back(cw);
    auto backbone = layoutRing();
    isg::NonMotorBackboneEdge owned; owned.id="owned"; owned.source_surface_ids={"owned"};
    owned.from="n0"; owned.to="n1"; owned.ring=true; owned.geometry.points={{0,0},{20,0}};
    backbone.edges.push_back(owned);
    for (size_t k = 0; k < 3; ++k) {
        isg::Boundary boundary; boundary.id=std::to_string(k);
        boundary.type = k == 0 ? isg::Boundary::Type::RoadEdge : isg::Boundary::Type::Other;
        boundary.geometry.points = k == 0 ? std::vector<isg::Vec3d>{{3,-3},{5,-3}} :
            k == 1 ? std::vector<isg::Vec3d>{{15,-3},{17,-3}} : std::vector<isg::Vec3d>{{10,-0.5},{10,0.5}};
        input.boundaries.push_back(boundary);
    }
    auto curves = isg::NonMotorConnectivityGenerator().generate(input,{},backbone);
    REQUIRE(curves.size() == 1); INFO(curves[0].violation.reason);
    CHECK(curves[0].status == isg::CurveStatus::OK);
    CHECK(crossedNonmotorBoundaries(curves[0].geometry,input.boundaries).empty());
    CHECK(curves[0].geometry.points.front() == input.lanes[0].geometry.points.back());
    CHECK(curves[0].geometry.points.back() == input.lanes[1].geometry.points.front());
    CHECK(std::abs(ordinateAtX(curves[0].geometry,10)) < 1.7);
    // 边缘下侧被障碍封住时必须选择上侧合法绕行，不能为避边缘新增穿障。
    isg::Obstacle obstacle; obstacle.geometry.outer={{8,-2},{12,-2},{12,-0.1},{8,-0.1}};
    input.obstacles.push_back(obstacle);
    curves = isg::NonMotorConnectivityGenerator().generate(input,{},backbone);
    INFO(curves[0].violation.reason);
    CHECK(curves[0].status == isg::CurveStatus::OK);
    CHECK(crossedNonmotorBoundaries(curves[0].geometry,input.boundaries).empty());
    CHECK(ordinateAtX(curves[0].geometry,10) > 0.5);
    input.obstacles.clear();
    // 圆角关闭仍能绕行且不越界，不以强制开启圆角改写配置语义。
    isg::NonmotorConifg sharp; sharp.corner_radius=0;
    curves = isg::NonMotorConnectivityGenerator(sharp).generate(input,{},backbone);
    CHECK(crossedNonmotorBoundaries(curves[0].geometry,input.boundaries).empty());
    // 固定形态是原始数据约束：不改写它，但穿越必须作为硬失败明确报告。
    input.connectivities[0].fixed_shape=true;
    input.connectivities[0].geometry.points={{4,-6,3},{4,0,4},{16,0,4},{16,-6,5}};
    curves = isg::NonMotorConnectivityGenerator().generate(input,{},backbone);
    CHECK(curves[0].geometry.points == input.connectivities[0].geometry.points);
    CHECK(curves[0].status == isg::CurveStatus::Infeasible);
    CHECK(curves[0].violation.reason.find("Boundary") != std::string::npos);
}

// 横跨整个fence的边界没有合法绕行：有界搜索必须结束并明确硬失败，不能降级为正常穿越。
TEST_CASE("nonmotor boundary blocking the complete fence is infeasible") {
    isg::IntersectionInput input; addBottomConnection(input,"blocked",4,16);
    input.area.geometry.outer={{0,-2},{20,-2},{20,22},{0,22}};
    isg::Boundary wall; wall.id="wall"; wall.type=isg::Boundary::Type::Other;
    wall.geometry.points={{10,-2},{10,22}}; input.boundaries.push_back(wall);
    auto backbone=layoutRing();
    const auto curves=isg::NonMotorConnectivityGenerator().generate(input,{},backbone);
    REQUIRE(curves.size()==1);
    CHECK(curves[0].status==isg::CurveStatus::Infeasible);
    CHECK(curves[0].violation.reason.find("Boundary")!=std::string::npos);
    CHECK(curves[0].geometry.points.front()==input.lanes[0].geometry.points.back());
    CHECK(curves[0].geometry.points.back()==input.lanes[1].geometry.points.front());
}

// 未相交但贴边也应触发绕行；所有类型使用配置净距，最终测量完整平滑线段。
TEST_CASE("nonmotor boundary avoidance clearance is configurable and detects near misses") {
    isg::IntersectionInput input; addBottomConnection(input,"clearance",4,16);
    input.lanes[0].geometry.points.back()={4,-6,3};
    input.lanes[1].geometry.points.front()={16,-6,5};
    isg::Boundary edge; edge.id="nearby-other"; edge.type=isg::Boundary::Type::Other;
    edge.geometry.points={{8,0.6},{12,0.6}}; input.boundaries.push_back(edge);
    auto backbone=layoutRing();
    CHECK(isg::NonmotorConifg().boundary_avoidance_clearance==1.0);
    for (double minimum:{0.4,1.0,1.4}) {
        isg::NonmotorConifg config; config.boundary_avoidance_clearance=minimum;
        const auto curves=isg::NonMotorConnectivityGenerator(config).generate(input,{},backbone);
        REQUIRE(curves.size()==1); INFO("configured minimum="<<minimum<<" violation="<<curves[0].violation.reason);
        CHECK(curves[0].status==isg::CurveStatus::OK);
        CHECK(nonmotorBoundaryMinimum(curves[0].geometry,input.boundaries)>minimum);
        CHECK(curves[0].geometry.points.front()==input.lanes[0].geometry.points.back());
        CHECK(curves[0].geometry.points.back()==input.lanes[1].geometry.points.front());
    }
    // RoadEdge原净距设置更大时，绕行不能用新参数较小值覆盖更严格的要求。
    auto road_input=input; road_input.boundaries[0].type=isg::Boundary::Type::RoadEdge;
    isg::NonmotorConifg road_config; road_config.boundary_avoidance_clearance=0.4; road_config.road_edge_clearance=1.6;
    const auto road_curves=isg::NonMotorConnectivityGenerator(road_config).generate(road_input,{},backbone);
    REQUIRE(road_curves.size()==1);
    CHECK(road_curves[0].status==isg::CurveStatus::OK);
    CHECK(nonmotorBoundaryMinimum(road_curves[0].geometry,road_input.boundaries)>1.6);
    // 固定端点距边缘不足时不能偷偷移动端点或按旧0.05米兜底。
    input.lanes[0].geometry.points.back()={4,0.2,3};
    input.boundaries[0].geometry.points={{3,0.6},{5,0.6}};
    const auto curves=isg::NonMotorConnectivityGenerator().generate(input,{},backbone);
    CHECK(curves[0].status==isg::CurveStatus::Infeasible);
    CHECK(curves[0].geometry.points.front()==input.lanes[0].geometry.points.back());
    CHECK(curves[0].violation.reason.find("Boundary avoidance clearance")!=std::string::npos);
}
