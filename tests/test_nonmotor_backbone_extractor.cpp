#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <algorithm>
#include <chrono>
#include <set>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>

#include "io/iodata_json.h"
#include "nonmotor/nonmotor_backbone_extractor.h"
#include "utils/filesystem.hpp"

namespace {

/// 构造由两个端点定义的三维车道中心线，供非机动车骨干提取测试使用。
/// @param id 车道 ID。
/// @param a 车道几何起点的平面坐标。
/// @param b 车道几何终点的平面坐标。
/// @return 高程为零、几何点为 a 到 b 的车道。
isg::Lane lane(const std::string& id, const isg::Vec2d& a, const isg::Vec2d& b) {
    isg::Lane result;
    result.id = id;
    result.geometry.points.push_back(isg::Vec3d(a, 0.0));
    result.geometry.points.push_back(isg::Vec3d(b, 0.0));
    return result;
}

/// 创建一个按逆时针顶点顺序描述的闭合矩形面。
/// @param xmin 矩形最小 X 坐标。
/// @param ymin 矩形最小 Y 坐标。
/// @param xmax 矩形最大 X 坐标。
/// @param ymax 矩形最大 Y 坐标。
/// @return 外环四个顶点高程均为零的矩形多边形。
isg::Polygon2d rectangle(double xmin, double ymin, double xmax, double ymax) {
    isg::Polygon2d result;
    result.outer = {isg::Vec3d(xmin, ymin, 0.0), isg::Vec3d(xmax, ymin, 0.0),
                    isg::Vec3d(xmax, ymax, 0.0), isg::Vec3d(xmin, ymax, 0.0)};
    return result;
}

/// 构造沿 y=x 方向布置的窄矩形，便于测试斜向横道骨干线。
/// @param center_x 矩形中心 X 坐标。
/// @param center_y 矩形中心 Y 坐标。
/// @return 长轴沿 y=x 的四边形面。
isg::Polygon2d diagonalRectangle(double center_x, double center_y) {
    isg::Polygon2d result;
    result.outer = {isg::Vec3d(center_x - 2.9, center_y - 3.3, 0.0),
                    isg::Vec3d(center_x + 3.1, center_y + 2.7, 0.0),
                    isg::Vec3d(center_x + 2.9, center_y + 3.3, 0.0),
                    isg::Vec3d(center_x - 3.1, center_y - 2.7, 0.0)};
    return result;
}

/// 构造给定中心、方向及半长的窄矩形，长轴用于产生测试横道骨干线。
isg::Polygon2d orientedRectangle(const isg::Vec2d& center, const isg::Vec2d& direction,
                                 double half_length) {
    const isg::Vec2d along = direction.normalized() * half_length;
    const isg::Vec2d across(-direction.normalized().y() * 0.2,
                             direction.normalized().x() * 0.2);
    isg::Polygon2d result;
    result.outer = {isg::Vec3d(center - along - across, 0),
                    isg::Vec3d(center + along - across, 0),
                    isg::Vec3d(center + along + across, 0),
                    isg::Vec3d(center - along + across, 0)};
    return result;
}

/// 向测试输入中添加一组进入/退出车道及其直行连通关系。
/// @param input 接收测试车道组、车道和连通关系的路口输入对象。
/// @param entry_id 进入车道 ID。
/// @param exit_id 退出车道 ID。
/// @param x 两条车道中心线共用的 X 坐标。
void addStraightPair(isg::IntersectionInput& input, const std::string& entry_id,
                     const std::string& exit_id, double x) {
    input.lanes.push_back(lane(entry_id, isg::Vec2d(x, -20), isg::Vec2d(x, -5)));
    input.lanes.push_back(lane(exit_id, isg::Vec2d(x, 5), isg::Vec2d(x, 20)));
    isg::LaneGroup entry;
    entry.id = entry_id + "_group";
    entry.role = isg::GroupRole::Entry;
    entry.lanes.push_back(entry_id);
    input.lane_groups.push_back(entry);
    isg::LaneGroup exit;
    exit.id = exit_id + "_group";
    exit.role = isg::GroupRole::Exit;
    exit.lanes.push_back(exit_id);
    input.lane_groups.push_back(exit);
    isg::Connectivity connection;
    connection.id = entry_id + "_straight";
    connection.entry_lane_id = entry_id;
    connection.exit_lane_id = exit_id;
    connection.turn_type = isg::ConnTurnType::Straight;
    input.connectivities.push_back(connection);
}

/// 向测试输入中添加一组沿 X 方向通行的直行进入/退出车道。
/// @param input 接收测试车道组、车道和连通关系的路口输入对象。
/// @param entry_id 进入车道 ID。
/// @param exit_id 退出车道 ID。
/// @param y 两条车道中心线共用的 Y 坐标。
void addHorizontalStraightPair(isg::IntersectionInput& input, const std::string& entry_id,
                               const std::string& exit_id, double y) {
    input.lanes.push_back(lane(entry_id, isg::Vec2d(-20, y), isg::Vec2d(-5, y)));
    input.lanes.push_back(lane(exit_id, isg::Vec2d(5, y), isg::Vec2d(20, y)));
    isg::LaneGroup entry;
    entry.id = entry_id + "_group";
    entry.role = isg::GroupRole::Entry;
    entry.lanes.push_back(entry_id);
    input.lane_groups.push_back(entry);
    isg::LaneGroup exit;
    exit.id = exit_id + "_group";
    exit.role = isg::GroupRole::Exit;
    exit.lanes.push_back(exit_id);
    input.lane_groups.push_back(exit);
    isg::Connectivity connection;
    connection.id = entry_id + "_straight";
    connection.entry_lane_id = entry_id;
    connection.exit_lane_id = exit_id;
    connection.turn_type = isg::ConnTurnType::Straight;
    input.connectivities.push_back(connection);
}

/// 添加一个指定进入/退出角色及固定北向切向的直行车道组。
void addRoleDirection(isg::IntersectionInput& input, const std::string& id,
                      isg::GroupRole role, double endpoint_y, double direction_sign = 1.0) {
    const isg::Vec2d direction(0, direction_sign);
    const isg::Vec2d endpoint(0, endpoint_y);
    const isg::Vec2d a = role == isg::GroupRole::Entry
                             ? endpoint - 10.0 * direction : endpoint;
    const isg::Vec2d b = role == isg::GroupRole::Entry
                             ? endpoint : endpoint + 10.0 * direction;
    input.lanes.push_back(lane(id, a, b));
    isg::LaneGroup group;
    group.id = id + "_group";
    group.role = role;
    group.lanes.push_back(id);
    input.lane_groups.push_back(group);
    isg::Connectivity connection;
    connection.id = id + "_straight";
    connection.entry_lane_id = id;
    connection.exit_lane_id = id;
    connection.turn_type = isg::ConnTurnType::Straight;
    input.connectivities.push_back(connection);
}

/// 返回线段首尾点按指定方向观察时的有符号侧向关系。
double orientedSide(const isg::NonMotorBackboneEdge& edge, const isg::Vec2d& group_direction) {
    const isg::Vec2d first_from_last = edge.geometry.points.front().xy() -
                                       edge.geometry.points.back().xy();
    const isg::Vec2d direction = group_direction.normalized();
    return direction.x() * first_from_last.y() - direction.y() * first_from_last.x();
}

bool hasEndpoint(const isg::NonMotorBackboneEdge& edge, const isg::Vec2d& point) {
    return (edge.geometry.points.front().xy() - point).norm() < 1e-7 ||
           (edge.geometry.points.back().xy() - point).norm() < 1e-7;
}

/// 按边 ID 中包含的来源骨干 ID 查找原始边或平滑合并后的 Composite 边。
const isg::NonMotorBackboneEdge* edgeContainingId(
    const isg::NonMotorBackbone& backbone, const std::string& id_part) {
    const auto found = std::find_if(backbone.edges.begin(), backbone.edges.end(),
        [&](const isg::NonMotorBackboneEdge& edge) {
            return edge.id.find(id_part) != std::string::npos;
        });
    return found == backbone.edges.end() ? nullptr : &*found;
}

/// 统计包含指定来源 ID 的原始或 Composite 骨干边数量。
size_t countEdgesContainingId(const isg::NonMotorBackbone& backbone,
                              const std::string& id_part) {
    return static_cast<size_t>(std::count_if(backbone.edges.begin(), backbone.edges.end(),
        [&](const isg::NonMotorBackboneEdge& edge) {
            return edge.id.find(id_part) != std::string::npos;
        }));
}

}  // namespace

// 验证横道骨干线根据附近直行车道射线确定中心，并裁剪到横道边界。
TEST_CASE("NonMotorBackboneExtractor builds a crosswalk spine from straight rays") {
    isg::IntersectionInput input;
    input.id = "crosswalk";
    addStraightPair(input, "e_left", "x_left", -1.0);
    addStraightPair(input, "e_right", "x_right", 1.0);
    isg::Crosswalk crosswalk;
    crosswalk.id = "cw";
    crosswalk.geometry = rectangle(-2.0, -1.0, 2.0, 1.0);
    // 故意与车道垂向冲突，确认骨干方向由直行车道射线中点确定。
    crosswalk.crossing_direction = isg::Vec2d(1, 0);
    input.crosswalks.push_back(crosswalk);

    const isg::NonMotorBackbone backbone = isg::NonMotorBackboneExtractor().run(input);
    REQUIRE(backbone.edges.size() == 1);
    REQUIRE(backbone.edges.front().type == isg::NonMotorBackboneEdgeType::CrosswalkSpine);
    REQUIRE(backbone.edges.front().geometry.points.size() == 2);
    CHECK(std::min(backbone.edges.front().geometry.points.front().x(),
                   backbone.edges.front().geometry.points.back().x()) == Catch::Approx(-2.0));
    CHECK(std::max(backbone.edges.front().geometry.points.front().x(),
                   backbone.edges.front().geometry.points.back().x()) == Catch::Approx(2.0));
    CHECK(backbone.edges.front().geometry.points.front().y() == Catch::Approx(0.0));
    CHECK(backbone.edges.front().geometry.points.back().y() == Catch::Approx(0.0));
}

// 验证进入组与退出组分别按矢量左侧、右侧约定统一横道骨干首尾点。
TEST_CASE("NonMotorBackboneExtractor orients crosswalk spines by lane group role") {
    for (isg::GroupRole role : {isg::GroupRole::Entry, isg::GroupRole::Exit}) {
        for (double direction_sign : {-1.0, 1.0}) {
            isg::IntersectionInput input;
            input.id = role == isg::GroupRole::Entry ? "entry_oriented_cw" : "exit_oriented_cw";
            isg::Crosswalk crosswalk;
            crosswalk.id = "cw";
            crosswalk.geometry = rectangle(-2, -1, 2, 1);
            input.crosswalks.push_back(crosswalk);
            const double endpoint_y = (role == isg::GroupRole::Entry ? -5.0 : 5.0) * direction_sign;
            addRoleDirection(input, "lane", role, endpoint_y, direction_sign);

            const auto backbone = isg::NonMotorBackboneExtractor().run(input);
            REQUIRE(backbone.edges.size() == 1);
            const double side = orientedSide(backbone.edges.front(), isg::Vec2d(0, direction_sign));
            if (role == isg::GroupRole::Entry) CHECK(side > 0.0);
            else CHECK(side < 0.0);
        }
    }
}

// 验证只有转弯车道时仍用组角色定向，且首尾不受横道外环顶点绕向影响。
TEST_CASE("NonMotorBackboneExtractor orients crosswalks without straight lanes") {
    for (isg::GroupRole role : {isg::GroupRole::Entry, isg::GroupRole::Exit}) {
        for (double direction_sign : {-1.0, 1.0}) {
            for (bool reverse_ring : {false, true}) {
                isg::IntersectionInput input;
                input.id = "turn_only_oriented_cw";
                isg::Crosswalk crosswalk;
                crosswalk.id = "cw";
                crosswalk.geometry = rectangle(-2, -1, 2, 1);
                if (reverse_ring)
                    std::reverse(crosswalk.geometry.outer.begin(), crosswalk.geometry.outer.end());
                input.crosswalks.push_back(crosswalk);
                const double endpoint_y =
                    (role == isg::GroupRole::Entry ? -5.0 : 5.0) * direction_sign;
                addRoleDirection(input, "lane", role, endpoint_y, direction_sign);
                input.connectivities.front().turn_type = isg::ConnTurnType::TurnLeft;

                const auto backbone = isg::NonMotorBackboneExtractor().run(input);
                REQUIRE(backbone.edges.size() == 1);
                const auto& spine = backbone.edges.front();
                const double side = orientedSide(spine, isg::Vec2d(0, direction_sign));
                if (role == isg::GroupRole::Entry) CHECK(side > 0.0);
                else CHECK(side < 0.0);
                CHECK(hasEndpoint(spine, isg::Vec2d(-2.0, 0.0)));
                CHECK(hasEndpoint(spine, isg::Vec2d(2.0, 0.0)));
            }
        }
    }
}

// 验证直行基准优先于转弯候选；回退时排除延长线不穿横道或不近似垂直的邻路车道。
TEST_CASE("NonMotorBackboneExtractor filters fallback orientation candidates") {
    isg::IntersectionInput input;
    input.id = "filtered_orientation_basis";
    isg::Crosswalk crosswalk;
    crosswalk.id = "cw";
    crosswalk.geometry = rectangle(-2, -1, 2, 1);
    input.crosswalks.push_back(crosswalk);
    addRoleDirection(input, "valid_exit", isg::GroupRole::Exit, 5.0);

    SECTION("valid straight group has priority") {
        addRoleDirection(input, "turning_entry", isg::GroupRole::Entry, -2.0);
        input.connectivities.back().turn_type = isg::ConnTurnType::TurnLeft;
    }
    SECTION("turn-only fallback ignores irrelevant nearby groups") {
        input.connectivities.front().turn_type = isg::ConnTurnType::TurnRight;
        input.lanes.push_back(lane("misses_surface", isg::Vec2d(3, -10), isg::Vec2d(3, 0)));
        isg::LaneGroup missed;
        missed.id = "missed_group";
        missed.role = isg::GroupRole::Entry;
        missed.lanes.push_back("misses_surface");
        input.lane_groups.push_back(missed);
        input.lanes.push_back(lane("along_spine", isg::Vec2d(-10, -1), isg::Vec2d(0, 0)));
        isg::LaneGroup parallel;
        parallel.id = "parallel_group";
        parallel.role = isg::GroupRole::Entry;
        parallel.lanes.push_back("along_spine");
        input.lane_groups.push_back(parallel);
    }

    const auto backbone = isg::NonMotorBackboneExtractor().run(input);
    REQUIRE(backbone.edges.size() == 1);
    CHECK(orientedSide(backbone.edges.front(), isg::Vec2d(0, 1)) < 0.0);
    CHECK(hasEndpoint(backbone.edges.front(), isg::Vec2d(-2.0, 0.0)));
    CHECK(hasEndpoint(backbone.edges.front(), isg::Vec2d(2.0, 0.0)));
}

// T 形路口横道110000601附近没有直行连通；使用进入/退出组统一方向，端点随后接入道路Gap公共点。
TEST_CASE("NonMotorBackboneExtractor orients crosswalk 110000601 in the T junction") {
    const std::string path = std::string(PROJECT_ROOT_DIR) +
                             "/datas/216798546_110003457.json";
    for (bool reverse_input : {false, true}) {
        isg::IntersectionInput input = isg::IntersectionIO::loadFromFile(path);
        if (reverse_input) {
            for (auto& crosswalk : input.crosswalks)
                std::reverse(crosswalk.geometry.outer.begin(), crosswalk.geometry.outer.end());
            std::reverse(input.lane_groups.begin(), input.lane_groups.end());
        }
        const auto backbone = isg::NonMotorBackboneExtractor().run(input);
        const auto* spine = edgeContainingId(backbone, "cw_110000601_spine");
        REQUIRE(spine != nullptr);
        REQUIRE(spine->geometry.points.size() == 2);
        CHECK(backbone.closed);
        CHECK(backbone.ring_edges.size() == 6);
        CHECK(backbone.branch_edges.empty());
        // 进入组向东，首点应在北侧；退出组向西，首点应在右侧（同为北侧）。
        CHECK(orientedSide(*spine, input.entryPtDir("43106113").second) > 0.0);
        CHECK(orientedSide(*spine, input.exitPtDir("43106024").second) < 0.0);
        // 三段道路Gap启用后，骨干端点应与0.9米内移线的延长公共点一致。
        CHECK((spine->geometry.points.front().xy() -
               isg::Vec2d(-21.6942217829, 14.2306013037)).norm() < 1e-7);
        CHECK((spine->geometry.points.back().xy() -
               isg::Vec2d(-32.7433199376, -17.6723321047)).norm() < 1e-7);
        CHECK(spine->geometry.points.front().z() == Catch::Approx(446.785581821));
        CHECK(spine->geometry.points.back().z() == Catch::Approx(446.671696357));
    }
}

// 110003457 的三条连续RoadEdge分别承担三段骨干间隙；不能因端点直线距离较长、
// 道路线中段远离两个端点，或骨干端点与边缘的合法公共接触而退化为gap_ring直线。
TEST_CASE("NonMotorBackboneExtractor uses all long road edge gaps in 110003457") {
    const auto input = isg::IntersectionIO::loadFromFile(
        std::string(PROJECT_ROOT_DIR) + "/datas/216798546_110003457.json");
    const auto backbone = isg::NonMotorBackboneExtractor(
        std::string(PROJECT_ROOT_DIR) + "/cmake-build-debug").run(input);
    std::vector<const isg::NonMotorBackboneEdge*> road_gaps;
    for (const auto& edge : backbone.edges)
        if (edge.id.find("gap_road_") == 0) road_gaps.push_back(&edge);
    REQUIRE(road_gaps.size() >= 3);
    CHECK(backbone.connected_components == 1);
    CHECK(backbone.diagnostics.empty());
    for (const auto* boundary_id : {"43102449", "288", "43247810"}) {
        const auto boundary = std::find_if(input.boundaries.begin(), input.boundaries.end(),
            [&](const isg::Boundary& value) { return value.id == boundary_id; });
        REQUIRE(boundary != input.boundaries.end());
        double nearest_gap = std::numeric_limits<double>::infinity();
        for (const auto* gap : road_gaps) {
            for (size_t i = 1; i < gap->geometry.points.size(); ++i) {
                const auto midpoint = 0.5 * (gap->geometry.points[i - 1] + gap->geometry.points[i]);
                for (size_t j = 1; j < boundary->geometry.points.size(); ++j)
                    nearest_gap = std::min(nearest_gap,
                        isg::pointToSegment(midpoint, boundary->geometry.points[j - 1],
                                            boundary->geometry.points[j]).first);
            }
        }
        INFO("road edge=" << boundary_id << " nearest inward gap distance=" << nearest_gap);
        CHECK(nearest_gap == Catch::Approx(0.9).margin(0.08));
    }
    CHECK(std::none_of(backbone.edges.begin(), backbone.edges.end(),
        [](const isg::NonMotorBackboneEdge& edge) { return edge.id.find("gap_ring_") == 0; }));
}

// 验证安全岛单方向组骨干同样按进入/退出组方向确定首尾。
TEST_CASE("NonMotorBackboneExtractor orients island spines by lane group role") {
    for (isg::GroupRole role : {isg::GroupRole::Entry, isg::GroupRole::Exit}) {
        for (double direction_sign : {-1.0, 1.0}) {
            isg::IntersectionInput input;
            input.id = role == isg::GroupRole::Entry ? "entry_oriented_island" : "exit_oriented_island";
            isg::SafetyIsland island;
            island.id = "island";
            island.geometry = rectangle(-2, -1, 2, 1);
            input.safety_islands.push_back(island);
            const double endpoint_y = (role == isg::GroupRole::Entry ? -2.0 : 2.0) * direction_sign;
            addRoleDirection(input, "lane", role, endpoint_y, direction_sign);

            const auto backbone = isg::NonMotorBackboneExtractor().run(input);
            const auto spine = std::find_if(backbone.edges.begin(), backbone.edges.end(),
                [](const isg::NonMotorBackboneEdge& edge) {
                    return edge.type == isg::NonMotorBackboneEdgeType::SafetyIslandSpine;
                });
            REQUIRE(spine != backbone.edges.end());
            const double side = orientedSide(*spine, isg::Vec2d(0, direction_sign));
            if (role == isg::GroupRole::Entry) CHECK(side > 0.0);
            else CHECK(side < 0.0);
        }
    }
}

// 验证连接岛线沿来源横道骨干向量排列；若环平滑合并，则检查合并路径仍保留该方向。
TEST_CASE("NonMotorBackboneExtractor orients connected island spines along their source spine") {
    isg::IntersectionInput input;
    input.id = "source_oriented_island";
    isg::SafetyIsland island;
    island.id = "island";
    island.geometry = rectangle(-1, -1, 1, 1);
    input.safety_islands.push_back(island);
    isg::Crosswalk crosswalk;
    crosswalk.id = "cw";
    crosswalk.geometry = orientedRectangle(isg::Vec2d(3.5, 0), isg::Vec2d(1, 0), 1.5);
    input.crosswalks.push_back(crosswalk);

    const auto backbone = isg::NonMotorBackboneExtractor().run(input);
    const isg::NonMotorBackboneEdge* crosswalk_spine = edgeContainingId(backbone, "cw_cw_spine");
    const isg::NonMotorBackboneEdge* island_spine = edgeContainingId(backbone, "si_island_spine");
    REQUIRE(crosswalk_spine != nullptr);
    REQUIRE(island_spine != nullptr);
    const isg::Vec2d source_direction = crosswalk_spine->geometry.points.back().xy() -
                                        crosswalk_spine->geometry.points.front().xy();
    const isg::Vec2d island_direction = island_spine->geometry.points.back().xy() -
                                        island_spine->geometry.points.front().xy();
    CHECK(island_direction.dot(source_direction) > 0.0);
}

// 验证多个骨干端点按逆时针顺序组织，并由 Gap 边闭合主环。
TEST_CASE("NonMotorBackboneExtractor orders multiple spines and closes the gap ring") {
    isg::IntersectionInput input;
    input.id = "ring";
    addStraightPair(input, "e_left", "x_left", -1.0);
    addStraightPair(input, "e_right", "x_right", 1.0);
    isg::Crosswalk lower;
    lower.id = "lower";
    lower.geometry = rectangle(-2.0, -1.0, 2.0, 1.0);
    lower.crossing_direction = isg::Vec2d(0, 1);
    input.crosswalks.push_back(lower);
    isg::Crosswalk upper;
    upper.id = "upper";
    upper.geometry = rectangle(-2.0, 7.0, 2.0, 9.0);
    upper.crossing_direction = isg::Vec2d(0, 1);
    input.crosswalks.push_back(upper);

    const isg::NonMotorBackbone backbone = isg::NonMotorBackboneExtractor().run(input);
    REQUIRE(backbone.closed);
    REQUIRE(backbone.ring_nodes.size() == 4);
    CHECK(backbone.ring_edges.size() == 4);
    CHECK(std::count_if(backbone.edges.begin(), backbone.edges.end(),
                        [](const isg::NonMotorBackboneEdge& edge) {
                            return edge.type == isg::NonMotorBackboneEdgeType::Gap;
                        }) == 2);
    REQUIRE(backbone.ring_edges.size() == backbone.ring_nodes.size());
    for (size_t i = 0; i < backbone.ring_edges.size(); ++i) {
        const auto edge = std::find_if(backbone.edges.begin(), backbone.edges.end(),
            [&](const isg::NonMotorBackboneEdge& candidate) {
                return candidate.id == backbone.ring_edges[i];
            });
        REQUIRE(edge != backbone.edges.end());
        CHECK(edge->from == backbone.ring_nodes[i]);
    }
    REQUIRE(backbone.ring_geometry.points.size() == 5);
    CHECK((backbone.ring_geometry.points.front().xy() -
           backbone.ring_geometry.points.back().xy()).norm() < 1e-9);
}

// 验证孤立端点延长交点落在来源对象面内时，用交点合并两个骨干端点。
TEST_CASE("NonMotorBackboneExtractor snaps a nearby gap to an in-surface spine intersection") {
    isg::IntersectionInput input;
    input.id = "nearby_gap_in_surface";
    isg::Crosswalk horizontal;
    horizontal.id = "horizontal";
    horizontal.geometry = rectangle(-3.0, -0.5, 0.0, 0.5);
    input.crosswalks.push_back(horizontal);
    isg::Crosswalk vertical;
    vertical.id = "vertical";
    vertical.geometry = rectangle(-0.5, 1.0, 0.5, 3.5);
    input.crosswalks.push_back(vertical);

    const isg::NonMotorBackbone backbone = isg::NonMotorBackboneExtractor().run(input);
    const auto horizontal_spine = std::find_if(backbone.edges.begin(), backbone.edges.end(),
        [](const isg::NonMotorBackboneEdge& edge) { return edge.id == "cw_horizontal_spine"; });
    const auto vertical_spine = std::find_if(backbone.edges.begin(), backbone.edges.end(),
        [](const isg::NonMotorBackboneEdge& edge) { return edge.id == "cw_vertical_spine"; });
    REQUIRE(horizontal_spine != backbone.edges.end());
    REQUIRE(vertical_spine != backbone.edges.end());

    std::string shared_node;
    for (const std::string& horizontal_node : {horizontal_spine->from, horizontal_spine->to}) {
        if (horizontal_node == vertical_spine->from || horizontal_node == vertical_spine->to)
            shared_node = horizontal_node;
    }
    REQUIRE_FALSE(shared_node.empty());
    const auto joint = std::find_if(backbone.nodes.begin(), backbone.nodes.end(),
        [&](const isg::NonMotorBackboneNode& node) { return node.id == shared_node; });
    REQUIRE(joint != backbone.nodes.end());
    CHECK((joint->point.xy() - isg::Vec2d(0.0, 0.0)).norm() < 1e-7);
    CHECK(backbone.ring_nodes.size() == 3);
}

// 验证长度超过2米且延长交点不在任一来源面内时用最近端点直接连接，并按骨干方向定向。
TEST_CASE("NonMotorBackboneExtractor keeps the directed endpoint gap when extensions miss both surfaces") {
    isg::IntersectionInput input;
    input.id = "nearby_gap_outside_surfaces";
    isg::Crosswalk horizontal;
    horizontal.id = "horizontal";
    horizontal.geometry = rectangle(-3.0, -0.5, 0.0, 0.5);
    input.crosswalks.push_back(horizontal);
    isg::Crosswalk vertical;
    vertical.id = "vertical";
    vertical.geometry = rectangle(1.2, 1.5, 1.8, 4.0);
    input.crosswalks.push_back(vertical);

    const isg::NonMotorBackbone backbone = isg::NonMotorBackboneExtractor().run(input);
    const auto gap = std::find_if(backbone.edges.begin(), backbone.edges.end(),
        [](const isg::NonMotorBackboneEdge& edge) {
            if (edge.type != isg::NonMotorBackboneEdgeType::Gap || edge.geometry.points.size() != 2)
                return false;
            if (edge.id.find("gap_near_") != 0) return false;
            const isg::Vec2d first = edge.geometry.points.front().xy();
            const isg::Vec2d last = edge.geometry.points.back().xy();
            return ((first - isg::Vec2d(0.0, 0.0)).norm() < 1e-7 &&
                    (last - isg::Vec2d(1.5, 1.5)).norm() < 1e-7) ||
                   ((last - isg::Vec2d(0.0, 0.0)).norm() < 1e-7 &&
                    (first - isg::Vec2d(1.5, 1.5)).norm() < 1e-7);
        });
    REQUIRE(gap != backbone.edges.end());
    CHECK((gap->geometry.points.front().xy() - isg::Vec2d(0.0, 0.0)).norm() < 1e-7);
    CHECK((gap->geometry.points.back().xy() - isg::Vec2d(1.5, 1.5)).norm() < 1e-7);
    CHECK(gap->from != gap->to);
}

// 验证端点连线严格小于1米时，直接用两条骨干端点延长线交点合并公共节点，
// 即使交点落在两个来源面外也不再生成过短的局部 Gap。
TEST_CASE("NonMotorBackboneExtractor snaps a submeter gap to the extension intersection") {
    isg::IntersectionInput input;
    input.id = "submeter_gap_extension_intersection";
    isg::Crosswalk horizontal;
    horizontal.id = "horizontal";
    horizontal.geometry = rectangle(-3.0, -0.5, 0.0, 0.5);
    input.crosswalks.push_back(horizontal);
    isg::Crosswalk vertical;
    vertical.id = "vertical";
    vertical.geometry = rectangle(0.2, 0.2, 0.8, 3.0);
    input.crosswalks.push_back(vertical);

    const isg::NonMotorBackbone backbone = isg::NonMotorBackboneExtractor().run(input);
    const auto horizontal_spine = std::find_if(backbone.edges.begin(), backbone.edges.end(),
        [](const isg::NonMotorBackboneEdge& edge) { return edge.id == "cw_horizontal_spine"; });
    const auto vertical_spine = std::find_if(backbone.edges.begin(), backbone.edges.end(),
        [](const isg::NonMotorBackboneEdge& edge) { return edge.id == "cw_vertical_spine"; });
    REQUIRE(horizontal_spine != backbone.edges.end());
    REQUIRE(vertical_spine != backbone.edges.end());
    CHECK(std::none_of(backbone.edges.begin(), backbone.edges.end(),
        [](const isg::NonMotorBackboneEdge& edge) { return edge.id.find("gap_near_") == 0; }));

    std::string shared_node;
    for (const std::string& horizontal_node : {horizontal_spine->from, horizontal_spine->to}) {
        if (horizontal_node == vertical_spine->from || horizontal_node == vertical_spine->to)
            shared_node = horizontal_node;
    }
    REQUIRE_FALSE(shared_node.empty());
    const auto joint = std::find_if(backbone.nodes.begin(), backbone.nodes.end(),
        [&](const isg::NonMotorBackboneNode& node) { return node.id == shared_node; });
    REQUIRE(joint != backbone.nodes.end());
    CHECK((joint->point.xy() - isg::Vec2d(0.5, 0.0)).norm() < 1e-7);
}

// 验证2米边界：略小于2米时合并面外交点，恰好2米时保留原有局部 Gap。
TEST_CASE("NonMotorBackboneExtractor applies a strict two meter gap threshold") {
    const double endpoint_y = GENERATE(1.59, 1.6);
    isg::IntersectionInput input;
    input.id = "two_meter_gap_extension_intersection";
    isg::Crosswalk horizontal;
    horizontal.id = "horizontal";
    horizontal.geometry = rectangle(-3.0, -0.5, 0.0, 0.5);
    input.crosswalks.push_back(horizontal);
    isg::Crosswalk vertical;
    vertical.id = "vertical";
    vertical.geometry = rectangle(0.9, endpoint_y, 1.5, 4.0);
    input.crosswalks.push_back(vertical);

    const isg::NonMotorBackbone backbone = isg::NonMotorBackboneExtractor().run(input);
    const bool has_gap = std::any_of(backbone.edges.begin(), backbone.edges.end(),
        [](const isg::NonMotorBackboneEdge& edge) { return edge.id.find("gap_near_") == 0; });
    CHECK(has_gap == (endpoint_y == 1.6));
    const auto* horizontal_spine = edgeContainingId(backbone, "cw_horizontal_spine");
    const auto* vertical_spine = edgeContainingId(backbone, "cw_vertical_spine");
    REQUIRE(horizontal_spine != nullptr);
    REQUIRE(vertical_spine != nullptr);
    if (endpoint_y < 1.6) {
        CHECK(horizontal_spine->to == vertical_spine->from);
        CHECK(hasEndpoint(*horizontal_spine, isg::Vec2d(1.2, 0.0)));
        CHECK(hasEndpoint(*vertical_spine, isg::Vec2d(1.2, 0.0)));
    }
}

// 验证短于2米时允许单侧或双侧反向交点，且高程沿两条原骨干线外推后取均值。
TEST_CASE("NonMotorBackboneExtractor allows reverse intersections for short gaps") {
    const double vertical_x = GENERATE(-0.8, 1.5);
    isg::IntersectionInput input;
    input.id = "short_gap_reverse_intersection";
    isg::Crosswalk horizontal;
    horizontal.id = "horizontal";
    horizontal.geometry = rectangle(-5.0, -0.1, 0.0, 0.1);
    for (auto& point : horizontal.geometry.outer) point.z() = 10.0 + point.x();
    input.crosswalks.push_back(horizontal);
    isg::Crosswalk vertical;
    vertical.id = "vertical";
    vertical.geometry = rectangle(vertical_x - 0.2, -0.4, vertical_x + 0.2, 3.0);
    for (auto& point : vertical.geometry.outer) point.z() = 20.0 + 2.0 * point.y();
    input.crosswalks.push_back(vertical);

    const auto backbone = isg::NonMotorBackboneExtractor().run(input);
    const auto* first = edgeContainingId(backbone, "cw_horizontal_spine");
    const auto* second = edgeContainingId(backbone, "cw_vertical_spine");
    REQUIRE(first != nullptr);
    REQUIRE(second != nullptr);
    CHECK(first->to == second->from);
    CHECK(hasEndpoint(*first, isg::Vec2d(vertical_x, 0.0)));
    CHECK(hasEndpoint(*second, isg::Vec2d(vertical_x, 0.0)));
    CHECK(first->geometry.points.back().z() == Catch::Approx(15.0 + 0.5 * vertical_x));
    CHECK(std::none_of(backbone.edges.begin(), backbone.edges.end(),
        [](const isg::NonMotorBackboneEdge& edge) { return edge.id.find("gap_near_") == 0; }));
}

// 验证超过短Gap阈值后，原有面内交点仍必须满足两端向外延长的方向限制。
TEST_CASE("NonMotorBackboneExtractor keeps long gaps with reverse intersections") {
    isg::IntersectionInput input;
    input.id = "long_gap_reverse_intersection";
    isg::Crosswalk horizontal;
    horizontal.id = "horizontal";
    horizontal.geometry = rectangle(-5.0, -0.1, 0.0, 0.1);
    input.crosswalks.push_back(horizontal);
    isg::Crosswalk vertical;
    vertical.id = "vertical";
    vertical.geometry = rectangle(1.9, -0.4, 2.3, 3.0);
    input.crosswalks.push_back(vertical);

    const auto backbone = isg::NonMotorBackboneExtractor().run(input);
    const auto* first = edgeContainingId(backbone, "cw_horizontal_spine");
    const auto* second = edgeContainingId(backbone, "cw_vertical_spine");
    REQUIRE(first != nullptr);
    REQUIRE(second != nullptr);
    CHECK(first->to != second->from);
    CHECK(hasEndpoint(*first, isg::Vec2d(0.0, 0.0)));
    CHECK(hasEndpoint(*second, isg::Vec2d(2.1, -0.4)));
    CHECK(std::any_of(backbone.edges.begin(), backbone.edges.end(),
        [](const isg::NonMotorBackboneEdge& edge) { return edge.id.find("gap_near_") == 0; }));
}

// 验证对侧端点仅空间接近但没有骨干端点连通路径时，仍按原逻辑生成局部Gap。
TEST_CASE("NonMotorBackboneExtractor keeps a gap when the opposite endpoint is only spatially nearby") {
    isg::IntersectionInput input;
    input.id = "skip_redundant_local_gap";
    isg::Crosswalk first;
    first.id = "first";
    first.geometry = rectangle(-2.0, -0.2, 0.0, 0.2);
    input.crosswalks.push_back(first);
    isg::Crosswalk second;
    second.id = "second";
    second.geometry = rectangle(0.0, 0.3, 2.0, 0.7);
    input.crosswalks.push_back(second);

    const isg::NonMotorBackbone backbone = isg::NonMotorBackboneExtractor().run(input);
    const auto first_spine = std::find_if(backbone.edges.begin(), backbone.edges.end(),
        [](const isg::NonMotorBackboneEdge& edge) { return edge.id == "cw_first_spine"; });
    const auto second_spine = std::find_if(backbone.edges.begin(), backbone.edges.end(),
        [](const isg::NonMotorBackboneEdge& edge) { return edge.id == "cw_second_spine"; });
    REQUIRE(first_spine != backbone.edges.end());
    REQUIRE(second_spine != backbone.edges.end());
    CHECK(std::any_of(backbone.edges.begin(), backbone.edges.end(),
        [](const isg::NonMotorBackboneEdge& edge) { return edge.id.find("gap_near_") == 0; }));
    CHECK((first_spine->geometry.points.back().xy() - isg::Vec2d(0.0, 0.0)).norm() < 1e-7);
    CHECK((second_spine->geometry.points.front().xy() - isg::Vec2d(0.0, 0.5)).norm() < 1e-7);
}

// 验证对侧端点可经多条骨干线和共用端点到达最近候选点时，跳过局部Gap并保留原端点。
TEST_CASE("NonMotorBackboneExtractor skips a gap reachable through multiple connected spines") {
    isg::IntersectionInput input;
    input.id = "skip_gap_through_spine_connectivity";
    isg::Crosswalk first;
    first.id = "first";
    first.geometry = rectangle(-2.0, -0.2, 0.0, 0.2);
    input.crosswalks.push_back(first);

    isg::Crosswalk first_bridge;
    first_bridge.id = "first_bridge";
    first_bridge.geometry = orientedRectangle(isg::Vec2d(-1.5, 0.5), isg::Vec2d(1.0, 1.0),
                                              std::sqrt(2.0) / 2.0);
    input.crosswalks.push_back(first_bridge);
    isg::Crosswalk second_bridge;
    second_bridge.id = "second_bridge";
    second_bridge.geometry = rectangle(-1.0, 0.8, 0.0, 1.2);
    input.crosswalks.push_back(second_bridge);

    const isg::NonMotorBackbone backbone = isg::NonMotorBackboneExtractor().run(input);
    const isg::NonMotorBackboneEdge* first_spine = edgeContainingId(backbone, "cw_first_spine");
    const isg::NonMotorBackboneEdge* first_bridge_spine =
        edgeContainingId(backbone, "cw_first_bridge_spine");
    const isg::NonMotorBackboneEdge* second_bridge_spine =
        edgeContainingId(backbone, "cw_second_bridge_spine");
    REQUIRE(first_spine != nullptr);
    REQUIRE(first_bridge_spine != nullptr);
    REQUIRE(second_bridge_spine != nullptr);
    CHECK(std::none_of(backbone.edges.begin(), backbone.edges.end(),
        [](const isg::NonMotorBackboneEdge& edge) { return edge.id.find("gap_near_") == 0; }));
    CHECK(hasEndpoint(*first_spine, isg::Vec2d(0.0, 0.0)));
    CHECK(hasEndpoint(*first_bridge_spine, isg::Vec2d(-1.0, 1.0)));
    CHECK(hasEndpoint(*second_bridge_spine, isg::Vec2d(0.0, 1.0)));
}

// 验证真实路口短Gap反向求交合并后保持闭环、端点一致与环/分支互斥，且不补出重复环Gap。
TEST_CASE("NonMotorBackboneExtractor does not duplicate connected endpoint ring gaps") {
    const std::string path = std::string(PROJECT_ROOT_DIR) +
                             "/datas/216784982_110004764.json";
    const isg::IntersectionInput input = isg::IntersectionIO::loadFromFile(path);
    const isg::NonMotorBackbone backbone = isg::NonMotorBackboneExtractor().run(input);
    CHECK(std::none_of(backbone.edges.begin(), backbone.edges.end(),
        [](const isg::NonMotorBackboneEdge& edge) {
            return edge.id.find("gap_near_") == 0;
        }));
    const auto* right_lower = edgeContainingId(backbone, "cw_110001585_spine");
    const auto* right_upper = edgeContainingId(backbone, "cw_110000726_spine");
    const auto* left = edgeContainingId(backbone, "cw_110003312_spine");
    const auto* lower = edgeContainingId(backbone, "cw_110002216_spine");
    REQUIRE(right_lower != nullptr);
    REQUIRE(right_upper != nullptr);
    REQUIRE(left != nullptr);
    REQUIRE(lower != nullptr);
    // 10度共点合并阈值和道路Gap可保留更多有效分段；不再绑定旧30度配置的Composite ID。
    REQUIRE(edgeContainingId(backbone, "si_110004263_spine_1") != nullptr);
    REQUIRE(edgeContainingId(backbone, "si_110004294_spine_2") != nullptr);
    CHECK(backbone.connected_components == 1);
    CHECK(backbone.diagnostics.empty());
    CHECK(std::any_of(backbone.edges.begin(), backbone.edges.end(), [](const isg::NonMotorBackboneEdge& e) {
        return e.id.find("gap_road_component_") == 0;
    }));
    CHECK(right_lower->to == right_upper->from);
    CHECK(hasEndpoint(*right_lower, isg::Vec2d(35.0731030380, 12.7233487584)));
    CHECK(hasEndpoint(*right_upper, isg::Vec2d(35.0731030380, 12.7233487584)));
    CHECK(backbone.ring_edges.size() >= 7);
    REQUIRE(backbone.ring_edges.size() == backbone.ring_nodes.size());
    REQUIRE(backbone.ring_geometry.points.size() >= backbone.ring_nodes.size() + 1);
    // 拓扑节点序与完整边几何分别核查，多点Gap不可被当作单段压平。
    for (size_t i = 0; i < backbone.ring_nodes.size(); ++i) {
        const auto edge = std::find_if(backbone.edges.begin(), backbone.edges.end(), [&](const isg::NonMotorBackboneEdge& e) {
            return e.id == backbone.ring_edges[i];
        });
        REQUIRE(edge != backbone.edges.end());
        const auto& from = backbone.ring_nodes[i];
        const auto& to = backbone.ring_nodes[(i + 1) % backbone.ring_nodes.size()];
        CHECK(((edge->from == from && edge->to == to) || (edge->from == to && edge->to == from)));
        CHECK_FALSE(edge->branch);
        for (const auto& p : edge->geometry.points)
            CHECK(std::any_of(backbone.ring_geometry.points.begin(), backbone.ring_geometry.points.end(),
                [&](const isg::Vec3d& q) { return isg::dist(p, q) < 1e-6; }));
    }
    double twice_area = 0;
    const auto& ring = backbone.ring_geometry.points;
    for (size_t i = 1; i < ring.size(); ++i) {
        twice_area += isg::cross2d(ring[i - 1].xy(), ring[i].xy());
        for (size_t j = i + 2; j < ring.size(); ++j) {
            if (i == 1 && j + 1 == ring.size()) continue;
            CHECK_FALSE(isg::segmentsIntersectStrict(ring[i - 1], ring[i], ring[j - 1], ring[j]));
        }
    }
    CHECK(twice_area > 0.0);
    CHECK(std::none_of(backbone.edges.begin(), backbone.edges.end(),
        [](const isg::NonMotorBackboneEdge& edge) {
            return edge.id.find("gap_ring_") == 0;
        }));
    CHECK(std::any_of(backbone.edges.begin(), backbone.edges.end(),
        [](const isg::NonMotorBackboneEdge& edge) {
            return edge.id == "cw_110001923_spine";
        }));
    CHECK(std::any_of(backbone.edges.begin(), backbone.edges.end(),
        [](const isg::NonMotorBackboneEdge& edge) {
            return edge.id == "cw_110003312_spine";
        }));
    CHECK(backbone.closed);
    CHECK(std::all_of(backbone.edges.begin(), backbone.edges.end(),
        [](const isg::NonMotorBackboneEdge& edge) {
            return edge.ring != edge.branch && edge.branch == !edge.ring;
        }));
    CHECK(std::all_of(backbone.branch_edges.begin(), backbone.branch_edges.end(),
        [&](const std::string& id) {
            return std::find(backbone.ring_edges.begin(), backbone.ring_edges.end(), id) ==
                   backbone.ring_edges.end();
        }));
    CHECK(std::none_of(backbone.ring_edges.begin(), backbone.ring_edges.end(),
        [](const std::string& id) {
            return id == "si_110004294_spine_2" || id == "cw_110001993_spine" ||
                   id == "si_110004294_spine_1" || id == "cw_110001923_spine" ||
                   id == "si_110004265_spine_0" || id == "cw_110000067_spine" ||
                   id == "si_110004263_spine_2" || id == "cw_110004685_spine";
        }));
    CHECK(backbone.ring_geometry.points.size() == backbone.ring_nodes.size() + 1);
}

// 验证孤立端点超过3米搜索半径时不生成局部连接线。
TEST_CASE("NonMotorBackboneExtractor skips isolated endpoints outside the three meter radius") {
    isg::IntersectionInput input;
    input.id = "isolated_gap_outside_radius";
    isg::Crosswalk horizontal;
    horizontal.id = "horizontal";
    horizontal.geometry = rectangle(-3.0, -0.5, 0.0, 0.5);
    input.crosswalks.push_back(horizontal);
    isg::Crosswalk vertical;
    vertical.id = "vertical";
    vertical.geometry = rectangle(-0.5, 3.1, 0.5, 5.1);
    input.crosswalks.push_back(vertical);

    const isg::NonMotorBackbone backbone = isg::NonMotorBackboneExtractor().run(input);
    const auto horizontal_spine = std::find_if(backbone.edges.begin(), backbone.edges.end(),
        [](const isg::NonMotorBackboneEdge& edge) { return edge.id == "cw_horizontal_spine"; });
    const auto vertical_spine = std::find_if(backbone.edges.begin(), backbone.edges.end(),
        [](const isg::NonMotorBackboneEdge& edge) { return edge.id == "cw_vertical_spine"; });
    REQUIRE(horizontal_spine != backbone.edges.end());
    REQUIRE(vertical_spine != backbone.edges.end());
    CHECK(std::none_of(backbone.edges.begin(), backbone.edges.end(),
        [](const isg::NonMotorBackboneEdge& edge) {
            return edge.id.find("gap_near_") == 0;
        }));
    CHECK(backbone.ring_nodes.size() == 4);
}

// 验证仅两条骨干共用且夹角小于30度时删除共点并合并线段。
TEST_CASE("NonMotorBackboneExtractor merges a shallow shared spine endpoint when boundary edges stay stable") {
    isg::IntersectionInput input;
    input.id = "smooth_shared_spine_endpoint";
    isg::Crosswalk left;
    left.id = "left";
    left.geometry = rectangle(-10.0, -0.5, 0.0, 0.5);
    input.crosswalks.push_back(left);
    isg::Crosswalk right;
    right.id = "right";
    right.geometry = rectangle(0.0, -0.5, 10.0, 0.5);
    input.crosswalks.push_back(right);
    isg::Crosswalk far;
    far.id = "far";
    far.geometry = rectangle(29.5, -10.0, 30.5, 10.0);
    input.crosswalks.push_back(far);

    const isg::NonMotorBackbone backbone = isg::NonMotorBackboneExtractor().run(input);
    CHECK(std::none_of(backbone.edges.begin(), backbone.edges.end(),
        [](const isg::NonMotorBackboneEdge& edge) {
            return edge.id == "cw_left_spine" || edge.id == "cw_right_spine";
        }));
    const auto merged = std::find_if(backbone.edges.begin(), backbone.edges.end(),
        [](const isg::NonMotorBackboneEdge& edge) {
            return edge.type == isg::NonMotorBackboneEdgeType::Composite &&
                   edge.id.find("cw_left_spine") != std::string::npos &&
                   edge.id.find("cw_right_spine") != std::string::npos;
        });
    REQUIRE(merged != backbone.edges.end());
    CHECK(backbone.ring_nodes.size() == 4);
    CHECK(std::none_of(backbone.nodes.begin(), backbone.nodes.end(),
        [](const isg::NonMotorBackboneNode& node) {
            return (node.point.xy() - isg::Vec2d(0.0, 0.0)).norm() < 1e-7;
        }));
}

// 验证只有两条骨干线共点且夹角小于30度时，即使合并线改变来源面边界交点也删除共点。
TEST_CASE("NonMotorBackboneExtractor merges a shallow shared endpoint regardless of boundary edges") {
    isg::IntersectionInput input;
    input.id = "preserve_shared_spine_endpoint";
    isg::Crosswalk left;
    left.id = "left";
    left.geometry = rectangle(-10.0, -0.5, 0.0, 0.5);
    input.crosswalks.push_back(left);
    isg::Crosswalk right;
    right.id = "right";
    const double angle = 10.0 * isg::DEG2RAD;
    right.geometry = orientedRectangle(
        isg::Vec2d(5.0 * std::cos(angle), 5.0 * std::sin(angle)),
        isg::Vec2d(std::cos(angle), std::sin(angle)), 5.0);
    input.crosswalks.push_back(right);
    isg::Crosswalk far;
    far.id = "far";
    far.geometry = rectangle(29.5, -10.0, 30.5, 10.0);
    input.crosswalks.push_back(far);

    const isg::NonMotorBackbone backbone = isg::NonMotorBackboneExtractor().run(input);
    CHECK(std::none_of(backbone.edges.begin(), backbone.edges.end(),
        [](const isg::NonMotorBackboneEdge& edge) {
            return edge.id == "cw_left_spine" || edge.id == "cw_right_spine";
        }));
    const auto merged = std::find_if(backbone.edges.begin(), backbone.edges.end(),
        [](const isg::NonMotorBackboneEdge& edge) {
            return edge.type == isg::NonMotorBackboneEdgeType::Composite &&
                   edge.id.find("cw_left_spine") != std::string::npos &&
                   edge.id.find("cw_right_spine") != std::string::npos;
        });
    REQUIRE(merged != backbone.edges.end());
    CHECK(backbone.nodes.end() == std::find_if(backbone.nodes.begin(), backbone.nodes.end(),
        [](const isg::NonMotorBackboneNode& node) {
            return (node.point.xy() - isg::Vec2d(0.0, 0.0)).norm() < 1e-7;
        }));
}

// 验证共端点夹角达到30度以上时保留原始两条骨干线，不递归合并。
TEST_CASE("NonMotorBackboneExtractor keeps a shared endpoint at thirty degrees or more") {
    isg::IntersectionInput input;
    input.id = "preserve_shared_spine_at_large_angle";
    isg::Crosswalk left;
    left.id = "left";
    left.geometry = rectangle(-10.0, -0.5, 0.0, 0.5);
    input.crosswalks.push_back(left);
    isg::Crosswalk right;
    right.id = "right";
    const double angle = 35.0 * isg::DEG2RAD;
    right.geometry = orientedRectangle(
        isg::Vec2d(5.0 * std::cos(angle), 5.0 * std::sin(angle)),
        isg::Vec2d(std::cos(angle), std::sin(angle)), 5.0);
    input.crosswalks.push_back(right);

    const isg::NonMotorBackbone backbone = isg::NonMotorBackboneExtractor().run(input);
    CHECK(std::any_of(backbone.edges.begin(), backbone.edges.end(),
        [](const isg::NonMotorBackboneEdge& edge) { return edge.id == "cw_left_spine"; }));
    CHECK(std::any_of(backbone.edges.begin(), backbone.edges.end(),
        [](const isg::NonMotorBackboneEdge& edge) { return edge.id == "cw_right_spine"; }));
    CHECK(std::none_of(backbone.edges.begin(), backbone.edges.end(),
        [](const isg::NonMotorBackboneEdge& edge) {
            return edge.type == isg::NonMotorBackboneEdgeType::Composite;
        }));
}

// 验证三条骨干线共用的主环端点不进入两线平滑，即使其中两线夹角小于30度。
TEST_CASE("NonMotorBackboneExtractor leaves a three-spine shared endpoint unchanged") {
    isg::IntersectionInput input;
    input.id = "three_spines_shared_endpoint";
    isg::Crosswalk left;
    left.id = "left";
    left.geometry = rectangle(-10.0, -0.5, 0.0, 0.5);
    input.crosswalks.push_back(left);
    isg::Crosswalk shallow;
    shallow.id = "shallow";
    const double angle = 10.0 * isg::DEG2RAD;
    shallow.geometry = orientedRectangle(
        isg::Vec2d(5.0 * std::cos(angle), 5.0 * std::sin(angle)),
        isg::Vec2d(std::cos(angle), std::sin(angle)), 5.0);
    input.crosswalks.push_back(shallow);
    isg::Crosswalk vertical;
    vertical.id = "vertical";
    vertical.geometry = orientedRectangle(isg::Vec2d(0.0, 5.0), isg::Vec2d(0.0, 1.0), 5.0);
    input.crosswalks.push_back(vertical);

    const isg::NonMotorBackbone backbone = isg::NonMotorBackboneExtractor().run(input);
    REQUIRE(std::count_if(backbone.edges.begin(), backbone.edges.end(),
        [](const isg::NonMotorBackboneEdge& edge) {
            return edge.type == isg::NonMotorBackboneEdgeType::CrosswalkSpine;
        }) == 3);
    CHECK(std::none_of(backbone.edges.begin(), backbone.edges.end(),
        [](const isg::NonMotorBackboneEdge& edge) {
            return edge.type == isg::NonMotorBackboneEdgeType::Composite;
        }));
    CHECK(std::count_if(backbone.nodes.begin(), backbone.nodes.end(),
        [](const isg::NonMotorBackboneNode& node) {
            return (node.point.xy() - isg::Vec2d(0.0, 0.0)).norm() < 1e-7;
        }) == 1);
}

// 验证安全岛优先连接相邻横道，并在没有横道时按附近直行方向提取骨干线。
TEST_CASE("NonMotorBackboneExtractor creates island branch and fallback spine") {
    SECTION("island connects to the nearest crosswalk endpoint") {
        isg::IntersectionInput input;
        input.id = "island_branch";
        isg::SafetyIsland island;
        island.id = "island";
        island.geometry = rectangle(-1.0, -1.0, 1.0, 1.0);
        input.safety_islands.push_back(island);
        isg::Crosswalk crosswalk;
        crosswalk.id = "cw";
        crosswalk.geometry = rectangle(-2.0, 2.0, 2.0, 4.0);
        crosswalk.crossing_direction = isg::Vec2d(0, 1);
        input.crosswalks.push_back(crosswalk);
        addStraightPair(input, "e", "x", 0.0);

        const isg::NonMotorBackbone backbone = isg::NonMotorBackboneExtractor().run(input);
        REQUIRE(backbone.edges.size() == 2);
        CHECK(std::count_if(backbone.edges.begin(), backbone.edges.end(),
                            [](const isg::NonMotorBackboneEdge& edge) {
                                return edge.type == isg::NonMotorBackboneEdgeType::SafetyIslandSpine && edge.branch;
                            }) == 1);
    }

    SECTION("island without crosswalk uses a perpendicular lane direction") {
        isg::IntersectionInput input;
        input.id = "island_fallback";
        isg::SafetyIsland island;
        island.id = "island";
        island.geometry = rectangle(-2.0, -1.0, 2.0, 1.0);
        input.safety_islands.push_back(island);
        addStraightPair(input, "e", "x", 8.0);

        const isg::NonMotorBackbone backbone = isg::NonMotorBackboneExtractor().run(input);
        REQUIRE(backbone.edges.size() == 1);
        CHECK(backbone.edges.front().type == isg::NonMotorBackboneEdgeType::SafetyIslandSpine);
        CHECK(std::min(backbone.edges.front().geometry.points.front().x(),
                       backbone.edges.front().geometry.points.back().x()) == Catch::Approx(-2.0));
        CHECK(std::max(backbone.edges.front().geometry.points.front().x(),
                       backbone.edges.front().geometry.points.back().x()) == Catch::Approx(2.0));
        CHECK(backbone.edges.front().geometry.points.front().y() == Catch::Approx(0.0));
        CHECK(backbone.edges.front().geometry.points.back().y() == Catch::Approx(0.0));
    }
}

// 验证单横道从两侧向安全岛延伸至远侧交点；递归平滑可继续合并共点的 Composite 路径。
TEST_CASE("NonMotorBackboneExtractor extends a single crosswalk through the island") {
    for (double side : {-1.0, 1.0}) {
        isg::IntersectionInput input;
        input.id = side < 0 ? "single_island_left" : "single_island_right";
        isg::SafetyIsland island;
        island.id = "island";
        island.geometry = rectangle(-1.0, -1.0, 1.0, 1.0);
        input.safety_islands.push_back(island);
        isg::Crosswalk crosswalk;
        crosswalk.id = "cw";
        crosswalk.geometry = orientedRectangle(isg::Vec2d(3.5 * side, 0), isg::Vec2d(1, 0), 1.5);
        input.crosswalks.push_back(crosswalk);

        const isg::NonMotorBackbone backbone = isg::NonMotorBackboneExtractor().run(input);
        const isg::NonMotorBackboneEdge* island_path = edgeContainingId(backbone, "si_island_spine");
        REQUIRE(island_path != nullptr);
        REQUIRE(island_path->geometry.points.size() == 2);
        const double min_x = std::min(island_path->geometry.points.front().x(),
                                      island_path->geometry.points.back().x());
        const double max_x = std::max(island_path->geometry.points.front().x(),
                                      island_path->geometry.points.back().x());
        CHECK(min_x == Catch::Approx(side < 0 ? -5.0 : -1.0));
        CHECK(max_x == Catch::Approx(side < 0 ? 1.0 : 5.0));
        CHECK(island_path->geometry.points.front().y() == Catch::Approx(0.0).margin(1e-9));
        CHECK(island_path->geometry.points.back().y() == Catch::Approx(0.0).margin(1e-9));
    }
}

// 验证双横道缺少有效交点中心时，仅在近似平行且最近端点连线接触安全岛时直接连接。
TEST_CASE("NonMotorBackboneExtractor connects parallel crosswalk endpoints through island") {
    isg::IntersectionInput input;
    input.id = "parallel_island";
    isg::SafetyIsland island;
    island.id = "island";
    island.geometry = rectangle(-1.0, -1.0, 1.0, 1.0);
    input.safety_islands.push_back(island);
    double y = 0.0;
    double slope = 0.0;
    bool expect_direct = true;
    SECTION("missing intersection area") {}
    SECTION("parallel spines have no intersection") {
        input.area.geometry = rectangle(-1, -1, 1, 1);
    }
    SECTION("the endpoint segment misses the island") {
        y = 2.0;
        expect_direct = false;
    }
    SECTION("near parallel intersection is outside the island") {
        input.area.geometry = rectangle(-1, -1, 1, 1);
        y = 2.0;
        slope = 0.4;
    }
    SECTION("angle exceeds the parallel threshold") {
        slope = 0.8;
        expect_direct = false;
    }
    isg::Crosswalk left;
    left.id = "left";
    left.geometry = orientedRectangle(isg::Vec2d(-3.5, slope == 0.0 ? y : 0.0),
                                      isg::Vec2d(1, 0), 1.5);
    input.crosswalks.push_back(left);
    isg::Crosswalk right;
    right.id = "right";
    right.geometry = orientedRectangle(isg::Vec2d(3.5, y), isg::Vec2d(1, slope), 1.5);
    input.crosswalks.push_back(right);

    const isg::NonMotorBackbone backbone = isg::NonMotorBackboneExtractor().run(input);
    std::vector<const isg::NonMotorBackboneEdge*> island_spines;
    for (const auto& edge : backbone.edges)
        if (edge.id.find("si_island_spine") != std::string::npos)
            island_spines.push_back(&edge);
    if (expect_direct) {
        REQUIRE_FALSE(island_spines.empty());
        const auto direct_path = std::find_if(island_spines.begin(), island_spines.end(),
            [](const isg::NonMotorBackboneEdge* edge) {
                const double min_x = std::min(edge->geometry.points.front().x(),
                                              edge->geometry.points.back().x());
                const double max_x = std::max(edge->geometry.points.front().x(),
                                              edge->geometry.points.back().x());
                return min_x <= -2.0 + 1e-7 && max_x >= 2.0 - 1e-7;
            });
        REQUIRE(direct_path != island_spines.end());
        // 直接岛路径的生成阈值为30度，共点合并阈值为10度；较大夹角允许左横道单独保留。
        CHECK((*direct_path)->type == isg::NonMotorBackboneEdgeType::Composite);
        CHECK(edgeContainingId(backbone, "cw_left_spine") != nullptr);
        CHECK(edgeContainingId(backbone, "cw_right_spine") != nullptr);
    } else {
        REQUIRE(island_spines.size() == 2);
        CHECK(island_spines.front()->branch);
        CHECK(island_spines.back()->branch);
    }
}

// 验证两块横道间隔小于 1 米时，外扩相交且距离安全岛更远的横道被过滤。
TEST_CASE("NonMotorBackboneExtractor filters farther crosswalk within one meter gap") {
    isg::IntersectionInput input;
    input.id = "filter_overlapping_candidates";
    isg::SafetyIsland island;
    island.id = "island";
    island.geometry = rectangle(-1.0, -1.0, 1.0, 1.0);
    input.safety_islands.push_back(island);

    isg::Crosswalk near_crosswalk;
    near_crosswalk.id = "near";
    near_crosswalk.geometry = rectangle(-2.1, -0.2, -1.5, 0.2);
    input.crosswalks.push_back(near_crosswalk);

    isg::Crosswalk farther_crosswalk;
    farther_crosswalk.id = "farther";
    farther_crosswalk.geometry = rectangle(-3.8, -0.2, -2.7, 0.2);
    input.crosswalks.push_back(farther_crosswalk);

    const auto backbone = isg::NonMotorBackboneExtractor().run(input);
    const isg::NonMotorBackboneEdge* island_path = edgeContainingId(backbone, "si_island_spine");
    REQUIRE(island_path != nullptr);
    const double min_x = std::min(island_path->geometry.points.front().x(),
                                  island_path->geometry.points.back().x());
    const double max_x = std::max(island_path->geometry.points.front().x(),
                                  island_path->geometry.points.back().x());
    CHECK(min_x == Catch::Approx(-2.1));
    CHECK(max_x == Catch::Approx(1.0));
}

// 验证两块横道间隔超过 1 米时，不因外扩相交规则而过滤较远横道。
TEST_CASE("NonMotorBackboneExtractor keeps farther crosswalk beyond one meter gap") {
    isg::IntersectionInput input;
    input.id = "keep_beyond_merge_gap";
    isg::SafetyIsland island;
    island.id = "island";
    island.geometry = rectangle(-1.0, -1.0, 1.0, 1.0);
    input.safety_islands.push_back(island);

    isg::Crosswalk near_crosswalk;
    near_crosswalk.id = "near";
    near_crosswalk.geometry = rectangle(-1.8, -0.2, -1.2, 0.2);
    input.crosswalks.push_back(near_crosswalk);
    isg::Crosswalk farther_crosswalk;
    farther_crosswalk.id = "farther";
    farther_crosswalk.geometry = rectangle(-0.2, 2.5, 0.2, 3.5);
    input.crosswalks.push_back(farther_crosswalk);

    const auto backbone = isg::NonMotorBackboneExtractor().run(input);
    CHECK(std::count_if(backbone.edges.begin(), backbone.edges.end(),
                        [](const isg::NonMotorBackboneEdge& edge) {
                            return edge.id.find("si_island_spine") != std::string::npos;
                        }) == 2);
}

// 验证相互分离的三个邻接横道不会被两两过滤，并分别保留含安全岛来源的骨干边。
TEST_CASE("NonMotorBackboneExtractor keeps separated adjacent crosswalk candidates") {
    isg::IntersectionInput input;
    input.id = "keep_separated_candidates";
    isg::SafetyIsland island;
    island.id = "island";
    island.geometry = rectangle(-1.0, -1.0, 1.0, 1.0);
    input.safety_islands.push_back(island);

    isg::Crosswalk left;
    left.id = "left";
    left.geometry = rectangle(-3.0, -0.2, -1.5, 0.2);
    input.crosswalks.push_back(left);
    isg::Crosswalk right;
    right.id = "right";
    right.geometry = rectangle(1.5, -0.2, 3.0, 0.2);
    input.crosswalks.push_back(right);
    isg::Crosswalk upper;
    upper.id = "upper";
    upper.geometry = rectangle(-0.2, 1.5, 0.2, 3.0);
    input.crosswalks.push_back(upper);

    const auto backbone = isg::NonMotorBackboneExtractor().run(input);
    CHECK(std::count_if(backbone.edges.begin(), backbone.edges.end(),
                        [](const isg::NonMotorBackboneEdge& edge) {
                            return edge.id.find("si_island_spine") != std::string::npos;
                        }) == 3);
}

// 验证两条近似平行横道的延长交点在岛内时，中心逻辑保留且浅角公共端点合并为 Composite。
TEST_CASE("NonMotorBackboneExtractor prefers a valid intersection center over a parallel link") {
    isg::IntersectionInput input;
    input.id = "valid_parallel_intersection";
    input.area.geometry = rectangle(-1, -1, 1, 1);
    isg::SafetyIsland island;
    island.id = "island";
    island.geometry = rectangle(-1, -1, 1, 1);
    input.safety_islands.push_back(island);
    isg::Crosswalk first;
    first.id = "first";
    first.geometry = orientedRectangle(isg::Vec2d(-3.5, 0), isg::Vec2d(1, 0), 1.5);
    input.crosswalks.push_back(first);
    isg::Crosswalk second;
    second.id = "second";
    second.geometry = orientedRectangle(isg::Vec2d(3.5, 0.35), isg::Vec2d(1, 0.1), 1.5);
    input.crosswalks.push_back(second);
    const auto backbone = isg::NonMotorBackboneExtractor().run(input);
    size_t count = 0;
    for (const auto& edge : backbone.edges) {
        if (edge.id.find("si_island_spine") == std::string::npos) continue;
        ++count;
        CHECK(edge.branch);
        CHECK(edge.type == isg::NonMotorBackboneEdgeType::Composite);
    }
    CHECK(count == 2);
}

// 验证缺少路口面时，相邻横道骨干线无法按路口中心筛选，安全岛中心回退逻辑仍生成 Composite。
TEST_CASE("NonMotorBackboneExtractor falls back to island centroid without intersection area") {
    isg::IntersectionInput input;
    isg::SafetyIsland island;
    island.id = "offset_centroid_island";
    island.geometry.outer = {isg::Vec3d(-2.0, -1.0, 0.0), isg::Vec3d(2.0, -1.0, 0.0),
                             isg::Vec3d(1.0, 1.0, 0.0), isg::Vec3d(-1.0, 1.0, 0.0)};
    input.safety_islands.push_back(island);

    isg::Crosswalk horizontal_crosswalk;
    horizontal_crosswalk.id = "horizontal_crosswalk";
    horizontal_crosswalk.geometry = rectangle(-3.0, -0.5, 3.0, 0.5);
    input.crosswalks.push_back(horizontal_crosswalk);

    isg::Crosswalk vertical_crosswalk;
    vertical_crosswalk.id = "vertical_crosswalk";
    vertical_crosswalk.geometry = rectangle(-0.5, -3.0, 0.5, 3.0);
    input.crosswalks.push_back(vertical_crosswalk);

    const isg::NonMotorBackbone backbone = isg::NonMotorBackboneExtractor().run(input);
    const auto island_spine = std::find_if(backbone.edges.begin(), backbone.edges.end(),
        [](const isg::NonMotorBackboneEdge& edge) {
            return edge.id.find("si_offset_centroid_island_spine") != std::string::npos;
        });
    REQUIRE(island_spine != backbone.edges.end());
    REQUIRE(island_spine->geometry.points.size() == 2);
    CHECK(island_spine->type == isg::NonMotorBackboneEdgeType::Composite);
}

// 验证从所有邻接横道中按路口面中心距离选取最近两条，并保留其交点中心对应的 Composite 来源。
TEST_CASE("NonMotorBackboneExtractor uses the two crosswalk spines nearest the intersection center") {
    isg::IntersectionInput input;
    input.area.geometry = rectangle(-0.5, -0.3, 0.5, 0.7);
    isg::SafetyIsland island;
    island.id = "nearest_intersection_island";
    island.geometry.outer = {isg::Vec3d(-2.0, -1.0, 0.0), isg::Vec3d(2.0, -1.0, 0.0),
                             isg::Vec3d(1.0, 1.0, 0.0), isg::Vec3d(-1.0, 1.0, 0.0)};
    input.safety_islands.push_back(island);

    isg::Crosswalk horizontal_crosswalk;
    horizontal_crosswalk.id = "nearest_horizontal_crosswalk";
    horizontal_crosswalk.geometry = rectangle(-3.0, -0.5, 3.0, 0.5);
    input.crosswalks.push_back(horizontal_crosswalk);

    isg::Crosswalk vertical_crosswalk;
    vertical_crosswalk.id = "nearest_vertical_crosswalk";
    vertical_crosswalk.geometry = rectangle(0.0, -3.0, 1.0, 3.0);
    input.crosswalks.push_back(vertical_crosswalk);

    isg::Crosswalk diagonal_crosswalk;
    diagonal_crosswalk.id = "nearest_diagonal_crosswalk";
    diagonal_crosswalk.geometry = diagonalRectangle(0.0, 0.3);
    input.crosswalks.push_back(diagonal_crosswalk);

    const isg::NonMotorBackbone backbone = isg::NonMotorBackboneExtractor().run(input);
    CHECK(countEdgesContainingId(backbone, "si_nearest_intersection_island_spine") == 3);
    CHECK(std::any_of(backbone.edges.begin(), backbone.edges.end(),
        [](const isg::NonMotorBackboneEdge& edge) {
            return edge.id.find("si_nearest_intersection_island_spine") != std::string::npos &&
                   edge.type == isg::NonMotorBackboneEdgeType::Composite;
        }));
}

// 验证过滤后仍有至少三条横道时，近似平行横道对的最近端点连线中点优先作为安全岛中心，公共端点按浅角规则合并。
TEST_CASE("NonMotorBackboneExtractor centers island on a near-parallel pair among three crosswalks") {
    isg::IntersectionInput input;
    input.id = "three_crosswalks_parallel_pair";
    input.area.geometry = rectangle(-0.1, -0.1, 0.1, 0.1);
    isg::SafetyIsland island;
    island.id = "island";
    island.geometry = rectangle(-1.0, -1.0, 1.0, 1.0);
    input.safety_islands.push_back(island);

    isg::Crosswalk first_parallel;
    first_parallel.id = "first_parallel";
    first_parallel.geometry = rectangle(-3.0, 0.1, -1.0, 0.3);
    input.crosswalks.push_back(first_parallel);
    isg::Crosswalk second_parallel;
    second_parallel.id = "second_parallel";
    second_parallel.geometry = rectangle(1.2, 0.3, 3.0, 0.5);
    input.crosswalks.push_back(second_parallel);
    isg::Crosswalk perpendicular;
    perpendicular.id = "perpendicular";
    perpendicular.geometry = rectangle(-0.1, -3.0, 0.1, -1.2);
    input.crosswalks.push_back(perpendicular);

    const auto backbone = isg::NonMotorBackboneExtractor().run(input);
    CHECK(countEdgesContainingId(backbone, "si_island_spine") == 3);
    CHECK(std::any_of(backbone.nodes.begin(), backbone.nodes.end(),
        [](const isg::NonMotorBackboneNode& node) {
            return (node.point.xy() - isg::Vec2d(0.1, 0.3)).norm() < 1e-7;
        }));
}

// 验证最近端点连线中点离安全岛边界不足0.5米时，改用岛中心在连线上的最近安全投影点，并保留中心节点。
TEST_CASE("NonMotorBackboneExtractor projects an unsafe parallel-pair midpoint onto the segment") {
    isg::IntersectionInput input;
    input.id = "parallel_pair_center_projection";
    isg::SafetyIsland island;
    island.id = "island";
    island.geometry = rectangle(-2.0, -2.0, 2.0, 2.0);
    input.safety_islands.push_back(island);

    isg::Crosswalk first_parallel;
    first_parallel.id = "first_parallel";
    first_parallel.geometry = rectangle(-1.5, 0.9, 0.2, 1.1);
    input.crosswalks.push_back(first_parallel);
    isg::Crosswalk second_parallel;
    second_parallel.id = "second_parallel";
    second_parallel.geometry = rectangle(0.2, 2.4, 1.5, 2.6);
    input.crosswalks.push_back(second_parallel);
    isg::Crosswalk perpendicular;
    perpendicular.id = "perpendicular";
    perpendicular.geometry = rectangle(-3.2, -1.0, -3.0, 1.0);
    input.crosswalks.push_back(perpendicular);

    const auto backbone = isg::NonMotorBackboneExtractor().run(input);
    size_t crosswalk_spine_count = 0;
    for (const auto& edge : backbone.edges) {
        if (edge.type == isg::NonMotorBackboneEdgeType::CrosswalkSpine)
            ++crosswalk_spine_count;
    }
    CHECK(crosswalk_spine_count == 3);
    CHECK(countEdgesContainingId(backbone, "si_island_spine") == 2);
    CHECK(std::any_of(backbone.nodes.begin(), backbone.nodes.end(),
        [](const isg::NonMotorBackboneNode& node) {
            return (node.point.xy() - isg::Vec2d(0.2, 1.0)).norm() < 1e-7;
        }));
}

// 验证连线中点落在安全岛外时，将岛中心投影到线段内部并作为安全岛中心，公共端点按浅角规则合并。
TEST_CASE("NonMotorBackboneExtractor projects an outside parallel-pair midpoint into the island") {
    isg::IntersectionInput input;
    input.id = "outside_parallel_pair_midpoint_projection";
    isg::SafetyIsland island;
    island.id = "island";
    island.geometry = rectangle(-1.0, -1.0, 1.0, 1.0);
    input.safety_islands.push_back(island);

    isg::Crosswalk left;
    left.id = "left_parallel";
    left.geometry = rectangle(-5.0, -0.1, -3.5, 0.1);
    input.crosswalks.push_back(left);
    isg::Crosswalk right;
    right.id = "right_parallel";
    right.geometry = rectangle(1.0, -0.1, 2.5, 0.1);
    input.crosswalks.push_back(right);
    isg::Crosswalk upper;
    upper.id = "upper_perpendicular";
    upper.geometry = rectangle(-0.1, -3.5, 0.1, -2.5);
    input.crosswalks.push_back(upper);

    const auto backbone = isg::NonMotorBackboneExtractor().run(input);
    CHECK(countEdgesContainingId(backbone, "si_island_spine") == 3);
    CHECK(std::any_of(backbone.nodes.begin(), backbone.nodes.end(),
        [](const isg::NonMotorBackboneNode& node) {
            return (node.point.xy() - isg::Vec2d(0.0, 0.0)).norm() < 1e-7;
        }));
}

// 验证三条及以上邻接横道的骨干夹角达到10度时不触发平行对中心规则，仍按来源保留 Composite。
TEST_CASE("NonMotorBackboneExtractor keeps the original center logic for a fifteen degree pair") {
    isg::IntersectionInput input;
    input.id = "fifteen_degree_pair_uses_original_center";
    input.area.geometry = rectangle(-0.1, -0.1, 0.1, 0.1);
    isg::SafetyIsland island;
    island.id = "island";
    island.geometry = rectangle(-1.0, -1.0, 1.0, 1.0);
    input.safety_islands.push_back(island);

    isg::Crosswalk first_parallel;
    first_parallel.id = "first_parallel";
    first_parallel.geometry = rectangle(-3.0, 0.1, -1.0, 0.3);
    input.crosswalks.push_back(first_parallel);
    isg::Crosswalk second_parallel;
    second_parallel.id = "second_parallel";
    second_parallel.geometry = orientedRectangle(
        isg::Vec2d(2.1, 0.4),
        isg::Vec2d(1.0, std::tan(15.0 * 3.14159265358979323846 / 180.0)), 0.9);
    input.crosswalks.push_back(second_parallel);
    isg::Crosswalk perpendicular;
    perpendicular.id = "perpendicular";
    perpendicular.geometry = rectangle(-0.1, -3.0, 0.1, -1.2);
    input.crosswalks.push_back(perpendicular);

    const auto backbone = isg::NonMotorBackboneExtractor().run(input);
    CHECK(countEdgesContainingId(backbone, "si_island_spine") == 3);
    CHECK(std::any_of(backbone.nodes.begin(), backbone.nodes.end(),
        [](const isg::NonMotorBackboneNode& node) {
            return (node.point.xy() - isg::Vec2d(0.0, 0.2)).norm() < 1e-7;
        }));
}

// 验证任一横道骨干延长线交点落在安全岛面外时，中心回退到安全岛面质心并合并浅角公共端点。
TEST_CASE("NonMotorBackboneExtractor falls back when a spine intersection is outside island") {
    isg::IntersectionInput input;
    input.area.geometry = rectangle(-1, -1, 1, 1);
    isg::SafetyIsland island;
    island.id = "outside_intersection_island";
    island.geometry.outer = {isg::Vec3d(-2.0, -1.0, 0.0), isg::Vec3d(2.0, -1.0, 0.0),
                             isg::Vec3d(1.0, 1.0, 0.0), isg::Vec3d(-1.0, 1.0, 0.0)};
    input.safety_islands.push_back(island);

    isg::Crosswalk horizontal_crosswalk;
    horizontal_crosswalk.id = "outside_horizontal_crosswalk";
    horizontal_crosswalk.geometry = rectangle(-3.0, -0.5, 3.0, 0.5);
    input.crosswalks.push_back(horizontal_crosswalk);

    isg::Crosswalk vertical_crosswalk;
    vertical_crosswalk.id = "outside_vertical_crosswalk";
    vertical_crosswalk.geometry = rectangle(2.5, -3.0, 3.5, 3.0);
    input.crosswalks.push_back(vertical_crosswalk);

    const isg::NonMotorBackbone backbone = isg::NonMotorBackboneExtractor().run(input);
    const auto island_spine = std::find_if(backbone.edges.begin(), backbone.edges.end(),
        [](const isg::NonMotorBackboneEdge& edge) {
            return edge.id.find("si_outside_intersection_island_spine") != std::string::npos;
        });
    REQUIRE(island_spine != backbone.edges.end());
    REQUIRE(island_spine->geometry.points.size() == 2);
    CHECK(island_spine->type == isg::NonMotorBackboneEdgeType::Composite);
}

// 验证横道和安全岛骨干线均写入带Z的调试CSV；不依赖控制台调试开关。
TEST_CASE("NonMotorBackboneExtractor logs crosswalk and island spine WKT") {
    isg::IntersectionInput input;
    addStraightPair(input, "entry", "exit", 0.0);
    isg::Crosswalk crosswalk;
    crosswalk.id = "debug_cw";
    crosswalk.geometry = rectangle(-2.0, 2.0, 2.0, 4.0);
    crosswalk.crossing_direction = isg::Vec2d(0, 1);
    input.crosswalks.push_back(crosswalk);
    isg::SafetyIsland island;
    island.id = "debug_island";
    island.geometry = rectangle(-1.0, -1.0, 1.0, 1.0);
    input.safety_islands.push_back(island);

    const std::string out_dir = std::string(PROJECT_ROOT_DIR) + "/cmake-build-debug/road_gap_csv_test";
    ghc::filesystem::create_directories(out_dir);
    const auto backbone = isg::NonMotorBackboneExtractor(out_dir).run(input);
    std::ifstream csv(out_dir + "/nonmotor_backbone.csv");
    std::ostringstream debug_output; debug_output << csv.rdbuf();
    REQUIRE_FALSE(backbone.edges.empty());
    CHECK(debug_output.str().find("人行横道骨干线, debug_cw, \"LINESTRING Z (") != std::string::npos);
    CHECK(debug_output.str().find("安全岛骨干线, debug_island, \"LINESTRING Z (") != std::string::npos);
}


// 连续RoadEdge可来自多条反向数字化的边；偏移后的折线和公共节点必须完整入图。
TEST_CASE("NonMotorBackboneExtractor preserves inward road gap geometry and shared endpoints") {
    isg::IntersectionInput input;
    input.id = "road_gap_polyline";
    input.area.geometry = rectangle(-5, -2, 15, 15);
    isg::Crosswalk left, right;
    left.id = "left"; left.geometry = rectangle(-0.2, 0, 0.2, 5);
    right.id = "right"; right.geometry = rectangle(9.8, 0, 10.2, 5);
    // 两侧骨干按环序取相反方向，顶部沿原Gap闭合；底部采用多点道路Gap。
    std::rotate(right.geometry.outer.begin(), right.geometry.outer.begin() + 2, right.geometry.outer.end());
    for (auto* crosswalk : {&left, &right})
        for (auto& p : crosswalk->geometry.outer) p.z() = 4.0;
    input.crosswalks = {left, right};
    isg::Boundary a, b;
    a.type = b.type = isg::Boundary::Type::RoadEdge;
    a.geometry.points = {isg::Vec3d(-1, 0, 4), isg::Vec3d(4, 0, 4), isg::Vec3d(5, -0.5, 4)};
    b.geometry.points = {isg::Vec3d(11, 0, 4), isg::Vec3d(6, 0, 4), isg::Vec3d(5, -0.5, 4)};
    input.boundaries = {a, b};
    SECTION("reversed fragments") {}
    SECTION("endpoint connects to a segment interior") {
        input.boundaries[0].geometry.points.back() = isg::Vec3d(6, -1, 4);
    }
    SECTION("fragment input order does not matter") { std::reverse(input.boundaries.begin(), input.boundaries.end()); }
    const auto backbone = isg::NonMotorBackboneExtractor().run(input);
    const auto gap = std::find_if(backbone.edges.begin(), backbone.edges.end(),
        [](const isg::NonMotorBackboneEdge& edge) { return edge.id.find("gap_road_") == 0; });
    REQUIRE(gap != backbone.edges.end());
    REQUIRE(gap->geometry.points.size() >= 5);
    CHECK(backbone.connected_components == 1);
    CHECK(backbone.diagnostics.empty());
    REQUIRE(backbone.closed);
    double signed_area = 0;
    for (size_t i = 1; i < backbone.ring_geometry.points.size(); ++i)
        signed_area += isg::cross2d(backbone.ring_geometry.points[i - 1].xy(), backbone.ring_geometry.points[i].xy());
    CHECK(signed_area > 0);
    CHECK(gap->geometry.points.front().y() == Catch::Approx(0.9));
    CHECK(gap->geometry.points.back().y() == Catch::Approx(0.9));
    for (const auto& p : gap->geometry.points) CHECK(p.z() == Catch::Approx(4.0));
    // 每个实际Gap直段的中点与来源RoadEdge保持0.9米法向间距，不能只平移首尾。
    for (size_t i = 1; i < gap->geometry.points.size(); ++i) {
        const auto midpoint = 0.5 * (gap->geometry.points[i - 1] + gap->geometry.points[i]);
        double distance = std::numeric_limits<double>::infinity();
        for (const auto& road : input.boundaries)
            for (size_t j = 1; j < road.geometry.points.size(); ++j)
                distance = std::min(distance, isg::pointToSegment(midpoint, road.geometry.points[j - 1], road.geometry.points[j]).first);
        CHECK(distance == Catch::Approx(0.9).margin(1e-6));
    }
    for (const auto& node : backbone.nodes) {
        if (node.id == gap->from) CHECK(isg::dist(node.point, gap->geometry.points.front()) < 1e-7);
        if (node.id == gap->to) CHECK(isg::dist(node.point, gap->geometry.points.back()) < 1e-7);
    }
    // 原骨干的公共端点随Gap交点移动，另一端仍保留，不允许变成跳接车道端点。
    for (const auto& edge : backbone.edges) if (edge.type == isg::NonMotorBackboneEdgeType::CrosswalkSpine) {
        CHECK(std::max(edge.geometry.points.front().y(), edge.geometry.points.back().y()) == Catch::Approx(5.0));
        CHECK(std::min(edge.geometry.points.front().y(), edge.geometry.points.back().y()) == Catch::Approx(0.9));
    }
    if (gap->ring) {
        REQUIRE(backbone.closed);
        CHECK(backbone.ring_geometry.points.size() > backbone.ring_nodes.size() + 1);
        for (const auto& p : gap->geometry.points)
            CHECK(std::any_of(backbone.ring_geometry.points.begin(), backbone.ring_geometry.points.end(),
                [&](const isg::Vec3d& q) { return isg::dist(p, q) < 1e-7; }));
    }
}

// 道路边缘连接可跨过普通20米端点搜索范围；最终仍沿道路有效段内移，并保持公共节点一致。
TEST_CASE("NonMotorBackboneExtractor uses an inward road edge for a long backbone gap") {
    isg::IntersectionInput input;
    input.id = "long_road_gap";
    input.area.geometry = rectangle(-5, -2, 30, 15);
    isg::Crosswalk left, right;
    left.id = "left"; left.geometry = rectangle(-0.2, 0, 0.2, 5);
    right.id = "right"; right.geometry = rectangle(24.8, 0, 25.2, 5);
    input.crosswalks = {left, right};
    isg::Boundary road;
    road.type = isg::Boundary::Type::RoadEdge;
    road.id = "long-road-edge";
    road.geometry.points = {isg::Vec3d(-1, 0, 2), isg::Vec3d(26, 0, 2)};
    input.boundaries = {road};

    const auto backbone = isg::NonMotorBackboneExtractor().run(input);
    const auto gap = std::find_if(backbone.edges.begin(), backbone.edges.end(),
        [](const isg::NonMotorBackboneEdge& edge) { return edge.id.find("gap_road_") == 0; });
    REQUIRE(gap != backbone.edges.end());
    REQUIRE(gap->geometry.points.size() >= 2);
    CHECK(backbone.connected_components == 1);
    CHECK(backbone.diagnostics.empty());
    for (const auto& point : gap->geometry.points)
        CHECK(point.y() == Catch::Approx(0.9).margin(1e-6));
    // 道路内部点沿RoadEdge保留高程；两端为与骨干求交后的公共点，沿骨干插值。
    for (size_t i = 1; i + 1 < gap->geometry.points.size(); ++i)
        CHECK(gap->geometry.points[i].z() == Catch::Approx(2.0).margin(1e-6));
    for (const auto& node : backbone.nodes) {
        if (node.id == gap->from) CHECK(isg::dist(node.point, gap->geometry.points.front()) < 1e-7);
        if (node.id == gap->to) CHECK(isg::dist(node.point, gap->geometry.points.back()) < 1e-7);
    }
}

// 断缝、外扩面不相交、平行无交点、搜索超限及障碍必须回退，不得伪造道路连续。
TEST_CASE("NonMotorBackboneExtractor rejects invalid road gap candidates") {
    isg::IntersectionInput input;
    input.area.geometry = rectangle(-5, -5, 15, 15);
    isg::Crosswalk left, right;
    left.id = "left"; left.geometry = rectangle(-0.2, 0, 0.2, 5);
    right.id = "right"; right.geometry = rectangle(9.8, 0, 10.2, 5);
    input.crosswalks = {left, right};
    isg::Boundary road; road.type = isg::Boundary::Type::RoadEdge;
    road.geometry.points = {isg::Vec3d(-1, 0), isg::Vec3d(11, 0)};
    input.boundaries = {road};
    isg::NonMotorGapOptions options;
    SECTION("disconnected digitization") {
        input.boundaries[0].geometry.points.back() = isg::Vec3d(4.99, 0);
        road.geometry.points.front() = isg::Vec3d(5.01, 0); input.boundaries.push_back(road);
    }
    SECTION("expanded surface misses road") {
        for (auto& p : input.boundaries[0].geometry.points) p.y() = -2;
    }
    SECTION("search radius") {
        options.search_radius = 9;
        // 普通端点搜索和道路边缘搜索是独立上限；本节同时收紧道路上限以验证超限回退。
        options.road_edge_search_radius = 9;
    }
    SECTION("endpoint adjustment is bounded") { options.endpoint_adjustment = 0.1; }
    SECTION("grade separated road fragments") {
        input.boundaries[0].geometry.points.back() = isg::Vec3d(5, 0);
        road.geometry.points.front() = isg::Vec3d(5, 0, 2);
        road.geometry.points.back().z() = 2;
        input.boundaries.push_back(road);
    }
    SECTION("parallel support") {
        input.boundaries[0].geometry.points = {isg::Vec3d(0, 0), isg::Vec3d(0, -1),
            isg::Vec3d(10, -1), isg::Vec3d(10, 0)};
    }
    SECTION("obstacle") {
        isg::Obstacle obstacle; obstacle.geometry = rectangle(4, 0.5, 6, 1.3);
        input.obstacles.push_back(obstacle);
    }
    SECTION("fence hole") { input.area.geometry.holes.push_back(rectangle(4, 0.7, 6, 1.1).outer); }
    SECTION("other boundary") {
        isg::Boundary wall; wall.geometry.points = {isg::Vec3d(5, 0.5), isg::Vec3d(5, 1.5)};
        input.boundaries.push_back(wall);
    }
    const auto backbone = isg::NonMotorBackboneExtractor("./", options).run(input);
    CHECK(std::none_of(backbone.edges.begin(), backbone.edges.end(),
        [](const isg::NonMotorBackboneEdge& edge) { return edge.id.find("gap_road_") == 0; }));
    CHECK_FALSE(backbone.edges.empty());
}

// 原数据留下两个连通分量：补连后应为一个；阻挡/搜索超限时保留明确诊断，不能穿障假连通。
TEST_CASE("NonMotorBackboneExtractor audits disconnected components without breaking the original ring") {
    auto input = isg::IntersectionIO::loadFromFile(std::string(PROJECT_ROOT_DIR) + "/datas/216784982_110004764.json");
    isg::NonMotorGapOptions options;
    bool expect_connected = true;
    SECTION("road priority joins the isolated crosswalk") {}
    SECTION("legacy straight gap fallback") { input.boundaries.clear(); }
    SECTION("boundary prevents any safe repair") {
        isg::Boundary wall;
        wall.geometry.points = {isg::Vec3d(-45, -100), isg::Vec3d(-45, 100)};
        input.boundaries.push_back(wall);
        expect_connected = false;
    }
    SECTION("disabled search does not fake connectivity") { options.search_radius = 0; expect_connected = false; }
    const auto started = std::chrono::steady_clock::now();
    const auto backbone = isg::NonMotorBackboneExtractor(std::string(PROJECT_ROOT_DIR) + "/cmake-build-debug", options).run(input);
    const double elapsed = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - started).count();
    std::cout << "Backbone road-gap generation: " << elapsed << " ms\n";
    CHECK(elapsed < 15000);
    REQUIRE(backbone.closed);
    CHECK(backbone.connected_components == (expect_connected ? 1 : 2));
    CHECK(backbone.diagnostics.empty() == expect_connected);
    // 用实际边图独立遍历核对计数，诊断字段不能掩盖残余孤立骨干。
    std::set<std::string> visited;
    size_t components = 0;
    for (const auto& node : backbone.nodes) {
        if (visited.count(node.id)) continue;
        ++components;
        std::vector<std::string> pending{node.id}; visited.insert(node.id);
        for (size_t cursor = 0; cursor < pending.size(); ++cursor)
            for (const auto& edge : backbone.edges) {
                const auto other = edge.from == pending[cursor] ? edge.to : edge.to == pending[cursor] ? edge.from : "";
                if (!other.empty() && visited.insert(other).second) pending.push_back(other);
            }
    }
    CHECK(components == backbone.connected_components);
}
