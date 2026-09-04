// 性能优化引入的缓存与包围盒剪枝的等价性回归。
//
// 本轮把三处热路径改成「先查缓存/先按包围盒拒绝，再做原判定」：
//   1. IntersectionInput::findLane 的 id→Lane* 两级缓存；
//   2. cachedBoundarySafetySegments 按内容指纹复用边界安全段；
//   3. curveBoundarySafety / curveRawIntersectsBoundariesImpl 的分段包围盒剪枝。
// 这些改动都以「跳过的工作结论必然相同」为前提，所以回归的重点是：缓存不会
// 跨输入串味、原地改写输入后不会返回过期结果、贴近包围盒边界的真实接触不会
// 被剪枝漏判。
#include <catch2/catch_test_macros.hpp>

#include <string>
#include <vector>

#include "constraints/boundary_safety.h"
#include "curve/bezier.h"
#include "types.h"
#include "utils.h"

using namespace isg;

namespace {

Lane makeLane(const std::string& id, const std::vector<Vec2d>& pts) {
    Lane lane;
    lane.id = id;
    lane.geometry.points = toVec3dArray(pts);
    return lane;
}

Boundary makeBoundary(const std::string& id, Boundary::Type type,
                      const std::vector<Vec2d>& pts) {
    Boundary bnd;
    bnd.id = id;
    bnd.type = type;
    bnd.geometry.points = toVec3dArray(pts);
    return bnd;
}

bool samePoint(const Vec2d& a, const Vec2d& b) {
    return a.x() == b.x() && a.y() == b.y();
}

BezierCurve straightCurve(const Vec2d& from, const Vec2d& to) {
    BezierSegment seg;
    seg.ctrl[0] = from;
    seg.ctrl[1] = from + (to - from) / 3.0;
    seg.ctrl[2] = from + 2.0 * (to - from) / 3.0;
    seg.ctrl[3] = to;
    BezierCurve curve;
    curve.segs.push_back(seg);
    return curve;
}

}  // namespace

TEST_CASE("findLane cache keeps per-input results separate", "[perf_cache]") {
    IntersectionInput a;
    a.lanes.push_back(makeLane("L1", {Vec2d(0, 0), Vec2d(10, 0)}));
    a.lanes.push_back(makeLane("L2", {Vec2d(0, 5), Vec2d(10, 5)}));

    IntersectionInput b;
    b.lanes.push_back(makeLane("L1", {Vec2d(0, 100), Vec2d(10, 100)}));

    // 两个输入交替查询同一个 id：缓存必须按输入对象隔离，不能互相串味。
    for (int i = 0; i < 4; ++i) {
        const Lane* la = a.findLane("L1");
        const Lane* lb = b.findLane("L1");
        REQUIRE(la != nullptr);
        REQUIRE(lb != nullptr);
        CHECK(la->geometry.points.front().y() == 0.0);
        CHECK(lb->geometry.points.front().y() == 100.0);
    }
    CHECK(a.findLane("L2") != nullptr);
    CHECK(a.findLane("missing") == nullptr);
    CHECK(b.findLane("L2") == nullptr);
}

TEST_CASE("findLane cache survives in-place lane rewrite", "[perf_cache]") {
    IntersectionInput input;
    input.lanes.push_back(makeLane("A", {Vec2d(0, 0), Vec2d(10, 0)}));
    input.lanes.push_back(makeLane("B", {Vec2d(0, 5), Vec2d(10, 5)}));
    REQUIRE(input.findLane("A") != nullptr);
    REQUIRE(input.findLane("B") != nullptr);

    // 元素数量不变、缓冲区地址不变，只把 id 换位：索引若不校验命中项的 id
    // 就会返回错误 Lane。
    input.lanes[0].id = "B";
    input.lanes[1].id = "A";
    const Lane* found_a = input.findLane("A");
    const Lane* found_b = input.findLane("B");
    REQUIRE(found_a != nullptr);
    REQUIRE(found_b != nullptr);
    CHECK(found_a->geometry.points.front().y() == 5.0);
    CHECK(found_b->geometry.points.front().y() == 0.0);

    // 几何原地改写必须立刻可见（缓存只存指针，不存坐标副本）。
    input.lanes[1].geometry.points.front() = Vec3d(Vec2d(0, -7));
    CHECK(input.findLane("A")->geometry.points.front().y() == -7.0);
}

TEST_CASE("findLane address cache validates string content", "[perf_cache]") {
    IntersectionInput input;
    input.lanes.push_back(makeLane("A", {Vec2d(0, 0), Vec2d(10, 0)}));
    input.lanes.push_back(makeLane("B", {Vec2d(0, 5), Vec2d(10, 5)}));

    // 循环里的局部字符串通常复用同一栈地址而内容不同，直接命中一级地址缓存
    // 的场景。校验必须落到字符串内容上。
    for (int i = 0; i < 8; ++i) {
        const std::string id = (i % 2 == 0) ? "A" : "B";
        const Lane* lane = input.findLane(id);
        REQUIRE(lane != nullptr);
        CHECK(lane->id == id);
    }
    for (int i = 0; i < 4; ++i) {
        const std::string id = "nope";
        CHECK(input.findLane(id) == nullptr);
    }
}

TEST_CASE("entry/exit endpoints follow in-place geometry edits", "[perf_cache]") {
    IntersectionInput input;
    input.lanes.push_back(makeLane("E", {Vec2d(0, 0), Vec2d(10, 0)}));
    const std::pair<Vec2d, Vec2d> before = input.entryPtDir("E");
    CHECK(before.first.x() == 10.0);

    input.lanes[0].geometry.points.back() = Vec3d(Vec2d(12, 0));
    const std::pair<Vec2d, Vec2d> after = input.entryPtDir("E");
    CHECK(after.first.x() == 12.0);
    CHECK(input.exitPtDir("E").first.x() == 0.0);
}

TEST_CASE("boundary safety segment cache matches direct rebuild",
          "[perf_cache]") {
    const Vec2d center(5, 0);
    std::vector<Boundary> set_a;
    set_a.push_back(makeBoundary("b1", Boundary::Type::RoadEdge,
                                 {Vec2d(0, 4), Vec2d(10, 4)}));
    set_a.push_back(makeBoundary("b2", Boundary::Type::RoadEdge,
                                 {Vec2d(0, -4), Vec2d(10, -4)}));
    std::vector<Boundary> set_b;
    set_b.push_back(makeBoundary("b3", Boundary::Type::RoadEdge,
                                 {Vec2d(0, 9), Vec2d(10, 9)}));

    const std::vector<BoundarySafetySegment> ref_a =
        buildBoundarySafetySegments(set_a, center, nullptr, nullptr);
    const std::vector<BoundarySafetySegment> ref_b =
        buildBoundarySafetySegments(set_b, center, nullptr, nullptr);

    // 交替查询两套边界，缓存不得互相覆盖，且必须与直接重建逐项一致。
    for (int i = 0; i < 3; ++i) {
        const std::vector<BoundarySafetySegment>& got_a =
            cachedBoundarySafetySegments(set_a, center, nullptr);
        const std::vector<BoundarySafetySegment>& got_b =
            cachedBoundarySafetySegments(set_b, center, nullptr);
        REQUIRE(got_a.size() == ref_a.size());
        REQUIRE(got_b.size() == ref_b.size());
        for (std::size_t k = 0; k < ref_a.size(); ++k) {
            CHECK(samePoint(got_a[k].a, ref_a[k].a));
            CHECK(samePoint(got_a[k].b, ref_a[k].b));
            CHECK(got_a[k].center_signed == ref_a[k].center_signed);
            CHECK(got_a[k].road_edge == ref_a[k].road_edge);
        }
        for (std::size_t k = 0; k < ref_b.size(); ++k) {
            CHECK(samePoint(got_b[k].a, ref_b[k].a));
            CHECK(got_b[k].center_signed == ref_b[k].center_signed);
        }
    }

    // 坐标改写会改变内容指纹，缓存必须重建而不是复用旧段。
    set_a[0].geometry.points[0] = Vec3d(Vec2d(0, 6));
    const std::vector<BoundarySafetySegment> ref_moved =
        buildBoundarySafetySegments(set_a, center, nullptr, nullptr);
    const std::vector<BoundarySafetySegment>& got_moved =
        cachedBoundarySafetySegments(set_a, center, nullptr);
    REQUIRE(got_moved.size() == ref_moved.size());
    for (std::size_t k = 0; k < ref_moved.size(); ++k)
        CHECK(samePoint(got_moved[k].a, ref_moved[k].a));
}

TEST_CASE("boundary contact survives bounding-box pruning", "[perf_cache]") {
    const Vec2d center(5, 0);
    // 曲线自下向上穿过横向 RoadEdge：剪枝后仍必须判为相交。
    std::vector<Boundary> crossing;
    crossing.push_back(makeBoundary("cross", Boundary::Type::RoadEdge,
                                    {Vec2d(-20, 1), Vec2d(20, 1)}));
    const BezierCurve up = straightCurve(Vec2d(5, -6), Vec2d(5, 6));
    const BoundarySafetyResult hit =
        curveBoundarySafety(up, crossing, center, 128, 0.15);
    CHECK(hit.intersects);

    // 同一条曲线整体位于边界一侧且远离：不得判为相交。
    std::vector<Boundary> apart;
    apart.push_back(makeBoundary("apart", Boundary::Type::RoadEdge,
                                 {Vec2d(-20, 30), Vec2d(20, 30)}));
    const BoundarySafetyResult miss =
        curveBoundarySafety(up, apart, center, 128, 0.15);
    CHECK_FALSE(miss.intersects);

    // 贴着采样段包围盒边缘的接触：包围盒剪枝一旦用了严格小于就会漏判。
    std::vector<Boundary> grazing;
    grazing.push_back(makeBoundary("graze", Boundary::Type::RoadEdge,
                                   {Vec2d(5, -3), Vec2d(9, -3)}));
    const BezierCurve across = straightCurve(Vec2d(7, -8), Vec2d(7, 8));
    const BoundarySafetyResult graze =
        curveBoundarySafety(across, grazing, center, 128, 0.15);
    CHECK(graze.intersects);
}
