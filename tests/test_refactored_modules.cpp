#include <catch2/catch_test_macros.hpp>

#include "constraints/constraint_evaluator.h"
#include "constraints/fence_check.h"
#include "constraints/uturn_envelope_constraint.h"
#include "curve/bezier.h"
#include "curve/curve_utils.h"
#include "utils.h"
#include "domain/scene_view.h"
#include "io/iodata_json.h"
#include "initialization/fixed_shape_initializer.h"
#include "initialization/ordinary_curve_initializer.h"
#include "initialization/uturn_curve_initializer.h"
#include "initialization/avoidance_candidate_generator.h"
#include "generation/candidate_selector.h"
#include "generation/bounded_uturn_candidate_search.h"
#include "generation/curve_sampling.h"
#include "generation/segmented_uturn_candidate_search.h"
#include "generation/uturn_multi_constraint_search.h"
#include "generation/uturn_shape.h"
#include "generation/curve_post_processor.h"
#include "initialization/curve_initializer_registry.h"
#include "preprocessing/uturn_family_builder.h"
#include "preprocessing/crosswalk_clearance_calculator.h"
#include "toolkits/toolkits.h"
#include "ordering/generation_planner.h"
#include "optimizer/sdf_field.h"
#include "optimizer/lbfgs_solver.h"
#include "generation/curve_optimizer.h"
#include "generation/optimization_result_processor.h"
#include "generation/physical_repair_coordinator.h"
#include "generation/final_curve_auditor.h"
#include "generation/repair_budget.h"
#include "generation/repair_impact_closure.h"
#include "generation/atomic_curve_patch.h"
#include "generation/final_pair_auditor.h"

namespace {

isg::Lane lane(const std::string& id, const isg::Vec3d& a, const isg::Vec3d& b) {
    isg::Lane result;
    result.id = id;
    result.geometry.points.push_back(a);
    result.geometry.points.push_back(b);
    return result;
}

isg::BezierCurve cubic(const isg::Vec2d& p0, const isg::Vec2d& p1,
                       const isg::Vec2d& p2, const isg::Vec2d& p3) {
    isg::BezierCurve curve;
    isg::BezierSegment segment;
    segment.ctrl = {{p0, p1, p2, p3}};
    curve.segs.push_back(segment);
    return curve;
}

}  // namespace

TEST_CASE("InputNormalizer fills group references without mutating source") {
    isg::IntersectionInput source;
    source.lanes.push_back(lane("entry", isg::Vec3d(0, -10), isg::Vec3d(0, 0)));
    source.lanes.push_back(lane("exit", isg::Vec3d(0, 10), isg::Vec3d(0, 20)));
    isg::LaneGroup entry;
    entry.id = "g-entry";
    entry.role = isg::GroupRole::Entry;
    entry.lanes.push_back("entry");
    isg::LaneGroup exit;
    exit.id = "g-exit";
    exit.role = isg::GroupRole::Exit;
    exit.lanes.push_back("exit");
    source.lane_groups.push_back(entry);
    source.lane_groups.push_back(exit);
    isg::Connectivity connection;
    connection.id = "c";
    connection.entry_lane_id = "entry";
    connection.exit_lane_id = "exit";
    source.connectivities.push_back(connection);

    const isg::IntersectionInput normalized = isg::InputNormalizer(source);
    REQUIRE(source.connectivities.front().enterGroupId.empty());
    REQUIRE(normalized.connectivities.front().enterGroupId == "g-entry");
    REQUIRE(normalized.connectivities.front().exitGroupId == "g-exit");
}

TEST_CASE("ConnectivityDirectionNormalizer reuses the indexed group direction") {
    isg::IntersectionInput input;
    input.lanes.push_back(lane("E1", isg::Vec3d(-10, 1), isg::Vec3d(0, 1)));
    // E2 与组方向只差 5.7°：同一进出口内的车道抖动，必须被统一到 +x。
    input.lanes.push_back(lane("E2", isg::Vec3d(-10, -1), isg::Vec3d(0, 0)));
    // E3 偏 90°：塞进同一 groupId 的另一条臂，跨臂保护必须原样保留它的端点切向。
    input.lanes.push_back(lane("E3", isg::Vec3d(0, -10), isg::Vec3d(0, -2)));
    input.lanes.push_back(lane("X", isg::Vec3d(10, 1), isg::Vec3d(20, 1)));
    isg::LaneGroup group;
    group.id = "entry";
    group.role = isg::GroupRole::Entry;
    group.lanes = {"E1", "E2", "E3"};
    input.lane_groups.push_back(group);
    isg::Connectivity straight;
    straight.id = "straight";
    straight.entry_lane_id = "E1";
    straight.exit_lane_id = "X";
    straight.turn_type = isg::ConnTurnType::Straight;
    input.connectivities.push_back(straight);

    const isg::Vec3d endpoint = input.lanes[1].geometry.points.back();
    const isg::Vec2d cross_arm_dir =
        input.entryPtDir("E3").second.normalized();
    isg::ConnectivityDirectionNormalizer(
        input, isg::ConnectivityDirectionConfig());

    REQUIRE(input.lanes[1].geometry.points.back().xy().isApprox(endpoint.xy()));
    REQUIRE(input.lanes[1].geometry.points.back().z() == endpoint.z());
    REQUIRE(input.entryPtDir("E2").second.normalized().dot(isg::Vec2d(1, 0)) > 0.99);
    REQUIRE(input.entryPtDir("E3").second.normalized().dot(cross_arm_dir) > 0.99);
}

TEST_CASE("FixedShapeInitializer preserves non-degenerate polyline segments") {
    isg::Connectivity connection;
    connection.geometry.points = {
        isg::Vec3d(0, 0, 4), isg::Vec3d(3, 0, 5),
        isg::Vec3d(3, 0, 6), isg::Vec3d(3, 4, 7)};

    const isg::FixedShapeInitializer initializer;
    REQUIRE(initializer.hasGeometry(connection));
    const isg::BezierCurve curve = initializer.build(connection);
    REQUIRE(curve.numSegments() == 2);
    REQUIRE(curve.segs.front().ctrl[0].isApprox(isg::Vec2d(0, 0)));
    REQUIRE(curve.segs.front().ctrl[3].isApprox(isg::Vec2d(3, 0)));
    REQUIRE(curve.segs.back().ctrl[0].isApprox(isg::Vec2d(3, 0)));
    REQUIRE(curve.segs.back().ctrl[3].isApprox(isg::Vec2d(3, 4)));
}

TEST_CASE("OrdinaryCurveInitializer keeps alpha order and tangent fallbacks") {
    const isg::OrdinaryCurveInitializer initializer;
    const std::vector<double> alphas = {0.50, 0.30, 0.12};
    const std::vector<isg::BezierCurve> candidates = initializer.buildAlphaCandidates(
        isg::Vec2d(0, 0), isg::Vec2d(0, 0),
        isg::Vec2d(10, 0), isg::Vec2d(0, 0), alphas);
    REQUIRE(candidates.size() == alphas.size());
    for (size_t i = 0; i < candidates.size(); ++i) {
        REQUIRE(candidates[i].numSegments() == 1);
        REQUIRE(candidates[i].segs.front().ctrl[0].isApprox(isg::Vec2d(0, 0)));
        REQUIRE(candidates[i].segs.front().ctrl[3].isApprox(isg::Vec2d(10, 0)));
    }
    REQUIRE(candidates[0].segs.front().ctrl[1].x() >
            candidates[1].segs.front().ctrl[1].x());
    REQUIRE(candidates[1].segs.front().ctrl[1].x() >
            candidates[2].segs.front().ctrl[1].x());
}

TEST_CASE("OrdinaryCurveInitializer uses quadratic elevation for ordinary turns") {
    const isg::OrdinaryCurveInitializer initializer;
    const isg::BezierCurve turn = initializer.buildPreferredSingleCubic(
        isg::Vec2d(0, 0), isg::Vec2d(1, 0),
        isg::Vec2d(10, 10), isg::Vec2d(0, 1));
    REQUIRE(turn.numSegments() == 1);
    CHECK(turn.segs.front().ctrl[1].isApprox(
        isg::Vec2d(20.0 / 3.0, 0.0), 1e-9));
    CHECK(turn.segs.front().ctrl[2].isApprox(
        isg::Vec2d(10.0, 10.0 / 3.0), 1e-9));

    const isg::BezierCurve straight = initializer.buildPreferredSingleCubic(
        isg::Vec2d(0, 0), isg::Vec2d(1, 0),
        isg::Vec2d(12, 0), isg::Vec2d(1, 0));
    REQUIRE(straight.numSegments() == 1);
    CHECK(straight.segs.front().ctrl[1].isApprox(isg::Vec2d(6, 0), 1e-9));
    CHECK(straight.segs.front().ctrl[2].isApprox(isg::Vec2d(6, 0), 1e-9));
    CHECK(straight.segs.front().ctrl[1].isApprox(
        straight.segs.front().ctrl[2], 1e-12));
    CHECK(std::abs(straight.arcLength() - 12.0) < 1e-8);

    const isg::BezierCurve skewed_straight =
        initializer.buildPreferredSingleCubic(
            isg::Vec2d(0, 0), isg::Vec2d(1, 0.1),
            isg::Vec2d(12, 0), isg::Vec2d(1, -0.1));
    REQUIRE(skewed_straight.numSegments() == 1);
    CHECK((skewed_straight.segs.front().ctrl[1] -
           skewed_straight.segs.front().ctrl[0]).norm() ==
          (skewed_straight.segs.front().ctrl[3] -
           skewed_straight.segs.front().ctrl[2]).norm());
    CHECK_FALSE(skewed_straight.segs.front().ctrl[1].isApprox(
        skewed_straight.segs.front().ctrl[2], 1e-9));
    CHECK(isg::ordinarySingleCubicControlsValid(
        skewed_straight, isg::Vec2d(0, 0), isg::Vec2d(1, 0.1),
        isg::Vec2d(12, 0), isg::Vec2d(1, -0.1), 1e-6, true));
}

TEST_CASE("UTurnCurveInitializer builds aligned and segmented shapes") {
    const isg::UTurnCurveInitializer initializer;
    const isg::Vec2d p0(0, 0);
    const isg::Vec2d p1(4, 0);
    const isg::Vec2d t0(0, 1);
    const isg::Vec2d t1(0, -1);

    const isg::BezierCurve aligned = initializer.buildAligned(
        p0, t0, p1, t1, isg::Vec2d(1, 0), 0.25, 1.0, 2.0, 2.0);
    REQUIRE(aligned.numSegments() == 1);
    REQUIRE(aligned.segs.front().ctrl[0].isApprox(p0));
    REQUIRE(aligned.segs.front().ctrl[3].isApprox(p1));

    const isg::BezierCurve segmented = initializer.buildSegmented(
        p0, t0, p1, t1, 2.0, 2.0);
    REQUIRE(segmented.numSegments() == 3);
    REQUIRE(segmented.segs.front().ctrl[0].isApprox(p0));
    REQUIRE(segmented.segs.back().ctrl[3].isApprox(p1));
    REQUIRE(segmented.segs.front().ctrl[3].y() > p0.y());
    REQUIRE(segmented.segs.back().ctrl[0].y() > p1.y());
}

TEST_CASE("UTurnCurveInitializer keeps a narrow U-turn shape usable after minimum leads") {
    // 0.24m 横向间隙的极窄调头：首尾必须先拉满 2m 直段，中弧只能在这条
    // 0.24m 走廊里完成反向。此时"圆整"不是硬门禁——候选搜索的圆整判定
    // 只在 chord_len >= 1.0 时启用(segmented_uturn_candidate_search.cpp)，
    // 因为 segmentedUTurnMiddleArcLooksRound 的凸出量下限含 0.25m 绝对项，
    // 而中弧凸出量是 0.75 * arc_alpha * gap，gap=0.24 时默认把手比例
    // 2/3 只能凸出 0.12m，先验地不可能达标。窄间隙下真正的硬门禁是
    // 段数、最小直段、自交与曲率上限。
    const isg::Vec2d p0(0, 0);
    const isg::Vec2d t0(0, 1);
    const isg::Vec2d p1(0.24, 0);
    const isg::Vec2d t1(0, -1);
    const isg::UTurnCurveInitializer initializer;
    const isg::BezierCurve curve = initializer.buildSegmented(p0, t0, p1, t1,
                                                              2.0, 2.0);

    REQUIRE(curve.numSegments() == 3);
    CHECK(isg::segmentedUTurnHasMinimumStraightLeads(curve, 2.0, 2.0));
    CHECK((curve.segs.front().ctrl[3] - curve.segs.front().ctrl[0]).norm() >= 2.0);
    CHECK((curve.segs.back().ctrl[3] - curve.segs.back().ctrl[0]).norm() >= 2.0);
    CHECK_FALSE(isg::curveSelfIntersectsBusiness(curve, 1.0));
    CHECK(curve.maxCurvature(40) < isg::segmentedUTurnMaxCurvatureLimit(
        p0, t0, p1, t1));

    // 圆整仍然是可达的，只是需要调用方把把手比例调大——候选搜索的
    // arc_alphas 一直枚举到 2.0，正是为窄间隙准备的档位。
    const isg::BezierCurve round_curve = initializer.buildSegmented(
        p0, t0, p1, t1, 2.0, 2.0, 2.0);
    REQUIRE(round_curve.numSegments() == 3);
    CHECK(isg::segmentedUTurnMiddleArcLooksRound(round_curve, isg::Vec2d(0, 1)));
    CHECK(isg::segmentedUTurnHasMinimumStraightLeads(round_curve, 2.0, 2.0));
}

TEST_CASE("UTurnCurveInitializer keeps minimum straight leads under lateral stagger") {
    // 首尾平齐点的横向错开沿 U 轴法线平移 q0/q1。切向与 U 轴不平行时，这个
    // 平移会在 -T 方向留下分量，把"最小直行段"的弦长读数压到 min_lead 之下，
    // 而轴向站位不变——于是候选搜索的最小直段门禁与家族平齐站位门禁互斥，
    // 三段式候选集被夹空，几何 U 型调头退化成跨不过人行横道的单段曲线。
    // buildSegmented 必须按错开量补偿 lead，使弦长读数重新达到 min_lead。
    const isg::UTurnCurveInitializer initializer;
    const isg::Vec2d p0(0.0, 0.0);
    const isg::Vec2d t0 = isg::Vec2d(1.0, -1.0).normalized();
    const isg::Vec2d p1(2.3, 2.3);
    // 出口切向与进入切向不严格反平行（真实路口的常态）：此时 U 轴是两者的
    // 角平分线，横向错开不再垂直于 T0，才会出现毫米级的弦长亏损。
    const isg::Vec2d t1 = isg::Vec2d(-1.05, 1.0).normalized();
    const double min_lead = 7.93;

    for (double stagger : {0.02, 0.10, 0.30}) {
        const isg::BezierCurve curve = initializer.buildSegmented(
            p0, t0, p1, t1, min_lead, min_lead, 2.0 / 3.0, stagger,
            0.0, 0.0, stagger, stagger);
        REQUIRE(curve.numSegments() == 3);
        // 严格容差下的最小直段门禁：与候选搜索使用的判定完全一致。
        CHECK(isg::segmentedUTurnHasMinimumStraightLeads(
            curve, min_lead, min_lead));
        // 补偿只沿切向加长，轴向站位仍须保持 q0/q1 平齐。
        isg::Vec2d axis = t0 - t1;
        axis.normalize();
        const isg::Vec2d q0 = curve.segs.front().ctrl[3];
        const isg::Vec2d q1 = curve.segs.back().ctrl[0];
        CHECK(std::abs(q0.dot(axis) - q1.dot(axis)) < 0.05);
        CHECK_FALSE(isg::curveSelfIntersectsBusiness(curve, 1.0));
    }
}

TEST_CASE("AvoidanceCandidateGenerator builds waypoint curves with endpoint G1") {
    const isg::AvoidanceCandidateGenerator generator;
    const std::vector<isg::Vec2d> points = {
        isg::Vec2d(0, 0), isg::Vec2d(3, 2),
        isg::Vec2d(7, 2), isg::Vec2d(10, 0)};
    const isg::BezierCurve curve = generator.buildWaypointCurve(
        points, isg::Vec2d(1, 0), isg::Vec2d(1, 0));

    REQUIRE(curve.numSegments() == 3);
    REQUIRE(curve.segs.front().ctrl[0].isApprox(points.front()));
    REQUIRE(curve.segs.back().ctrl[3].isApprox(points.back()));
    REQUIRE(curve.startTan().normalized().isApprox(isg::Vec2d(1, 0)));
    REQUIRE(curve.endTan().normalized().isApprox(isg::Vec2d(1, 0)));
    REQUIRE(generator.buildWaypointCurve(
        std::vector<isg::Vec2d>(1, isg::Vec2d(0, 0)),
        isg::Vec2d(1, 0), isg::Vec2d(1, 0)).empty());
}

TEST_CASE("CandidateSelector rejects stale and hard-violating reports stably") {
    std::vector<isg::CurveCandidate> candidates(4);
    candidates[0].quality_score = 1.0;

    candidates[1].quality_score = 3.0;
    candidates[1].acceptReport(isg::ConstraintReport());

    candidates[2].quality_score = 0.5;
    isg::ConstraintReport hard_report;
    isg::ConstraintResult violation;
    violation.severity = isg::ConstraintSeverity::Hard;
    violation.state = isg::ConstraintState::Violated;
    hard_report.results.push_back(violation);
    candidates[2].acceptReport(hard_report);

    candidates[3].quality_score = 3.0;
    candidates[3].acceptReport(isg::ConstraintReport());

    const isg::CandidateSelector selector;
    REQUIRE(selector.selectBest(candidates) == 1);
    candidates[1].replaceCurve(isg::BezierCurve());
    REQUIRE(selector.selectBest(candidates) == 3);
    candidates[3].replaceCurve(isg::BezierCurve());
    REQUIRE(selector.selectBest(candidates) == isg::CandidateSelector::npos);
}

TEST_CASE("CandidateSelector follows structured business lexicographic order") {
    std::vector<isg::CurveCandidate> candidates(4);
    for (auto& candidate : candidates)
        candidate.acceptReport(isg::ConstraintReport());

    candidates[0].physical_violation = 0.1;
    candidates[0].quality_score = 0.0;

    candidates[1].new_cluster_crosses = 1;
    candidates[1].quality_score = 0.0;

    candidates[2].fixed_shape_loss = 1;
    candidates[2].quality_score = 1.0;

    candidates[3].fixed_shape_loss = 1;
    candidates[3].quality_score = 1.0;

    const isg::CandidateSelector selector;
    REQUIRE(selector.selectBest(candidates) == 2);
    candidates[2].shape_hard_violations = 1;
    REQUIRE(selector.selectBest(candidates) == 3);
}

TEST_CASE("CurveInitializerRegistry applies U-turn fixed and ordinary priority") {
    isg::Connectivity connection;
    connection.id = "candidate";
    isg::CurveGenerationContext context;
    context.connectivity = &connection;
    context.entry = std::make_pair(isg::Vec2d(0, 0), isg::Vec2d(0, 1));
    context.exit = std::make_pair(isg::Vec2d(4, 0), isg::Vec2d(0, -1));
    isg::CurveInitializationOptions options;
    options.uturn_min_lead0 = 2.0;
    options.uturn_min_lead1 = 2.0;

    const isg::CurveInitializerRegistry registry;
    std::vector<isg::CurveCandidate> candidates = registry.build(context, options);
    REQUIRE(candidates.size() == 1);
    REQUIRE(candidates.front().origin == isg::CandidateOrigin::SegmentedUTurn);
    REQUIRE(candidates.front().curve.numSegments() == 3);

    connection.fixed_shape = true;
    connection.geometry.points = {isg::Vec3d(0, 0), isg::Vec3d(4, 4)};
    context.exit.second = isg::Vec2d(1, 0);
    candidates = registry.build(context, options);
    REQUIRE(candidates.size() == 1);
    REQUIRE(candidates.front().origin == isg::CandidateOrigin::FixedShape);
    REQUIRE(candidates.front().preserves_fixed_shape);

    options.allow_fixed_shape = false;
    options.ordinary_alphas = {0.4, 0.3};
    candidates = registry.build(context, options);
    REQUIRE(candidates.size() == 3);
    REQUIRE(candidates[0].origin == isg::CandidateOrigin::NaturalSingleCubic);
    REQUIRE(candidates[1].quality_score > candidates[0].quality_score);
    REQUIRE_FALSE(candidates[0].hasCurrentReport());
}

TEST_CASE("UTurnFamilyBuilder computes shared axis, radius and aligned station") {
    isg::IntersectionInput input;
    input.lanes.push_back(lane("entry", isg::Vec3d(0, -10), isg::Vec3d(0, 0)));
    input.lanes.push_back(lane("exit", isg::Vec3d(4, 0), isg::Vec3d(4, -10)));
    isg::Connectivity connection;
    connection.id = "uturn-family";
    connection.entry_lane_id = "entry";
    connection.exit_lane_id = "exit";
    input.connectivities.push_back(connection);

    const isg::UTurnFamilyBuilder builder;
    REQUIRE(builder.isGeometricUTurn(connection, input));
    REQUIRE(builder.alignmentAxis(isg::Vec2d(0, 1), isg::Vec2d(0, -1))
                .isApprox(isg::Vec2d(0, 1)));
    REQUIRE(std::abs(builder.radiusKey(connection, input) - 4.0) < 1e-9);
    REQUIRE(std::abs(builder.requiredAlignedStation(
                isg::Vec2d(0, 0), isg::Vec2d(0, 1),
                isg::Vec2d(4, 0), isg::Vec2d(0, -1), 2.0, 3.0) - 3.0) < 1e-9);
}

TEST_CASE("UTurnFamilyBuilder closes transitive lane-group families and ranks radius") {
    isg::IntersectionInput input;
    input.lanes.push_back(lane("entry-a", isg::Vec3d(0, -10), isg::Vec3d(0, 0)));
    input.lanes.push_back(lane("entry-b", isg::Vec3d(2, -10), isg::Vec3d(2, 0)));
    input.lanes.push_back(lane("exit-a", isg::Vec3d(4, 0), isg::Vec3d(4, -10)));
    input.lanes.push_back(lane("exit-b", isg::Vec3d(8, 0), isg::Vec3d(8, -10)));

    isg::LaneGroup entry_group;
    entry_group.id = "entry-family";
    entry_group.role = isg::GroupRole::Entry;
    entry_group.lanes = {"entry-a", "entry-b"};
    isg::LaneGroup exit_group;
    exit_group.id = "exit-family";
    exit_group.role = isg::GroupRole::Exit;
    exit_group.lanes = {"exit-a", "exit-b"};
    input.lane_groups = {entry_group, exit_group};

    isg::Connectivity a;
    a.id = "a";
    a.entry_lane_id = "entry-a";
    a.exit_lane_id = "exit-a";
    isg::Connectivity b = a;
    b.id = "b";
    b.exit_lane_id = "exit-b";
    isg::Connectivity c = b;
    c.id = "c";
    c.entry_lane_id = "entry-b";
    input.connectivities = {a, b, c};

    const isg::UTurnFamilyBuilder builder;
    const std::vector<const isg::Connectivity*> component =
        builder.alignmentComponent(
            input.connectivities.front(), input,
            isg::UTurnAlignmentScope::LaneGroup);
    REQUIRE(component.size() == 3);
    const isg::UTurnFamilyRank a_rank = builder.radiusRank(
        input.connectivities[0], input, isg::UTurnAlignmentScope::LaneGroup);
    const isg::UTurnFamilyRank b_rank = builder.radiusRank(
        input.connectivities[1], input, isg::UTurnAlignmentScope::LaneGroup);
    REQUIRE(a_rank.family_size == 3);
    REQUIRE(a_rank.reverse_radius_rank > b_rank.reverse_radius_rank);
}

TEST_CASE("CrosswalkClearanceCalculator selects nearest crossing and far projection") {
    isg::IntersectionInput input;
    isg::Crosswalk near_crosswalk;
    near_crosswalk.id = "near";
    near_crosswalk.geometry.outer = {
        isg::Vec3d(-1, 3), isg::Vec3d(1, 3),
        isg::Vec3d(1, 5), isg::Vec3d(-1, 5)};
    isg::Crosswalk far_crosswalk;
    far_crosswalk.id = "far";
    far_crosswalk.geometry.outer = {
        isg::Vec3d(-1, 8), isg::Vec3d(1, 8),
        isg::Vec3d(1, 10), isg::Vec3d(-1, 10)};
    input.crosswalks = {far_crosswalk, near_crosswalk};

    const isg::CrosswalkClearanceCalculator calculator;
    const isg::CrosswalkClearanceResult ahead = calculator.ahead(
        isg::Vec2d(0, 0), isg::Vec2d(0, 1), input);
    REQUIRE(ahead.found);
    REQUIRE(ahead.crosswalk_id == "near");
    REQUIRE(std::abs(ahead.near - 3.0) < 1e-9);
    REQUIRE(std::abs(ahead.far - 5.0) < 1e-9);
    REQUIRE(std::abs(ahead.clearance - 5.3) < 1e-9);

    const isg::CrosswalkClearanceResult behind = calculator.behind(
        isg::Vec2d(0, 10), isg::Vec2d(0, 1), input);
    REQUIRE(behind.found);
    REQUIRE(behind.crosswalk_id == "far");
}

TEST_CASE("UTurnFamilyBuilder derives one-sided and fallback lead floors") {
    isg::IntersectionInput input;
    isg::Crosswalk crosswalk;
    crosswalk.id = "entry-cw";
    crosswalk.geometry.outer = {
        isg::Vec3d(-1, 3), isg::Vec3d(1, 3),
        isg::Vec3d(1, 5), isg::Vec3d(-1, 5)};
    input.crosswalks.push_back(crosswalk);

    const isg::UTurnFamilyBuilder builder;
    const isg::UTurnLeadFloors one_sided = builder.leadFloors(
        isg::Vec2d(0, 0), isg::Vec2d(0, 1),
        isg::Vec2d(10, 0), isg::Vec2d(0, -1), input, 2.0);
    REQUIRE(one_sided.crosswalk0);
    REQUIRE_FALSE(one_sided.crosswalk1);
    REQUIRE(one_sided.lead0 > 5.0);
    REQUIRE(one_sided.lead1 == 0.0);

    input.crosswalks.clear();
    const isg::UTurnLeadFloors fallback = builder.leadFloors(
        isg::Vec2d(0, 0), isg::Vec2d(0, 1),
        isg::Vec2d(10, 0), isg::Vec2d(0, -1), input, 2.0);
    REQUIRE_FALSE(fallback.crosswalk0);
    REQUIRE_FALSE(fallback.crosswalk1);
    REQUIRE(fallback.lead0 == 2.0);
    REQUIRE(fallback.lead1 == 2.0);
}

TEST_CASE("UTurnFamilyBuilder shares the deepest station across a family") {
    isg::IntersectionInput input;
    input.lanes.push_back(lane("entry", isg::Vec3d(0, -10), isg::Vec3d(0, 0)));
    input.lanes.push_back(lane("exit-near", isg::Vec3d(4, 0), isg::Vec3d(4, -10)));
    input.lanes.push_back(lane("exit-far", isg::Vec3d(8, 0), isg::Vec3d(8, -10)));
    isg::Connectivity near_connection;
    near_connection.id = "near";
    near_connection.entry_lane_id = "entry";
    near_connection.exit_lane_id = "exit-near";
    isg::Connectivity far_connection = near_connection;
    far_connection.id = "far";
    far_connection.exit_lane_id = "exit-far";
    input.connectivities = {near_connection, far_connection};

    const isg::UTurnFamilyBuilder builder;
    const double station = builder.alignmentComponentStation(
        input.connectivities.front(), input,
        isg::UTurnAlignmentScope::LaneEndpoint);
    REQUIRE(std::abs(station - 2.0) < 1e-9);

    double lead0 = 0.0;
    double lead1 = 0.0;
    builder.enforceAlignmentStation(
        isg::Vec2d(0, 0), isg::Vec2d(0, 1),
        isg::Vec2d(4, 0), isg::Vec2d(0, -1),
        station, lead0, lead1);
    REQUIRE(std::abs(lead0 - 2.0) < 1e-9);
    REQUIRE(std::abs(lead1 - 2.0) < 1e-9);
}

TEST_CASE("UTurnFamilyBuilder assembles a reusable initialization snapshot") {
    isg::IntersectionInput input;
    input.lanes.push_back(lane("entry", isg::Vec3d(0, -10), isg::Vec3d(0, 0)));
    input.lanes.push_back(lane("exit-near", isg::Vec3d(4, 0), isg::Vec3d(4, -10)));
    input.lanes.push_back(lane("exit-far", isg::Vec3d(8, 0), isg::Vec3d(8, -10)));
    isg::Crosswalk crosswalk;
    crosswalk.id = "entry-crosswalk";
    crosswalk.geometry.outer = {
        isg::Vec3d(-1, 3), isg::Vec3d(1, 3),
        isg::Vec3d(1, 5), isg::Vec3d(-1, 5)};
    input.crosswalks.push_back(crosswalk);

    isg::Connectivity near_connection;
    near_connection.id = "near";
    near_connection.entry_lane_id = "entry";
    near_connection.exit_lane_id = "exit-near";
    isg::Connectivity far_connection = near_connection;
    far_connection.id = "far";
    far_connection.exit_lane_id = "exit-far";
    input.connectivities = {near_connection, far_connection};

    const isg::UTurnFamilyInfo info = isg::UTurnFamilyBuilder().build(
        input.connectivities.front(), input,
        isg::UTurnAlignmentScope::LaneEndpoint);
    REQUIRE(info.geometric_uturn);
    REQUIRE(info.axis.isApprox(isg::Vec2d(0, 1)));
    REQUIRE(std::abs(info.radius - 4.0) < 1e-9);
    REQUIRE(std::abs(info.aligned_station - 5.3) < 1e-9);
    REQUIRE(std::abs(info.lead0 - 5.3) < 1e-9);
    REQUIRE(std::abs(info.lead1 - 5.3) < 1e-9);
    REQUIRE(info.rank.family_size == 2);
    REQUIRE(info.rank.reverse_radius_rank == 2);
    REQUIRE(info.clearance_crosswalks.size() == 1);

    isg::CurveInitializationOptions options;
    options.uturn_arc_alpha = 0.75;
    options.applyUTurnFamily(info);
    REQUIRE(std::abs(options.uturn_min_lead0 - info.lead0) < 1e-9);
    REQUIRE(std::abs(options.uturn_min_lead1 - info.lead1) < 1e-9);
    REQUIRE(std::abs(options.uturn_arc_alpha - 0.75) < 1e-9);
}

TEST_CASE("CurveGenerationContextBuilder binds scene frames and U-turn family") {
    isg::IntersectionInput input;
    input.lanes.push_back(lane("entry", isg::Vec3d(0, -10), isg::Vec3d(0, 0)));
    input.lanes.push_back(lane("exit", isg::Vec3d(4, 0), isg::Vec3d(4, -10)));
    isg::Connectivity connection;
    connection.id = "context-uturn";
    connection.entry_lane_id = "entry";
    connection.exit_lane_id = "exit";
    connection.turn_type = isg::ConnTurnType::UTurnLeft;
    input.connectivities.push_back(connection);

    isg::SceneContext scene(input);
    isg::ClusterOrderSolver topology;
    topology.build(input.connectivities, input.lanes, input.lane_groups,
                   input.crosswalks);
    const isg::CurveGenerationContext context =
        isg::CurveGenerationContextBuilder().build(
            scene, input.connectivities.front(),
            isg::UTurnAlignmentScope::LaneEndpoint, &topology);

    REQUIRE(context.scene == &scene);
    REQUIRE(context.connectivity == &input.connectivities.front());
    REQUIRE(context.entry.first.isApprox(isg::Vec2d(0, 0)));
    REQUIRE(context.exit.first.isApprox(isg::Vec2d(4, 0)));
    REQUIRE(context.turn == isg::ConnTurnType::UTurnLeft);
    REQUIRE(context.uturn_family.geometric_uturn);
    REQUIRE(std::abs(context.uturn_family.radius - 4.0) < 1e-9);
    REQUIRE(context.uturn_family.rank.family_size == 1);
}

TEST_CASE("SceneView indexes domain objects and endpoint frames") {
    isg::IntersectionInput input;
    input.lanes.push_back(lane("entry", isg::Vec3d(0, -10), isg::Vec3d(0, 0)));
    isg::SceneView view(input);
    REQUIRE(view.lane("entry") != nullptr);
    REQUIRE(view.lane("missing") == nullptr);
    REQUIRE(view.entryFrame("entry").first.isApprox(isg::Vec2d(0, 0)));
}

TEST_CASE("SceneContext exposes bound fields and turn metadata") {
    isg::IntersectionInput input;
    input.lanes.push_back(lane("entry", isg::Vec3d(0, -10), isg::Vec3d(0, 0)));
    input.lanes.push_back(lane("exit", isg::Vec3d(0, 10), isg::Vec3d(0, 20)));
    isg::Connectivity connection;
    connection.id = "connection";
    connection.entry_lane_id = "entry";
    connection.exit_lane_id = "exit";
    connection.turn_type = isg::ConnTurnType::Straight;
    input.connectivities.push_back(connection);
    isg::SDFField fine;
    isg::SDFField coarse;
    isg::ClusterTopology topology;
    isg::SceneContext context(input);
    context.bindSdf(&fine, &coarse);
    context.bindClusterTopology(&topology);

    REQUIRE(context.sdf == &fine);
    REQUIRE(context.coarse_sdf == &coarse);
    REQUIRE(context.cluster_topology == &topology);
    REQUIRE(context.turn_metadata.at("connection").turn == isg::ConnTurnType::Straight);
    REQUIRE(context.turn_metadata.at("connection").entry_point.isApprox(isg::Vec2d(0, 0)));
}

TEST_CASE("ConstraintEvaluator audits obstacle boundary shape and continuity") {
    isg::IntersectionInput input;
    isg::Obstacle obstacle;
    obstacle.id = "obstacle";
    obstacle.geometry.outer = {isg::Vec3d(-1, -1), isg::Vec3d(1, -1),
                               isg::Vec3d(1, 1), isg::Vec3d(-1, 1)};
    input.obstacles.push_back(obstacle);
    isg::Boundary boundary;
    boundary.id = "median";
    boundary.type = isg::Boundary::Type::MedianStrip;
    boundary.geometry.points = {isg::Vec3d(-3, 3), isg::Vec3d(3, 3)};
    input.boundaries.push_back(boundary);
    isg::SceneContext context(input);
    isg::ConstraintProfile profile;
    profile.enforce_fence = false;

    const isg::ConstraintReport obstacle_report = isg::ConstraintEvaluator().evaluate(
        cubic(isg::Vec2d(-2, 0), isg::Vec2d(-1, 0), isg::Vec2d(1, 0), isg::Vec2d(2, 0)),
        context, profile);
    REQUIRE(obstacle_report.hasHardViolation());

    profile.check_obstacle = false;
    const isg::ConstraintReport boundary_report = isg::ConstraintEvaluator().evaluate(
        cubic(isg::Vec2d(0, 2), isg::Vec2d(0, 2.5), isg::Vec2d(0, 3.5), isg::Vec2d(0, 4)),
        context, profile);
    REQUIRE(boundary_report.hasHardViolation());

    input.boundaries.front().type = isg::Boundary::Type::RoadEdge;
    isg::SceneContext road_edge_context(input);
    profile.road_edge_clearance = 1.0;
    const isg::ConstraintReport clearance_report = isg::ConstraintEvaluator().evaluate(
        cubic(isg::Vec2d(-2, 2.5), isg::Vec2d(-1, 2.5),
              isg::Vec2d(1, 2.5), isg::Vec2d(2, 2.5)), road_edge_context, profile);
    REQUIRE(clearance_report.hasHardViolation());

    profile.check_boundary = false;
    profile.check_curvature = true;
    profile.max_curvature = 0.01;
    const isg::ConstraintReport curvature_report = isg::ConstraintEvaluator().evaluate(
        cubic(isg::Vec2d(0, 0), isg::Vec2d(0, 2), isg::Vec2d(2, 2), isg::Vec2d(2, 0)),
        context, profile);
    REQUIRE(curvature_report.hasHardViolation());

    profile.check_curvature = false;
    profile.check_g1 = true;
    isg::BezierCurve broken_g1;
    broken_g1.segs.push_back(cubic(isg::Vec2d(0, 0), isg::Vec2d(1, 0),
                                    isg::Vec2d(2, 0), isg::Vec2d(3, 0)).segs.front());
    broken_g1.segs.push_back(cubic(isg::Vec2d(3, 0), isg::Vec2d(3, 1),
                                    isg::Vec2d(3, 2), isg::Vec2d(3, 3)).segs.front());
    REQUIRE(isg::ConstraintEvaluator().evaluate(broken_g1, context, profile).hasHardViolation());

    profile.check_g1 = false;
    isg::BezierCurve self_crossing;
    self_crossing.segs.push_back(cubic(isg::Vec2d(-2, -2), isg::Vec2d(-1, -1),
                                       isg::Vec2d(1, 1), isg::Vec2d(2, 2)).segs.front());
    self_crossing.segs.push_back(cubic(isg::Vec2d(2, 2), isg::Vec2d(1, 2),
                                       isg::Vec2d(-1, 2), isg::Vec2d(-2, 2)).segs.front());
    self_crossing.segs.push_back(cubic(isg::Vec2d(-2, 2), isg::Vec2d(-1, 1),
                                       isg::Vec2d(1, -1), isg::Vec2d(2, -2)).segs.front());
    REQUIRE(isg::ConstraintEvaluator().evaluate(self_crossing, context, profile).hasHardViolation());
}

TEST_CASE("GenerationPlanner keeps straight turns ahead of turns and U-turns") {
    isg::IntersectionInput input;
    input.lanes.push_back(lane("entry", isg::Vec3d(0, -10), isg::Vec3d(0, 0)));
    input.lanes.push_back(lane("straight", isg::Vec3d(0, 10), isg::Vec3d(0, 20)));
    input.lanes.push_back(lane("left", isg::Vec3d(-10, 0), isg::Vec3d(-20, 0)));
    input.lanes.push_back(lane("uturn", isg::Vec3d(1, -1), isg::Vec3d(1, -10)));
    const char* ids[] = {"straight-conn", "left-conn", "uturn-conn"};
    const char* exits[] = {"straight", "left", "uturn"};
    for (int i = 0; i < 3; ++i) {
        isg::Connectivity conn;
        conn.id = ids[i];
        conn.entry_lane_id = "entry";
        conn.exit_lane_id = exits[i];
        input.connectivities.push_back(conn);
    }
    const isg::GenerationPlan plan = isg::GenerationPlanner().build(input);
    REQUIRE(plan.batches.size() == 3);
    REQUIRE(plan.batches[0].conn_ids.size() == 1);
    REQUIRE(plan.batches[1].conn_ids.size() == 1);
    REQUIRE(plan.batches[2].conn_ids.size() == 1);
}

TEST_CASE("Unified evaluator checks ordinary shape and fixed-shape policy") {
    isg::IntersectionInput input;
    input.lanes.push_back(lane("entry", isg::Vec3d(0, -10), isg::Vec3d(0, 0)));
    input.lanes.push_back(lane("exit", isg::Vec3d(0, 10), isg::Vec3d(0, 20)));
    isg::Connectivity conn;
    conn.id = "ordinary";
    conn.entry_lane_id = "entry";
    conn.exit_lane_id = "exit";
    input.connectivities.push_back(conn);
    isg::SceneContext scene(input);
    isg::CurveGenerationContext context;
    context.scene = &scene;
    context.connectivity = &input.connectivities.front();
    context.entry = input.entryPtDir("entry");
    context.exit = input.exitPtDir("exit");
    context.profile.enforce_fence = false;
    context.profile.check_ordinary_shape = true;
    isg::GenerationState state;

    const isg::BezierCurve ordinary = cubic(isg::Vec2d(0, 0), isg::Vec2d(0, 3),
                                             isg::Vec2d(0, 7), isg::Vec2d(0, 10));
    REQUIRE_FALSE(isg::ConstraintEvaluator().evaluate(ordinary, context, state).hasHardViolation());

    const isg::BezierCurve multi = isg::makeCurveFromKnots(
        {isg::Vec2d(0, 0), isg::Vec2d(0, 5), isg::Vec2d(0, 10)},
        {isg::Vec2d(0, 1), isg::Vec2d(0, 1), isg::Vec2d(0, 1)});
    REQUIRE(isg::ConstraintEvaluator().evaluate(multi, context, state).hasHardViolation());

    isg::Connectivity fixed_conn = input.connectivities.front();
    fixed_conn.id = "fixed";
    fixed_conn.fixed_shape = true;
    input.connectivities.push_back(fixed_conn);
    isg::CurveGenerationContext fixed_context = context;
    fixed_context.connectivity = &input.connectivities.back();
    fixed_context.profile.check_ordinary_shape = false;
    fixed_context.profile.require_fixed_shape_preservation = true;
    isg::CurveCandidate candidate;
    candidate.curve = ordinary;
    candidate.preserves_fixed_shape = false;
    REQUIRE(isg::ConstraintEvaluator().evaluate(candidate, fixed_context, state).hasHardViolation());
    candidate.preserves_fixed_shape = true;
    REQUIRE_FALSE(isg::ConstraintEvaluator().evaluate(candidate, fixed_context, state).hasHardViolation());

    isg::ConstraintEvaluator evaluator;
    evaluator.audit(candidate, fixed_context, state);
    REQUIRE(candidate.hasCurrentReport());
    candidate.replaceCurve(multi);
    REQUIRE_FALSE(candidate.hasCurrentReport());
    REQUIRE(candidate.report.results.empty());
}

TEST_CASE("Unified evaluator checks segmented U-turn, Crosswalk and cluster pairs") {
    isg::IntersectionInput input;
    input.lanes.push_back(lane("entry", isg::Vec3d(0, -12), isg::Vec3d(0, 0)));
    input.lanes.push_back(lane("exit", isg::Vec3d(4, 0), isg::Vec3d(4, -12)));
    input.crosswalks.push_back(isg::Crosswalk());
    input.crosswalks.front().id = "cw";
    input.crosswalks.front().geometry.outer = {
        isg::Vec3d(-2, 2), isg::Vec3d(2, 2), isg::Vec3d(2, 4), isg::Vec3d(-2, 4)};
    isg::Connectivity conn;
    conn.id = "uturn";
    conn.entry_lane_id = "entry";
    conn.exit_lane_id = "exit";
    input.connectivities.push_back(conn);
    isg::ClusterTopology topology;
    isg::SceneContext scene(input);
    scene.bindClusterTopology(&topology);
    isg::CurveGenerationContext context;
    context.scene = &scene;
    context.connectivity = &input.connectivities.front();
    context.entry = input.entryPtDir("entry");
    context.exit = input.exitPtDir("exit");
    context.profile.enforce_fence = false;
    context.profile.check_uturn_shape = true;
    context.profile.check_crosswalk = true;
    isg::GenerationState state;

    isg::BezierCurve uturn;
    uturn.segs.push_back(cubic(isg::Vec2d(0, 0), isg::Vec2d(0, 1.7),
                                isg::Vec2d(0, 3.3), isg::Vec2d(0, 5)).segs.front());
    uturn.segs.push_back(cubic(isg::Vec2d(0, 5), isg::Vec2d(0, 7.7),
                                isg::Vec2d(4, 7.7), isg::Vec2d(4, 5)).segs.front());
    uturn.segs.push_back(cubic(isg::Vec2d(4, 5), isg::Vec2d(4, 3.3),
                                isg::Vec2d(4, 1.7), isg::Vec2d(4, 0)).segs.front());
    REQUIRE_FALSE(isg::ConstraintEvaluator().evaluate(uturn, context, state).hasHardViolation());

    isg::BezierCurve crosswalk_intrusion = uturn;
    crosswalk_intrusion.segs[1] = cubic(isg::Vec2d(0, 3), isg::Vec2d(3, 3),
                                        isg::Vec2d(3, 3), isg::Vec2d(0, 3)).segs.front();
    REQUIRE(isg::ConstraintEvaluator().evaluate(crosswalk_intrusion, context, state).hasHardViolation());

    isg::ClusterTopology crossing_topology;
    isg::CurvePair pair;
    pair.id_a = "uturn";
    pair.id_b = "sibling";
    crossing_topology.pairs.push_back(pair);
    scene.bindClusterTopology(&crossing_topology);
    isg::BezierCurve sibling = cubic(isg::Vec2d(-2, 8), isg::Vec2d(-1, 8),
                                      isg::Vec2d(1, 8), isg::Vec2d(2, 8));
    state.accepted_curves["sibling"] = sibling;
    context.profile.check_uturn_shape = false;
    context.profile.check_crosswalk = false;
    context.profile.check_cluster = true;
    const isg::BezierCurve crossing = cubic(isg::Vec2d(-2, 6), isg::Vec2d(-1, 7),
                                             isg::Vec2d(1, 9), isg::Vec2d(2, 10));
    REQUIRE(isg::ConstraintEvaluator().evaluate(crossing, context, state).hasHardViolation());

    crossing_topology.pairs.front().exempt = isg::CrossExemption::StructuralCross;
    context.profile.allow_structural_cross = true;
    REQUIRE_FALSE(isg::ConstraintEvaluator().evaluate(crossing, context, state).hasHardViolation());
    const isg::BezierCurve overlap = cubic(isg::Vec2d(-2, 8.1), isg::Vec2d(-1, 8.1),
                                            isg::Vec2d(1, 8.1), isg::Vec2d(2, 8.1));
    REQUIRE(isg::ConstraintEvaluator().evaluate(overlap, context, state).hasHardViolation());
}

TEST_CASE("JsonCodec accepts legacy 2D points and preserves elevation") {
    const isg::Vec3d legacy = isg::vec3dFromJson(
        nlohmann::json::parse("{\"x\":1.5,\"y\":-2.0}"));
    REQUIRE(legacy.x() == 1.5);
    REQUIRE(legacy.y() == -2.0);
    REQUIRE(legacy.z() == 0.0);

    const isg::Vec3d elevated(3.0, 4.0, 12.5);
    const isg::Vec3d decoded = isg::vec3dFromJson(isg::vec3dToJson(elevated));
    REQUIRE(decoded.x() == elevated.x());
    REQUIRE(decoded.y() == elevated.y());
    REQUIRE(decoded.z() == elevated.z());
}

TEST_CASE("CurvePostProcessor restores exact endpoints and endpoint G1") {
    isg::LBFGSSolver solver;
    const isg::CurvePostProcessor processor(solver);
    isg::SDFField sdf;
    const isg::Polygon2d fence;
    const isg::Vec2d entry(0, 0);
    const isg::Vec2d exit(10, 0);
    const isg::Vec2d entry_tangent(1, 0);
    const isg::Vec2d exit_tangent(1, 0);
    const isg::BezierCurve source = cubic(
        isg::Vec2d(0.1, 0.2), isg::Vec2d(2, 1),
        isg::Vec2d(8, -1), isg::Vec2d(9.9, -0.2));

    const isg::BezierCurve processed = processor.process(
        source, sdf, fence, 10.0, entry_tangent, exit_tangent,
        true, &entry, &exit);
    REQUIRE(processed.numSegments() == 1);
    REQUIRE(processed.startPt().isApprox(entry));
    REQUIRE(processed.endPt().isApprox(exit));
    REQUIRE(processed.startTan().normalized().isApprox(entry_tangent));
    REQUIRE(processed.endTan().normalized().isApprox(exit_tangent));
}

TEST_CASE("CurveOptimizer preserves endpoint contract for an ordinary cubic") {
    isg::IntersectionInput input;
    input.lanes.push_back(lane(
        "entry", isg::Vec3d(0, 0), isg::Vec3d(0, 1)));
    input.lanes.push_back(lane(
        "exit", isg::Vec3d(10, 0), isg::Vec3d(10, 1)));
    isg::SDFField sdf;
    isg::LBFGSSolver solver;
    const isg::CurveOptimizer optimizer(solver);
    const isg::BezierCurve initial = cubic(
        isg::Vec2d(0, 0), isg::Vec2d(2, 0.5),
        isg::Vec2d(8, -0.5), isg::Vec2d(10, 0));
    isg::CurveOptimizationOptions options;
    options.outer_iterations = 1;
    options.enforce_fence = false;
    options.constrain_single_cubic_axes = true;

    const isg::BezierCurve optimized = optimizer.optimize(
        initial, input, sdf, std::vector<isg::SiblingCurve>(),
        isg::Vec2d(1, 0), isg::Vec2d(1, 0), options);
    REQUIRE_FALSE(optimized.empty());
    REQUIRE(optimized.startPt().isApprox(initial.startPt(), 1e-6));
    REQUIRE(optimized.endPt().isApprox(initial.endPt(), 1e-6));
    REQUIRE(optimized.startTan().normalized().dot(isg::Vec2d(1, 0)) > 0.99);
    REQUIRE(optimized.endTan().normalized().dot(isg::Vec2d(1, 0)) > 0.99);
}

TEST_CASE("OptimizationResultProcessor falls back from divergent optimization") {
    isg::LBFGSSolver solver;
    const isg::OptimizationResultProcessor processor(solver);
    isg::SDFField sdf;
    const isg::Polygon2d fence;
    const isg::BezierCurve initial = cubic(
        isg::Vec2d(0, 0), isg::Vec2d(2, 0),
        isg::Vec2d(8, 0), isg::Vec2d(10, 0));
    const isg::BezierCurve divergent = cubic(
        isg::Vec2d(0, 0), isg::Vec2d(1000, 1000),
        isg::Vec2d(-1000, -1000), isg::Vec2d(10, 0));

    const isg::BezierCurve result = processor.process(
        initial, divergent, sdf, fence,
        isg::Vec2d(0, 0), isg::Vec2d(1, 0),
        isg::Vec2d(10, 0), isg::Vec2d(1, 0));
    REQUIRE(result.startPt().isApprox(initial.startPt()));
    REQUIRE(result.endPt().isApprox(initial.endPt()));
    REQUIRE(result.arcLength() < divergent.arcLength());
}

TEST_CASE("UTurnEnvelopeConstraint detects expansion and collapse") {
    const isg::UTurnEnvelopeConstraint constraint;
    const isg::Vec2d entry(0, 0);
    const isg::Vec2d exit(10, 0);
    const isg::BezierCurve reference = cubic(
        entry, isg::Vec2d(0, 8), isg::Vec2d(10, 8), exit);
    const isg::BezierCurve collapsed = cubic(
        entry, isg::Vec2d(3, 0.1), isg::Vec2d(7, 0.1), exit);
    const isg::BezierCurve expanded = cubic(
        entry, isg::Vec2d(200, 200), isg::Vec2d(-200, 200), exit);

    REQUIRE_FALSE(constraint.exceeds(reference, reference, entry, exit));
    REQUIRE_FALSE(constraint.collapses(reference, reference, entry, exit));
    REQUIRE(constraint.collapses(collapsed, reference, entry, exit));
    REQUIRE(constraint.exceeds(expanded, reference, entry, exit));
}

TEST_CASE("BoundedUTurnCandidateSearch delegates audit and selects safely") {
    const isg::Vec2d entry(0, 0);
    const isg::Vec2d exit(4, 0);
    const isg::Vec2d entry_tangent(0, 1);
    const isg::Vec2d exit_tangent(0, -1);
    const isg::BezierCurve reference =
        isg::UTurnCurveInitializer().buildAligned(
            entry, entry_tangent, exit, exit_tangent,
            entry_tangent, 0.0, 1.0, 2.0, 2.0);
    const isg::BoundedUTurnCandidateSearch search;
    isg::BezierCurve rejected;
    const bool has_rejected = search.search(
        entry, entry_tangent, exit, exit_tangent, reference,
        [](const isg::BezierCurve&) {
            isg::BoundedUTurnCandidateAudit audit;
            audit.physical_violation = true;
            return audit;
        },
        rejected, 2.0, 2.0);
    REQUIRE_FALSE(has_rejected);

    int audited = 0;
    isg::BezierCurve selected;
    const bool found = search.search(
        entry, entry_tangent, exit, exit_tangent, reference,
        [&](const isg::BezierCurve&) {
            ++audited;
            isg::BoundedUTurnCandidateAudit audit;
            audit.sibling_crosses = 0;
            return audit;
        },
        selected, 2.0, 2.0);
    REQUIRE(found);
    REQUIRE(audited > 0);
    REQUIRE(selected.startPt().isApprox(entry));
    REQUIRE(selected.endPt().isApprox(exit));
    REQUIRE(selected.startTan().normalized().dot(entry_tangent) > 0.99);
    REQUIRE(selected.endTan().normalized().dot(exit_tangent) > 0.99);
}

TEST_CASE("PhysicalRepairCoordinator preserves repair order and U-turn fallback") {
    const isg::PhysicalRepairCoordinator coordinator;
    const isg::BezierCurve initial = cubic(
        isg::Vec2d(0, 0), isg::Vec2d(1, 0),
        isg::Vec2d(3, 0), isg::Vec2d(4, 0));
    const isg::BezierCurve current = cubic(
        isg::Vec2d(0, 0), isg::Vec2d(1, 3),
        isg::Vec2d(3, 3), isg::Vec2d(4, 0));
    std::vector<std::string> calls;
    const isg::PhysicalRiskAuditor auditor =
        [&](const isg::BezierCurve& curve) {
            calls.push_back("audit");
            isg::PhysicalRiskSnapshot risk;
            const double height = curve.segs.front().ctrl[1].y();
            risk.boundary = height > 2.5;
            risk.obstacle = height > 0.5 && height <= 2.5;
            return risk;
        };
    const isg::PhysicalRepairAttempt boundary_repair =
        [&](const isg::BezierCurve&, isg::BezierCurve& repaired) {
            calls.push_back("boundary");
            repaired = cubic(
                isg::Vec2d(0, 0), isg::Vec2d(1, 2),
                isg::Vec2d(3, 2), isg::Vec2d(4, 0));
            return true;
        };
    const isg::PhysicalRepairAttempt obstacle_repair =
        [&](const isg::BezierCurve&, isg::BezierCurve& repaired) {
            calls.push_back("obstacle");
            repaired = initial;
            return true;
        };

    const isg::PhysicalRepairResult ordinary = coordinator.repair(
        initial, current, false, auditor,
        boundary_repair, obstacle_repair);
    REQUIRE(ordinary.boundary_repaired);
    REQUIRE_FALSE(ordinary.risk.physical());
    REQUIRE(ordinary.curve.segs.front().ctrl[1].y() == 0.0);
    REQUIRE(calls == std::vector<std::string>{
        "audit", "boundary", "audit", "obstacle", "audit"});

    calls.clear();
    const isg::PhysicalRepairResult uturn = coordinator.repair(
        initial, current, true, auditor,
        boundary_repair, obstacle_repair);
    REQUIRE_FALSE(uturn.boundary_repaired);
    REQUIRE_FALSE(uturn.risk.physical());
    REQUIRE(uturn.curve.segs.front().ctrl[1].y() == 0.0);
    REQUIRE(calls == std::vector<std::string>{"audit", "audit"});
}

TEST_CASE("FinalCurveAuditor preserves status and reason precedence") {
    const isg::FinalCurveAuditor auditor;
    isg::ConnectivityCurve curve;
    curve.violation.max_fence_overflow = 3.0;

    isg::FinalCurveAuditSnapshot clean;
    auditor.apply(clean, curve);
    REQUIRE(curve.status == isg::CurveStatus::OK);
    REQUIRE(curve.violation.reason.empty());
    REQUIRE(curve.violation.max_fence_overflow == 3.0);

    curve.violation.exempt_crosses.push_back(isg::Vec2d(1, 2));
    auditor.apply(clean, curve);
    REQUIRE(curve.status == isg::CurveStatus::WarnA2);

    isg::FinalCurveAuditSnapshot violations;
    violations.obstacle_intersection = true;
    violations.self_intersection = true;
    violations.boundary_intersection = true;
    violations.uturn = true;
    violations.road_edge_clearance_violation = true;
    violations.ordinary_single_axis_violation = true;
    violations.update_fence_overflow = true;
    violations.max_fence_overflow = 0.25;
    auditor.apply(violations, curve);
    REQUIRE(curve.status == isg::CurveStatus::Degraded);
    REQUIRE(curve.violation.max_obstacle_penetration == 0.01);
    REQUIRE(curve.violation.max_fence_overflow == 0.25);
    REQUIRE(curve.violation.reason ==
            "ordinary single cubic control points leave endpoint direction axes");

    violations.ordinary_single_axis_violation = false;
    auditor.apply(violations, curve);
    REQUIRE(curve.violation.reason == "curve violates RoadEdge clearance");
}

TEST_CASE("RepairBudget preserves bounded global repair defaults") {
    const isg::RepairBudget budget;
    REQUIRE(budget.mode2_uturn_pair_passes == 2);
    REQUIRE(budget.mode2_uturn_pairs_per_pass == 4);
    REQUIRE(budget.shared_endpoint_pair_repairs == 20);
    REQUIRE(budget.early_straight_uturn_repairs == 4);
    REQUIRE(budget.final_straight_uturn_repairs == 6);
    REQUIRE(budget.pure_topology_passes == 3);
    REQUIRE(budget.pure_topology_repairs_per_pass == 6);
    REQUIRE(budget.shape_restores == 8);
    REQUIRE(budget.base_expression_passes == 1);
    REQUIRE(budget.pair_safe_passes == 3);
    REQUIRE(budget.pair_safe_phases == 4);
    REQUIRE(budget.regeneration_passes == 1);
    REQUIRE(budget.regenerations_per_pass == 8);
}

TEST_CASE("RepairImpactClosureBuilder limits audit to constrained neighbors") {
    isg::CurvePair ab;
    ab.id_a = "a";
    ab.id_b = "b";
    isg::CurvePair bc;
    bc.id_a = "b";
    bc.id_b = "c";
    isg::CurvePair ad;
    ad.id_a = "a";
    ad.id_b = "d";
    ad.exempt = isg::CrossExemption::StructuralCross;
    const std::vector<isg::CurvePair> pairs{ab, bc, ad};

    const isg::RepairImpactClosure closure =
        isg::RepairImpactClosureBuilder().build(
            pairs, std::unordered_set<isg::ConnId>{"a"});
    REQUIRE(closure.pairs.size() == 1);
    REQUIRE(closure.curve_ids.size() == 2);
    REQUIRE(closure.curve_ids.count("a") == 1);
    REQUIRE(closure.curve_ids.count("b") == 1);
    REQUIRE(closure.curve_ids.count("c") == 0);
    REQUIRE(closure.curve_ids.count("d") == 0);
    REQUIRE(closure.neighbors.at("a") == std::vector<isg::ConnId>{"b"});
    REQUIRE(closure.neighbors.at("b") == std::vector<isg::ConnId>{"a"});

    const auto neighbors = isg::RepairImpactClosureBuilder().neighbors(
        pairs, std::unordered_set<isg::ConnId>{"a", "b", "c"});
    REQUIRE(neighbors.at("a") == std::vector<isg::ConnId>{"b"});
    REQUIRE(neighbors.at("b") ==
            std::vector<isg::ConnId>{"a", "c"});
    REQUIRE(neighbors.at("c") == std::vector<isg::ConnId>{"b"});
    REQUIRE(neighbors.count("d") == 0);
}

TEST_CASE("AtomicCurvePatch commits all prepared curves or none") {
    std::vector<isg::ConnectivityCurve> results(2);
    results[0].id = "a";
    results[1].id = "b";
    const isg::BezierCurve original_a = cubic(
        isg::Vec2d(0, 0), isg::Vec2d(1, 0),
        isg::Vec2d(2, 0), isg::Vec2d(3, 0));
    const isg::BezierCurve original_b = cubic(
        isg::Vec2d(0, 1), isg::Vec2d(1, 1),
        isg::Vec2d(2, 1), isg::Vec2d(3, 1));
    results[0].curve.reset(new isg::BezierCurve(original_a));
    results[1].curve.reset(new isg::BezierCurve(original_b));

    const isg::BezierCurve replacement_a = cubic(
        isg::Vec2d(0, 0), isg::Vec2d(1, 2),
        isg::Vec2d(2, 2), isg::Vec2d(3, 0));
    const isg::BezierCurve replacement_b = cubic(
        isg::Vec2d(0, 1), isg::Vec2d(1, 3),
        isg::Vec2d(2, 3), isg::Vec2d(3, 1));
    int prepared = 0;
    const isg::CurvePatchPreparer prepare =
        [&](isg::ConnectivityCurve& result, const isg::BezierCurve& curve) {
            ++prepared;
            result.curve.reset(new isg::BezierCurve(curve));
            result.status = isg::CurveStatus::WarnA2;
        };

    const bool committed = isg::AtomicCurvePatch().apply(
        std::vector<isg::CurvePatchEntry>{
            isg::CurvePatchEntry{"a", replacement_a},
            isg::CurvePatchEntry{"b", replacement_b}},
        results, prepare);
    REQUIRE(committed);
    REQUIRE(prepared == 2);
    REQUIRE(results[0].curve->segs.front().ctrl[1].y() == 2.0);
    REQUIRE(results[1].curve->segs.front().ctrl[1].y() == 3.0);
    REQUIRE(results[0].status == isg::CurveStatus::WarnA2);
    REQUIRE(results[1].status == isg::CurveStatus::WarnA2);

    prepared = 0;
    const bool rejected = isg::AtomicCurvePatch().apply(
        std::vector<isg::CurvePatchEntry>{
            isg::CurvePatchEntry{"a", original_a},
            isg::CurvePatchEntry{"missing", original_b}},
        results, prepare);
    REQUIRE_FALSE(rejected);
    REQUIRE(prepared == 0);
    REQUIRE(results[0].curve->segs.front().ctrl[1].y() == 2.0);
    REQUIRE(results[1].curve->segs.front().ctrl[1].y() == 3.0);
}

TEST_CASE("FinalPairAuditor audits every applicable pair and maps severity") {
    std::vector<isg::ConnectivityCurve> results(3);
    results[0].id = "a";
    results[1].id = "b";
    results[2].id = "c";
    const isg::BezierCurve curve = cubic(
        isg::Vec2d(0, 0), isg::Vec2d(1, 0),
        isg::Vec2d(2, 0), isg::Vec2d(3, 0));
    for (auto& result : results)
        result.curve.reset(new isg::BezierCurve(curve));

    isg::CurvePair ab;
    ab.id_a = "a";
    ab.id_b = "b";
    ab.exempt = isg::CrossExemption::ObstacleCross;
    isg::CurvePair bc;
    bc.id_a = "b";
    bc.id_b = "c";
    isg::CurvePair ac;
    ac.id_a = "a";
    ac.id_b = "c";
    ac.exempt = isg::CrossExemption::StructuralCross;
    int audited = 0;
    const isg::FinalPairViolationDetector detector =
        [&](const isg::CurvePair&, const isg::BezierCurve&,
            const isg::BezierCurve&, isg::Vec2d& location) {
            ++audited;
            location = isg::Vec2d(audited, -audited);
            return true;
        };

    const isg::FinalPairAuditor auditor;
    const std::vector<isg::FinalPairFinding> findings = auditor.audit(
        std::vector<isg::CurvePair>{ab, bc, ac}, results, detector);
    REQUIRE(audited == 2);
    REQUIRE(findings.size() == 2);
    REQUIRE(findings[0].id_a == "a");
    REQUIRE(findings[1].id_b == "c");

    auditor.apply(findings, results);
    REQUIRE(results[0].status == isg::CurveStatus::WarnA2);
    REQUIRE(results[1].status == isg::CurveStatus::Degraded);
    REQUIRE(results[2].status == isg::CurveStatus::Degraded);
    REQUIRE(results[0].violation.exempt_crosses.size() == 1);
    REQUIRE(results[1].violation.exempt_crosses.size() == 2);
    REQUIRE(results[2].violation.exempt_crosses.size() == 1);
}

TEST_CASE("IntersectionIO string roundtrip preserves domain references") {
    isg::IntersectionInput input;
    input.id = "roundtrip";
    input.mode = 2;
    input.lanes.push_back(lane("entry", isg::Vec3d(0, -10, 2), isg::Vec3d(0, 0, 3)));
    isg::Connectivity connection;
    connection.id = "connection";
    connection.entry_lane_id = "entry";
    connection.enterGroupId = "group-entry";
    input.connectivities.push_back(connection);

    const isg::IntersectionInput decoded = isg::IntersectionIO::fromJsonString(
        isg::IntersectionIO::toJsonString(input));
    REQUIRE(decoded.id == input.id);
    REQUIRE(decoded.mode == input.mode);
    REQUIRE(decoded.lanes.front().geometry.points.back().z() == 3.0);
    REQUIRE(decoded.connectivities.front().enterGroupId == "group-entry");
}

TEST_CASE("CurveSampling preserves sibling metadata and exempt sampling policy") {
    isg::SiblingCurve regular;
    regular.id = "regular";
    regular.curve = cubic(
        isg::Vec2d(0, 0), isg::Vec2d(1, 0),
        isg::Vec2d(2, 0), isg::Vec2d(3, 0));
    regular.expected_side = 1;
    regular.ref_perp = isg::Vec2d(0, 1);
    regular.shared_endpoint = true;
    isg::SiblingCurve exempt = regular;
    exempt.id = "exempt";
    exempt.exempt_a1 = true;

    const auto sampled = isg::sampleSiblingsForIntersections(
        std::vector<isg::SiblingCurve>{regular, exempt}, 8, false);
    REQUIRE(sampled.size() == 2);
    REQUIRE_FALSE(sampled[0].sampled.pts.empty());
    REQUIRE(sampled[0].expected_side == 1);
    REQUIRE(sampled[0].ref_perp.isApprox(isg::Vec2d(0, 1)));
    REQUIRE(sampled[0].shared_endpoint);
    REQUIRE(sampled[1].sampled.pts.empty());
    REQUIRE(sampled[1].exempt_a1);
}

TEST_CASE("UTurnMultiConstraintSearch delegates audits and preserves endpoint G1") {
    const isg::Vec2d p0(0, 0);
    const isg::Vec2d t0(1, 0);
    const isg::Vec2d p1(0, 4);
    const isg::Vec2d t1(-1, 0);
    const isg::BezierCurve reference = isg::UTurnCurveInitializer().buildAligned(
        p0, t0, p1, t1, isg::Vec2d(0, 1), 0.0);
    int samples = 0;
    int audits = 0;
    isg::UTurnSearchBackend backend;
    backend.sample = [&](const isg::BezierCurve& curve) {
        ++samples;
        return isg::sampleCurveForIntersections(curve);
    };
    backend.sibling_cross_count = [&](const isg::BezierCurve&) {
        ++audits;
        return 0;
    };
    backend.obstacle_penalty = [](const isg::BezierCurve&) { return 0.0; };
    backend.boundary_penalty = [](const isg::BezierCurve&) { return 0.0; };
    backend.endpoint_side_violation = [](const isg::SampledCurve&) { return 0.0; };

    isg::IntersectionInput input;
    input.area.is_rough = true;
    const isg::UTurnSolverResult result = isg::UTurnMultiConstraintSearch().search(
        p0, t0, p1, t1, input, {}, reference, backend, 0.0, 0.0);
    REQUIRE(std::isfinite(result.cost));
    REQUIRE_FALSE(result.curve.empty());
    REQUIRE(samples > 0);
    REQUIRE(audits == samples);
    REQUIRE(result.curve.startTan().normalized().dot(t0) > 0.99);
    REQUIRE(result.curve.endTan().normalized().dot(t1) > 0.99);
}

TEST_CASE("UTurnShape recognizes segmented leads and crosswalk clearance") {
    const isg::Vec2d p0(0, 0);
    const isg::Vec2d t0(1, 0);
    const isg::Vec2d p1(0, 4);
    const isg::Vec2d t1(-1, 0);
    const isg::BezierCurve curve = isg::UTurnCurveInitializer().buildSegmented(
        p0, t0, p1, t1, 1.0, 1.0, 2.0 / 3.0);
    REQUIRE(curve.numSegments() == 3);
    REQUIRE(isg::segmentLooksStraight(curve.segs.front()));
    REQUIRE(isg::segmentLooksStraight(curve.segs.back()));
    REQUIRE(isg::curveLooksUTurnForClusterExemption(curve));
    REQUIRE(isg::segmentedUTurnMiddleArcLooksRound(curve, t0));
    REQUIRE(isg::segmentedUTurnMiddleArcClearsCrosswalks(curve, {}));
}

TEST_CASE("SegmentedUTurnCandidateSearch delegates hard audits") {
    const isg::Vec2d p0(0, 0);
    const isg::Vec2d t0(1, 0);
    const isg::Vec2d p1(0, 4);
    const isg::Vec2d t1(-1, 0);
    isg::BezierCurve curve = isg::UTurnCurveInitializer().buildAligned(
        p0, t0, p1, t1, isg::Vec2d(0, 1), 0.0);
    isg::IntersectionInput input;
    input.area.is_rough = true;
    int audits = 0;
    const isg::SegmentedUTurnAuditor auditor =
        [&](const isg::BezierCurve&, bool) {
            ++audits;
            return isg::SegmentedUTurnAudit();
        };
    const bool found = isg::SegmentedUTurnCandidateSearch().search(
        p0, t0, p1, t1, input, {}, auditor, false, curve, 1.0, 1.0);
    REQUIRE(found);
    REQUIRE(audits > 1);
    REQUIRE(curve.numSegments() == 3);
    REQUIRE(curve.startTan().normalized().dot(t0) > 0.99);
    REQUIRE(curve.endTan().normalized().dot(t1) > 0.99);
}

TEST_CASE("singleCubicSignedEndCurvature orders a shared-endpoint fan-out") {
    // 共享端点与共享端点切向的一族单段三次曲线，其端点邻域左右次序完全由
    // 有符号端点曲率决定；本例验证符号语义、闭式解和"h0 单调不足以定序"。
    const isg::Vec2d p0(0, 0);
    const isg::Vec2d t0(1, 0);
    const isg::Vec2d p1(20, 10);       // 远端在进入切向左侧
    const isg::Vec2d t1 = isg::Vec2d(1, 1).normalized();
    auto make = [&](double h0, double h1) {
        isg::BezierSegment seg;
        seg.ctrl[0] = p0;
        seg.ctrl[1] = p0 + t0 * h0;
        seg.ctrl[2] = p1 - t1 * h1;
        seg.ctrl[3] = p1;
        isg::BezierCurve c;
        c.segs.push_back(seg);
        return c;
    };
    // 左偏的曲线转入率为正，直线为零。
    REQUIRE(isg::singleCubicSignedEndCurvature(make(6.0, 6.0), true) > 0.0);
    const isg::Vec2d chord = p1 - p0;
    isg::BezierSegment line;
    line.ctrl[0] = p0;
    line.ctrl[1] = p0 + 0.25 * chord;
    line.ctrl[2] = p0 + 0.75 * chord;
    line.ctrl[3] = p1;
    isg::BezierCurve straight;
    straight.segs.push_back(line);
    REQUIRE(std::abs(isg::singleCubicSignedEndCurvature(straight, true)) < 1e-12);

    // 闭式解 kappa0 = (2/3)*(cross(T0,chord) - h1*cross(T0,T1)) / h0^2。
    const double h0 = 5.0, h1 = 7.0;
    const double expected = (2.0 / 3.0) *
        (isg::cross2d(t0, chord) - h1 * isg::cross2d(t0, t1)) / (h0 * h0);
    REQUIRE(std::abs(isg::singleCubicSignedEndCurvature(make(h0, h1), true) -
                     expected) < 1e-9);

    // 尾把手相同时，共享侧把手 h0 单调即可定序（h0 越短转入越急）。
    const double ka = isg::singleCubicSignedEndCurvature(make(5.0, 2.0), true);
    const double kb = isg::singleCubicSignedEndCurvature(make(6.0, 2.0), true);
    REQUIRE(ka > kb);
    // 但尾把手独立选取即可把这一次序翻转：h0 仍是 5 < 6，转入率却反过来。
    // 这正是逐对比较共享侧把手投影无法表达扇出总序的原因。
    const double k_short_h0 = isg::singleCubicSignedEndCurvature(
        make(5.0, 6.0), true);
    const double k_long_h0 = isg::singleCubicSignedEndCurvature(
        make(6.0, 0.5), true);
    REQUIRE(k_short_h0 < k_long_h0);

    // 尾端点侧的语义与首端点一致：值越大越偏向该端点切向的左侧。
    const isg::BezierCurve c = make(6.0, 6.0);
    REQUIRE(isg::singleCubicSignedEndCurvature(c, false) < 0.0);
    REQUIRE(std::abs(isg::singleCubicSignedEndCurvature(isg::BezierCurve(), true))
            < 1e-12);
}

// ── 围栏（粗糙路口面）越界判定 ────────────────────────────────────────────────
// 粗糙路口面的外环是从各进出口车道端头连出来的，所以外环恰好穿过每条曲线的连接点。
// 由此产生两类不可归因于生成器的"越界"，必须与真实外溢分开（见 fence_check.h）。
namespace {

isg::Polygon2d makeSquareFence() {
    isg::Polygon2d fence;
    fence.outer = {isg::Vec2d(0, 0), isg::Vec2d(100, 0), isg::Vec2d(100, 100),
                   isg::Vec2d(0, 100)};
    return fence;
}

// 底边中段被切掉一块（y=0 抬到 y=5），用来构造"面本身不含两连接点之间直线弦"的情形，
// 即粗糙面把转角切掉的等价几何。
isg::Polygon2d makeNotchedFence() {
    isg::Polygon2d fence;
    fence.outer = {isg::Vec2d(0, 0),   isg::Vec2d(30, 0),   isg::Vec2d(30, 5),
                   isg::Vec2d(70, 5),  isg::Vec2d(70, 0),   isg::Vec2d(100, 0),
                   isg::Vec2d(100, 100), isg::Vec2d(0, 100)};
    return fence;
}

isg::BezierCurve makeCubic(const isg::Vec2d& p0, const isg::Vec2d& c1,
                           const isg::Vec2d& c2, const isg::Vec2d& p1) {
    isg::BezierCurve curve;
    isg::BezierSegment seg;
    seg.ctrl = {p0, c1, c2, p1};
    curve.segs.push_back(seg);
    return curve;
}

}  // namespace

TEST_CASE("fence exemption separates connection-point rounding from real overflow",
          "[fence][constraint]") {
    const isg::Polygon2d fence = makeSquareFence();
    const isg::Vec2d p0(10, -0.0005);   // 连接点本身，毫米级落在外环外侧
    const isg::Vec2d p1(90, 3.0);

    // 一、连接点本身：生成器自由度为零，允许量到 kFenceEndpointOutsideTol。
    REQUIRE(isg::fenceOutsideIsConnectionPointRounding(fence, p0, p0, p1));
    // 但"车道端点被数据放在面外很远"仍要报出来。
    const isg::Vec2d far_ep(10, -0.20);
    REQUIRE_FALSE(isg::fenceOutsideIsConnectionPointRounding(fence, far_ep, far_ep, p1));

    // 二、连接点近侧贴外环的一小段：允许量收紧到 kFenceConnectionRoundingTol。
    //     实测 sps 加密后会在此处多打出 0.23mm 一点（intersection_cross 61）。
    REQUIRE(isg::fenceOutsideIsConnectionPointRounding(
        fence, isg::Vec2d(10.04, -0.00023), p0, p1));
    // 同一位置但深 10mm，已超出坐标舍入尺度，不豁免。
    REQUIRE_FALSE(isg::fenceOutsideIsConnectionPointRounding(
        fence, isg::Vec2d(10.04, -0.010), p0, p1));
    // 跨度之外的真实中段外溢不豁免，哪怕只有 3.7mm（intersection_cross 65 实测值）。
    REQUIRE_FALSE(isg::fenceOutsideIsConnectionPointRounding(
        fence, isg::Vec2d(11.9, -0.0037), p0, p1));
}

TEST_CASE("curveInsideFence tolerates endpoint rounding but not interior overflow",
          "[fence][constraint]") {
    const isg::Polygon2d fence = makeSquareFence();
    // 首尾连接点各外溢 0.5mm，中段整体在面内：判为面内。
    const isg::BezierCurve clean = makeCubic(
        isg::Vec2d(10, -0.0005), isg::Vec2d(30, 10), isg::Vec2d(70, 10),
        isg::Vec2d(90, -0.0005));
    CHECK(isg::curveInsideFence(clean, fence, 25));
    CHECK(isg::curveInsideFence(clean, fence, 64));   // 审计用密度
    CHECK(isg::curveInsideFence(clean, fence, 160));  // 加密也不应翻结论
    CHECK(isg::curveFenceOverflow(clean, fence, 160) == 0.0);

    // 同样的端点，但中段压到外环之外：必须报出。
    const isg::BezierCurve dipping = makeCubic(
        isg::Vec2d(10, -0.0005), isg::Vec2d(30, -0.02), isg::Vec2d(70, -0.02),
        isg::Vec2d(90, -0.0005));
    CHECK_FALSE(isg::curveInsideFence(dipping, fence, 64));
    CHECK(isg::curveFenceOverflow(dipping, fence, 160) > 0.005);
}

TEST_CASE("fence overflow forced by the face is separated from shape defects",
          "[fence][constraint]") {
    const isg::Polygon2d fence = makeNotchedFence();
    const isg::Vec2d p0(20, 1), p1(80, 1);
    // 面本身不含两连接点之间的直线弦：缺口底面在 y=5，弦却走 y=1。
    const double chord = isg::fenceChordOverflow(fence, p0, p1, 160);
    REQUIRE(chord > 3.9);

    // 绕到缺口上方：完全在面内，不产生任何外溢。
    const isg::BezierCurve over = makeCubic(p0, isg::Vec2d(30, 20),
                                            isg::Vec2d(70, 20), p1);
    CHECK(isg::curveInsideFence(over, fence, 64));

    // 贴着弦走：外溢不超过弦强制的下限，归因于面画小了。
    const isg::BezierCurve along = makeCubic(p0, isg::Vec2d(35, 1),
                                             isg::Vec2d(65, 1), p1);
    const double along_out = isg::curveFenceOverflow(along, fence, 160);
    CHECK_FALSE(isg::curveInsideFence(along, fence, 64));
    CHECK(isg::fenceOverflowForcedByFace(along_out, chord));

    // 比直连还差：越过底边继续向外，超出弦强制下限，属形态缺陷。
    const isg::BezierCurve below = makeCubic(p0, isg::Vec2d(35, -6),
                                             isg::Vec2d(65, -6), p1);
    const double below_out = isg::curveFenceOverflow(below, fence, 160);
    CHECK(below_out > chord + 1.0);
    CHECK_FALSE(isg::fenceOverflowForcedByFace(below_out, chord));
}
