#include "constraints/constraint_evaluator.h"
#include "constraints/boundary_safety.h"
#include "constraints/cluster_constraint.h"
#include "constraints/fence_check.h"
#include "constraints/shape_constraint.h"
#include "geometry/predicates.h"
#include "optimizer/sdf_field.h"
#include "curve/curve_utils.h"
#include "utils.h"

#include <cmath>

namespace isg {

ConstraintReport ConstraintEvaluator::evaluate(const BezierCurve& curve,
                                               const SceneContext& context,
                                               const ConstraintProfile& profile) const {
    ConstraintReport report;
    if (curve.empty()) {
        ConstraintResult result;
        result.id = "geometry.empty";
        result.severity = ConstraintSeverity::Hard;
        result.state = ConstraintState::Violated;
        result.reason = "curve has no segments";
        report.results.push_back(result);
        return report;
    }
    if (profile.require_single_segment && curve.numSegments() != 1) {
        ConstraintResult result;
        result.id = "shape.single_segment";
        result.severity = ConstraintSeverity::Hard;
        result.state = ConstraintState::Violated;
        result.violation = static_cast<double>(curve.numSegments() - 1);
        result.reason = "ordinary curve must remain a single cubic segment";
        report.results.push_back(result);
    }
    if (profile.enforce_fence && !context.view.input().area.geometry.outer.empty() &&
        !curveInsideFence(curve, context.view.input().area.geometry, profile.samples)) {
        // 越界分两类：面本身不含"两个固定连接点之间的直线弦"时，任何曲线都出不来，
        // 外溢是输入面画小了（右转被切角尤甚），记 Exempt 并带上定量证据；反之才是
        // 生成器把曲线画到了比直连还差的位置，记 Violated。判据见 fence_check.h。
        const Polygon2d& fence = context.view.input().area.geometry;
        const int n = std::max(32, profile.samples * 2);
        const double curve_overflow = curveFenceOverflow(curve, fence, n);
        const double chord_overflow =
            fenceChordOverflow(fence, curve.startPt(), curve.endPt(), n);
        const bool forced = fenceOverflowForcedByFace(curve_overflow, chord_overflow);
        ConstraintResult result;
        result.id = "fence.containment";
        result.severity = ConstraintSeverity::Hard;
        result.state = forced ? ConstraintState::Exempt : ConstraintState::Violated;
        result.violation = forced ? 0.0 : curve_overflow - chord_overflow;
        result.reason = forced
            ? "curve leaves coarse intersection fence, but the fence itself excludes "
              "the straight chord between the two fixed connection points"
            : "curve leaves coarse intersection fence beyond the chord-forced floor";
        report.results.push_back(result);
    }

    if (profile.check_self_intersection && curveSelfIntersectsForAudit(curve)) {
        ConstraintResult result;
        result.id = "geometry.self_intersection";
        result.severity = ConstraintSeverity::Hard;
        result.state = ConstraintState::Violated;
        result.reason = "curve self-intersects away from endpoints";
        report.results.push_back(result);
    }

    if (profile.check_obstacle && !context.view.input().obstacles.empty()) {
        Vec2d location(0, 0);
        bool hit = curveIntersectsObstaclesForAudit(curve, context.view.input().obstacles, &location);
        double min_sdf = 1e18;
        if (context.sdf && context.sdf->valid()) {
            const std::vector<Vec2d> points = sampleCurveForAudit(curve, profile.samples);
            for (const auto& point : points) {
                const double distance = context.sdf->queryWithGrad(point).first;
                if (distance < min_sdf) {
                    min_sdf = distance;
                    location = point;
                }
            }
            hit = hit || min_sdf < profile.obstacle_clearance;
        }
        if (hit) {
            ConstraintResult result;
            result.id = "physical.obstacle";
            result.severity = ConstraintSeverity::Hard;
            result.state = ConstraintState::Violated;
            result.violation = min_sdf < 1e17 ? std::max(0.0, profile.obstacle_clearance - min_sdf) : 1.0;
            result.locations.push_back(location);
            result.reason = "curve intersects obstacle or violates obstacle clearance";
            report.results.push_back(result);
        }
    }

    if (profile.check_boundary && !context.view.input().boundaries.empty()) {
        // 与生成侧同口径：只允许曲线真实首/尾连接点的浮点误差；Boundary
        // 端点及端点后的共线贴合、重叠和穿越都参与判违。
        const BoundarySafetyResult safety = curveBoundarySafetyForInput(
            curve, context.view.input(), std::max(32, profile.samples * 2),
            0.75, 0.10, 0.05);
        if (safety.intersects || safety.outside_road_edge) {
            ConstraintResult result;
            result.id = "physical.boundary";
            result.severity = ConstraintSeverity::Hard;
            result.state = ConstraintState::Violated;
            result.violation = safety.outside_penalty;
            result.reason = safety.outside_road_edge ? "curve leaves road-edge center side" :
                                                       "curve intersects boundary away from endpoints";
            report.results.push_back(result);
        }
        if (profile.road_edge_clearance > 0.0) {
            // 连接点位置由输入固定、端点切向被 G1 锁定在车道朝向上：连接点自身就落在净距
            // 之内时（corpus 里 110000741-u 出口连接点距 RoadEdge 仅 0.9521m，要求 1.0m），
            // 端点邻域的亏欠不可能由曲线形状消除。逐点判据见 road_edge_clearance.h。
            const RoadEdgeClearanceMeasure measure =
                measureCurveRoadEdgeClearanceForAudit(
                    curve, context.view.input().boundaries, Boundary::Type::RoadEdge,
                    profile.road_edge_clearance, std::max(32, profile.samples * 2),
                    kConnectionPointTolerance);
            const double deficit = roadEdgeClearanceDeficit(
                measure, profile.road_edge_clearance);
            const bool forced = roadEdgeClearanceForcedByEndpoint(
                measure, profile.road_edge_clearance);
            if (deficit > kRoadEdgeClearanceRoundingTol || forced) {
                ConstraintResult result;
                result.id = "physical.road_edge_clearance";
                result.severity = ConstraintSeverity::Hard;
                result.state = forced ? ConstraintState::Exempt
                                      : ConstraintState::Violated;
                result.violation = forced ? 0.0 : deficit;
                result.locations.push_back(
                    forced ? measure.location : measure.deficit_location);
                result.reason = forced
                    ? "road-edge clearance deficit forced by connection point"
                    : "curve violates non-endpoint road-edge clearance";
                report.results.push_back(result);
            }
        }
    }

    if (profile.check_curvature && curve.maxCurvature(profile.samples) > profile.max_curvature) {
        ConstraintResult result;
        result.id = "shape.curvature";
        result.severity = ConstraintSeverity::Hard;
        result.state = ConstraintState::Violated;
        result.violation = curve.maxCurvature(profile.samples) - profile.max_curvature;
        result.reason = "curve exceeds maximum curvature";
        report.results.push_back(result);
    }

    if (profile.check_g1 && curve.numSegments() > 1) {
        const double angle_limit = profile.g1_angle_deg * DEG2RAD;
        for (int i = 0; i + 1 < curve.numSegments(); ++i) {
            const Vec2d a = curve.segs[i].evalDeriv1(1.0);
            const Vec2d b = curve.segs[i + 1].evalDeriv1(0.0);
            if (a.norm() < 1e-9 || b.norm() < 1e-9 || angleBetween(a, b) > angle_limit) {
                ConstraintResult result;
                result.id = "shape.g1";
                result.severity = ConstraintSeverity::Hard;
                result.state = ConstraintState::Violated;
                result.reason = "adjacent cubic segments are not G1 continuous";
                report.results.push_back(result);
                break;
            }
        }
    }

    if (report.results.empty()) {
        ConstraintResult result;
        result.id = "common.audit";
        result.severity = ConstraintSeverity::Advisory;
        result.state = ConstraintState::Satisfied;
        report.results.push_back(result);
    }
    return report;
}

ConstraintReport ConstraintEvaluator::evaluate(const BezierCurve& curve,
                                               const CurveGenerationContext& context,
                                               const GenerationState& state) const {
    ConstraintReport report;
    if (context.scene == nullptr) {
        ConstraintResult result;
        result.id = "context.missing_scene";
        result.severity = ConstraintSeverity::Hard;
        result.state = ConstraintState::Violated;
        result.reason = "curve generation context has no scene";
        report.results.push_back(result);
        return report;
    }
    report = evaluate(curve, *context.scene, context.profile);
    if (context.profile.check_ordinary_shape)
        report.results.push_back(evaluateOrdinaryShape(curve, context));
    if (context.profile.check_uturn_shape)
        report.results.push_back(evaluateUTurnShape(curve, context));
    if (context.profile.check_crosswalk)
        report.results.push_back(evaluateUTurnCrosswalk(curve, context));
    if (context.profile.check_cluster) {
        const std::vector<ConstraintResult> cluster =
            evaluateClusterIntersections(curve, context, state);
        report.results.insert(report.results.end(), cluster.begin(), cluster.end());
    }
    return report;
}

ConstraintReport ConstraintEvaluator::evaluate(const CurveCandidate& candidate,
                                               const CurveGenerationContext& context,
                                               const GenerationState& state) const {
    ConstraintReport report = evaluate(candidate.curve, context, state);
    if (context.profile.require_fixed_shape_preservation && context.connectivity &&
        context.connectivity->fixed_shape && !candidate.preserves_fixed_shape) {
        ConstraintResult result;
        result.id = "fixed_shape.preservation";
        result.severity = ConstraintSeverity::Hard;
        result.state = ConstraintState::Violated;
        result.reason = "fixed connectivity candidate does not preserve declared geometry";
        report.results.push_back(result);
    }
    return report;
}

void ConstraintEvaluator::audit(CurveCandidate& candidate,
                                const CurveGenerationContext& context,
                                const GenerationState& state) const {
    candidate.acceptReport(evaluate(candidate, context, state));
}

}  // 命名空间 isg
