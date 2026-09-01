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
        ConstraintResult result;
        result.id = "fence.containment";
        result.severity = ConstraintSeverity::Hard;
        result.state = ConstraintState::Violated;
        result.reason = "curve leaves coarse intersection fence";
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
        // 与生成侧同口径：剔除端点贴合 RoadEdge 的共线擦碰（曲线端点与折线端点
        // 重合导致的厘米级数值穿越），否则生成侧接受的单段 cubic 会在审计侧
        // 重新报 physical.boundary。发夹形鼻端折回的另一条腿仍参与判违。
        const BoundarySafetyResult safety = curveBoundarySafetyIgnoringEndpointGraze(
            curve, context.view.input().boundaries, boundarySafetyCenter(context.view.input()),
            std::max(32, profile.samples * 2), 0.75, 0.10, 0.05);
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
            Vec2d location(0, 0);
            const double distance = minimumCurveBoundaryDistanceForAudit(
                curve, context.view.input().boundaries, Boundary::Type::RoadEdge,
                std::max(32, profile.samples * 2), 0.75, &location);
            if (distance < profile.road_edge_clearance) {
                ConstraintResult result;
                result.id = "physical.road_edge_clearance";
                result.severity = ConstraintSeverity::Hard;
                result.state = ConstraintState::Violated;
                result.violation = profile.road_edge_clearance - distance;
                result.locations.push_back(location);
                result.reason = "curve violates non-endpoint road-edge clearance";
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
