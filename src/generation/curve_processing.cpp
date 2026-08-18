#include "generation/curve_processing.h"

#include "constraints/boundary_safety.h"
#include "curve/bezier.h"
#include "curve/curve_utils.h"
#include "optimizer/lbfgs_solver.h"
#include "optimizer/penalty_cost.h"
#include "utils.h"

#include <algorithm>

namespace isg {

// 连续优化：只组装既有代价项和端点切向约束，不改变求解器权重。
CurveOptimizationOptions::CurveOptimizationOptions()
    : enforce_fence(false),
      constrain_single_cubic_axes(true),
      obstacle_clearance(0.0),
      outer_iterations(3) {}

BezierCurve CurveOptimizer::optimize(
    const BezierCurve& initial, const IntersectionInput& input,
    const SDFField& sdf, const std::vector<SiblingCurve>& siblings,
    const Vec2d& entry_tangent, const Vec2d& exit_tangent,
    const CurveOptimizationOptions& options) const {
    PenaltyCost cost;
    cost.proto = initial;
    cost.sdf = &sdf;
    cost.boundaries = input.boundaries;
    cost.road_center = boundarySafetyCenter(input);
    if (options.enforce_fence && !input.area.is_rough)
        cost.fence = input.area.geometry;
    cost.siblings = siblings;
    cost.obstacle_clearance = options.obstacle_clearance;
    cost.road_edge_clearance = roadEdgeAvoidanceClearanceForMode(input.mode);
    cost.constrain_single_cubic_axes = options.constrain_single_cubic_axes;
    cost.start_tan_dir = entry_tangent.norm() > 1e-8
        ? entry_tangent.normalized() : Vec2d(1, 0);
    cost.end_tan_dir = exit_tangent.norm() > 1e-8
        ? exit_tangent.normalized() : Vec2d(1, 0);
    cost.full_param_mode = initial.numSegments() > 1;
    return optimiseCurve(cost, solver_, initial, options.outer_iterations);
}

// 后处理：自适应细分后恢复端点与 G1，并按配置决定是否执行弹性带平滑。
BezierCurve CurvePostProcessor::process(
    const BezierCurve& curve, const SDFField& sdf, const Polygon2d& fence,
    double max_curvature, const Vec2d& entry_tangent,
    const Vec2d& exit_tangent, bool skip_elastic_band,
    const Vec2d* exact_entry, const Vec2d* exact_exit) const {
    const AdaptiveRefineResult refined = adaptiveRefine(curve, sdf, max_curvature);
    BezierCurve current = refined.curve;
    if (refined.was_split) {
        PenaltyCost cost;
        cost.proto = current;
        cost.sdf = &sdf;
        cost.fence = fence;
        cost.start_tan_dir = entry_tangent.norm() > 1e-8
            ? entry_tangent.normalized() : curve.startTan();
        cost.end_tan_dir = exit_tangent.norm() > 1e-8
            ? exit_tangent.normalized() : curve.endTan();
        cost.obstacle_clearance = 0.0;
        cost.full_param_mode = current.numSegments() > 1;
        cost.buildCache();
        current = optimiseCurve(cost, solver_, current, 2);
    }

    const Vec2d start_tangent = entry_tangent.norm() > 1e-8
        ? entry_tangent.normalized() : current.startTan();
    const Vec2d end_tangent = exit_tangent.norm() > 1e-8
        ? exit_tangent.normalized() : current.endTan();
    const Vec2d entry = exact_entry ? *exact_entry : current.startPt();
    const Vec2d exit = exact_exit ? *exact_exit : current.endPt();

    if (skip_elastic_band) {
        if (!current.segs.empty()) {
            current.segs.front().ctrl[0] = entry;
            current.segs.back().ctrl[3] = exit;
            if (current.numSegments() == 1) {
                constrainOrdinarySingleCubicControls(
                    current, entry, start_tangent, exit, end_tangent, true);
                return current;
            }
            Vec2d& first_control = current.segs.front().ctrl[1];
            const double start_length =
                (current.segs.front().ctrl[3] - entry).norm();
            const double start_scale = std::max(
                (first_control - entry).dot(start_tangent), start_length * 0.1);
            first_control = entry + start_scale * start_tangent;
            Vec2d& last_control = current.segs.back().ctrl[2];
            const double end_length =
                (exit - current.segs.back().ctrl[0]).norm();
            const double end_scale = std::max(
                (exit - last_control).dot(end_tangent), end_length * 0.1);
            last_control = exit - end_scale * end_tangent;
        }
        return current;
    }

    const double arc_length = current.arcLength();
    const int sample_count = std::max(
        40, std::min(100, static_cast<int>(arc_length / 0.15)));
    std::vector<Vec2d> points = current.sampleByArcLength(sample_count);
    if (points.size() >= 3) {
        points.front() = entry;
        points.back() = exit;
        const double curve_max_curvature = current.maxCurvature(20);
        const double move_step = std::min(
            0.15, std::max(0.02, curve_max_curvature * 0.3));
        std::vector<Vec2d> smoothed = elasticBandSmooth(
            points, sdf, fence, max_curvature, move_step, 30, 0.1);
        smoothed.front() = entry;
        smoothed.back() = exit;
        current = rebuildFromSmoothedPts(smoothed, start_tangent, end_tangent);
        if (!current.segs.empty()) {
            current.segs.front().ctrl[0] = entry;
            current.segs.back().ctrl[3] = exit;
        }
    }
    return current;
}

bool isCurveSeverelyDivergent(
    const BezierCurve& curve, const Vec2d& entry, const Vec2d& exit,
    double chord_length) {
    if (curve.empty() || chord_length < 1e-6)
        return false;
    if (curveSelfIntersectsBusiness(curve, 1.0))
        return true;
    const double arc_length = curve.arcLength();
    if (arc_length > std::max(3.0 * chord_length + 20.0, chord_length * 5.0))
        return true;
    if (curve.maxCurvature(20) > 5.0)
        return true;

    const double far_threshold = 3.0 * chord_length + 10.0;
    const double far_threshold_squared = far_threshold * far_threshold;
    for (const auto& segment : curve.segs) {
        for (int control = 0; control < 4; ++control) {
            const Vec2d& point = segment.ctrl[control];
            if ((point - entry).squaredNorm() > far_threshold_squared &&
                (point - exit).squaredNorm() > far_threshold_squared) {
                return true;
            }
        }
    }
    return false;
}

// 结果回退：拒绝严重发散或破坏 U-turn 曲率上限的优化结果。
OptimizationResultOptions::OptimizationResultOptions()
    : geometric_uturn(false),
      skip_elastic_band(true),
      postprocess_max_curvature(0.25),
      uturn_fallback_curvature(3.0) {}

BezierCurve OptimizationResultProcessor::process(
    const BezierCurve& initial, const BezierCurve& optimized,
    const SDFField& sdf, const Polygon2d& fence,
    const Vec2d& entry, const Vec2d& entry_tangent,
    const Vec2d& exit, const Vec2d& exit_tangent,
    const OptimizationResultOptions& options) const {
    const double chord_length = (exit - entry).norm();
    BezierCurve accepted = isCurveSeverelyDivergent(
        optimized, entry, exit, chord_length) ? initial : optimized;
    BezierCurve result = options.geometric_uturn
        ? accepted
        : CurvePostProcessor(solver_).process(
              accepted, sdf, fence, options.postprocess_max_curvature,
              entry_tangent, exit_tangent, options.skip_elastic_band,
              &entry, &exit);
    if (isCurveSeverelyDivergent(result, entry, exit, chord_length))
        result = initial;
    if (options.geometric_uturn &&
        result.maxCurvature(40) > options.uturn_fallback_curvature &&
        initial.maxCurvature(40) <= options.uturn_fallback_curvature) {
        result = initial;
    }
    return result;
}

}  // 命名空间 isg
