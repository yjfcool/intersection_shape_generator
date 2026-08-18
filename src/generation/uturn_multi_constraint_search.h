#pragma once

#include "generation/curve_sampling.h"

#include <functional>
#include <limits>

namespace isg {

struct UTurnSolverResult {
    BezierCurve curve;
    double cost = std::numeric_limits<double>::infinity();
    int sibling_crosses = std::numeric_limits<int>::max();
    double obst_pen = 0.0;
    double boundary_pen = 0.0;
    double fence_overflow = 0.0;
    double lateral_bias_used = 0.0;
    double lead0_used = 0.0;
};

struct UTurnSearchBackend {
    std::function<SampledCurve(const BezierCurve&)> sample;
    std::function<int(const BezierCurve&)> sibling_cross_count;
    std::function<double(const BezierCurve&)> obstacle_penalty;
    std::function<double(const BezierCurve&)> boundary_penalty;
    std::function<double(const SampledCurve&)> endpoint_side_violation;
};

/// 固定网格 U-turn 多约束搜索；场景专用审计由后端回调提供。
class UTurnMultiConstraintSearch {
public:
    UTurnSolverResult search(
        const Vec2d& entry, const Vec2d& entry_tangent,
        const Vec2d& exit, const Vec2d& exit_tangent,
        const IntersectionInput& input,
        const std::vector<SampledSiblingCurve>& sampled_siblings,
        const BezierCurve& reference, const UTurnSearchBackend& backend,
        double min_lead0_floor, double min_lead1_floor,
        int depth_preference = 0, int lead1_depth_preference = 0,
        double lateral_preference = 0.0,
        const std::vector<Crosswalk>* crosswalks_for_clearance = nullptr) const;
};

}  // 命名空间 isg
