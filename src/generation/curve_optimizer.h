#pragma once

#include "optimizer/penalty_cost.h"
#include "types.h"

#include <vector>

namespace isg {

class SDFField;

struct CurveOptimizationOptions {
    bool enforce_fence;
    bool constrain_single_cubic_axes;
    double obstacle_clearance;
    int outer_iterations;

    CurveOptimizationOptions();
};

/// PenaltyCost 与 LBFGS 的单曲线适配层；离散录取仍由调用方约束门禁负责。
class CurveOptimizer {
public:
    explicit CurveOptimizer(LBFGSSolver& solver) : solver_(solver) {}

    BezierCurve optimize(
        const BezierCurve& initial,
        const IntersectionInput& input,
        const SDFField& sdf,
        const std::vector<SiblingCurve>& siblings,
        const Vec2d& entry_tangent,
        const Vec2d& exit_tangent,
        const CurveOptimizationOptions& options =
            CurveOptimizationOptions()) const;

private:
    LBFGSSolver& solver_;
};

}  // 命名空间 isg
