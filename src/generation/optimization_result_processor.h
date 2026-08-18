#pragma once

#include "types.h"

namespace isg {

class LBFGSSolver;
class SDFField;

bool isCurveSeverelyDivergent(
    const BezierCurve& curve, const Vec2d& entry, const Vec2d& exit,
    double chord_length);

struct OptimizationResultOptions {
    bool geometric_uturn;
    bool skip_elastic_band;
    double postprocess_max_curvature;
    double uturn_fallback_curvature;

    OptimizationResultOptions();
};

/// 对连续优化结果执行通用后处理和确定性回退，不负责业务候选搜索。
class OptimizationResultProcessor {
public:
    explicit OptimizationResultProcessor(LBFGSSolver& solver)
        : solver_(solver) {}

    BezierCurve process(
        const BezierCurve& initial,
        const BezierCurve& optimized,
        const SDFField& sdf,
        const Polygon2d& fence,
        const Vec2d& entry,
        const Vec2d& entry_tangent,
        const Vec2d& exit,
        const Vec2d& exit_tangent,
        const OptimizationResultOptions& options =
            OptimizationResultOptions()) const;

private:
    LBFGSSolver& solver_;
};

}  // 命名空间 isg
