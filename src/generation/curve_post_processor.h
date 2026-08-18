#pragma once

#include "types.h"

namespace isg {

class LBFGSSolver;
class SDFField;

/// 自适应细分、必要的二次优化、弹性带平滑和端点 G1 恢复。
class CurvePostProcessor {
public:
    explicit CurvePostProcessor(LBFGSSolver& solver) : solver_(solver) {}

    BezierCurve process(const BezierCurve& curve,
                        const SDFField& sdf,
                        const Polygon2d& fence,
                        double max_curvature,
                        const Vec2d& entry_tangent,
                        const Vec2d& exit_tangent,
                        bool skip_elastic_band,
                        const Vec2d* exact_entry = nullptr,
                        const Vec2d* exact_exit = nullptr) const;

private:
    LBFGSSolver& solver_;
};

}  // 命名空间 isg
