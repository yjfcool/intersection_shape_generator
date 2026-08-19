#pragma once

#include "types.h"

namespace isg {

/// 普通曲线的无状态单段候选构造器。
/// 物理、簇和形态约束由 ConstraintEvaluator/选择器负责。
class OrdinaryCurveInitializer {
public:
    /// 无避让普通曲线的首选基础表达：直行的两个中间控制点重合于
    /// 等弦中点 M；
    /// 左右转在前向方向交点稳定且不过远时使用交点距离的 2/3。
    BezierCurve buildPreferredSingleCubic(
        const Vec2d& entry_point, const Vec2d& entry_tangent,
        const Vec2d& exit_point, const Vec2d& exit_tangent) const;

    BezierCurve buildSingleCubic(const Vec2d& entry_point,
                                 const Vec2d& entry_tangent,
                                 const Vec2d& exit_point,
                                 const Vec2d& exit_tangent,
                                 double alpha = 0.4) const;

    std::vector<BezierCurve> buildAlphaCandidates(
        const Vec2d& entry_point, const Vec2d& entry_tangent,
        const Vec2d& exit_point, const Vec2d& exit_tangent,
        const std::vector<double>& alphas) const;
};

}  // 命名空间 isg
