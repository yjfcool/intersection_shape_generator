#pragma once

#include "types.h"

namespace isg {

/// U-turn 的纯几何初始化器。
/// 只负责按给定轴向、首尾直行长度和错开量表达曲线，不执行约束筛选。
class UTurnCurveInitializer {
public:
    BezierCurve buildAligned(const Vec2d& entry_point,
                             const Vec2d& entry_tangent,
                             const Vec2d& exit_point,
                             const Vec2d& exit_tangent,
                             const Vec2d& offset_direction,
                             double offset_m,
                             double handle_scale = 1.0,
                             double min_lead0 = 0.0,
                             double min_lead1 = 0.0) const;

    BezierCurve buildSegmented(const Vec2d& entry_point,
                               const Vec2d& entry_tangent,
                               const Vec2d& exit_point,
                               const Vec2d& exit_tangent,
                               double min_lead0 = 0.0,
                               double min_lead1 = 0.0,
                               double arc_alpha = 2.0 / 3.0,
                               double aligned_point_stagger = 0.0,
                               double lead0_extra_after_align = 0.0,
                               double lead1_extra_after_align = 0.0,
                               double aligned_entry_stagger = -1.0,
                               double aligned_exit_stagger = -1.0) const;
};

}  // 命名空间 isg
