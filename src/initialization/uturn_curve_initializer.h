#pragma once

#include "types.h"

#include <limits>

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

    /// 三段式(直-弧-直)表达。
    /// entry_lateral_bias / exit_lateral_bias 是有符号的横向偏置(米)，沿 U 轴
    /// 法线把首/尾平齐点 q0/q1 各自推离轴线。它与 aligned_*_stagger 的相向内缩
    /// 不同：内缩用于打散同家族调头的直行段重叠，偏置用于把某一侧直行段推离
    /// 共享出入车道中心线，从而与同簇的直行/转向曲线彻底分侧。
    /// 偏置沿轴法线，因此不改变 q0/q1 的轴向站位；lead 会按偏置反解补偿，
    /// 保证最小直行段的弦长读数不被偏置吃掉。代价是首尾直段方向相对车道切向
    /// 出现 atan(bias/lead) 的偏差，调用方需自行限制偏置幅度以满足 G1 门禁。
    ///
    /// aligned_entry_stagger / aligned_exit_stagger 是入口侧、出口侧各自的分档量，
    /// 两侧独立生效、不合并成一个共享标量：掉头家族是一条交替由"共享入口端点"和
    /// "共享出口端点"连接起来的链，而内缩在入口侧沿 +side*lateral、出口侧沿
    /// -side*lateral 施加，合并后压小一侧会连带压小另一侧，成员相对共享另一端点
    /// 的邻居随即反号。两侧的分档由 UTurnFamilyBuilder::familyLateralLadder 分别
    /// 裁定。允许为负：表示首尾平齐点向远离出口的一侧拓宽（相背而不是相向），
    /// 几何上同样合法，且不受走廊宽度上限约束，只受 G1 余量约束。
    /// 未指定时传 NaN（默认值），此时回退到 aligned_point_stagger。
    ///
    /// family_stagger_step 是同入/同出掉头家族的分档步长(米)，即相邻半径名次
    /// 的横向站位之差，也就是家族的横向"槽位"宽度。传入非零值后，与分档方向
    /// 相反的偏置会被钳制在 0.4 个槽位内（仅对该侧分档确实完整生效时钳制），
    /// 使家族左右次序与生成顺序无关地保持成立；传 0 表示调用方不掌握家族分档，
    /// 不做钳制。
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
                               double aligned_entry_stagger =
                                   std::numeric_limits<double>::quiet_NaN(),
                               double aligned_exit_stagger =
                                   std::numeric_limits<double>::quiet_NaN(),
                               double entry_lateral_bias = 0.0,
                               double exit_lateral_bias = 0.0,
                               double family_stagger_step = 0.0) const;
};

}  // 命名空间 isg
