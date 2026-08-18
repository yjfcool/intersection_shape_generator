#pragma once
#include "types.h"

namespace isg {

struct SDFField;

// ── BezierSegment 与 BezierCurve 结构体定义在 types.h 中 ──

/// 构造G1连续的三次Bezier段，端点切向由t0/t1给出，alpha 控制把手长度比例（0.33=标准，较小=更紧凑，较大=更平缓）
BezierSegment makeCubicG1(const Vec2d& p0, const Vec2d& t0, const Vec2d& p1, const Vec2d& t1, double alpha = 0.33);

/// 构造对齐U型调头的单段三次Bezier：
///   ctrl[0] = p0                              (进入端点)
///   ctrl[1] = p0 + (lead0 + arc_handle) * T0  (进入延长 + 拱弧切向把手)
///   ctrl[2] = p1 - (lead1 + arc_handle) * T1  (退出反向延长 + 拱弧切向把手)
///   ctrl[3] = p1                              (退出端点)
///
/// - lead0/lead1 使拱弧顶点(q0=p0+lead0*T0, q1=p1+lead1*(-T1))沿U型轴线对齐
/// - arc_handle 物理基础为 (2/3)*turn_gap (经典4/3*r半圆近似, r=turn_gap/2)
/// - min_lead0 强制最小进入延长(米)，用于"调头需跨越进入侧人行横道后开始拱弧"
/// - min_lead1 强制最小退出延长(米)，用于"调头拱弧需跨越退出侧人行横道后结束"
/// - lateral_bias(米,有符号) 沿U型横向方向对称平移两个内部控制点，
///   保持两端G1连续，用于扇开共享进入车道的多个U型调头
BezierSegment makeAlignedUTurnCubic(
    const Vec2d& p0, const Vec2d& t0, const Vec2d& p1, const Vec2d& t1,
    double handle_scale = 1.0, double handle_bias = 0.0, double lateral_bias = 0.0,
    double min_lead0 = 0.0, double min_lead1 = 0.0);

/// 由节点点列与切向序列构造拼接Bezier曲线
BezierCurve makeCurveFromKnots(const std::vector<Vec2d>& pts, const std::vector<Vec2d>& tans, double alpha = 0.33);

/// 用首尾切向对点列进行Bezier拟合
BezierCurve fitBezierWithEndTangents(const std::vector<Vec2d>& pts, const Vec2d& start_tan, const Vec2d& end_tan);

/// 自适应细分：当曲线曲率超限或穿透障碍物时细分曲线段
AdaptiveRefineResult adaptiveRefine(const BezierCurve& c, const SDFField& sdf, double kappa_max = 0.25, double pen_tol = 0.05);

/// 曲线 → 参数向量（仅含段间连接点）
VecXd curveToParams(const BezierCurve& c);

/// 参数向量 → 曲线（按原型恢复段结构）
BezierCurve curveFromParams(const VecXd& params, const BezierCurve& proto);

/// 曲线 → 参数向量（所有控制点含连接点均为优化变量）
VecXd curveToParamsFull(const BezierCurve& c);

/// 参数向量 → 曲线（所有控制点含连接点）
BezierCurve curveFromParamsFull(const VecXd& params, const BezierCurve& proto);

}
