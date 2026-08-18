#pragma once
#include "bezier.h"

namespace isg {

struct SDFField;

/// 初始曲线构造（方案A：几何直接构造）
/// sibling_polys: 已生成同簇曲线的采样点列，Level-1用于避免选择与对向拱交叉的绕行侧
BezierCurve buildInitialCurve(
    const Vec2d&, const Vec2d&, const Vec2d&, const Vec2d&,
    const SDFField&, const Polygon2d&, const std::vector<std::vector<Vec2d>>& sibling_polys = {});

/// 两段式U型调头构造
BezierCurve buildTwoSegmentUTurn(
    const Vec2d&, const Vec2d&, const Vec2d&, const Vec2d&, const SDFField&, const Polygon2d&);

}
