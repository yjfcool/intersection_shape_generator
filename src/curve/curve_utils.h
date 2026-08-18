#pragma once
#include "bezier.h"

namespace isg {

struct SDFField;

/// 三点局部曲率（基于外接圆半径倒数）
double localCurvature(const Vec2d& a, const Vec2d& b, const Vec2d& c);

/// 三点外接圆圆心
Vec2d circumcenter(const Vec2d& a, const Vec2d& b, const Vec2d& c);

/// 弹性带平滑：将点列视为弹性带，在曲率超限处沿外接圆径向移动
/// 顶点以降低曲率，受SDF和围栏约束
std::vector<Vec2d> elasticBandSmooth(const std::vector<Vec2d>&, const SDFField&,
                                     const Polygon2d&, double kappa_max = 0.25, double move_step = 0.03,
                                     int max_iter = 50, double min_sdf = 0.1);

/// 自适应版本弹性带平滑：根据曲线长度与曲率变化自适应调整采样数
std::vector<Vec2d> elasticBandSmoothAdaptive(const std::vector<Vec2d>&, const SDFField&,
                                     const Polygon2d&, double kappa_max = 0.25, double move_step = 0.03,
                                     int max_iter = 50, double min_sdf = 0.1);

/// 抽稀点列至目标数量
std::vector<Vec2d> downsamplePoints(const std::vector<Vec2d>&, int target_count);

/// 由平滑后点列重建Bezier曲线
BezierCurve rebuildFromSmoothedPts(const std::vector<Vec2d>&, const Vec2d&, const Vec2d&);

/// 沿两条曲线弧长方向均匀采样中点
std::vector<Vec2d> midlineSampleByArcLength(const BezierCurve&, const BezierCurve&, int n = 20);

/// 点到直线的有符号距离
double signedDistToLine(const Vec2d& pt, const Vec2d& p0, const Vec2d& p1);

/// 两线段是否相交，out返回交点
bool segmentsIntersect(const Vec2d& a, const Vec2d& b,
                       const Vec2d& c, const Vec2d& d, Vec2d* out = nullptr);

/// 曲线沿采样点的最小SDF值
double minSDFAlongCurve(const BezierCurve&, const SDFField&, int sps = 20);

/// 交点到两条曲线所有端点的最短距离
double distToAllEndpoints(const Vec2d&, const BezierCurve&, const BezierCurve&);

/// 返回两条曲线的所有交点（容差tol）
std::vector<Vec2d> curveCrossings(const BezierCurve&, const BezierCurve&, double tol = 0.01);

/// 两条曲线包围盒是否重叠
bool bboxOverlap(const BezierCurve&, const BezierCurve&);

/// 业务相交判定：两条曲线是否存在非端点附近相交（ep=端点排除半径，米）
bool curvesIntersectBusiness(const BezierCurve&, const BezierCurve&, double ep = 0.01);

/// 曲线自身相交判定：排除首尾相接的情况
bool curveSelfIntersectsBusiness(const BezierCurve&, double ep = 0.01);

/// 将普通单段三次Bezier的两个内部控制点限制在端点切向轴的有效范围内。
/// U型调头和物理避让候选由调用方显式跳过本约束。
void constrainOrdinarySingleCubicControls(
    BezierCurve& curve, const Vec2d& p0, const Vec2d& start_tan,
    const Vec2d& p1, const Vec2d& end_tan,
    bool cap_at_direction_intersection = true);

/// 检查普通单段三次Bezier是否满足端点方向轴和有效范围约束。
bool ordinarySingleCubicControlsValid(
    const BezierCurve& curve, const Vec2d& p0, const Vec2d& start_tan,
    const Vec2d& p1, const Vec2d& end_tan, double tol = 1e-6,
    bool cap_at_direction_intersection = true);

}
