#pragma once

#include "types.h"
#include "constraints/road_edge_clearance.h"

namespace isg {

/// 统一的曲线采样，供约束审计和诊断使用。
std::vector<Vec2d> sampleCurveForAudit(const BezierCurve& curve, int samples = 64);

/// 曲线是否存在排除端点后的自交。
bool curveSelfIntersectsForAudit(const BezierCurve& curve, double endpoint_tol = 0.01);

/// 曲线是否穿过多边形（包含内部点和边界交叉）。
bool curveIntersectsPolygonForAudit(const BezierCurve& curve, const Polygon2d& polygon,
                                    int samples = 64, Vec2d* location = nullptr);

/// 曲线是否穿过任一障碍物；硬碰撞优先使用原始 geometry，与旧生成器语义一致。
bool curveIntersectsObstaclesForAudit(const BezierCurve& curve,
                                      const std::vector<Obstacle>& obstacles,
                                      Vec2d* location = nullptr);

/// 两条曲线是否存在非端点的平行贴合/重叠。
bool curvesHaveForbiddenAdherenceForAudit(const BezierCurve& a,
                                          const BezierCurve& b,
                                          double endpoint_tol = 0.15,
                                          double adherence_tol = 0.18);

/// 曲线非端点采样到指定类型边界折线的最小距离。
double minimumCurveBoundaryDistanceForAudit(const BezierCurve& curve,
                                            const std::vector<Boundary>& boundaries,
                                            Boundary::Type type,
                                            int samples = 64,
                                            double endpoint_exclusion = 1e-4,
                                            Vec2d* location = nullptr);

/// 测量曲线对指定类型边界的净距，并把"端点邻域"与"内部"分开统计，供
/// `constraints/road_edge_clearance.h` 的统一口径判定使用。clearance 参与端点邻域
/// 判定下限的计算，故必须与调用方最终使用的净距要求一致。
RoadEdgeClearanceMeasure measureCurveRoadEdgeClearanceForAudit(
    const BezierCurve& curve, const std::vector<Boundary>& boundaries,
    Boundary::Type type, double clearance, int samples = 64,
    double endpoint_exclusion = 1e-4);

}  // 命名空间 isg
