#pragma once

#include "types.h"

namespace isg {

/// 最终曲线的兼容审计快照；几何检测由统一约束谓词提供。
struct FinalCurveAuditSnapshot {
    double max_obstacle_penetration;
    double max_fence_overflow;
    bool update_fence_overflow;
    bool obstacle_intersection;
    bool self_intersection;
    bool boundary_intersection;
    bool uturn;
    bool road_edge_clearance_violation;
    bool ordinary_single_axis_violation;

    FinalCurveAuditSnapshot()
        : max_obstacle_penetration(0.0), max_fence_overflow(0.0),
          update_fence_overflow(false), obstacle_intersection(false),
          self_intersection(false), boundary_intersection(false),
          uturn(false), road_edge_clearance_violation(false),
          ordinary_single_axis_violation(false) {}
};

/// 将完整审计快照确定性映射到兼容 ViolationInfo 和 CurveStatus。
class FinalCurveAuditor {
public:
    void apply(const FinalCurveAuditSnapshot& snapshot,
               ConnectivityCurve& curve) const;
};

}  // 命名空间 isg
