#pragma once

#include "types.h"

namespace isg {

/// 物理避让候选的无状态几何构造器。
/// 场景分析、约束审计和候选评分由生成会话负责。
class AvoidanceCandidateGenerator {
public:
    BezierCurve buildWaypointCurve(const std::vector<Vec2d>& points,
                                   const Vec2d& start_tangent,
                                   const Vec2d& end_tangent) const;
};

}  // 命名空间 isg
