#pragma once

#include "types.h"

#include <limits>
#include <string>

namespace isg {

struct CrosswalkClearanceResult {
    bool found = false;
    std::string crosswalk_id;
    double near = std::numeric_limits<double>::infinity();
    double far = 0.0;
    double clearance = 0.0;
};

/// 沿给定射线选择最近 Crosswalk 并计算完整远边投影。
class CrosswalkClearanceCalculator {
public:
    CrosswalkClearanceResult alongRay(
            const Vec2d& origin, const Vec2d& direction, const IntersectionInput& input) const;

    CrosswalkClearanceResult ahead(
            const Vec2d& origin, const Vec2d& tangent, const IntersectionInput& input) const;

    CrosswalkClearanceResult behind(
            const Vec2d& endpoint, const Vec2d& tangent, const IntersectionInput& input) const;
};

}  // 命名空间 isg
