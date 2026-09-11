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

/// 判断射线命中的横道是否属于当前 U 型掉头走廊。
///
/// mode=2 的几何数据中，单条射线可能擦到相邻道路横道。端点弦距离是
/// 一个必要的快速判据，但对长掉头走廊会误杀真正横穿入口/出口两条车道的
/// 同一横道；当两侧射线命中同一横道且近端投影对齐时，保留该横道。
bool crosswalkRelevantToUTurn(
        const Crosswalk& crosswalk,
        const Vec2d& turn_entry, const Vec2d& turn_exit,
        const CrosswalkClearanceResult& candidate,
        const CrosswalkClearanceResult* paired_candidate,
        double chord_distance = 4.0,
        double paired_near_delta = 1.0);

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
