#pragma once

#include "types.h"
#include "utils.h"

#include <algorithm>
#include <limits>
#include <string>

namespace isg {

// ─────────────────────────────────────────────────────────────────────────────
// 安全岛辅助
// ─────────────────────────────────────────────────────────────────────────────

/// 安全岛边界（各边段）到指定点的最短距离（Q4：边界到端点，非顶点到端点）
inline double distIslandBoundaryToPoint(const SafetyIsland& island, const Vec2d& pt) {
    double min_dist = std::numeric_limits<double>::infinity();
    const auto& outer = island.geometry.outer;
    for (int i = 0; i + 1 < (int)outer.size(); ++i) {
        double d = pointToSegment(pt, outer[i], outer[i + 1]).first;
        min_dist = std::min(min_dist, d);
    }
    // 处理最后一条闭合边
    if (outer.size() >= 2) {
        double d = pointToSegment(pt, outer.back(), outer.front()).first;
        min_dist = std::min(min_dist, d);
    }
    return min_dist;
}

/// 查找与指定车道端点最近的安全岛（Q6：距进入车道端点最近原则）
/// 返回 nullptr 表示无安全岛
inline const SafetyIsland* findNearestSafetyIsland(
    const Vec2d& lane_endpoint,
    const std::vector<SafetyIsland>& safety_islands)
{
    const SafetyIsland* nearest = nullptr;
    double min_dist = std::numeric_limits<double>::infinity();
    for (const auto& si : safety_islands) {
        double d = distIslandBoundaryToPoint(si, lane_endpoint);
        if (d < min_dist) {
            min_dist = d;
            nearest = &si;
        }
    }
    return nearest;
}

// ─────────────────────────────────────────────────────────────────────────────
// 人行横道辅助
// ─────────────────────────────────────────────────────────────────────────────

/// 按 crosswalk_id 查找人行横道
inline const Crosswalk* findCrosswalk(
    const std::string& cw_id, const IntersectionInput& input)
{
    for (const auto& cw : input.crosswalks) {
        if (cw.id == cw_id) return &cw;
    }
    return nullptr;
}

// ─────────────────────────────────────────────────────────────────────────────
// 车道辅助
// ─────────────────────────────────────────────────────────────────────────────

/// 判断车道是否为非机动车道（motoway == false）
inline bool isNonMotorLane(const LaneId& lane_id, const IntersectionInput& input) {
    const Lane* lane = input.findLane(lane_id);
    return lane && !lane->motoway;
}

/// 统计组内非机动车道数量
inline int countNonMotorLanesInGroup(
    const LaneGroupId& gid, const IntersectionInput& input)
{
    const LaneGroup* grp = input.findGroup(gid);
    if (!grp) return 0;
    int cnt = 0;
    for (const auto& lid : grp->lanes) {
        const Lane* lane = input.findLane(lid);
        if (lane && !lane->motoway) ++cnt;
    }
    return cnt;
}

/// 判断连通关系是否为提右车道（CON_NUM 属性存在于进入车道）
inline bool isYieldRightLane(const Connectivity& conn, const IntersectionInput& input) {
    const Lane* entry_lane = input.findLane(conn.entry_lane_id);
    if (!entry_lane) return false;
    auto it = entry_lane->attrs.find("CON_NUM");
    return it != entry_lane->attrs.end();
}

/// 统计路口内机动车道的往返总数（所有车道组的机动车道数之和）
inline int countMotorLanesTotalIntersection(const IntersectionInput& input) {
    int cnt = 0;
    for (const auto& lane : input.lanes) {
        if (lane.motoway) ++cnt;
    }
    return cnt;
}

/// 获取组内 laneOrder 最大的车道（最外侧/最右侧车道）
inline const Lane* getRightmostLane(
    const LaneGroupId& gid, const IntersectionInput& input)
{
    const LaneGroup* grp = input.findGroup(gid);
    if (!grp || grp->lanes.empty()) return nullptr;
    const Lane* result = nullptr;
    int max_order = -1;
    for (const auto& lid : grp->lanes) {
        const Lane* lane = input.findLane(lid);
        if (lane && lane->laneOrder > max_order) {
            max_order = lane->laneOrder;
            result = lane;
        }
    }
    return result;
}

// ─────────────────────────────────────────────────────────────────────────────
// "右侧"判定（Q1确认）：进入组与退出组最右侧车道端点连线的右侧
// ─────────────────────────────────────────────────────────────────────────────

/// 计算进入/退出组最右侧连线的右侧法向量
/// 右侧 = 连线方向（进入端点→退出端点）顺时针旋转 90°
inline Vec2d computeRightSideNormal(
    const LaneGroupId& entry_gid,
    const LaneGroupId& exit_gid,
    const IntersectionInput& input)
{
    const Lane* entry_right = getRightmostLane(entry_gid, input);
    const Lane* exit_right  = getRightmostLane(exit_gid,  input);
    if (!entry_right || !exit_right) return Vec2d(0, -1);

    Vec2d entry_pt = input.entryPtDir(entry_right->id).first;
    Vec2d exit_pt  = input.exitPtDir(exit_right->id).first;
    Vec2d line_dir = exit_pt - entry_pt;
    if (line_dir.norm() < 1e-10) return Vec2d(0, -1);
    line_dir.normalize();
    // 顺时针旋转 90° = (dy, -dx)
    return Vec2d(line_dir.y(), -line_dir.x());
}

/// 统计某侧（进入组或退出组）在"右侧连线右侧"的同侧组数量（不含自身所在组）
inline int countRightSideGroups(
    const Connectivity& conn,
    GroupRole role,
    const IntersectionInput& input)
{
    Vec2d right_normal = computeRightSideNormal(
        conn.enterGroupId, conn.exitGroupId, input);
    // 参考原点取进入组最右侧车道端点
    const Lane* entry_right = getRightmostLane(conn.enterGroupId, input);
    Vec2d origin = entry_right
        ? input.entryPtDir(entry_right->id).first
        : Vec2d(0, 0);

    const LaneGroupId& self_gid =
        (role == GroupRole::Entry) ? conn.enterGroupId : conn.exitGroupId;

    int count = 0;
    for (const auto& grp : input.lane_groups) {
        if (grp.role != role) continue;
        if (grp.id == self_gid) continue;
        // 取组内所有端点重心
        Vec2d centroid = Vec2d::Zero();
        int n = 0;
        for (const auto& lid : grp.lanes) {
            const Lane* lane = input.findLane(lid);
            if (!lane) continue;
            Vec2d pt = (role == GroupRole::Entry)
                ? input.entryPtDir(lid).first
                : input.exitPtDir(lid).first;
            centroid += pt;
            ++n;
        }
        if (n == 0) continue;
        centroid /= (double)n;
        if ((centroid - origin).dot(right_normal) > 0.0)
            ++count;
    }
    return count;
}

// ─────────────────────────────────────────────────────────────────────────────
// 两线段交点（无穷延长线）
// ─────────────────────────────────────────────────────────────────────────────

/// 计算两条直线（由点+方向定义）的交点。
/// 若近乎平行则返回两起点的中点作兜底。
inline Vec2d lineIntersection(
    const Vec2d& p1, const Vec2d& d1,
    const Vec2d& p2, const Vec2d& d2)
{
    double denom = cross2d(d1, d2);
    if (std::abs(denom) < 1e-10) {
        // 近乎平行，返回中点
        return 0.5 * (p1 + p2);
    }
    double t = cross2d(p2 - p1, d2) / denom;
    return p1 + t * d1;
}

// ─────────────────────────────────────────────────────────────────────────────
// 微偏移辅助（Q5：消除折线重叠）
// ─────────────────────────────────────────────────────────────────────────────

/// 检查两条 polyline 是否在中间段（非端点）存在重叠（近乎平行且距离 < tol）
inline bool polylinesHaveOverlap(
    const LineString2d& a,
    const LineString2d& b,
    double tol)
{
    if (a.points.size() < 2 || b.points.size() < 2) return false;
    const auto& ap = a.points;
    const auto& bp = b.points;
    // 简化检查：取各线段中点到另一折线的距离
    for (size_t i = 1; i < ap.size(); ++i) {
        Vec2d mid = 0.5 * (ap[i - 1].xy() + ap[i].xy());
        double d = pointToPolyline(mid, toVec2dArray(bp));
        if (d < tol) return true;
    }
    return false;
}

/// 对 polyline 的中间点（非首尾）施加向左侧法向的偏移（消重叠）
inline void applyOverlapOffsetToMiddle(
    LineString2d& path,
    double offset)
{
    auto& pts = path.points;
    if (pts.size() < 3) return;
    // 只偏移中间点，保持首尾端点不变
    for (size_t i = 1; i + 1 < pts.size(); ++i) {
        Vec2d d = pts[i + 1].xy() - pts[i - 1].xy();
        if (d.norm() < 1e-10) continue;
        d.normalize();
        Vec2d left_normal{-d.y(), d.x()};  // 左侧法向
        pts[i] = offsetXY(pts[i], left_normal * offset);
    }
}

}  // namespace isg
