/// InfeasibilityDetector —— 无解预检测
///   - BFS可达性 (C类: 拓扑死锁)
///   - 通道宽度检测 (A类: 过窄)
///   - 围栏-障碍夹断检测 (B类)
#pragma once
#include "types.h"

namespace isg {

/// 预检测结果
struct PreCheckResult {
    bool topological_block = false;  ///< 拓扑死锁
    bool narrow_passage = false;     ///< 通道过窄
    bool fence_sandwich = false;     ///< 围栏-障碍夹断
    double min_corridor_width = 0.0;
    ViolationInfo::InfeasibilityType type = ViolationInfo::InfeasibilityType::None;
};

/// 综合预检测：通道宽度+围栏+拓扑可达性
PreCheckResult preCheck(const SDFField& sdf, const Polygon2d& fence,
                        const Vec2d& entry_pt, const Vec2d& exit_pt,
                        double required_width, const std::vector<Boundary>& boundaries);

/// BFS可达性：从entry_pt到exit_pt在SDF+围栏限定的可行域内是否可达
bool bfsReachable(const SDFField&, const Polygon2d&, const Vec2d&, const Vec2d&, double mc = 0.3);

/// 检测围栏与障碍物是否形成夹断（sandwich）
bool detectSandwich(const SDFField&, const Polygon2d&, double cl = 0.5);

/// 围栏放松结果
struct FenceRelaxResult {
    bool success = false;
    Polygon2d relaxed_fence;
    double max_expansion = 0;

    FenceRelaxResult(bool success, const Polygon2d& relaxed_fence, double max_expansion)
        : success(success), relaxed_fence(relaxed_fence) {}
};

/// 尝试放松围栏以避开夹断
FenceRelaxResult tryRelaxFence(const Polygon2d&, const std::vector<Boundary>&, const Vec2d&, double mx = 1.5);

/// 构造回退曲线（当主流程无解时使用）
ConnectivityCurve makeFallbackCurve(const PreCheckResult&, const Connectivity&, const Vec2d&, const Vec2d&);

/// 多边形是否简单（无自相交）
bool isSimplePolygon(const Polygon2d&);

/// 两线段是否相交（与utils.h中segmentsIntersect独立实现，用于无依赖场景）
bool segmentsIntersect2(const Vec2d&, const Vec2d&, const Vec2d&, const Vec2d&);

/// 重载：默认空边界列表
inline PreCheckResult preCheck(
    const SDFField& sdf, const Polygon2d& fence, const Vec2d& ep, const Vec2d& xp, double rw) {
    return preCheck(sdf, fence, ep, xp, rw, {});
}

}
