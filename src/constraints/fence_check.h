#pragma once
#include "types.h"
#include "curve/curve_utils.h"

namespace isg {

/// 点是否在多边形内（射线法，支持洞）
static bool polygonContains(const Polygon2d& poly, const Vec2d& pt) {
    auto& ring = poly.outer;
    int n = (int)ring.size();
    if (n < 3) return false;
    bool inside = false;
    for (int i = 0, j = n - 1; i < n; j = i++) {
        double xi = ring[i][0], yi = ring[i][1], xj = ring[j][0], yj = ring[j][1];
        if (((yi > pt.y()) != (yj > pt.y())) && (pt[0] < (xj - xi) * (pt[1] - yi) / (yj - yi + 1e-20) + xi))
            inside = !inside;
    }
    if (inside)
        for (auto& h : poly.holes) {
            Polygon2d hp;
            hp.outer = h;
            if (polygonContains(hp, pt))
                return false;
        }
    return inside;
}

/// 点到多边形外环的最短距离
static double pointToPolygonDist(const Vec2d& pt, const Polygon2d& poly) {
    auto& ring = poly.outer;
    int n = (int)ring.size();
    if (n < 2) return 1e18;
    double m = 1e18;
    for (int i = 0; i < n; ++i) {
        Vec2d a = ring[i], b = ring[(i + 1) % n], ab = b - a, ap = pt - a;
        double t = std::max(0.0, std::min(1.0, ab.dot(ap) / std::max(1e-20, ab.squaredNorm())));
        m = std::min(m, (pt - (a + t * ab)).norm());
    }
    return m;
}

/// 首尾连接点附近落在粗糙路口面外环外侧时的判定容差。
///
/// 连接点由输入车道端点决定，生成器无法移动它；而粗糙路口面的外环恰好是从各进出口
/// 车道端头连出来的，所以外环在连接点处**穿过曲线本身**。于是有两类不可归因于生成器
/// 的"越界"：
///
///   一、连接点采样点本身。"在环内还是环外"完全由坐标舍入决定，生成器在该点自由度为零，
///       报它必然是误伤。允许量取 `kFenceEndpointOutsideTol`，只用来防住"车道端点被数据
///       放在路口面外很远"的真实数据错误。
///   二、连接点近侧一小段。曲线在连接点的切向也由输入车道方向固定，若它几乎与外环该段
///       平行，则离开后数厘米内的采样点会成片读出亚毫米级外溢，且采样越密读出越多
///       （sps=25 只命中端点，sps=64 会多命中 t=0.9844 处 0.23mm 一点）。这一段的允许量
///       必须远小于任何真实外溢：取 `kFenceConnectionRoundingTol`。
///
/// corpus 定标（sps=64 全量参数采样，967 个界外采样点）：
///   * 全部界外点都落在连接点 0.20m 内的曲线共 18 条，最大外溢 `0.001665m`；
///   * 存在 0.20m 以外界外点的曲线共 26 条，最大外溢 `0.0038~8.7518m`。
/// 两类在"最大外溢"上的分界是 `0.001665m` 对 `0.003837m`，取几何中点 `0.0025m`，
/// 两侧各留 1.5 倍余量。跨度 `0.20m` 与 §6.8.7 离开楔形的 `anchor_tol` 一致。
///
/// 安全性：豁免是逐采样点的。真实外溢即便从连接点附近起步（实测最小起步值 `0.000244m`），
/// 其增长段的采样点 `ep_d` 必然远超 0.20m，仍会被报出——实测 26 条真实外溢曲线的
/// 最大外溢点 `ep_d` 均 ≥ 2.0m。只有"整段外溢都挤在连接点 0.20m 内且始终不超 2.5mm"
/// 才会被整条放过，而这个量级已低于输入坐标本身的精度。
constexpr double kFenceConnectionPointRadius = 1e-3;
constexpr double kFenceEndpointOutsideTol = 0.05;
constexpr double kFenceConnectionSpan = 0.20;
constexpr double kFenceConnectionRoundingTol = 0.0025;

/// 该越界采样点是否只是连接点附近的坐标舍入。
static bool fenceOutsideIsConnectionPointRounding(
    const Polygon2d& fence, const Vec2d& pt, const Vec2d& p0, const Vec2d& p1) {
    const double ep = std::min((pt - p0).norm(), (pt - p1).norm());
    if (ep > kFenceConnectionSpan)
        return false;
    const double d = pointToPolygonDist(pt, fence);
    if (ep <= kFenceConnectionPointRadius)
        return d <= kFenceEndpointOutsideTol;  // 连接点本身：生成器自由度为零
    return d <= kFenceConnectionRoundingTol;   // 近侧贴外环的一小段
}

/// 曲线是否完全位于围栏内
static bool curveInsideFence(const BezierCurve& c, const Polygon2d& fence, int sps = 25) {
    if (fence.outer.empty()) return true;
    if (c.empty()) return true;
    const Vec2d p0 = c.startPt(), p1 = c.endPt();
    for (auto& seg : c.segs)
        for (int i = 0; i <= sps; ++i) {
            const Vec2d pt = seg.evaluate((double)i / sps);
            if (polygonContains(fence, pt)) continue;
            if (fenceOutsideIsConnectionPointRounding(fence, pt, p0, p1)) continue;
            return false;
        }
    return true;
}

/// 曲线相对围栏的最大有效外溢量（连接点舍入豁免同 curveInsideFence）。
/// 采样按弧长均匀，`samples` 为采样点数上限的建议值。
static double curveFenceOverflow(const BezierCurve& c, const Polygon2d& fence, int samples) {
    if (fence.outer.empty() || c.empty()) return 0.0;
    const Vec2d p0 = c.startPt(), p1 = c.endPt();
    double overflow = 0.0;
    for (const auto& pt : c.sampleByArcLength(samples)) {
        if (polygonContains(fence, pt)) continue;
        if (fenceOutsideIsConnectionPointRounding(fence, pt, p0, p1)) continue;
        overflow = std::max(overflow, pointToPolygonDist(pt, fence));
    }
    return overflow;
}

/// 连接两端连接点的直线弦相对围栏的最大外溢（连接点舍入豁免同上）。
///
/// 两个连接点由输入车道端点固定，生成器不能移动它们。若连它们的**直线弦**都在粗糙
/// 路口面之外，说明该面在这两点之间是内凹的——典型情形是粗糙面把各进出口"喉部"用
/// 直线弦连起来，于是把转角整块切掉。此时任何贴弦而行的曲线都不可能待在面内，右转
/// 更是只能从被切掉的那一侧过去：外溢的成因是面画小了，不是曲线画歪了。
///
/// 该量给出"几何强制的外溢下限"，供审计把不可归因于生成器的外溢与真正的形态缺陷
/// 分开。corpus 实测（见架构文档 6.8.8）：8 条右转的曲线外溢 `0.14~1.96m`，而其弦
/// 外溢 `2.91~4.39m`——曲线已经比直连好得多，是优化器在向内挤压的结果。
static double fenceChordOverflow(const Polygon2d& fence, const Vec2d& p0,
                                 const Vec2d& p1, int samples) {
    if (fence.outer.empty() || samples < 1)
        return 0.0;
    double overflow = 0.0;
    for (int i = 0; i <= samples; ++i) {
        const Vec2d pt = p0 + (p1 - p0) * ((double)i / samples);
        if (polygonContains(fence, pt)) continue;
        if (fenceOutsideIsConnectionPointRounding(fence, pt, p0, p1)) continue;
        overflow = std::max(overflow, pointToPolygonDist(pt, fence));
    }
    return overflow;
}

/// 曲线外溢是否已被"面本身不含该弦"这一几何事实强制。
///
/// 判据：曲线外溢不超过弦外溢（容差取坐标舍入尺度 `kFenceConnectionRoundingTol`，
/// 低于该尺度的差值没有物理意义，且两者采样点集不同本身就带这个量级的误差）。
/// 实测分界：右转与直行 10 例全部满足（曲线 ≤ 弦）；真正超出的 3 例分别超出
/// `0.0018m`（低于舍入尺度，判为强制）、`0.2386m`、`5.2154m`（判为形态缺陷）。
static bool fenceOverflowForcedByFace(double curve_overflow, double chord_overflow) {
    return curve_overflow <= chord_overflow + kFenceConnectionRoundingTol;
}

}
