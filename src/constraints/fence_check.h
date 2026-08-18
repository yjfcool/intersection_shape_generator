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

/// 曲线是否完全位于围栏内
static bool curveInsideFence(const BezierCurve& c, const Polygon2d& fence, int sps = 25) {
    if (fence.outer.empty()) return true;
    for (auto& seg : c.segs)
        for (int i = 0; i <= sps; ++i)
            if (!polygonContains(fence, seg.evaluate((double)i / sps)))
                return false;
    return true;
}

}
