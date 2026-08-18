#pragma once
#include "types.h"
#include <string>
#include <vector>
#include <array>
#include <map>
#include <unordered_map>
#include <functional>
#include <cmath>
#include <algorithm>
#include <limits>

namespace isg {

static constexpr double PI = 3.14159265358979323846;
static constexpr double DEG2RAD = PI / 180.0;
static constexpr double RAD2DEG = 180.0 / PI;
static constexpr double EPS = 1e-9;

/// 障碍避让清距规则: mode=2 需要保留1m绕障距离,其他模式只禁止穿行。
inline double obstacleAvoidanceClearanceForMode(int mode) {
    return mode == 2 ? 1.0 : 0.0;
}

/// RoadEdge避让清距规则: mode=2 对道路边缘也保留1m非端点净距。
inline double roadEdgeAvoidanceClearanceForMode(int mode) {
    return mode == 2 ? 1.0 : 0.0;
}

/// 路口面生成缩扩距离: mode=2时缩扩0.05米,其他模式0.5米
inline double intersectionAreaExtendDistance(int mode) {
    return mode == 2 ? 0.05 : 0.5;
}

// 内联数学工具
inline Vec2d xyOf(const Vec2d& p) { return p; }
inline Vec2d xyOf(const Vec3d& p) { return p.xy(); }
inline Vec3d xyzOf(const Vec2d& p) { return Vec3d(p, 0.0); }
inline Vec3d xyzOf(const Vec3d& p) { return p; }

inline double dot(const Vec2d& a, const Vec2d& b) { return a[0] * b[0] + a[1] * b[1]; }
inline double cross2d(const Vec2d& a, const Vec2d& b) { return a[0] * b[1] - a[1] * b[0]; }
inline double dist(const Vec2d& a, const Vec2d& b) { return (a - b).norm(); }
inline double dist(const Vec3d& a, const Vec3d& b) { return (a.xy() - b.xy()).norm(); }
inline double dist(const Vec3d& a, const Vec2d& b) { return (a.xy() - b).norm(); }
inline double dist(const Vec2d& a, const Vec3d& b) { return (a - b.xy()).norm(); }

inline Vec3d offsetXY(const Vec3d& p, const Vec2d& delta) {
    return Vec3d(p[0] + delta[0], p[1] + delta[1], p[2]);
}

/// 两向量夹角（弧度）
inline double angleBetween(const Vec2d& a, const Vec2d& b) {
    double c = a.normalized().dot(b.normalized());
    return std::acos(std::max(-1.0, std::min(1.0, c)));
}

// 向左逆时针旋转90°
inline Vec2d rotLeft(Vec2d p) { return Vec2d{-p[1], p[0]}; }
// 向右顺时针旋转90°
inline Vec2d rotRight(Vec2d p) { return Vec2d{p[1], -p[0]}; }
// 点到线段的最短距离（和最近点参数t∈[0,1]）
template <typename P, typename A, typename B>
inline std::pair<double, double> pointToSegment(
    const P& p, const A& a, const B& b) {
    Vec2d pp = xyOf(p);
    Vec2d aa = xyOf(a);
    Vec2d bb = xyOf(b);
    Vec2d ab = bb - aa;
    double len2 = ab[0] * ab[0] + ab[1] * ab[1];
    if (len2 < EPS * EPS) {
        return {dist(pp, aa), 0.0};
    }
    double t = std::max(0.0, std::min(1.0, (pp - aa).dot(ab) / len2));
    Vec2d closest = aa + ab * t;
    return {dist(pp, closest), t};
}

// 点到折线的最短距离
template <typename P, typename Q>
inline double pointToPolyline(const P& p, const std::vector<Q>& poly) {
    if (poly.empty())
        return std::numeric_limits<double>::max();
    if (poly.size() == 1)
        return dist(p, poly[0]);
    double minD = std::numeric_limits<double>::max();
    for (size_t i = 0; i + 1 < poly.size(); ++i) {
        std::pair<double, double> _seg = pointToSegment(p, poly[i], poly[i + 1]);
        minD = std::min(minD, _seg.first);
    }
    return minD;
}

// 两线段是否相交（不含共端点）
// 返回 true 表示严格相交（内部相交）
template <typename A, typename B, typename C, typename D>
inline bool segmentsIntersectStrict(
    const A& a_in, const B& b_in, const C& c_in, const D& d_in) {
    Vec2d a = xyOf(a_in);
    Vec2d b = xyOf(b_in);
    Vec2d c = xyOf(c_in);
    Vec2d d = xyOf(d_in);
    // 使用叉积方向判断
    auto sign = [](double v) -> int {
        if (v > EPS)
            return 1;
        if (v < -EPS)
            return -1;
        return 0;
    };
    Vec2d ab = b - a, ac = c - a, ad = d - a;
    Vec2d cd = d - c, ca = a - c, cb = b - c;
    int d1 = sign(cross2d(ab, ac));
    int d2 = sign(cross2d(ab, ad));
    int d3 = sign(cross2d(cd, ca));
    int d4 = sign(cross2d(cd, cb));

    if (((d1 > 0 && d2 < 0) || (d1 < 0 && d2 > 0)) && ((d3 > 0 && d4 < 0) || (d3 < 0 && d4 > 0)))
        return true;
    return false;
}

// 两折线是否有内部相交（排除端点）
template <typename PA, typename PB>
inline bool polylinesIntersect(const std::vector<PA>& A, const std::vector<PB>& B) {
    if (A.size() < 2 || B.size() < 2)
        return false;
    for (size_t i = 0; i + 1 < A.size(); ++i) {
        for (size_t j = 0; j + 1 < B.size(); ++j) {
            if (segmentsIntersectStrict(A[i], A[i + 1], B[j], B[j + 1])) {
                // 排除真正的端点（首点或尾点相交）,检查交点是否恰好是折线的首/尾端点,如果两段的端点相同则不算相交
                if (dist(A[i], B[j]) < 1e-4 ||
                    dist(A[i], B[j + 1]) < 1e-4 ||
                    dist(A[i + 1], B[j]) < 1e-4 ||
                    dist(A[i + 1], B[j + 1]) < 1e-4)
                    continue;
                return true;
            }
        }
    }
    return false;
}

// 检查两条折线间排除首尾连接点后的相交
// endPtsA: A的允许相交端点集; endPtsB: B的
template <typename PA, typename PB>
inline bool polylinesIntersectExcludeEndpoints(
    const std::vector<PA>& A, const std::vector<PB>& B, double endPtTol = 1e-4) {
    if (A.size() < 2 || B.size() < 2)
        return false;
    for (size_t i = 0; i + 1 < A.size(); ++i) {
        for (size_t j = 0; j + 1 < B.size(); ++j) {
            if (!segmentsIntersectStrict(A[i], A[i + 1], B[j], B[j + 1]))
                continue;
            // 判断是否是端点附近的共点（允许的连接点相交）
            bool nearEndA = (i == 0 || i + 1 == A.size() - 1);
            bool nearEndB = (j == 0 || j + 1 == B.size() - 1);
            // 如果交点附近存在公共端点 → 跳过
            auto checkEndpts = [&]() {
                const Vec2d* ptsA[2] = {&A[i], &A[i + 1]};
                const Vec2d* ptsB[2] = {&B[j], &B[j + 1]};
                // 仅当两折线的真正首/尾端点重合时才允许
                bool aStart = (dist(A.front(), A[i]) < endPtTol || dist(A.front(), A[i + 1]) < endPtTol);
                bool aEnd = (dist(A.back(), A[i]) < endPtTol || dist(A.back(), A[i + 1]) < endPtTol);
                bool bStart = (dist(B.front(), B[j]) < endPtTol || dist(B.front(), B[j + 1]) < endPtTol);
                bool bEnd = (dist(B.back(), B[j]) < endPtTol || dist(B.back(), B[j + 1]) < endPtTol);
                (void)ptsA;
                (void)ptsB;
                if ((aStart || aEnd) && (bStart || bEnd)) {
                    // 检查端点是否重合
                    if (dist(A.front(), B.front()) < endPtTol || dist(A.front(), B.back()) < endPtTol
                        || dist(A.back(), B.front()) < endPtTol || dist(A.back(), B.back()) < endPtTol)
                        return true; // 端点相交，允许
                }
                return false;
            };
            (void)nearEndA;
            (void)nearEndB;
            if (checkEndpts())
                continue;
            return true; // 真正的内部相交
        }
    }
    return false;
}

/// 折线按弧长重采样，固定间距spacing
template <typename P>
inline std::vector<P> resampleBySpacing(const std::vector<P>& src, double spacing) {
    if (src.size() < 2 || spacing < EPS)
        return src;
    std::vector<P> out;
    out.push_back(src.front());
    double traveled = 0;
    double nextSample = spacing;
    for (size_t i = 1; i < src.size(); ++i) {
        double seg = dist(src[i - 1], src[i]);
        while (traveled + seg >= nextSample) {
            double t = (nextSample - traveled) / std::max(seg, EPS);
            P p = src[i - 1] * (1 - t) + src[i] * t;
            out.push_back(p);
            nextSample += spacing;
        }
        traveled += seg;
    }
    if (dist(out.back(), src.back()) > EPS * 10)
        out.push_back(src.back());
    return out;
}

// 点是否在多边形内（射线法）
template <typename P, typename Q>
inline bool pointInPolygon(const P& pt_in, const std::vector<Q>& poly) {
    if (poly.size() < 3)
        return false;
    Vec2d pt = xyOf(pt_in);
    int cnt = 0;
    size_t n = poly.size();
    for (size_t i = 0, j = n - 1; i < n; j = i++) {
        auto& pi = poly[i];
        auto& pj = poly[j];
        if (((pi[1] > pt[1]) != (pj[1] > pt[1])) &&
            (pt[0] < (pj[0] - pi[0]) * (pt[1] - pi[1]) / (pj[1] - pi[1]) + pi[0]))
            cnt++;
    }
    return cnt % 2 == 1;
}

// 多边形面积（Shoelace公式，有符号）
template <typename P>
inline double polygonSignedArea(const std::vector<P>& poly) {
    double A = 0;
    size_t n = poly.size();
    for (size_t i = 0; i < n; ++i) {
        const Vec2d& a = poly[i];
        const Vec2d& b = poly[(i + 1) % n];
        A += a[0] * b[1] - b[0] * a[1];
    }
    return A * 0.5;
}

template <typename P>
inline double polygonArea(const std::vector<P>& p) { return std::abs(polygonSignedArea(p)); }

// ── 车道几何辅助 ─────────────────────────────────────────────
/// 进入车道：geometry 从外向路口，末点为路口端点
/// 退出车道：geometry 从路口向外，首点为路口端点
template <typename P>
inline Vec2d getConnPoint(const std::vector<P>& points, bool is_entryline) {
    if (points.size() < 1)
        return Vec2d(0, 0);
    return xyOf(is_entryline ? points.back() : points.front());
}

/// 获取车道端点三维坐标；XY 参与平面几何，Z 用于路口面/输出高程保留。
template <typename P>
inline Vec3d getConnPoint3d(const std::vector<P>& points, bool is_entryline) {
    if (points.size() < 1)
        return Vec3d(0, 0, 0);
    return xyzOf(is_entryline ? points.back() : points.front());
}

/// 获取车道端点切向（指向路口内/外）
template <typename P>
inline Vec2d getConnTangent(const std::vector<P>& points, bool is_entryline) {
    int n = points.size();
    if (n < 2)
        return Vec2d(1, 0);
    Vec2d d = is_entryline ? xyOf(points[n - 1]) - xyOf(points[n - 2]) : xyOf(points[1]) - xyOf(points[0]);
    return d.norm() > 1e-10 ? d.normalized() : Vec2d(1, 0);
}

/// 同时获取车道端点与切向
template <typename P>
inline std::pair<Vec2d, Vec2d> getConnInfo(const std::vector<P>& points, bool is_entryline) {
    return {getConnPoint(points, is_entryline), getConnTangent(points, is_entryline)};
}

/// 进入车道末点
template <typename P>
inline Vec2d entryLinePoint(const std::vector<P>& points) {
    if (points.size() < 1)
        return Vec2d(0, 0);
    return xyOf(points.back());
}

/// 进入车道末点切向
template <typename P>
inline Vec2d entryLineTangent(const std::vector<P>& points) {
    int n = (int)points.size();
    if (n < 2)
        return Vec2d(1, 0);
    Vec2d d = xyOf(points[n - 1]) - xyOf(points[n - 2]);
    return d.norm() > 1e-10 ? d.normalized() : Vec2d(1, 0);
}

/// 退出车道首点
template <typename P>
inline Vec2d exitLinePoint(const std::vector<P>& points) {
    if (points.size() < 1)
        return Vec2d(0, 0);
    return xyOf(points.front());
}

/// 退出车道首点切向
template <typename P>
inline Vec2d exitLineTangent(const std::vector<P>& points) {
    int n = (int)points.size();
    if (n < 2)
        return Vec2d(1, 0);
    Vec2d d = xyOf(points[1]) - xyOf(points[0]);
    return d.norm() > 1e-10 ? d.normalized() : Vec2d(1, 0);
}

}
