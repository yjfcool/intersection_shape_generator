#include "infeasibility_detector.h"
#include "fence_check.h"
#include "optimizer/sdf_field.h"
#include "utils.h"
#include "utils/clipper.hpp"
#include <queue>
#include <unordered_set>
#include <cmath>
#include <algorithm>

namespace isg {

bool segmentsIntersect2(const Vec2d& a, const Vec2d& b, const Vec2d& c, const Vec2d& d) {
    Vec2d r = b - a, s = d - c;
    double den = cross2d(r, s);
    if (std::abs(den) < 1e-12)
        return false;
    Vec2d ac = c - a;
    double t = cross2d(ac, s) / den, u = cross2d(ac, r) / den;
    return t >= 0 && t <= 1 && u >= 0 && u <= 1;
}

bool isSimplePolygon(const Polygon2d& poly) {
    auto ring = poly.outer;
    if (ring.front() == ring.back())
        ring.pop_back();
    int n = (int)ring.size();
    if (n < 4)
        return true;
    for (int i = 0; i < n; ++i)
        for (int j = i + 2; j < n; ++j) {
            if (i == 0 && j == n - 1)
                continue;
            if (segmentsIntersect2(ring[i], ring[(i + 1) % n], ring[j], ring[(j + 1) % n]))
                return false;
        }
    return true;
}

bool bfsReachable(const SDFField& sdf, const Polygon2d& fence,
                  const Vec2d& start, const Vec2d& goal, double mc) {
    if (!sdf.valid())
        return true;
    std::pair<int, int> _sc = sdf.worldToCell(start);
    int sr = _sc.first;
    int sc_ = _sc.second;
    std::pair<int, int> _gc = sdf.worldToCell(goal);
    int gr = _gc.first;
    int gc_ = _gc.second;
    int rows = sdf.rows(), cols = sdf.cols();
    auto key = [&](int r, int c) { return r * cols + c; };
    std::unordered_set<int> vis;
    std::queue<std::pair<int, int>> q;
    auto push = [&](int r, int c) {
        if (r < 0 || r >= rows || c < 0 || c >= cols)
            return;
        if (vis.count(key(r, c)))
            return;
        Vec2d w = sdf.cellToWorld(r, c);
        std::pair<double, Vec2d> _q = sdf.queryWithGrad(w);
        if (_q.first < mc)
            return;
        if (!fence.outer.empty() && !polygonContains(fence, w))
            return;
        vis.insert(key(r, c));
        q.push({r, c});
    };
    push(sr, sc_);
    int dx[] = {0, 0, 1, -1, 1, 1, -1, -1}, dy[] = {1, -1, 0, 0, 1, -1, 1, -1};
    while (!q.empty()) {
        std::pair<int, int> _fr = q.front();
        int r = _fr.first;
        int c = _fr.second;
        q.pop();
        if (r == gr && c == gc_)
            return true;
        for (int d = 0; d < 8; ++d)
            push(r + dx[d], c + dy[d]);
    }
    return false;
}

bool detectSandwich(const SDFField& sdf, const Polygon2d& fence, double cl) {
    if (!sdf.valid() || fence.outer.empty())
        return false;
    auto& ring = fence.outer;
    int n = (int)ring.size();
    for (int i = 0; i < n; ++i) {
        Vec2d a = ring[i], b = ring[(i + 1) % n], mid = 0.5 * (a + b);
        Vec2d edge = (b - a).normalized();
        Vec2d inward{edge[1], -edge[0]};
        Vec2d inner = mid + cl * inward;
        std::pair<double, Vec2d> _q2 = sdf.queryWithGrad(inner);
        if (_q2.first < cl)
            return true;
    }
    return false;
}

FenceRelaxResult tryRelaxFence(const Polygon2d& orig, const std::vector<Boundary>& edges,
                               const Vec2d&, double mx) {
    if (orig.outer.empty())
        return {false, {}, 0};
    Polygon2d relaxed;
    auto pths = ClipperUtil::InflatePaths({toArray(orig.outer)}, mx, ClipperLib::jtMiter, ClipperLib::etClosedPolygon);
    relaxed.outer = toArray(pths[0]);
    for (auto& bnd : edges) {
        if (bnd.type != Boundary::Type::RoadEdge)
            continue;
        auto& pts = bnd.geometry.points;
        auto& ring = relaxed.outer;
        int nr = (int)ring.size(), nb = (int)pts.size();
        for (int i = 0; i < nr; ++i)
            for (int j = 0; j + 1 < nb; ++j)
                if (segmentsIntersect2(ring[i], ring[(i + 1) % nr], pts[j], pts[j + 1]))
                    return {false, {}, 0};
    }
    return {true, relaxed, mx};
}

PreCheckResult preCheck(const SDFField& sdf, const Polygon2d& fence,
                        const Vec2d& ep, const Vec2d& xp, double rw, const std::vector<Boundary>&) {
    PreCheckResult r;
    r.topological_block = !bfsReachable(sdf, fence, ep, xp, 0.3);
    double minW = 1e18;
    for (int i = 1; i < 20; ++i) {
        double t = (double)i / 20;
        Vec2d pt = (1 - t) * ep + t * xp;
        if (!fence.outer.empty() && !polygonContains(fence, pt))
            continue;
        std::pair<double, Vec2d> _q3 = sdf.queryWithGrad(pt);
        minW = std::min(minW, _q3.first * 2.0);
    }
    r.min_corridor_width = minW;
    r.narrow_passage = (minW < rw * 1.2);
    r.fence_sandwich = detectSandwich(sdf, fence, 0.5);
    if (r.topological_block)
        r.type = ViolationInfo::InfeasibilityType::TopologicalBlock;
    else if (r.fence_sandwich)
        r.type = ViolationInfo::InfeasibilityType::Sandwich;
    else if (r.narrow_passage)
        r.type = ViolationInfo::InfeasibilityType::NarrowPassage;
    return r;
}

ConnectivityCurve makeFallbackCurve(const PreCheckResult& pre, const Connectivity& conn,
                                    const Vec2d& p0, const Vec2d& p1) {
    ConnectivityCurve out;
    out.id = conn.id;
    out.entry_lane_id = conn.entry_lane_id;
    out.exit_lane_id = conn.exit_lane_id;
    out.turn_type = conn.turn_type;
    out.violation.type = pre.type;
    out.fixed_shape = conn.fixed_shape;
    if (pre.type == ViolationInfo::InfeasibilityType::TopologicalBlock) {
        out.curve = nullptr;
        out.status = CurveStatus::Infeasible;
        out.violation.reason = "No obstacle-free path (BFS unreachable)";
    } else {
        Vec2d dir = (p1 - p0);
        double d = dir.norm();
        if (d > 1e-6)
            dir.normalize();
        else
            dir = Vec2d(1, 0);
        BezierSegment seg;
        seg.ctrl[0] = p0;
        seg.ctrl[1] = p0 + dir * (d / 3);
        seg.ctrl[2] = p1 - dir * (d / 3);
        seg.ctrl[3] = p1;
        BezierCurve c;
        c.segs.push_back(seg);
        out.curve = std::make_shared<BezierCurve>(c);
        out.geometry.points = {p0, p1};
        out.status = CurveStatus::Degraded;
        out.violation.reason = "Sandwich or narrow; degraded to straight";
    }
    return out;
}

}
