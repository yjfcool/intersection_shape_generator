#include "curve_utils.h"
#include "optimizer/sdf_field.h"
#include "constraints/fence_check.h"
#include "utils.h"
#include <cmath>
#include <algorithm>

namespace isg {

double localCurvature(const Vec2d& a, const Vec2d& b, const Vec2d& c) {
    double ab = (b - a).norm(), bc = (c - b).norm(), ac = (c - a).norm();
    double area = 0.5 * std::abs(cross2d(b - a, c - a));
    double den = ab * bc * ac;
    return den > 1e-12 ? (2 * area / den) : 0.0;
}

Vec2d circumcenter(const Vec2d& a, const Vec2d& b, const Vec2d& c) {
    double ax = a[0], ay = a[1], bx = b[0], by = b[1], cx = c[0], cy = c[1];
    double D = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by));
    if (std::abs(D) < 1e-12)return 0.5 * (a + c);
    double ux = ((ax * ax + ay * ay) * (by - cy) + (bx * bx + by * by) * (cy - ay) + (cx * cx + cy * cy) * (ay - by)) / D;
    double uy = ((ax * ax + ay * ay) * (cx - bx) + (bx * bx + by * by) * (ax - cx) + (cx * cx + cy * cy) * (bx - ax)) / D;
    return Vec2d(ux, uy);
}

std::vector<Vec2d> elasticBandSmooth(const std::vector<Vec2d>& pts, const SDFField& sdf,
                                     const Polygon2d& fence, double km, double ms, int mi, double msc) {
    auto res = pts;
    int N = (int)res.size();
    if (N < 3) return res;

    for (int iter = 0; iter < mi; ++iter) {
        bool ch = false;
        for (int i = 1; i < N - 1; ++i) {
            double k = localCurvature(res[i - 1], res[i], res[i + 1]);
            if (k <= km) continue;
            Vec2d cen = circumcenter(res[i - 1], res[i], res[i + 1]);
            Vec2d dir = res[i] - cen;
            if (dir.norm() < 1e-10) continue;
            dir.normalize();
            Vec2d cand = res[i] + dir * ms;
            std::pair<double, Vec2d> _q = sdf.queryWithGrad(cand);
            if (_q.first >= msc && polygonContains(fence, cand)) {
                res[i] = cand;
                ch = true;
            }
        }
        if (!ch) break;
    }
    return res;
}

// 自适应版本: 根据曲线特征调整采样数
std::vector<Vec2d> elasticBandSmoothAdaptive(const std::vector<Vec2d>& pts, const SDFField& sdf,
                                             const Polygon2d& fence, double km, double ms, int mi, double msc) {
    if (pts.size() < 3) return pts;

    // 计算曲线长度确定初始采样数
    double curve_length = 0.0;
    for (size_t i = 1; i < pts.size(); ++i) {
        curve_length += (pts[i] - pts[i-1]).norm();
    }

    // 基于长度的基准采样数,带上下界
    int initial_samples = std::max(10, std::min(200, static_cast<int>(curve_length / 0.5)));

    // 估计曲率变化以指导自适应采样
    double total_curvature = 0.0;
    int curvature_samples = std::min(static_cast<int>(pts.size()), 20); // 限制采样数以提高效率
    for (int i = 1; i < curvature_samples - 1; ++i) {
        double k = localCurvature(pts[std::min(i-1, static_cast<int>(pts.size()-1))],
                                  pts[i],
                                  pts[std::min(i+1, static_cast<int>(pts.size()-1))]);
        total_curvature += k;
    }

    // 根据曲率调整采样数
    double avg_curvature = total_curvature / std::max(1, curvature_samples - 2);
    int adaptive_samples = std::max(10, static_cast<int>(initial_samples * (1.0 + avg_curvature * 0.5)));

    // 点过多时抽稀
    std::vector<Vec2d> current_pts = pts;
    if (current_pts.size() > adaptive_samples) {
        current_pts = downsamplePoints(current_pts, adaptive_samples);
    }

    int N = (int)current_pts.size();
    if (N < 3) return current_pts;

    // 基于改进阈值的早停
    double improvement_threshold = 1e-4;
    std::vector<Vec2d> prev_pts = current_pts;

    for (int iter = 0; iter < mi; ++iter) {
        bool ch = false;
        for (int i = 1; i < N - 1; ++i) {
            double k = localCurvature(current_pts[i - 1], current_pts[i], current_pts[i + 1]);
            if (k <= km) continue;
            Vec2d cen = circumcenter(current_pts[i - 1], current_pts[i], current_pts[i + 1]);
            Vec2d dir = current_pts[i] - cen;
            if (dir.norm() < 1e-10) continue;
            dir.normalize();
            Vec2d cand = current_pts[i] + dir * ms;
            std::pair<double, Vec2d> _q = sdf.queryWithGrad(cand);
            if (_q.first >= msc && polygonContains(fence, cand)) {
                current_pts[i] = cand;
                ch = true;
            }
        }

        // 没有发生有效移动时提前终止。
        if (!ch) break;

        // 若本轮最大位移低于阈值，则认为已收敛。
        double max_displacement = 0.0;
        for (int i = 0; i < N; ++i) {
            double displacement = (current_pts[i] - prev_pts[i]).norm();
            max_displacement = std::max(max_displacement, displacement);
        }

        if (max_displacement < improvement_threshold) {
            break;
        }

        prev_pts = current_pts;
    }

    return current_pts;
}

// 将点列按索引均匀抽稀到目标数量。
std::vector<Vec2d> downsamplePoints(const std::vector<Vec2d>& points, int target_count) {
    if (points.size() <= target_count) {
        return points;
    }

    std::vector<Vec2d> result;
    result.reserve(target_count);

    double step = static_cast<double>(points.size() - 1) / (target_count - 1);
    for (int i = 0; i < target_count; i++) {
        int idx = std::min(static_cast<int>(i * step), static_cast<int>(points.size() - 1));
        result.push_back(points[idx]);
    }

    return result;
}

BezierCurve rebuildFromSmoothedPts(const std::vector<Vec2d>& s, const Vec2d& st, const Vec2d& et) {
    return fitBezierWithEndTangents(s, st, et);
}

std::vector<Vec2d> midlineSampleByArcLength(const BezierCurve& a, const BezierCurve& b, int n) {
    auto sa = a.sampleByArcLength(n), sb = b.sampleByArcLength(n);
    std::vector<Vec2d> m(n);
    for (int i = 0; i < n; ++i)m[i] = 0.5 * (sa[i] + sb[i]);
    return m;
}

double signedDistToLine(const Vec2d& pt, const Vec2d& p0, const Vec2d& p1) {
    Vec2d d = p1 - p0;
    double len = d.norm();
    return len < 1e-12 ? (pt - p0).norm() : cross2d(d, pt - p0) / len;
}

bool segmentsIntersect(const Vec2d& a, const Vec2d& b, const Vec2d& c, const Vec2d& d, Vec2d* out) {
    Vec2d r = b - a, s = d - c;
    double den = cross2d(r, s);
    if (std::abs(den) < 1e-12) {
        if (std::abs(cross2d(c - a, r)) > 1e-9)
            return false;
        double rr = r.squaredNorm();
        double ss = s.squaredNorm();
        if (rr < 1e-18 || ss < 1e-18)
            return false;
        double t0 = (c - a).dot(r) / rr;
        double t1 = (d - a).dot(r) / rr;
        double lo = std::max(0.0, std::min(t0, t1));
        double hi = std::min(1.0, std::max(t0, t1));
        if (lo > hi + 1e-9)
            return false;
        if (out)
            *out = a + 0.5 * (lo + hi) * r;
        return true;
    }
    Vec2d ac = c - a;
    double t = cross2d(ac, s) / den, u = cross2d(ac, r) / den;
    if (t >= 0 && t <= 1 && u >= 0 && u <= 1) {
        if (out)*out = a + t * r;
        return true;
    }
    return false;
}

double minSDFAlongCurve(const BezierCurve& c, const SDFField& sdf, int sps) {
    double m = 1e18;
    for (auto& seg : c.segs)
        for (int i = 0; i <= sps; ++i) {
            std::pair<double, Vec2d> _q = sdf.queryWithGrad(seg.evaluate((double)i / sps));
            m = std::min(m, _q.first);
        }
    return m;
}

double distToAllEndpoints(const Vec2d& pt, const BezierCurve& a, const BezierCurve& b) {
    double d = 1e18;
    if (!a.empty()) {
        d = std::min(d, (pt - a.startPt()).norm());
        d = std::min(d, (pt - a.endPt()).norm());
    }
    if (!b.empty()) {
        d = std::min(d, (pt - b.startPt()).norm());
        d = std::min(d, (pt - b.endPt()).norm());
    }
    return d;
}


bool bboxOverlap(const BezierCurve& a, const BezierCurve& b) {
    return a.bbox().intersects(b.bbox());
}

static void segPairCross(const BezierSegment& a, const BezierSegment& b,
                         std::vector<Vec2d>& out, double tol, int depth = 0) {
    if (!a.bbox().intersects(b.bbox()))
        return;
    if (depth > 10) {
        out.push_back(a.evaluate(0.5));
        return;
    }
    if (a.arcLength(4) < tol && b.arcLength(4) < tol) {
        out.push_back(a.evaluate(0.5));
        return;
    }
    std::pair<BezierSegment, BezierSegment> _sa = a.splitAt(0.5);
    std::pair<BezierSegment, BezierSegment> _sb = b.splitAt(0.5);
    BezierSegment al = _sa.first, ar = _sa.second;
    BezierSegment bl = _sb.first, br = _sb.second;
    segPairCross(al, bl, out, tol, depth + 1);
    segPairCross(al, br, out, tol, depth + 1);
    segPairCross(ar, bl, out, tol, depth + 1);
    segPairCross(ar, br, out, tol, depth + 1);
}

std::vector<Vec2d> curveCrossings(const BezierCurve& a, const BezierCurve& b, double tol) {
    std::vector<Vec2d> raw;
    if (!bboxOverlap(a, b))
        return raw;
    for (auto& sa : a.segs)
        for (auto& sb : b.segs)
            segPairCross(sa, sb, raw, tol);
    std::vector<Vec2d> out;
    for (auto& p : raw) {
        bool dup = false;
        for (auto& q : out)
            if ((p - q).norm() < tol) {
                dup = true;
                break;
            }
        if (!dup)
            out.push_back(p);
    }
    return out;
}

bool curvesIntersectBusiness(const BezierCurve& a, const BezierCurve& b, double ep) {
    if (!bboxOverlap(a, b))
        return false;
    auto adaptiveSample = [](const BezierCurve& c) {
        int n = std::max(48, (int)std::ceil(c.arcLength() / 0.20) + 1);
        n = std::min(n, 240);
        return c.sampleByArcLength(n);
    };
    auto pa = adaptiveSample(a);
    auto pb = adaptiveSample(b);
    for (int i = 0; i + 1 < (int)pa.size(); ++i) {
        for (int j = 0; j + 1 < (int)pb.size(); ++j) {
            Vec2d isect;
            if (!segmentsIntersect(pa[i], pa[i + 1], pb[j], pb[j + 1], &isect))
                continue;
            if (distToAllEndpoints(isect, a, b) > ep)
                return true;
        }
    }
    return false;
}

bool curveSelfIntersectsBusiness(const BezierCurve& c, double ep) {
    // 性能优化: 采样从80降至40 (O(N²) → 减少4倍)
    auto pts = c.sample(40);
    if (pts.size() < 4)
        return false;
    for (int i = 0; i + 1 < (int)pts.size(); ++i) {
        for (int j = i + 2; j + 1 < (int)pts.size(); ++j) {
            if (i == 0 && j + 1 == (int)pts.size() - 1)
                continue;
            Vec2d isect;
            if (!segmentsIntersect(pts[i], pts[i + 1], pts[j], pts[j + 1], &isect))
                continue;
            double d0 = std::min((isect - pts.front()).norm(), (isect - pts.back()).norm());
            if (d0 <= ep)
                continue;
            return true;
        }
    }
    return false;
}

namespace {

struct SingleCubicHandleBounds {
    Vec2d start_dir{1, 0};
    Vec2d end_dir{1, 0};
    double start_max = 0.0;
    double end_max = 0.0;
};

static SingleCubicHandleBounds singleCubicHandleBounds(
    const Vec2d& p0, const Vec2d& start_tan,
    const Vec2d& p1, const Vec2d& end_tan,
    bool cap_at_direction_intersection) {
    SingleCubicHandleBounds bounds;
    Vec2d chord = p1 - p0;
    double chord_len = chord.norm();
    bounds.start_dir = start_tan.norm() > 1e-8
        ? start_tan.normalized()
        : (chord_len > 1e-8 ? chord.normalized() : Vec2d(1, 0));
    bounds.end_dir = end_tan.norm() > 1e-8
        ? end_tan.normalized()
        : (chord_len > 1e-8 ? chord.normalized() : bounds.start_dir);
    bounds.start_max = std::max(0.05, chord_len);
    bounds.end_max = std::max(0.05, chord_len);
    if (!cap_at_direction_intersection || chord_len < 1e-8)
        return bounds;

    // 直行不使用方向交点上限；普通左右转使用首/尾方向射线的前向交点
    // 作为控制点有效范围，避免把手越过交点后形成反向拱弧。
    Vec2d exit_back = -bounds.end_dir;
    bool straight_like =
        bounds.start_dir.dot(bounds.end_dir) > 0.90 &&
        std::abs(cross2d(bounds.start_dir, bounds.end_dir)) < 0.25 &&
        std::abs(cross2d(bounds.start_dir, chord.normalized())) < 0.25;
    if (std::abs(cross2d(bounds.start_dir, exit_back)) < 1e-8 || straight_like)
        return bounds;

    Vec2d delta = p1 - p0;
    double den = cross2d(bounds.start_dir, exit_back);
    double start_station = cross2d(delta, exit_back) / den;
    double end_station = cross2d(delta, bounds.start_dir) / den;
    if (start_station > 0.05 && end_station > 0.05) {
        bounds.start_max = std::min(bounds.start_max, start_station);
        bounds.end_max = std::min(bounds.end_max, end_station);
    }
    return bounds;
}

static double startHandleStation(
    const Vec2d& control, const Vec2d& endpoint, const Vec2d& direction) {
    return (control - endpoint).dot(direction);
}

} // namespace

void constrainOrdinarySingleCubicControls(
    BezierCurve& curve, const Vec2d& p0, const Vec2d& start_tan,
    const Vec2d& p1, const Vec2d& end_tan,
    bool cap_at_direction_intersection) {
    if (curve.numSegments() != 1)
        return;

    const SingleCubicHandleBounds bounds = singleCubicHandleBounds(
        p0, start_tan, p1, end_tan, cap_at_direction_intersection);
    BezierSegment& seg = curve.segs.front();
    seg.ctrl[0] = p0;
    seg.ctrl[3] = p1;

    const double lambda = std::max(
        0.05, std::min(bounds.start_max,
                       startHandleStation(seg.ctrl[1], p0, bounds.start_dir)));
    const double mu = std::max(
        0.05, std::min(bounds.end_max,
                       (p1 - seg.ctrl[2]).dot(bounds.end_dir)));
    seg.ctrl[1] = p0 + lambda * bounds.start_dir;
    seg.ctrl[2] = p1 - mu * bounds.end_dir;
}

bool ordinarySingleCubicControlsValid(
    const BezierCurve& curve, const Vec2d& p0, const Vec2d& start_tan,
    const Vec2d& p1, const Vec2d& end_tan, double tol,
    bool cap_at_direction_intersection) {
    if (curve.numSegments() != 1)
        return true;

    const SingleCubicHandleBounds bounds = singleCubicHandleBounds(
        p0, start_tan, p1, end_tan, cap_at_direction_intersection);
    const BezierSegment& seg = curve.segs.front();
    if ((seg.ctrl[0] - p0).norm() > tol || (seg.ctrl[3] - p1).norm() > tol)
        return false;

    const double lambda = startHandleStation(seg.ctrl[1], p0, bounds.start_dir);
    const double mu = (p1 - seg.ctrl[2]).dot(bounds.end_dir);
    const double start_lateral =
        std::abs(cross2d(bounds.start_dir, seg.ctrl[1] - p0));
    const double end_lateral =
        std::abs(cross2d(bounds.end_dir, p1 - seg.ctrl[2]));
    return start_lateral <= tol && end_lateral <= tol &&
           lambda >= 0.05 - tol && lambda <= bounds.start_max + tol &&
           mu >= 0.05 - tol && mu <= bounds.end_max + tol;
}

}
