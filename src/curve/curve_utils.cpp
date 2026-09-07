#include "curve_utils.h"
#include "optimizer/sdf_field.h"
#include "constraints/fence_check.h"
#include "utils.h"
#include <cmath>
#include <algorithm>
#include <cstdint>
#include <cstring>
#include <limits>
#include <vector>

namespace isg {

namespace {

////////////////////////////////////////////////////////////
// 成对相交判定的采样索引缓存
////////////////////////////////////////////////////////////
// curvesIntersectBusinessOutsideBalls 是全流程的最热函数（病态数据上占
// 自身耗时的八成以上）。原实现的两处浪费：
//   1. 每次调用都按弧长重新采样两条曲线（各最多 240 点），而生成/修复阶段
//      会用同一批候选反复两两配对，同一条曲线在一轮里被重复采样几十次；
//   2. 内层是 O(Na×Nb) 的裸双层循环，240×240 约 5.7 万次配对全部逐个判断。
//
// 采样只依赖控制点，是纯函数，因此按控制点精确相等复用采样结果，得到的
// 折线与每次重算逐位一致。分块包围盒只用于整段跳过“分段包围盒必然分离”
// 的区间——这些配对在逐段判断里也会被同一条件拒掉，且 i、j 仍按升序访问，
// 首个命中的交点不变，故判定结果与原实现完全一致。
constexpr int kPairBlockSize = 16;
constexpr std::size_t kPairSampleSlots = 256;  // 2 的幂，直接映射

struct PairSampleIndex {
    std::vector<Vec2d> pts;
    std::vector<BoundingBox2d> seg_box;
    std::vector<BoundingBox2d> block_box;
    BoundingBox2d whole;
};

void fillPairSampleIndex(const BezierCurve& c, PairSampleIndex& idx) {
    int n = std::max(48, (int)std::ceil(c.arcLength() / 0.20) + 1);
    n = std::min(n, 240);
    idx.pts = c.sampleByArcLength(n);
    idx.seg_box.clear();
    idx.block_box.clear();
    idx.whole = BoundingBox2d();
    const int nseg = (int)idx.pts.size() - 1;
    if (nseg <= 0)
        return;
    idx.seg_box.resize(nseg);
    idx.block_box.resize((nseg + kPairBlockSize - 1) / kPairBlockSize);
    for (int i = 0; i < nseg; ++i) {
        idx.seg_box[i].expand(idx.pts[i]);
        idx.seg_box[i].expand(idx.pts[i + 1]);
        BoundingBox2d& block = idx.block_box[i / kPairBlockSize];
        block.expand(idx.pts[i]);
        block.expand(idx.pts[i + 1]);
        idx.whole.expand(idx.pts[i]);
    }
    idx.whole.expand(idx.pts[nseg]);
}

bool sameControlPoints(const std::vector<BezierSegment>& x,
                       const std::vector<BezierSegment>& y) {
    if (x.size() != y.size())
        return false;
    for (std::size_t s = 0; s < x.size(); ++s)
        for (int i = 0; i < 4; ++i)
            if (x[s].ctrl[i][0] != y[s].ctrl[i][0] ||
                x[s].ctrl[i][1] != y[s].ctrl[i][1])
                return false;
    return true;
}

std::size_t controlPointSlot(const std::vector<BezierSegment>& segs) {
    std::uint64_t h = 1469598103934665603ULL;
    for (const auto& s : segs)
        for (int i = 0; i < 4; ++i)
            for (int k = 0; k < 2; ++k) {
                std::uint64_t bits = 0;
                const double v = s.ctrl[i][k];
                std::memcpy(&bits, &v, sizeof(bits));
                h = (h ^ bits) * 1099511628211ULL;
            }
    h ^= h >> 29;
    return static_cast<std::size_t>(h) & (kPairSampleSlots - 1);
}

struct PairSampleSlot {
    std::vector<BezierSegment> key;
    BoundingBox2d bbox;    // 与 BezierCurve::bbox() 逐位一致的缓存值
    PairSampleIndex idx;
    bool valid = false;    // key/bbox 已就绪
    bool sampled = false;  // idx 已就绪
};

// 取（必要时建立）该曲线的缓存槽位，只保证 key 与 bbox 就绪；采样按需惰性填充，
// 这样包围盒判否的配对不必付出采样代价。
PairSampleSlot& acquirePairSlot(std::vector<PairSampleSlot>& cache,
                                std::size_t slot, const BezierCurve& c) {
    PairSampleSlot& e = cache[slot];
    if (e.valid && sameControlPoints(e.key, c.segs))
        return e;
    e.key = c.segs;
    e.bbox = c.bbox();
    e.sampled = false;
    e.valid = true;
    return e;
}

const PairSampleIndex& ensurePairSamples(PairSampleSlot& e, const BezierCurve& c) {
    if (!e.sampled) {
        fillPairSampleIndex(c, e.idx);
        e.sampled = true;
    }
    return e.idx;
}

}  // 匿名命名空间

double localCurvature(const Vec2d& a, const Vec2d& b, const Vec2d& c) {
    double ab = (b - a).norm(), bc = (c - b).norm(), ac = (c - a).norm();
    double area = 0.5 * std::abs(cross2d(b - a, c - a));
    double den = ab * bc * ac;
    // 外接圆半径 R = |ab|·|bc|·|ac| / (4·S)，故曲率 κ = 1/R = 4·S / den。
    // 早期实现写作 2·S/den，只有真实曲率的一半，使
    // `elasticBandSmooth` 的 kappa_max 门槛在实际尺度上被放大一倍。
    return den > 1e-12 ? (4 * area / den) : 0.0;
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
    return curvesIntersectBusinessOutsideBalls(a, b, ep, std::vector<Vec2d>(), 0.0);
}

bool curvesIntersectBusinessOutsideBalls(const BezierCurve& a,
                                        const BezierCurve& b, double ep,
                                        const std::vector<Vec2d>& centers,
                                        double radius) {
    if (a.empty() || b.empty())
        return false;
    static thread_local std::vector<PairSampleSlot> cache(kPairSampleSlots);
    const std::size_t slot_a = controlPointSlot(a.segs);
    const std::size_t slot_b = controlPointSlot(b.segs);
    PairSampleSlot& ea = acquirePairSlot(cache, slot_a, a);
    // b 与 a 争用同一槽位且控制点不同时改用局部槽位，避免覆盖 a 的缓存。
    PairSampleSlot local_b;
    PairSampleSlot* eb = nullptr;
    if (slot_b == slot_a) {
        if (sameControlPoints(ea.key, b.segs)) {
            eb = &ea;
        } else {
            local_b.key = b.segs;
            local_b.bbox = b.bbox();
            local_b.valid = true;
            eb = &local_b;
        }
    } else {
        eb = &acquirePairSlot(cache, slot_b, b);
    }
    // 等价于 bboxOverlap(a, b)，只是包围盒来自缓存。
    if (!ea.bbox.intersects(eb->bbox))
        return false;
    const PairSampleIndex& ia = ensurePairSamples(ea, a);
    const PairSampleIndex& ib = ensurePairSamples(*eb, b);
    const std::vector<Vec2d>& pa = ia.pts;
    const std::vector<Vec2d>& pb = ib.pts;
    const int na = (int)pa.size() - 1;
    const int nb = (int)pb.size() - 1;
    if (na <= 0 || nb <= 0)
        return false;
    const int nblock = (int)ib.block_box.size();
    for (int i = 0; i < na; ++i) {
        // a 侧整块与 b 的总包围盒分离时，跳过该块的全部 a 段。
        if ((i % kPairBlockSize) == 0 &&
            !ia.block_box[i / kPairBlockSize].intersects(ib.whole)) {
            i += kPairBlockSize - 1;
            continue;
        }
        const BoundingBox2d& abox = ia.seg_box[i];
        if (!abox.intersects(ib.whole))
            continue;
        for (int blk = 0; blk < nblock; ++blk) {
            if (!abox.intersects(ib.block_box[blk]))
                continue;
            const int lo = blk * kPairBlockSize;
            const int hi = std::min(nb, lo + kPairBlockSize);
            for (int j = lo; j < hi; ++j) {
                if (!abox.intersects(ib.seg_box[j]))
                    continue;
                Vec2d isect;
                if (!segmentsIntersect(pa[i], pa[i + 1], pb[j], pb[j + 1], &isect))
                    continue;
                if (distToAllEndpoints(isect, a, b) <= ep)
                    continue;
                bool in_ball = false;
                for (const Vec2d& c : centers) {
                    if ((isect - c).norm() <= radius) {
                        in_ball = true;
                        break;
                    }
                }
                if (!in_ball)
                    return true;
            }
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

constexpr double kMinOrdinaryHandle = 0.05;

// 弦方向联合预算的判定松量：与 makeCubicG1 的裁剪共用同一比例常量，
// 保证「构造时裁剪到的形态」必然通过「闸门与审计的判定」。
// 松量的取值依据见 bezier.h 的 kChordProjectionBudgetSlack。
//
// 注意这是闸门的松量；修复路径 applyChordProjectionBudget 仍按等号裁剪，
// 即修复结果严格单调，比闸门更严。
constexpr double kChordBudgetSlackFraction = kChordProjectionBudgetSlack;

/// 弦方向联合预算的系数：`a = T0·u`、`b = T1·u`（u 为弦方向）。
/// 返回 false 表示不适用（弦退化，或某侧方向轴不朝弦前方）。
bool chordProjectionCoeffs(const Vec2d& p0, const Vec2d& start_dir,
                           const Vec2d& p1, const Vec2d& end_dir,
                           double* a, double* b, double* chord_len) {
    const Vec2d chord = p1 - p0;
    const double len = chord.norm();
    if (len < 1e-6) return false;
    const Vec2d u = chord / len;
    const double pa = start_dir.dot(u);
    const double pb = end_dir.dot(u);
    if (pa <= 1e-6 || pb <= 1e-6) return false;
    *a = pa;
    *b = pb;
    *chord_len = len;
    return true;
}

/// 超支时在最小把手之上等比缩放，保持两侧把手比例（拱顶位置）不变，且两侧
/// 仍不小于最小把手。缩放后恒有 `λ·a + μ·b == chord_len`。
void applyChordProjectionBudget(double a, double b, double chord_len,
                                double* lambda, double* mu) {
    const double used = *lambda * a + *mu * b;
    if (used <= chord_len) return;
    const double floor_used = kMinOrdinaryHandle * (a + b);
    if (floor_used >= chord_len) {
        *lambda = kMinOrdinaryHandle;
        *mu = kMinOrdinaryHandle;
        return;
    }
    const double scale = (chord_len - floor_used) / (used - floor_used);
    *lambda = kMinOrdinaryHandle + (*lambda - kMinOrdinaryHandle) * scale;
    *mu = kMinOrdinaryHandle + (*mu - kMinOrdinaryHandle) * scale;
}

OrdinarySingleCubicHandleBounds makeOrdinarySingleCubicHandleBounds(
    const Vec2d& p0, const Vec2d& start_tan,
    const Vec2d& p1, const Vec2d& end_tan,
    bool cap_at_direction_intersection) {
    OrdinarySingleCubicHandleBounds bounds;
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
    bounds.has_chord_budget = chordProjectionCoeffs(
        p0, bounds.start_dir, p1, bounds.end_dir, &bounds.start_chord_proj,
        &bounds.end_chord_proj, &bounds.chord_len);
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
        bounds.has_direction_intersection = true;
        bounds.direction_intersection_start = start_station;
        bounds.direction_intersection_end = end_station;
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

OrdinarySingleCubicHandleBounds ordinarySingleCubicHandleBounds(
    const Vec2d& p0, const Vec2d& start_tan,
    const Vec2d& p1, const Vec2d& end_tan,
    bool cap_at_direction_intersection) {
    return makeOrdinarySingleCubicHandleBounds(
        p0, start_tan, p1, end_tan, cap_at_direction_intersection);
}

bool cubicControlPolygonMonotone(const BezierSegment& segment, double tol) {
    const Vec2d start_handle = segment.ctrl[1] - segment.ctrl[0];
    const Vec2d end_handle = segment.ctrl[3] - segment.ctrl[2];
    if (start_handle.norm() < 1e-9 || end_handle.norm() < 1e-9)
        return true;
    double a = 0.0;
    double b = 0.0;
    double chord_len = 0.0;
    if (!chordProjectionCoeffs(segment.ctrl[0], start_handle.normalized(),
                               segment.ctrl[3], end_handle.normalized(),
                               &a, &b, &chord_len))
        return true;
    return start_handle.norm() * a + end_handle.norm() * b <=
           chord_len * (1.0 + kChordBudgetSlackFraction) + tol;
}

double curveChordBudgetOvershoot(const BezierCurve& curve) {
    double worst = 0.0;
    for (const BezierSegment& segment : curve.segs) {
        const Vec2d start_handle = segment.ctrl[1] - segment.ctrl[0];
        const Vec2d end_handle = segment.ctrl[3] - segment.ctrl[2];
        if (start_handle.norm() < 1e-9 || end_handle.norm() < 1e-9)
            continue;
        double a = 0.0;
        double b = 0.0;
        double chord_len = 0.0;
        if (!chordProjectionCoeffs(segment.ctrl[0], start_handle.normalized(),
                                   segment.ctrl[3], end_handle.normalized(),
                                   &a, &b, &chord_len))
            continue;
        const double used = start_handle.norm() * a + end_handle.norm() * b;
        const double allowed = chord_len * (1.0 + kChordBudgetSlackFraction);
        if (used > allowed)
            worst = std::max(worst, (used - allowed) / chord_len);
    }
    return worst;
}

void constrainOrdinarySingleCubicControls(
    BezierCurve& curve, const Vec2d& p0, const Vec2d& start_tan,
    const Vec2d& p1, const Vec2d& end_tan,
    bool cap_at_direction_intersection) {
    if (curve.numSegments() != 1)
        return;

    const OrdinarySingleCubicHandleBounds bounds = ordinarySingleCubicHandleBounds(
        p0, start_tan, p1, end_tan, cap_at_direction_intersection);
    BezierSegment& seg = curve.segs.front();
    seg.ctrl[0] = p0;
    seg.ctrl[3] = p1;

    double lambda = std::max(
        0.05, std::min(bounds.start_max,
                       startHandleStation(seg.ctrl[1], p0, bounds.start_dir)));
    double mu = std::max(
        0.05, std::min(bounds.end_max,
                       (p1 - seg.ctrl[2]).dot(bounds.end_dir)));
    if (bounds.has_chord_budget)
        applyChordProjectionBudget(bounds.start_chord_proj,
                                   bounds.end_chord_proj, bounds.chord_len,
                                   &lambda, &mu);
    seg.ctrl[1] = p0 + lambda * bounds.start_dir;
    seg.ctrl[2] = p1 - mu * bounds.end_dir;
}

bool ordinarySingleCubicControlsValid(
    const BezierCurve& curve, const Vec2d& p0, const Vec2d& start_tan,
    const Vec2d& p1, const Vec2d& end_tan, double tol,
    bool cap_at_direction_intersection) {
    if (curve.numSegments() != 1)
        return true;

    const OrdinarySingleCubicHandleBounds bounds = ordinarySingleCubicHandleBounds(
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
    // 弦方向联合预算：控制多边形必须沿弦单调（允许 kChordBudgetSlackFraction
    // 的比例松量），否则曲线被折成 Z / L 形。
    if (bounds.has_chord_budget &&
        lambda * bounds.start_chord_proj + mu * bounds.end_chord_proj >
            bounds.chord_len * (1.0 + kChordBudgetSlackFraction) + tol)
        return false;
    return start_lateral <= tol && end_lateral <= tol &&
           lambda >= 0.05 - tol && lambda <= bounds.start_max + tol &&
           mu >= 0.05 - tol && mu <= bounds.end_max + tol;
}

bool sharedEndpointMidControlSegmentsCross(
    const BezierCurve& a, const BezierCurve& b) {
    if (a.numSegments() != 1 || b.numSegments() != 1)
        return false;
    const std::array<Vec2d, 4>& ga = a.segs.front().ctrl;
    const std::array<Vec2d, 4>& gb = b.segs.front().ctrl;
    return segmentsIntersect(ga[1], ga[2], gb[1], gb[2]);
}

bool sharedEndpointControlPolylinesCross(
    const BezierCurve& a, const BezierCurve& b, double ctrl_tol) {
    if (a.numSegments() != 1 || b.numSegments() != 1)
        return false;
    const std::array<Vec2d, 4>& ga = a.segs.front().ctrl;
    const std::array<Vec2d, 4>& gb = b.segs.front().ctrl;
    const double tol = std::max(1e-9, ctrl_tol);
    for (int i = 0; i + 1 < 4; ++i) {
        for (int j = 0; j + 1 < 4; ++j) {
            Vec2d r = ga[i + 1] - ga[i];
            Vec2d s = gb[j + 1] - gb[j];
            const double rn = r.norm(), sn = s.norm();
            // 退化把手（控制点重合）不构成折线段。
            if (rn < tol || sn < tol)
                continue;
            // 共线/平行视为重合接触，不算违约：共享端点侧的首把手沿同一
            // 切向天然共线重叠。
            if (std::abs(cross2d(r, s)) <= 1e-9 * rn * sn)
                continue;
            Vec2d ipt;
            if (!segmentsIntersect(ga[i], ga[i + 1], gb[j], gb[j + 1], &ipt))
                continue;
            // 落在任一控制点邻域内的相接是合法的（共享端点即属此类）。
            bool at_control_point = false;
            for (int k = 0; k < 4 && !at_control_point; ++k) {
                if ((ipt - ga[k]).norm() <= tol || (ipt - gb[k]).norm() <= tol)
                    at_control_point = true;
            }
            if (at_control_point)
                continue;
            return true;
        }
    }
    return false;
}

double singleCubicSignedEndCurvature(const BezierCurve& c, bool at_start) {
    if (c.empty())
        return 0.0;
    const std::array<Vec2d, 4>& g = at_start
        ? c.segs.front().ctrl : c.segs.back().ctrl;
    // 端点处的有符号曲率 kappa = 2/3 * cross(d1, d2) / |d1|^3，
    // 其中 d1 是端点方向的把手向量，d2 是相邻的控制点差向量。
    // 尾端点侧把方向整体反向（沿曲线倒序看），叉积随之取反，
    // 保证"值越大越偏向该端点切向的左侧"这一语义在两端一致。
    const Vec2d d1 = at_start ? (g[1] - g[0]) : (g[2] - g[3]);
    const Vec2d d2 = at_start ? (g[2] - g[1]) : (g[1] - g[2]);
    const double h = d1.norm();
    if (h < 1e-9)
        return 0.0;
    return (2.0 / 3.0) * cross2d(d1, d2) / (h * h * h);
}

namespace {

// 采样折线上点到另一条折线的最近距离。
double distToPolyline(const Vec2d& p, const std::vector<Vec2d>& poly) {
    double best = std::numeric_limits<double>::infinity();
    for (std::size_t i = 0; i + 1 < poly.size(); ++i) {
        const Vec2d ab = poly[i + 1] - poly[i];
        const double len2 = ab.squaredNorm();
        double t = len2 > 1e-18 ? (p - poly[i]).dot(ab) / len2 : 0.0;
        t = std::max(0.0, std::min(1.0, t));
        best = std::min(best, (p - (poly[i] + t * ab)).norm());
    }
    return best;
}

// 自 `poly` 的一端（`from_start` 决定哪一端）沿弧长前进，累计"到 other 的距离
// 始终 < band"的长度；一旦超出 band 立即停止并返回已累计长度（截断到 cap）。
double funnelRunFrom(const std::vector<Vec2d>& poly,
                     const std::vector<Vec2d>& other,
                     bool from_start, double band, double cap) {
    double run = 0.0;
    const std::size_t n = poly.size();
    for (std::size_t k = 0; k + 1 < n; ++k) {
        const std::size_t i = from_start ? k : n - 1 - k;
        const std::size_t j = from_start ? k + 1 : n - 2 - k;
        if (distToPolyline(poly[j], other) >= band)
            break;
        run += (poly[j] - poly[i]).norm();
        if (run >= cap)
            return cap;
    }
    return run;
}

}  // namespace

std::vector<Vec2d> sharedEndpointsOf(const BezierCurve& a, const BezierCurve& b,
                                     double endpoint_tol) {
    std::vector<Vec2d> shared;
    if (a.empty() || b.empty())
        return shared;
    const Vec2d ends_a[2] = {a.startPt(), a.endPt()};
    const Vec2d ends_b[2] = {b.startPt(), b.endPt()};
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            if ((ends_a[i] - ends_b[j]).norm() > endpoint_tol)
                continue;
            const Vec2d mid = 0.5 * (ends_a[i] + ends_b[j]);
            bool dup = false;
            for (const Vec2d& s : shared)
                if ((s - mid).norm() <= endpoint_tol)
                    dup = true;
            if (!dup)
                shared.push_back(mid);
        }
    }
    return shared;
}

double sharedEndpointMergeFunnelRadius(const BezierCurve& a,
                                       const BezierCurve& b,
                                       double band, double cap,
                                       double endpoint_tol) {
    if (a.empty() || b.empty() || band <= 0.0 || cap <= 0.0)
        return 0.0;
    const int n = 128;
    const std::vector<Vec2d> pa = a.sample(n);
    const std::vector<Vec2d> pb = b.sample(n);
    if (pa.size() < 2 || pb.size() < 2)
        return 0.0;
    double radius = 0.0;
    for (int ia = 0; ia < 2; ++ia) {
        const bool a_start = ia == 0;
        const Vec2d ea = a_start ? a.startPt() : a.endPt();
        for (int ib = 0; ib < 2; ++ib) {
            const bool b_start = ib == 0;
            const Vec2d eb = b_start ? b.startPt() : b.endPt();
            if ((ea - eb).norm() > endpoint_tol)
                continue;
            // 两条曲线各自的收敛段取小：任一条先离开贴近带，收敛段即结束。
            const double run = std::min(
                funnelRunFrom(pa, pb, a_start, band, cap),
                funnelRunFrom(pb, pa, b_start, band, cap));
            radius = std::max(radius, run);
        }
    }
    return radius;
}

double curveTurningSpan(const BezierCurve& curve, int samples_per_seg) {
    if (curve.empty() || samples_per_seg < 1)
        return 0.0;
    std::vector<Vec2d> tans;
    tans.reserve(curve.segs.size() * (size_t)(samples_per_seg + 1));
    for (const auto& seg : curve.segs)
        for (int i = 0; i <= samples_per_seg; ++i) {
            Vec2d d = seg.evalDeriv1((double)i / (double)samples_per_seg);
            if (d.norm() > 1e-9)
                tans.push_back(d.normalized());
        }
    if (tans.size() < 2)
        return 0.0;
    // running 有符号累加，跨度取 max-min：同向连续绕转会把两端拉开，
    // S 形的两段反向弯则互相抵消，只剩较大的单侧弯角。
    double running = 0.0;
    double lo = 0.0;
    double hi = 0.0;
    for (size_t i = 1; i < tans.size(); ++i) {
        running += std::atan2(
            cross2d(tans[i - 1], tans[i]), tans[i - 1].dot(tans[i]));
        lo = std::min(lo, running);
        hi = std::max(hi, running);
    }
    return hi - lo;
}

namespace {

// arc/chord 下限的公共口径：documented 是需求按 90° 折算出来的业务常数，
// 它只作上界；真正的阈值是"同转角圆弧比例的固定百分比"，两者取小。
double arcChordFloorWithCap(double turn_angle_rad, double chord_len,
                            double documented) {
    constexpr double kAbsoluteFloor = 1.005;
    const double theta = std::abs(turn_angle_rad);
    if (!std::isfinite(theta) || theta < 1e-6)
        return kAbsoluteFloor;
    // 同转角圆弧的 arc/chord = θ / (2 sin(θ/2))。θ >= π 时（几何掉头）该式
    // 退化，掉头由 evaluateUTurnShape 单独裁定，这里按 π 截断即可。
    const double clamped = std::min(theta, M_PI - 1e-9);
    const double ideal = clamped / (2.0 * std::sin(0.5 * clamped));
    // 长弦取圆弧的 97%、短弦取 93%：θ=90° 时分别还原 1.08 / 1.03 这两个
    // 需求写明的"长距离/短距离"数字。
    const double fraction = chord_len < 12.0 ? 0.93 : 0.97;
    return std::max(kAbsoluteFloor, std::min(documented, fraction * ideal));
}

}  // namespace

double ordinaryTurnArcChordFloor(double turn_angle_rad, double chord_len) {
    return arcChordFloorWithCap(turn_angle_rad, chord_len,
                                chord_len < 12.0 ? 1.02 : 1.06);
}

double ordinaryTurnArcChordRestoreFloor(double turn_angle_rad,
                                        double chord_len) {
    return arcChordFloorWithCap(turn_angle_rad, chord_len,
                                chord_len < 12.0 ? 1.03 : 1.08);
}

double curveEndpointTurnAngle(const BezierCurve& curve) {
    if (curve.empty())
        return 0.0;
    const Vec2d t0 = curve.startTan();
    const Vec2d t1 = curve.endTan();
    if (t0.norm() < 1e-9 || t1.norm() < 1e-9)
        return 0.0;
    const Vec2d a = t0.normalized();
    const Vec2d b = t1.normalized();
    return std::atan2(std::abs(cross2d(a, b)), a.dot(b));
}

}
