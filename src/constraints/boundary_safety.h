#pragma once

#include "types.h"
#include "curve/curve_utils.h"
#include "utils.h"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <vector>

namespace isg {

// 临时对照开关：设置 ISG_LEGACY_SHAPE=1 时，端点擦碰豁免与转角超量门禁全部退回
// 旧行为，用于在同一个二进制里对比回归基线。
inline bool isgLegacyShapeGates() {
    static const bool legacy = std::getenv("ISG_LEGACY_SHAPE") != nullptr;
    return legacy;
}

// 临时归因开关：单独退回三项改动之一，用于定位回归来源。
inline bool isgLegacyGraze() {
    static const bool legacy = std::getenv("ISG_LEGACY_GRAZE") != nullptr;
    return legacy || isgLegacyShapeGates();
}
inline bool isgLegacyExcess() {
    static const bool legacy = std::getenv("ISG_LEGACY_EXCESS") != nullptr;
    return legacy || isgLegacyShapeGates();
}
inline bool isgLegacySegTie() {
    static const bool legacy = std::getenv("ISG_LEGACY_SEGTIE") != nullptr;
    return legacy || isgLegacyShapeGates();
}

struct BoundarySafetySegment {
    Vec2d a{0, 0};
    Vec2d b{0, 0};
    double len = 0.0;
    double center_signed = 0.0;
    bool road_edge = false;
    bool center_side_reliable = true;
    bool a_is_boundary_endpoint = false;
    bool b_is_boundary_endpoint = false;
};

struct BoundarySafetyResult {
    bool intersects = false;
    bool outside_road_edge = false;
    double outside_penalty = 0.0;
};

template <typename SegmentT>
inline bool roadEdgeAdherentContactAllowed(
        const Vec2d& curve_a, const Vec2d& curve_b,
        const SegmentT& boundary_seg, double min_abs_dot = 0.92,
        double adherent_tol = 0.12) {
    if (!boundary_seg.road_edge)
        return false;
    Vec2d cd = curve_b - curve_a;
    Vec2d bd = boundary_seg.b - boundary_seg.a;
    if (cd.norm() < 1e-8 || bd.norm() < 1e-8)
        return false;
    if (std::abs(cd.normalized().dot(bd.normalized())) < min_abs_dot)
        return false;
    return pointToSegment(curve_a, boundary_seg.a, boundary_seg.b).first <= adherent_tol &&
           pointToSegment(curve_b, boundary_seg.a, boundary_seg.b).first <= adherent_tol;
}

template <typename SegmentT>
inline double signedDistanceToBoundarySegmentLine(
        const Vec2d& pt, const SegmentT& boundary_seg) {
    Vec2d bd = boundary_seg.b - boundary_seg.a;
    double len = bd.norm();
    if (len < 1e-8)
        return 0.0;
    return cross2d(bd, pt - boundary_seg.a) / len;
}

template <typename SegmentT>
inline bool roadEdgeVertexTangentialContactAllowed(
        const std::vector<Vec2d>& pts, int seg_index, const Vec2d& isect,
        const SegmentT& boundary_seg, double min_abs_dot = 0.90,
        double curve_cut_tol = 0.12, double adherent_tol = 0.12,
        double side_tol = 1e-3) {
    if (!boundary_seg.road_edge)
        return false;
    Vec2d bd = boundary_seg.b - boundary_seg.a;
    if (bd.norm() < 1e-8)
        return false;

    int vertex = -1;
    if ((isect - pts[seg_index]).norm() <= curve_cut_tol)
        vertex = seg_index;
    if ((isect - pts[seg_index + 1]).norm() <= curve_cut_tol)
        vertex = seg_index + 1;
    if (vertex < 0)
        return false;

    bool prev_adherent = vertex > 0 &&
        roadEdgeAdherentContactAllowed(
            pts[vertex - 1], pts[vertex], boundary_seg, min_abs_dot, adherent_tol);
    bool next_adherent = vertex + 1 < (int)pts.size() &&
        roadEdgeAdherentContactAllowed(
            pts[vertex], pts[vertex + 1], boundary_seg, min_abs_dot, adherent_tol);

    int prev_idx = -1;
    int next_idx = -1;
    double prev_signed = 0.0;
    double next_signed = 0.0;
    for (int k = vertex - 1; k >= 0; --k) {
        double s = signedDistanceToBoundarySegmentLine(pts[k], boundary_seg);
        if (std::abs(s) > side_tol) {
            prev_idx = k;
            prev_signed = s;
            break;
        }
    }
    for (int k = vertex + 1; k < (int)pts.size(); ++k) {
        double s = signedDistanceToBoundarySegmentLine(pts[k], boundary_seg);
        if (std::abs(s) > side_tol) {
            next_idx = k;
            next_signed = s;
            break;
        }
    }

    bool has_parallel_leg = false;
    if (vertex > 0) {
        Vec2d d = pts[vertex] - pts[vertex - 1];
        has_parallel_leg = has_parallel_leg ||
            (d.norm() > 1e-8 && std::abs(d.normalized().dot(bd.normalized())) >= min_abs_dot);
    }
    if (vertex + 1 < (int)pts.size()) {
        Vec2d d = pts[vertex + 1] - pts[vertex];
        has_parallel_leg = has_parallel_leg ||
            (d.norm() > 1e-8 && std::abs(d.normalized().dot(bd.normalized())) >= min_abs_dot);
    }
    if (!has_parallel_leg && !prev_adherent && !next_adherent)
        return false;

    if (prev_idx >= 0 && next_idx >= 0 &&
        prev_signed * next_signed < -side_tol * side_tol) {
        bool adherent_join = prev_adherent || next_adherent;
        if (!adherent_join)
            return false;
    }

    if (prev_idx < 0 || next_idx < 0)
        return prev_adherent || next_adherent;

    return pointToSegment(isect, boundary_seg.a, boundary_seg.b).first <=
           adherent_tol;
}

inline Vec2d boundarySafetyCenter(const IntersectionInput& input) {
    Vec2d center(0, 0);
    int count = 0;
    for (const auto& conn : input.connectivities) {
        auto entry = input.entryPtDir(conn.entry_lane_id);
        auto exit_ = input.exitPtDir(conn.exit_lane_id);
        center += entry.first;
        center += exit_.first;
        count += 2;
    }
    if (count > 0)
        return center / (double)count;

    for (const auto& lane : input.lanes) {
        for (const auto& pt : lane.geometry.points) {
            center += xyOf(pt);
            ++count;
        }
    }
    if (count > 0)
        return center / (double)count;

    BoundingBox2d box;
    for (const auto& bnd : input.boundaries)
        for (const auto& pt : bnd.geometry.points)
            box.expand(pt);
    if (!box.empty())
        return 0.5 * (box.min_pt + box.max_pt);
    return center;
}

inline std::vector<BoundarySafetySegment> buildBoundarySafetySegments(
        const std::vector<Boundary>& boundaries, const Vec2d& center,
        const std::vector<std::vector<bool>>* exempt_masks = nullptr) {
    std::vector<BoundarySafetySegment> segments;
    for (size_t bi = 0; bi < boundaries.size(); ++bi) {
        const auto& bnd = boundaries[bi];
        const auto& pts = bnd.geometry.points;
        const std::vector<bool>* mask =
            (exempt_masks && bi < exempt_masks->size() && !(*exempt_masks)[bi].empty())
                ? &(*exempt_masks)[bi] : nullptr;
        std::vector<BoundarySafetySegment> boundary_segments;
        bool has_pos_center_side = false;
        bool has_neg_center_side = false;
        for (int i = 0; i + 1 < (int)pts.size(); ++i) {
            Vec2d d = pts[i + 1] - pts[i];
            double len = d.norm();
            if (len < 1e-8)
                continue;
            BoundarySafetySegment seg;
            seg.a = pts[i];
            seg.b = pts[i + 1];
            seg.len = len;
            seg.center_signed = cross2d(d, center - pts[i]) / len;
            seg.road_edge = (bnd.type == Boundary::Type::RoadEdge);
            // 中心侧可靠性必须在完整折线上统计，豁免段只是不参与判违，
            // 不能改变该 boundary 的中心侧一致性结论。
            if (seg.road_edge && seg.center_signed > 1e-6)
                has_pos_center_side = true;
            if (seg.road_edge && seg.center_signed < -1e-6)
                has_neg_center_side = true;
            if (mask && i < (int)mask->size() && (*mask)[i])
                continue;
            seg.a_is_boundary_endpoint = (i == 0);
            seg.b_is_boundary_endpoint = (i + 1 == (int)pts.size() - 1);
            boundary_segments.push_back(seg);
        }
        bool reliable = !(bnd.type == Boundary::Type::RoadEdge &&
                          has_pos_center_side && has_neg_center_side);
        for (auto& seg : boundary_segments) {
            seg.center_side_reliable = reliable;
            segments.push_back(seg);
        }
    }
    return segments;
}

// ── 端点贴合 RoadEdge 的共线擦碰豁免 ─────────────────────────────────────
// 背景：路口臂的 RoadEdge 折线端点常常正好落在车道端点上（中央分隔带鼻端、
// 渠化岛边缘、非机动车道外侧路缘）。连接曲线必须从该点出发或到达，于是在与
// 折线共线的那一小段上以厘米级横向摆动跨过折线一到两次。这是共享端点造成的
// 数值退化，不是切进路缘。
//
// 豁免判据（三条必须同时成立）：
//   1. boundary 折线的首点或末点落在曲线端点 pos_tol 内；
//   2. 从该端点沿折线走出的连续"共线段区间"覆盖了违规接触。区间的走法是：
//      逐段比较段方向与该端点侧车道切向，|dot| >= dir_min_dot 则纳入，遇到
//      第一个不共线的段立即停止。发夹形鼻端折回的另一条腿因此留在区间外，
//      仍按穿越判定——这正是整条折线首尾弦方向判据（首尾弦近乎垂直于两条腿）
//      失效的地方。
//   3. 曲线从该端点到区间内最远交点之间，相对区间折线的最大横向偏移不超过
//      graze_tol。真正切进路缘会偏离数米，远超该阈值。
inline bool boundarySegmentAlignsWithTangent(
        const Vec2d& a, const Vec2d& b, const Vec2d& tangent, double min_dot) {
    Vec2d d = b - a;
    if (d.norm() < 1e-9 || tangent.norm() < 1e-9)
        return false;
    return std::abs(d.normalized().dot(tangent.normalized())) >= min_dot;
}

// 从折线的首端（tip_at_front）或末端向内标记连续共线段，遇首个不共线段停止。
inline void markBoundaryTipAlignedRun(
        const std::vector<Vec2d>& bpts, bool tip_at_front, const Vec2d& tangent,
        double min_dot, std::vector<bool>& mask) {
    const int nseg = (int)bpts.size() - 1;
    if (nseg <= 0 || (int)mask.size() < nseg)
        return;
    if (tip_at_front) {
        for (int s = 0; s < nseg; ++s) {
            if (!boundarySegmentAlignsWithTangent(bpts[s], bpts[s + 1], tangent, min_dot))
                break;
            mask[s] = true;
        }
    } else {
        for (int s = nseg - 1; s >= 0; --s) {
            if (!boundarySegmentAlignsWithTangent(bpts[s], bpts[s + 1], tangent, min_dot))
                break;
            mask[s] = true;
        }
    }
}

inline double roadEdgeOutsidePenalty(
        const Vec2d& pt, const std::vector<BoundarySafetySegment>& segments,
        double tol = 0.05, double max_influence_dist = 3.0) {
    const BoundarySafetySegment* best = nullptr;
    double best_dist = std::numeric_limits<double>::infinity();
    double best_signed = 0.0;
    for (const auto& seg : segments) {
        if (!seg.road_edge || !seg.center_side_reliable ||
            std::abs(seg.center_signed) < 1e-6)
            continue;
        Vec2d d = seg.b - seg.a;
        double len2 = d.squaredNorm();
        if (len2 < 1e-12)
            continue;
        double t = std::max(0.0, std::min(1.0, (pt - seg.a).dot(d) / len2));
        if (t <= 0.02 || t >= 0.98)
            continue;
        Vec2d closest = seg.a + t * d;
        double dist2 = (pt - closest).squaredNorm();
        if (dist2 < best_dist) {
            best_dist = dist2;
            best = &seg;
            best_signed = cross2d(d, pt - seg.a) / seg.len;
        }
    }
    if (!best || best_dist > max_influence_dist * max_influence_dist)
        return 0.0;
    double center_sign = best->center_signed > 0.0 ? 1.0 : -1.0;
    return std::max(0.0, -center_sign * best_signed - tol);
}

inline BoundarySafetyResult curveBoundarySafety(
        const BezierCurve& curve, const std::vector<BoundarySafetySegment>& segments,
        int min_samples = 64, double curve_endpoint_tol = 0.75,
        double boundary_endpoint_tol = 0.10, double outside_tol = 0.05) {
    BoundarySafetyResult result;
    if (segments.empty() || curve.empty())
        return result;

    std::vector<Vec2d> pts;
    pts.reserve(std::max(min_samples, 2));
    for (int si = 0; si < (int)curve.segs.size(); ++si) {
        const auto& curve_seg = curve.segs[si];
        int seg_samples = std::max(
            2, (int)std::ceil(curve_seg.arcLength(12) / 0.18) + 1);
        seg_samples = std::min(seg_samples, 80);
        for (int i = 0; i < seg_samples; ++i) {
            if (si > 0 && i == 0)
                continue;
            double t = seg_samples == 1 ? 0.0 : (double)i / (double)(seg_samples - 1);
            pts.push_back(curve_seg.evaluate(t));
        }
    }
    if (pts.size() < 2)
        return result;

    BoundingBox2d curve_box;
    for (const auto& pt : pts)
        curve_box.expand(pt);

    for (const auto& seg : segments) {
        BoundingBox2d seg_box;
        seg_box.expand(seg.a);
        seg_box.expand(seg.b);
        if (!curve_box.intersects(seg_box))
            continue;
        for (int i = 0; i + 1 < (int)pts.size(); ++i) {
            Vec2d isect;
            if (!segmentsIntersect(pts[i], pts[i + 1], seg.a, seg.b, &isect))
                continue;
            if ((isect - pts.front()).norm() <= curve_endpoint_tol ||
                (isect - pts.back()).norm() <= curve_endpoint_tol)
                continue;
            if ((seg.a_is_boundary_endpoint && (isect - seg.a).norm() <= boundary_endpoint_tol) ||
                (seg.b_is_boundary_endpoint && (isect - seg.b).norm() <= boundary_endpoint_tol))
                continue;
            result.intersects = true;
            break;
        }
        if (result.intersects)
            break;
    }

    for (int i = 1; i + 1 < (int)pts.size(); ++i) {
        if ((pts[i] - pts.front()).norm() <= curve_endpoint_tol ||
            (pts[i] - pts.back()).norm() <= curve_endpoint_tol)
            continue;
        double p = roadEdgeOutsidePenalty(pts[i], segments, outside_tol);
        if (p > 1e-3) {
            result.outside_road_edge = true;
            result.intersects = true;
            result.outside_penalty += p * p;
        }
    }
    return result;
}

inline BoundarySafetyResult curveBoundarySafety(
        const BezierCurve& curve, const std::vector<Boundary>& boundaries,
        const Vec2d& center, int min_samples = 64,
        double curve_endpoint_tol = 0.75,
        double boundary_endpoint_tol = 0.10,
        double outside_tol = 0.05) {
    auto all_segments = buildBoundarySafetySegments(boundaries, center);
    // 候选检查通常只涉及密集 Boundary 集合中的局部区域。
    // 保留 3 米缓冲带，因为 roadEdgeOutsidePenalty 会使用该范围；
    // 先过滤可避免每个采样点都扫描完整集合，同时保持结果不变。
    std::vector<BoundarySafetySegment> segments;
    if (!curve.empty() && !all_segments.empty()) {
        BoundingBox2d query_box = curve.bbox();
        query_box.min_pt -= Vec2d(3.05, 3.05);
        query_box.max_pt += Vec2d(3.05, 3.05);
        segments.reserve(all_segments.size());
        for (const auto& seg : all_segments) {
            BoundingBox2d seg_box;
            seg_box.expand(seg.a);
            seg_box.expand(seg.b);
            if (query_box.intersects(seg_box))
                segments.push_back(seg);
        }
    } else {
        segments = std::move(all_segments);
    }
    return curveBoundarySafety(
        curve, segments, min_samples, curve_endpoint_tol,
        boundary_endpoint_tol, outside_tol);
}

// 曲线相对折线中被标记区间的最大横向偏移（"擦碰深度"）。
// 采样范围限定在贴合端点与区间内最远交点之间；区间内无交点时返回 0。
inline double endpointGrazeDepth(
        const std::vector<Vec2d>& pts, bool from_start,
        const std::vector<Vec2d>& bpts, const std::vector<bool>& mask) {
    const int nseg = (int)bpts.size() - 1;
    int far = -1;
    for (int i = 0; i + 1 < (int)pts.size(); ++i)
        for (int j = 0; j < nseg; ++j) {
            if (j >= (int)mask.size() || !mask[j])
                continue;
            if (!segmentsIntersect(pts[i], pts[i + 1], bpts[j], bpts[j + 1], nullptr))
                continue;
            const int idx = from_start ? i + 1 : i;
            if (far < 0 || (from_start ? idx > far : idx < far))
                far = idx;
        }
    if (far < 0)
        return 0.0;
    const int lo = from_start ? 0 : far;
    const int hi = from_start ? far : (int)pts.size() - 1;
    double depth = 0.0;
    for (int i = lo; i <= hi; ++i) {
        double d = std::numeric_limits<double>::infinity();
        for (int j = 0; j < nseg; ++j)
            if (j < (int)mask.size() && mask[j])
                d = std::min(d, pointToSegment(pts[i], bpts[j], bpts[j + 1]).first);
        if (std::isfinite(d))
            depth = std::max(depth, d);
    }
    return depth;
}

// 为每条 boundary 计算可豁免的段掩码；未命中豁免判据的 boundary 掩码为空。
// 切向直接取自曲线端点，因此对已修复的多段曲线同样成立（端点 G1 保持不变）。
// 曲线采样是延迟的：只有确实存在"端点贴合"的 boundary 时才做，避免在没有贴合
// 关系的绝大多数调用上白付采样代价。
inline std::vector<std::vector<bool>> endpointGrazeExemptSegmentMasks(
        const BezierCurve& curve, const std::vector<Boundary>& boundaries,
        double pos_tol = 0.60, double dir_min_dot = 0.70, double graze_tol = 0.25) {
    std::vector<std::vector<bool>> masks(boundaries.size());
    if (curve.empty() || boundaries.empty())
        return masks;
    Vec2d t0 = curve.startTan();
    Vec2d t1 = curve.endTan();
    if (t0.norm() < 1e-9 || t1.norm() < 1e-9)
        return masks;
    t0.normalize();
    t1.normalize();
    const Vec2d p0 = curve.startPt();
    const Vec2d p1 = curve.endPt();
    std::vector<Vec2d> pts;

    for (size_t bi = 0; bi < boundaries.size(); ++bi) {
        const auto& bgeom = boundaries[bi].geometry.points;
        if (bgeom.size() < 2)
            continue;
        const Vec2d bfront = xyOf(bgeom.front());
        const Vec2d bback = xyOf(bgeom.back());
        if ((bfront - p0).norm() > pos_tol && (bback - p0).norm() > pos_tol &&
            (bfront - p1).norm() > pos_tol && (bback - p1).norm() > pos_tol)
            continue;
        if (pts.empty()) {
            pts = curve.sampleByArcLength(std::max(
                64, std::min(240, (int)std::ceil(curve.arcLength() / 0.18) + 1)));
            if (pts.size() < 2)
                return masks;
        }
        const std::vector<Vec2d> bpts = toVec2dArray(bgeom);
        const int nseg = (int)bpts.size() - 1;
        for (int side = 0; side < 2; ++side) {
            const Vec2d& cp = side == 0 ? p0 : p1;
            const Vec2d& ct = side == 0 ? t0 : t1;
            const bool tip_front = (bfront - cp).norm() <= pos_tol;
            const bool tip_back = (bback - cp).norm() <= pos_tol;
            if (!tip_front && !tip_back)
                continue;
            std::vector<bool> run(nseg, false);
            if (tip_front)
                markBoundaryTipAlignedRun(bpts, true, ct, dir_min_dot, run);
            if (tip_back)
                markBoundaryTipAlignedRun(bpts, false, ct, dir_min_dot, run);
            if (std::find(run.begin(), run.end(), true) == run.end())
                continue;
            if (endpointGrazeDepth(pts, side == 0, bpts, run) > graze_tol)
                continue;
            if (masks[bi].empty())
                masks[bi].assign(nseg, false);
            for (int s = 0; s < nseg; ++s)
                if (run[s])
                    masks[bi][s] = true;
        }
    }
    return masks;
}

// 掩码集合是否为空（没有任何段被豁免）。空掩码意味着豁免判定不会改变结论，
// 调用方可以直接沿用严格判定结果。
inline bool boundaryExemptMasksEmpty(const std::vector<std::vector<bool>>& masks) {
    for (const auto& m : masks)
        if (!m.empty())
            return false;
    return true;
}

// 忽略端点贴合共线擦碰后重新判定边界安全性。豁免段完全退出判违集合，
// 其余段（含发夹形鼻端折回的另一条腿）仍按原规则参与穿越与路缘外侧判定。
//
// 性能：先跑一次严格判定；只有严格判定报违时才计算豁免掩码并复判。绝大多数
// 候选曲线本来就不碰边界，因此这条快路径让豁免逻辑几乎不产生额外开销。
inline BoundarySafetyResult curveBoundarySafetyIgnoringEndpointGraze(
        const BezierCurve& curve, const std::vector<Boundary>& boundaries,
        const Vec2d& center, int min_samples = 64,
        double curve_endpoint_tol = 0.75, double boundary_endpoint_tol = 0.10,
        double outside_tol = 0.05, double pos_tol = 0.60,
        double dir_min_dot = 0.70, double graze_tol = 0.25) {
    const BoundarySafetyResult strict = curveBoundarySafety(
        curve, boundaries, center, min_samples, curve_endpoint_tol,
        boundary_endpoint_tol, outside_tol);
    if (isgLegacyGraze() || (!strict.intersects && !strict.outside_road_edge))
        return strict;
    const std::vector<std::vector<bool>> masks = endpointGrazeExemptSegmentMasks(
        curve, boundaries, pos_tol, dir_min_dot, graze_tol);
    if (boundaryExemptMasksEmpty(masks))
        return strict;
    std::vector<BoundarySafetySegment> segments =
        buildBoundarySafetySegments(boundaries, center, &masks);
    return curveBoundarySafety(curve, segments, min_samples, curve_endpoint_tol,
                              boundary_endpoint_tol, outside_tol);
}

} // 命名空间 isg
