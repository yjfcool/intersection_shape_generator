#pragma once

#include "types.h"
#include "curve/curve_utils.h"
#include "utils.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace isg {

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
        const std::vector<Boundary>& boundaries, const Vec2d& center) {
    std::vector<BoundarySafetySegment> segments;
    for (const auto& bnd : boundaries) {
        const auto& pts = bnd.geometry.points;
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
            if (seg.road_edge && seg.center_signed > 1e-6)
                has_pos_center_side = true;
            if (seg.road_edge && seg.center_signed < -1e-6)
                has_neg_center_side = true;
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

} // 命名空间 isg
