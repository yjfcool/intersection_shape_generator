#pragma once

#include "types.h"
#include "curve/curve_utils.h"
#include "utils.h"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <vector>

namespace isg {

// 避让约束唯一允许的端点误差。该误差只用于识别曲线自身的真实首点/尾点，
// 不用于放宽 Boundary 端点，也不允许把端点后的贴合、重叠或穿越隐藏掉。
static constexpr double kConnectionPointTolerance = 1e-4;

// 临时对照开关：设置 ISG_LEGACY_SHAPE=1 时，端点擦碰豁免与转角超量门禁全部退回
// 旧行为，用于在同一个二进制里对比回归基线。
inline bool isgLegacyShapeGates() {
    static const bool legacy = std::getenv("ISG_LEGACY_SHAPE") != nullptr;
    return legacy;
}

// 临时归因开关：单独退回三项改动之一，用于定位回归来源。
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
    BoundingBox2d bbox;
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

// 判断曲线采样线段与 Boundary 线段的接触是否超出曲线真实首/尾连接点。
// 非平行线段只有一个交点；共线时必须检查完整重叠区间，否则仅用
// segmentsIntersect 返回的中点可能把“从首点开始的贴合/重叠”错误当成合法连接。
inline bool segmentHasForbiddenBoundaryContact(
        const Vec2d& curve_a, const Vec2d& curve_b,
        const Vec2d& boundary_a, const Vec2d& boundary_b,
        const Vec2d& curve_start, const Vec2d& curve_end,
        double endpoint_tol = kConnectionPointTolerance) {
    const double tol = std::min(
        std::max(0.0, endpoint_tol), kConnectionPointTolerance);
    Vec2d curve_delta = curve_b - curve_a;
    Vec2d boundary_delta = boundary_b - boundary_a;
    if (curve_delta.norm() < 1e-10 || boundary_delta.norm() < 1e-10)
        return false;
    if (!segmentsIntersect(
            curve_a, curve_b, boundary_a, boundary_b, nullptr))
        return false;

    auto outside_endpoint_tolerance = [&](const Vec2d& point) {
        return (point - curve_start).norm() > tol &&
               (point - curve_end).norm() > tol;
    };

    const double denominator = cross2d(curve_delta, boundary_delta);
    if (std::abs(denominator) > 1e-12) {
        Vec2d intersection;
        if (!segmentsIntersect(
                curve_a, curve_b, boundary_a, boundary_b, &intersection))
            return false;
        return outside_endpoint_tolerance(intersection);
    }

    const double curve_len2 = curve_delta.squaredNorm();
    const double t0 = (boundary_a - curve_a).dot(curve_delta) / curve_len2;
    const double t1 = (boundary_b - curve_a).dot(curve_delta) / curve_len2;
    const double lo = std::max(0.0, std::min(t0, t1));
    const double hi = std::min(1.0, std::max(t0, t1));
    if (lo > hi + 1e-9)
        return false;
    // Five probes cover both ends and the interior of the overlap interval.
    // The interval is at most one sampled curve segment, so this also catches
    // a contact that starts at a curve endpoint and continues beyond it.
    for (int i = 0; i <= 4; ++i) {
        const double t = lo + (hi - lo) * static_cast<double>(i) / 4.0;
        if (outside_endpoint_tolerance(curve_a + t * curve_delta))
            return true;
    }
    return false;
}

/// RoadEdge 若整条折线都贴着精细 area.outer，则它只是围栏轮廓的重复报告，
/// 不能再以路口中心侧把有限边外推成无限半平面。逐段采样避免仅端点贴近围栏、
/// 中间却横穿路口的长边被误分类；该判据只影响中心侧，真实相交和净距仍照常检查。
inline bool roadEdgeDuplicatesFenceOutline(
        const Boundary& boundary, const Polygon2d* fence_outline,
        double max_distance = 0.60, double sample_spacing = 0.25) {
    if (boundary.type != Boundary::Type::RoadEdge || !fence_outline ||
        fence_outline->outer.size() < 2 || boundary.geometry.points.size() < 2)
        return false;

    const auto& edge = boundary.geometry.points;
    const auto& ring = fence_outline->outer;
    const double spacing = std::max(0.05, sample_spacing);
    for (int i = 0; i + 1 < (int)edge.size(); ++i) {
        const Vec2d a = edge[i];
        const Vec2d b = edge[i + 1];
        const double length = (b - a).norm();
        if (length < 1e-8)
            continue;
        const int samples = std::max(1, (int)std::ceil(length / spacing));
        for (int s = 0; s <= samples; ++s) {
            const Vec2d point = a + (double)s / (double)samples * (b - a);
            double distance = std::numeric_limits<double>::infinity();
            for (int r = 0; r < (int)ring.size(); ++r) {
                const Vec2d ra = ring[r];
                const Vec2d rb = ring[(r + 1) % ring.size()];
                distance = std::min(distance, pointToSegment(point, ra, rb).first);
            }
            if (distance > max_distance)
                return false;
        }
    }
    return true;
}

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
        const std::vector<std::vector<bool>>* exempt_masks = nullptr,
        const Polygon2d* fence_outline = nullptr) {
    auto endpointTouchesOtherBoundary = [&](size_t boundary_index,
                                             const Vec2d& point) {
        constexpr double kJunctionTolerance = 1e-3;
        for (size_t other_index = 0; other_index < boundaries.size(); ++other_index) {
            if (other_index == boundary_index)
                continue;
            const auto& other_points = boundaries[other_index].geometry.points;
            for (int i = 0; i + 1 < (int)other_points.size(); ++i) {
                if (pointToSegment(point, other_points[i], other_points[i + 1]).first <=
                    kJunctionTolerance)
                    return true;
            }
        }
        return false;
    };

    std::vector<BoundarySafetySegment> segments;
    for (size_t bi = 0; bi < boundaries.size(); ++bi) {
        const auto& bnd = boundaries[bi];
        const auto& pts = bnd.geometry.points;
        const bool duplicate_fence_outline =
            roadEdgeDuplicatesFenceOutline(bnd, fence_outline);
        (void)exempt_masks;
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
            seg.bbox.expand(seg.a);
            seg.bbox.expand(seg.b);
            seg.len = len;
            seg.center_signed = cross2d(d, center - pts[i]) / len;
            seg.road_edge = (bnd.type == Boundary::Type::RoadEdge);
            // RoadEdge 与 Connector/其它 Boundary 端点相接后已经是有限的
            // 连续链，单段的中心侧不能再向无穷半平面外推。这里仅关闭“外侧”
            // 推断；真实相交、共线贴合、重叠和 mode=2 净距仍由线段检查执行。
            const bool chain_junction = seg.road_edge &&
                (endpointTouchesOtherBoundary(bi, seg.a) ||
                 endpointTouchesOtherBoundary(bi, seg.b));
            seg.center_side_reliable = !duplicate_fence_outline && !chain_junction;
            // 中心侧可靠性必须在完整折线上统计，豁免段只是不参与判违，
            // 不能改变该 boundary 的中心侧一致性结论。
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
            seg.center_side_reliable = seg.center_side_reliable && reliable;
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
        // 很短的 RoadEdge 通常是车道端点/鼻端的连接片段，有限线段附近
        // 不足以支持稳定的“无限中心侧半平面”推断；真实接触和 mode=2
        // 净距仍由独立检查负责。仅对长度至少 2m 的边段启用外侧判定。
        if (seg.len < 2.0)
            continue;
        const double dx = pt.x() < seg.bbox.min_pt.x()
            ? seg.bbox.min_pt.x() - pt.x()
            : (pt.x() > seg.bbox.max_pt.x()
                ? pt.x() - seg.bbox.max_pt.x() : 0.0);
        const double dy = pt.y() < seg.bbox.min_pt.y()
            ? seg.bbox.min_pt.y() - pt.y()
            : (pt.y() > seg.bbox.max_pt.y()
                ? pt.y() - seg.bbox.max_pt.y() : 0.0);
        if (dx * dx + dy * dy > max_influence_dist * max_influence_dist)
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
        int min_samples = 64, double curve_endpoint_tol = kConnectionPointTolerance,
        double boundary_endpoint_tol = 0.10, double outside_tol = 0.05,
        double sample_spacing = 0.18) {
    (void)boundary_endpoint_tol;
    const double endpoint_tol = std::min(
        std::max(0.0, curve_endpoint_tol), kConnectionPointTolerance);
    BoundarySafetyResult result;
    if (segments.empty() || curve.empty())
        return result;

    std::vector<Vec2d> pts;
    pts.reserve(std::max(min_samples, 2));
    for (int si = 0; si < (int)curve.segs.size(); ++si) {
        const auto& curve_seg = curve.segs[si];
        int seg_samples = std::max(
            2, (int)std::ceil(curve_seg.arcLength(12) /
                              std::max(0.05, sample_spacing)) + 1);
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

    // 采样折线 × 边界安全段同样是 O(Ns×Nseg) 的双层循环。先算好每个采样段的
    // 包围盒，逐段拒绝掉与 seg.bbox 不相交的配对：
    // segmentHasForbiddenBoundaryContact 的首个判据是 segmentsIntersect，
    // 包围盒不相交时必然返回 false，因此剪枝不改变判定结果。
    // 1e-9 的余量对应 segmentsIntersect 共线重叠判据里的参数容差。
    const double box_slack = 1e-9;
    const std::size_t sample_segs = pts.size() - 1;
    std::vector<BoundingBox2d> sample_box(sample_segs);
    for (std::size_t i = 0; i < sample_segs; ++i) {
        sample_box[i].expand(pts[i]);
        sample_box[i].expand(pts[i + 1]);
    }
    const auto boxesApart = [box_slack](const BoundingBox2d& x,
                                       const BoundingBox2d& y) {
        return x.max_pt[0] < y.min_pt[0] - box_slack ||
               y.max_pt[0] < x.min_pt[0] - box_slack ||
               x.max_pt[1] < y.min_pt[1] - box_slack ||
               y.max_pt[1] < x.min_pt[1] - box_slack;
    };

    for (const auto& seg : segments) {
        if (!curve_box.intersects(seg.bbox))
            continue;
        for (std::size_t i = 0; i < sample_segs; ++i) {
            if (boxesApart(sample_box[i], seg.bbox))
                continue;
            if (segmentHasForbiddenBoundaryContact(
                    pts[i], pts[i + 1], seg.a, seg.b,
                    pts.front(), pts.back(), endpoint_tol)) {
                result.intersects = true;
                break;
            }
        }
        if (result.intersects)
            break;
    }

    for (int i = 1; i + 1 < (int)pts.size(); ++i) {
        if ((pts[i] - pts.front()).norm() <= endpoint_tol ||
            (pts[i] - pts.back()).norm() <= endpoint_tol)
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

// ── buildBoundarySafetySegments 结果缓存 ─────────────────────────────────
// 该函数的输出只取决于 boundaries 内容、center 与 fence_outline
// （exempt_masks 未参与计算）。它为每个 RoadEdge 段的两个端点扫描其它所有
// boundary 的全部线段，复杂度 O(S²)：110004764 有 91 条 boundary，单次重建
// 就要做数百万次点到线段距离，而候选评估会对同一批 boundary 反复调用它，
// 实测占该数据集 62% 的自身耗时。
//
// 这里按内容指纹缓存结果。指纹覆盖全部 boundary 类型与坐标、center 和围栏
// 外环，代价与点数成正比，相对 O(S²) 重建可以忽略；只要内容不变就复用，
// 返回的线段集合与每次重建逐位一致。
inline std::uint64_t boundarySafetyFingerprint(
        const std::vector<Boundary>& boundaries, const Vec2d& center,
        const Polygon2d* fence_outline) {
    std::uint64_t h = 1469598103934665603ULL;
    const auto mix = [&h](double v) {
        std::uint64_t bits = 0;
        std::memcpy(&bits, &v, sizeof(bits));
        h = (h ^ bits) * 1099511628211ULL;
    };
    const auto mixSize = [&h](std::size_t v) {
        h = (h ^ static_cast<std::uint64_t>(v)) * 1099511628211ULL;
    };
    mixSize(boundaries.size());
    for (const auto& bnd : boundaries) {
        mixSize(static_cast<std::size_t>(bnd.type));
        mixSize(bnd.geometry.points.size());
        for (const auto& pt : bnd.geometry.points) {
            mix(pt[0]);
            mix(pt[1]);
        }
    }
    mix(center[0]);
    mix(center[1]);
    if (fence_outline) {
        mixSize(fence_outline->outer.size() + 1);
        for (const auto& pt : fence_outline->outer) {
            mix(pt[0]);
            mix(pt[1]);
        }
    } else {
        mixSize(0);
    }
    return h;
}

inline const std::vector<BoundarySafetySegment>& cachedBoundarySafetySegments(
        const std::vector<Boundary>& boundaries, const Vec2d& center,
        const Polygon2d* fence_outline) {
    struct Entry {
        std::uint64_t key = 0;
        bool valid = false;
        std::vector<BoundarySafetySegment> segments;
    };
    // 同一轮生成里只会交替使用少量 (center, fence) 组合，8 槽轮转足够。
    static thread_local std::vector<Entry> cache(8);
    static thread_local std::size_t next_slot = 0;
    const std::uint64_t key =
        boundarySafetyFingerprint(boundaries, center, fence_outline);
    for (const Entry& entry : cache) {
        if (entry.valid && entry.key == key)
            return entry.segments;
    }
    Entry& slot = cache[next_slot];
    next_slot = (next_slot + 1) % cache.size();
    slot.segments =
        buildBoundarySafetySegments(boundaries, center, nullptr, fence_outline);
    slot.key = key;
    slot.valid = true;
    return slot.segments;
}

inline BoundarySafetyResult curveBoundarySafety(
        const BezierCurve& curve, const std::vector<Boundary>& boundaries,
        const Vec2d& center, int min_samples = 64,
        double curve_endpoint_tol = kConnectionPointTolerance,
        double boundary_endpoint_tol = 0.10,
        double outside_tol = 0.05,
        const Polygon2d* fence_outline = nullptr,
        double sample_spacing = 0.18) {
    const std::vector<BoundarySafetySegment>& all_segments =
        cachedBoundarySafetySegments(boundaries, center, fence_outline);
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
            if (query_box.intersects(seg.bbox))
                segments.push_back(seg);
        }
    } else {
        segments = all_segments;
    }
    return curveBoundarySafety(
        curve, segments, min_samples, curve_endpoint_tol,
        boundary_endpoint_tol, outside_tol, sample_spacing);
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
        double curve_endpoint_tol = kConnectionPointTolerance, double boundary_endpoint_tol = 0.10,
        double outside_tol = 0.05, double pos_tol = 0.60,
        double dir_min_dot = 0.70, double graze_tol = 0.25,
        const Polygon2d* fence_outline = nullptr) {
    // 保留旧函数名以兼容已有调用方，但“忽略端点擦碰”已不再是产品规则：
    // 只有曲线真实首/尾连接点允许浮点误差，端点后的贴合、重叠和穿越一律保留。
    (void)pos_tol;
    (void)dir_min_dot;
    (void)graze_tol;
    return curveBoundarySafety(
        curve, boundaries, center, min_samples, curve_endpoint_tol,
        boundary_endpoint_tol, outside_tol, fence_outline);
}

/// 生成、优化后审计与公共约束验证共用的输入级边界口径。只有精细路口面可作为
/// RoadEdge 重复轮廓证据；粗糙包络不参与该判断，避免放宽真实道路边界。
inline BoundarySafetyResult curveBoundarySafetyForInput(
        const BezierCurve& curve, const IntersectionInput& input,
        int min_samples = 64, double curve_endpoint_tol = kConnectionPointTolerance,
        double boundary_endpoint_tol = 0.10, double outside_tol = 0.05) {
    const Polygon2d* fence_outline =
        !input.area.is_rough && !input.area.geometry.outer.empty()
            ? &input.area.geometry : nullptr;
    return curveBoundarySafetyIgnoringEndpointGraze(
        curve, input.boundaries, boundarySafetyCenter(input), min_samples,
        curve_endpoint_tol, boundary_endpoint_tol, outside_tol,
        0.60, 0.70, 0.25, fence_outline);
}

} // 命名空间 isg
