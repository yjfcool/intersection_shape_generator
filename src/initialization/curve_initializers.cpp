#include "initialization/avoidance_candidate_generator.h"
#include "initialization/curve_initializer_registry.h"
#include "initialization/fixed_shape_initializer.h"
#include "initialization/ordinary_curve_initializer.h"
#include "initialization/uturn_curve_initializer.h"

#include "curve/bezier.h"
#include "curve/curve_utils.h"
#include "preprocessing/uturn_family_builder.h"
#include "utils.h"

#include <algorithm>

namespace isg {

// 固有形态初始化：逐段转换输入折线，不引入额外平滑或避让。
bool FixedShapeInitializer::hasGeometry(
    const Connectivity& connectivity) const {
    return connectivity.geometry.points.size() >= 2;
}

BezierCurve FixedShapeInitializer::build(
    const Connectivity& connectivity) const {
    BezierCurve curve;
    if (!hasGeometry(connectivity))
        return curve;
    const std::vector<Vec2d> points =
        toVec2dArray(connectivity.geometry.points);
    for (std::size_t i = 0; i + 1 < points.size(); ++i) {
        const Vec2d chord = points[i + 1] - points[i];
        if (chord.norm() < 1e-8)
            continue;
        const Vec2d direction = chord.normalized();
        curve.segs.push_back(makeCubicG1(
            points[i], direction, points[i + 1], direction, 1.0 / 3.0));
    }
    return curve;
}

// 普通曲线初始化：严格沿端点切向构造单段三次 Bezier 候选。
BezierCurve OrdinaryCurveInitializer::buildPreferredSingleCubic(
    const Vec2d& entry_point, const Vec2d& entry_tangent,
    const Vec2d& exit_point, const Vec2d& exit_tangent) const {
    BezierCurve curve;
    const Vec2d chord = exit_point - entry_point;
    const double chord_len = chord.norm();
    if (chord_len < 1e-8)
        return curve;

    const OrdinarySingleCubicHandleBounds bounds =
        ordinarySingleCubicHandleBounds(
            entry_point, entry_tangent, exit_point, exit_tangent, true);
    const Vec2d chord_dir = chord / chord_len;
    const bool straight_like =
        bounds.start_dir.dot(bounds.end_dir) > 0.90 &&
        std::abs(cross2d(bounds.start_dir, bounds.end_dir)) < 0.25 &&
        std::abs(cross2d(bounds.start_dir, chord_dir)) < 0.25;

    // 直行把手采用等弦中点对应的统一长度；控制点仍分别落在
    // 首、尾切线轴上，因此切线不一致时不强制 C1=C2。
    double start_handle = chord_len / 2.0;
    double end_handle = chord_len / 2.0;
    const bool stable_turn_intersection =
        !straight_like && bounds.has_direction_intersection &&
        bounds.direction_intersection_start <= 1.5 * chord_len &&
        bounds.direction_intersection_end <= 1.5 * chord_len;
    if (stable_turn_intersection) {
        start_handle = (2.0 / 3.0) *
            bounds.direction_intersection_start;
        end_handle = (2.0 / 3.0) *
            bounds.direction_intersection_end;
    } else if (!straight_like) {
        // 平行、反向或过远交点不适合作为形态基准，回退到有界自然弧。
        start_handle = 0.4 * chord_len;
        end_handle = 0.4 * chord_len;
    }

    BezierSegment segment;
    segment.ctrl[0] = entry_point;
    segment.ctrl[1] = entry_point + bounds.start_dir * start_handle;
    segment.ctrl[2] = exit_point - bounds.end_dir * end_handle;
    segment.ctrl[3] = exit_point;
    curve.segs.push_back(segment);
    constrainOrdinarySingleCubicControls(
        curve, entry_point, bounds.start_dir,
        exit_point, bounds.end_dir, true);

    // 2/3 交点升阶是首选形态而不是无条件硬编码。短急弯若直接采用
    // 该比例可能产生过高曲率，回退到有界自然弧，后续仍可由候选搜索
    // 在不破坏同簇和物理约束的前提下微调。
    if (!straight_like) {
        const double arc_chord = curve.arcLength() / chord_len;
        const double max_curvature = curve.maxCurvature(40);
        const double max_allowed_curvature = chord_len <= 10.0 ? 6.0 : 2.5;
        if (arc_chord < 1.005 || arc_chord > 1.35 ||
            max_curvature > max_allowed_curvature) {
            curve.segs.front() = makeCubicG1(
                entry_point, bounds.start_dir,
                exit_point, bounds.end_dir, 0.4);
            constrainOrdinarySingleCubicControls(
                curve, entry_point, bounds.start_dir,
                exit_point, bounds.end_dir, true);
        }
    }
    return curve;
}

BezierCurve OrdinaryCurveInitializer::buildSingleCubic(
    const Vec2d& entry_point, const Vec2d& entry_tangent,
    const Vec2d& exit_point, const Vec2d& exit_tangent, double alpha) const {
    BezierCurve curve;
    const Vec2d fallback = exit_point - entry_point;
    const Vec2d start = entry_tangent.norm() > 1e-8
        ? entry_tangent.normalized()
        : (fallback.norm() > 1e-8 ? fallback.normalized() : Vec2d(1, 0));
    const Vec2d end = exit_tangent.norm() > 1e-8
        ? exit_tangent.normalized() : start;
    curve.segs.push_back(
        makeCubicG1(entry_point, start, exit_point, end, alpha));
    return curve;
}

std::vector<BezierCurve> OrdinaryCurveInitializer::buildAlphaCandidates(
    const Vec2d& entry_point, const Vec2d& entry_tangent,
    const Vec2d& exit_point, const Vec2d& exit_tangent,
    const std::vector<double>& alphas) const {
    std::vector<BezierCurve> candidates;
    candidates.reserve(alphas.size());
    for (double alpha : alphas) {
        candidates.push_back(buildSingleCubic(
            entry_point, entry_tangent, exit_point, exit_tangent, alpha));
    }
    return candidates;
}

// 物理避让初始化：将有序路点转换为保持首尾切向的分段曲线。
BezierCurve AvoidanceCandidateGenerator::buildWaypointCurve(
    const std::vector<Vec2d>& points, const Vec2d& start_tangent,
    const Vec2d& end_tangent) const {
    if (points.size() < 2)
        return BezierCurve();
    std::vector<Vec2d> tangents(points.size(), Vec2d(1, 0));
    tangents.front() = start_tangent.norm() > 1e-8
        ? start_tangent.normalized() : (points[1] - points[0]).normalized();
    tangents.back() = end_tangent.norm() > 1e-8
        ? end_tangent.normalized()
        : (points.back() - points[points.size() - 2]).normalized();
    for (std::size_t i = 1; i + 1 < points.size(); ++i) {
        const Vec2d direction = points[i + 1] - points[i - 1];
        tangents[i] = direction.norm() > 1e-8
            ? direction.normalized() : tangents[i - 1];
    }
    return makeCurveFromKnots(points, tangents, 0.34);
}

// U-turn 初始化：构造轴向平齐的单段或三段几何表达。

BezierCurve UTurnCurveInitializer::buildAligned(
    const Vec2d& p0, const Vec2d& t0, const Vec2d& p1, const Vec2d& t1,
    const Vec2d& offset_dir, double offset_m, double handle_scale,
    double min_lead0, double min_lead1) const {
    Vec2d T0 = t0.norm() > 1e-8 ? t0.normalized() : Vec2d(1, 0);
    Vec2d T1 = t1.norm() > 1e-8 ? t1.normalized() : -T0;
    Vec2d axis = T0 - T1;
    if (axis.norm() < 1e-8)
        axis = T0;
    axis.normalize();
    if (axis.dot(T0) < 0.0)
        axis = -axis;

    double fwd_bias = 0.0;
    double lat_bias = 0.0;
    if (offset_m > 0.0 && offset_dir.norm() > 1e-8) {
        Vec2d dir = offset_dir.normalized();
        fwd_bias = offset_m * dir.dot(axis);
        Vec2d lat_dir{-axis.y(), axis.x()};
        lat_bias = offset_m * dir.dot(lat_dir);
    }

    BezierCurve curve;
    curve.segs.push_back(makeAlignedUTurnCubic(
        p0, T0, p1, T1, handle_scale, fwd_bias, lat_bias,
        min_lead0, min_lead1));
    return curve;
}

BezierCurve UTurnCurveInitializer::buildSegmented(
    const Vec2d& p0, const Vec2d& t0, const Vec2d& p1, const Vec2d& t1,
    double min_lead0, double min_lead1, double arc_alpha,
    double aligned_point_stagger, double lead0_extra_after_align,
    double lead1_extra_after_align, double aligned_entry_stagger,
    double aligned_exit_stagger) const {
    Vec2d T0 = t0.norm() > 1e-8 ? t0.normalized() : Vec2d(1, 0);
    Vec2d T1 = t1.norm() > 1e-8 ? t1.normalized() : -T0;
    Vec2d exit_back = -T1;
    Vec2d axis = T0 + exit_back;
    if (axis.norm() < 1e-8)
        axis = T0;
    axis.normalize();
    if (axis.dot(T0) < 0.0)
        axis = -axis;

    double s0 = p0.dot(axis);
    double s1 = p1.dot(axis);
    double common_s = std::max(s0, s1);
    double c0 = std::max(0.2, T0.dot(axis));
    double c1 = std::max(0.2, exit_back.dot(axis));
    double lead0 = std::max(0.0, (common_s - s0) / c0);
    double lead1 = std::max(0.0, (common_s - s1) / c1);
    lead0 = std::max(lead0, min_lead0);
    lead1 = std::max(lead1, min_lead1);

    double s0_new = p0.dot(axis) + lead0 * T0.dot(axis);
    double s1_new = p1.dot(axis) + lead1 * exit_back.dot(axis);
    double common_s_new = std::max(s0_new, s1_new);
    if (s0_new < common_s_new - 1e-6)
        lead0 += (common_s_new - s0_new) / c0;
    if (s1_new < common_s_new - 1e-6)
        lead1 += (common_s_new - s1_new) / c1;

    double extra_axis_advance = std::max(
        std::max(0.0, lead0_extra_after_align) * T0.dot(axis),
        std::max(0.0, lead1_extra_after_align) * exit_back.dot(axis));
    if (extra_axis_advance > 0.0) {
        lead0 += extra_axis_advance / c0;
        lead1 += extra_axis_advance / c1;
    }

    Vec2d q0 = p0 + lead0 * T0;
    Vec2d q1 = p1 + lead1 * exit_back;
    double entry_stagger = aligned_entry_stagger >= 0.0
        ? aligned_entry_stagger : aligned_point_stagger;
    double exit_stagger = aligned_exit_stagger >= 0.0
        ? aligned_exit_stagger : aligned_point_stagger;
    const double shared_stagger = std::max(
        std::max(0.0, entry_stagger), std::max(0.0, exit_stagger));
    entry_stagger = shared_stagger;
    exit_stagger = shared_stagger;
    if (entry_stagger > 0.0 || exit_stagger > 0.0) {
        Vec2d chord = q1 - q0;
        double gap = chord.norm();
        if (gap > 1e-6) {
            Vec2d lateral{-axis.y(), axis.x()};
            double lateral_gap = std::abs(chord.dot(lateral));
            if (lateral_gap > 1e-6) {
                double lead0_len = (q0 - p0).norm();
                double lead1_len = (p1 - q1).norm();
                double g1_safe_inset = 0.14 * std::min(lead0_len, lead1_len);
                double common_inset = std::min(
                    shared_stagger, std::min(0.25 * lateral_gap, g1_safe_inset));
                double side = chord.dot(lateral) >= 0.0 ? 1.0 : -1.0;
                q0 += side * common_inset * lateral;
                q1 -= side * common_inset * lateral;
            }
        }
    }

    BezierCurve curve;
    const double lead_eps = 0.20;
    Vec2d lead_dir0 = (q0 - p0).norm() > 1e-8 ? (q0 - p0).normalized() : T0;
    Vec2d lead_dir1 = (p1 - q1).norm() > 1e-8 ? (p1 - q1).normalized() : T1;
    if (lead0 > lead_eps)
        curve.segs.push_back(makeCubicG1(p0, lead_dir0, q0, lead_dir0, 1.0 / 3.0));
    if ((q1 - q0).norm() > 1e-6)
        curve.segs.push_back(makeCubicG1(q0, lead_dir0, q1, lead_dir1, arc_alpha));
    if (lead1 > lead_eps)
        curve.segs.push_back(makeCubicG1(q1, lead_dir1, p1, lead_dir1, 1.0 / 3.0));
    if (curve.empty())
        curve.segs.push_back(makeAlignedUTurnCubic(
            p0, T0, p1, T1, 1.0, 0.0, 0.0, min_lead0, min_lead1));
    return curve;
}

// 注册表保持 U-turn、fixed shape、普通曲线的既有优先级和候选顺序。
CurveInitializationOptions::CurveInitializationOptions()
    : ordinary_alphas(1, 0.4),
      allow_fixed_shape(true),
      uturn_min_lead0(0.0),
      uturn_min_lead1(0.0),
      uturn_arc_alpha(2.0 / 3.0),
      uturn_aligned_point_stagger(0.0),
      uturn_lead0_extra_after_align(0.0),
      uturn_lead1_extra_after_align(0.0),
      uturn_aligned_entry_stagger(-1.0),
      uturn_aligned_exit_stagger(-1.0) {}

void CurveInitializationOptions::applyUTurnFamily(
    const UTurnFamilyInfo& family) {
    if (!family.geometric_uturn)
        return;
    uturn_min_lead0 = family.lead0;
    uturn_min_lead1 = family.lead1;
}

std::vector<CurveCandidate> CurveInitializerRegistry::build(
    const CurveGenerationContext& context,
    const CurveInitializationOptions& options) const {
    std::vector<CurveCandidate> candidates;
    if (!context.connectivity)
        return candidates;

    const Vec2d& p0 = context.entry.first;
    const Vec2d& t0 = context.entry.second;
    const Vec2d& p1 = context.exit.first;
    const Vec2d& t1 = context.exit.second;
    const bool geometric_uturn =
        t0.norm() > 1e-8 && t1.norm() > 1e-8 &&
        t0.normalized().dot(t1.normalized()) < -0.5;

    if (geometric_uturn) {
        CurveCandidate candidate;
        candidate.origin = CandidateOrigin::SegmentedUTurn;
        candidate.curve = UTurnCurveInitializer().buildSegmented(
            p0, t0, p1, t1,
            options.uturn_min_lead0, options.uturn_min_lead1,
            options.uturn_arc_alpha, options.uturn_aligned_point_stagger,
            options.uturn_lead0_extra_after_align,
            options.uturn_lead1_extra_after_align,
            options.uturn_aligned_entry_stagger,
            options.uturn_aligned_exit_stagger);
        candidates.push_back(candidate);
        return candidates;
    }

    const FixedShapeInitializer fixed_initializer;
    if (options.allow_fixed_shape && context.connectivity->fixed_shape &&
        fixed_initializer.hasGeometry(*context.connectivity)) {
        CurveCandidate candidate;
        candidate.origin = CandidateOrigin::FixedShape;
        candidate.curve = fixed_initializer.build(*context.connectivity);
        candidate.preserves_fixed_shape = true;
        candidates.push_back(candidate);
        return candidates;
    }

    const OrdinaryCurveInitializer ordinary_initializer;
    std::vector<BezierCurve> ordinary;
    ordinary.push_back(ordinary_initializer.buildPreferredSingleCubic(
        p0, t0, p1, t1));
    const std::vector<BezierCurve> alpha_candidates =
        ordinary_initializer.buildAlphaCandidates(
            p0, t0, p1, t1, options.ordinary_alphas);
    ordinary.insert(ordinary.end(), alpha_candidates.begin(), alpha_candidates.end());
    candidates.reserve(ordinary.size());
    for (std::size_t i = 0; i < ordinary.size(); ++i) {
        CurveCandidate candidate;
        candidate.origin = CandidateOrigin::NaturalSingleCubic;
        candidate.curve = ordinary[i];
        candidate.quality_score = static_cast<double>(i);
        candidates.push_back(candidate);
    }
    return candidates;
}

}  // 命名空间 isg
