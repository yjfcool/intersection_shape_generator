#include "generation/uturn_multi_constraint_search.h"

#include "constraints/uturn_envelope_constraint.h"
#include "constraints/fence_check.h"
#include "curve/curve_utils.h"
#include "curve/hermite_init.h"
#include "utils.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>

namespace isg {

static bool debugUTurnSearch() {
    static const bool enabled = std::getenv("ISG_DEBUG_UTURN") != nullptr;
    return enabled;
}

UTurnSolverResult UTurnMultiConstraintSearch::search(
    const Vec2d& p0, const Vec2d& t0, const Vec2d& p1, const Vec2d& t1,
    const IntersectionInput& input,
    const std::vector<SampledSiblingCurve>& sampled_siblings,
    const BezierCurve& reference, const UTurnSearchBackend& backend,
    double min_lead0_floor, double min_lead1_floor,
    int depth_preference, int lead1_depth_preference,
    double lateral_preference,
    const std::vector<Crosswalk>* crosswalks_for_clearance) const {

    UTurnSolverResult best;

    bool debug = debugUTurnSearch();
    if (debug) {
        fprintf(stderr, "[UTURN-SOLVER] p0=(%.2f,%.2f) p1=(%.2f,%.2f) siblings=%zu (exempt_a1 count: ",
                p0.x(), p0.y(), p1.x(), p1.y(), sampled_siblings.size());
        int ex_cnt = 0;
        for (auto& s : sampled_siblings) if (s.exempt_a1) ++ex_cnt;
        fprintf(stderr, "%d/%zu)\n", ex_cnt, sampled_siblings.size());
        fprintf(stderr, "[UTURN-SOLVER] siblings:");
        for (const auto& s : sampled_siblings)
            fprintf(stderr, " %s(ex=%d)", s.curve.empty() ? "?" : s.id.c_str(),
                    s.exempt_a1 ? 1 : 0);
        fprintf(stderr, "\n");
    }

    Vec2d T0 = t0.norm() > 1e-8 ? t0.normalized() : Vec2d(1, 0);
    Vec2d T1 = t1.norm() > 1e-8 ? t1.normalized() : -T0;
    Vec2d axis = T0 - T1;
    if (axis.norm() < 1e-8) axis = T0;
    axis.normalize();

    Vec2d lat_dir{-axis[1], axis[0]};
    double turn_gap = std::abs((p1 - p0).dot(lat_dir));
    double chord_len = (p1 - p0).norm();
    double ref_len = std::max(reference.arcLength(), chord_len);
    double min_arch_ratio = chord_len >= 10.0 ? 1.40 : 1.0;
    double max_arch_curvature = chord_len >= 10.0 ? 0.50 : std::numeric_limits<double>::infinity();
    double endpoint_box_margin = std::max(chord_len, turn_gap) * 1.10 + 1.0;
    auto escapes_endpoint_box = [&](const BezierCurve& curve) {
        BoundingBox2d box = curve.bbox();
        double min_x = std::min(p0.x(), p1.x()) - endpoint_box_margin;
        double max_x = std::max(p0.x(), p1.x()) + endpoint_box_margin;
        double min_y = std::min(p0.y(), p1.y()) - endpoint_box_margin;
        double max_y = std::max(p0.y(), p1.y()) + endpoint_box_margin;
        return box.min_pt.x() < min_x || box.max_pt.x() > max_x ||
               box.min_pt.y() < min_y || box.max_pt.y() > max_y;
    };

    // 按U型调头横向间距给出物理曲率上限；小半径天然曲率更高。
    double maxk_phys_bound = (turn_gap < 1.0)
        ? std::max(8.0, 4.0 / std::max(0.1, turn_gap))
        : 1.5;

    // 同簇交叉权重最高；一次非豁免交叉比形态细节更严重。
    const double W_CROSS  = 1000.0;
    const double W_OBST   =  200.0;   // 每米障碍物侵入惩罚。
    const double W_BOUNDARY = 500.0;  // RoadEdge穿越或mode=2非端点1m净距惩罚。
    const double W_FENCE  =  200.0;   // 每米围栏越界惩罚。
    const double W_G1     =   50.0;   // G1 硬门槛之后的排序项。
    const double W_MAXK   =   20.0;   // 曲率超过1.0后的惩罚。
    const double W_ARCH   =    5.0;   // 偏离目标弧长比的惩罚。
    const double W_COMPACT =   0.5;   // 横向/把手/延长量紧凑度惩罚。
    const double W_LEN    =    0.1;   // 偏离参考弧长的惩罚。
    const double W_LATPREF = 200.0;

    constexpr double G1_HARD_MIN = 0.99;
    auto physical_bad = [](double obst_pen, double boundary_pen, double fence_ovf) {
        return obst_pen > 0.05 || boundary_pen > 0.05 || fence_ovf > 0.05;
    };
    auto better_uturn_candidate = [&](int sib_x, double obst_pen, double boundary_pen,
                                      double fence_ovf,
                                      double cost, const UTurnSolverResult& incumbent) {
        if (!std::isfinite(incumbent.cost))
            return true;
        // 避让约束优先于同簇交叉排序；mode=2 RoadEdge 还必须保留
        // 非端点1m净距，不能只靠后处理对普通转向修复。
        bool bad = physical_bad(obst_pen, boundary_pen, fence_ovf);
        bool incumbent_bad = physical_bad(
            incumbent.obst_pen, incumbent.boundary_pen, incumbent.fence_overflow);
        if (bad != incumbent_bad)
            return !bad;
        if (std::abs(obst_pen - incumbent.obst_pen) > 0.02)
            return obst_pen < incumbent.obst_pen;
        if (std::abs(boundary_pen - incumbent.boundary_pen) > 0.02)
            return boundary_pen < incumbent.boundary_pen;
        if (std::abs(fence_ovf - incumbent.fence_overflow) > 0.02)
            return fence_ovf < incumbent.fence_overflow;
        if (sib_x != incumbent.sibling_crosses)
            return sib_x < incumbent.sibling_crosses;
        return cost < incumbent.cost;
    };

    // 候选网格: 尺度、横向偏置、首尾直行延长和把手偏置。
    // 小半径候选少、约束紧；大半径扩大横向和深度搜索以分离同簇U型调头。
    std::vector<double> scales;
    std::vector<double> lat_biases;
    std::vector<double> lead0_extras;
    std::vector<double> handle_biases;

    scales = turn_gap >= 10.0
        ? std::vector<double>{1.0, 1.50, 2.20, 3.20, 4.50}
        : (turn_gap >= 4.0
            ? std::vector<double>{0.85, 1.0, 1.25}
            : std::vector<double>{0.7, 1.0});
    if (turn_gap < 1.5) {
        lat_biases = {0.0};
        lead0_extras = {0.0};
        handle_biases = {0.0};
    } else if (turn_gap < 4.0) {
        lat_biases = {0.0, 1.5, -1.5, 3.0, -3.0};
        lead0_extras = {0.0, 2.0};
        handle_biases = {0.0};
    } else {
        if (lateral_preference > 0.5)
            lat_biases = {0.0, 1.5, 3.0, 5.0, -1.5, -3.0};
        else if (lateral_preference < -0.5)
            lat_biases = {0.0, -1.5, -3.0, -5.0, 1.5, 3.0};
        else
            lat_biases = {0.0, 1.5, -1.5, 3.0, -3.0};

        if (depth_preference > 0)
            lead0_extras = {0.0, 2.0, 4.0};
        else if (depth_preference < 0)
            lead0_extras = {0.0, 1.0};
        else
            lead0_extras = {0.0, 2.0, 4.0};
        handle_biases = turn_gap >= 10.0
            ? std::vector<double>{0.0, 1.5, -1.5}
            : std::vector<double>{0.0};
    }

    // 首尾深度偏好只参与软排序，零交叉候选可覆盖偏好方向。
    std::vector<double> lead1_extras = {0.0};  // 默认不变化
    if (turn_gap >= 1.5) {
        if (lead1_depth_preference > 0)
            lead1_extras = {0.0, 2.0, 4.0};
        else if (lead1_depth_preference < 0)
            lead1_extras = {0.0, 1.0};
        else
            lead1_extras = {0.0, 2.0, 4.0};
    }

    if (input.obstacles.empty()) {
        scales = turn_gap >= 10.0
            ? std::vector<double>{1.0, 1.80, 2.80, 4.00}
            : std::vector<double>{0.7, 1.0};
        if (turn_gap >= 1.5) {
            // 含 Boundary 且无障碍物的场景仍需搜索 U-turn 首尾深度；
            // 过大的横向扫描通常只会增加 Boundary/Fence 风险，因此保留局部横向带，
            // 同簇顺序交给首尾深度和最终曲线对修复处理。
            bool complex_boundary_grid = input.boundaries.size() >= 20;
            lat_biases = complex_boundary_grid
                ? std::vector<double>{0.0, 3.0, -3.0}
                : std::vector<double>{0.0, 3.0, -3.0, 6.0, -6.0, 10.0, -10.0};
            lead0_extras = turn_gap < 4.0
                ? std::vector<double>{0.0, 1.5, 3.0}
                : std::vector<double>{0.0, 2.0, 4.0};
            lead1_extras = {0.0, 2.5, 5.0};
            handle_biases = {0.0};
        }
    }

    auto enforce_fence = !input.area.is_rough && !input.area.geometry.outer.empty();
    const Polygon2d& fence = input.area.geometry;
    const std::vector<Crosswalk>& clearance_crosswalks =
        crosswalks_for_clearance ? *crosswalks_for_clearance : input.crosswalks;
    auto fence_overflow = [&](const BezierCurve& curve) {
        double overflow = 0.0;
        if (!enforce_fence)
            return overflow;
        const auto points = curve.sample(24);
        for (int i = 1; i + 1 < static_cast<int>(points.size()); ++i)
            if (!polygonContains(fence, points[i]))
                overflow = std::max(overflow, pointToPolygonDist(points[i], fence));
        return overflow;
    };

    for (double scale : scales) {
        for (double lat_bias : lat_biases) {
            for (double lead0_extra : lead0_extras) {
                for (double lead1_extra : lead1_extras) {
                for (double handle_bias : handle_biases) {
                    double eff_min_lead0 = std::max(min_lead0_floor, min_lead0_floor + lead0_extra);
                    double eff_min_lead1 = std::max(min_lead1_floor, min_lead1_floor + lead1_extra);
                    BezierCurve cand;
                    cand.segs.push_back(makeAlignedUTurnCubic(
                        p0, T0, p1, T1, scale, handle_bias, lat_bias,
                        eff_min_lead0, eff_min_lead1));
                    if (cand.empty()) continue;

                    // 硬过滤: 自交、包络过大或包络塌缩。
                    if (curveSelfIntersectsBusiness(cand, 1.0)) continue;
                    if (UTurnEnvelopeConstraint().exceeds(
                            cand, reference, p0, p1)) continue;
                    if (UTurnEnvelopeConstraint().collapses(
                            cand, reference, p0, p1)) continue;
                    if (escapes_endpoint_box(cand)) continue;

                    double maxk = cand.maxCurvature(20);
                    if (maxk > maxk_phys_bound) continue;
                    double arc = cand.arcLength();
                    double arc_chord = chord_len > 1e-6 ? arc / chord_len : 1.0;
                    if (arc_chord <= min_arch_ratio || maxk >= max_arch_curvature)
                        continue;

                    // 端点 G1 方向一致性。
                    Vec2d st = cand.startTan().norm() > 1e-8 ? cand.startTan().normalized() : Vec2d(1, 0);
                    Vec2d et = cand.endTan().norm()   > 1e-8 ? cand.endTan().normalized()   : Vec2d(1, 0);
                    double g1_0 = st.dot(T0);
                    double g1_1 = et.dot(T1);
                    double g1_min = std::min(g1_0, g1_1);
                    if (g1_min < G1_HARD_MIN)
                        continue;

                    // 候选采样同时用于同簇交叉与端点侧向评分，避免每个
                    // U型候选重复按弧长采样。
                    SampledCurve cand_sample = backend.sample(cand);
                    int sib_x = backend.sibling_cross_count(cand);

                    double obst_pen = backend.obstacle_penalty(cand);
                    double boundary_pen = backend.boundary_penalty(cand);
                    double fence_ovf = fence_overflow(cand);

                    // 只惩罚顶点附近的中弧侵入；首尾直行段允许跨越 Crosswalk。
                    Vec2d ut_axis_cand = (T0 - T1).normalized();
                    if (ut_axis_cand.dot(T0) < 0) ut_axis_cand = -ut_axis_cand;
                    auto pts_c = cand.sample(40);
                    Vec2d apex_c = p0;
                    double max_proj_c = -1e18;
                    for (auto& pt : pts_c) {
                        double proj = (pt - p0).dot(ut_axis_cand);
                        if (proj > max_proj_c) { max_proj_c = proj; apex_c = pt; }
                    }

                    double xwalk_pen = 0.0;
                    if (!clearance_crosswalks.empty()) {
                        for (auto& pt : pts_c) {
                            if ((pt - apex_c).norm() > 3.0) continue;
                            for (const auto& cw : clearance_crosswalks) {
                                if (polygonContains(cw.geometry, pt)) {
                                    xwalk_pen += pointToPolygonDist(pt, cw.geometry);
                                }
                            }
                        }
                    }

                    // 非豁免共端点兄弟按 expected_side 保持顶点侧向顺序。
                    double order_violation = 0.0;
                    for (auto& sib : sampled_siblings) {
                        if (sib.exempt_a1) continue;
                        if (sib.expected_side == 0) continue;
                        if (sib.ref_perp.norm() < 1e-9) continue;
                        if (sib.sampled.pts.size() < 3) continue;

                        double cur_apex_lat = apex_c.dot(sib.ref_perp);
                        Vec2d sib_apex = sib.sampled.pts[0];
                        double sib_max_proj = -1e18;
                        for (auto& pt : sib.sampled.pts) {
                            double proj = (pt - p0).dot(ut_axis_cand);
                            if (proj > sib_max_proj) { sib_max_proj = proj; sib_apex = pt; }
                        }
                        double sib_apex_lat = sib_apex.dot(sib.ref_perp);

                        double diff = cur_apex_lat - sib_apex_lat;
                        double ORDER_MARGIN = 0.5;  // 米 — 最小分离量
                        double viol = 0;
                        if (sib.expected_side == +1) {
                            viol = std::max(0.0, diff + ORDER_MARGIN);
                        } else {
                            viol = std::max(0.0, -diff + ORDER_MARGIN);
                        }
                        order_violation += viol;
                    }
                    double endpoint_side_violation = backend.endpoint_side_violation(cand_sample);

                    if (debug && sib_x == 0) {
                        fprintf(stderr, "[UTURN-SOLVER] ZERO-X: scale=%.2f lat=%.1f lead0=%.2f hand=%.1f maxk=%.3f g1=%.3f arc/c=%.3f xwalk=%.2f ord_viol=%.2f\n",
                                scale, lat_bias, eff_min_lead0, handle_bias, maxk, g1_min, arc_chord, xwalk_pen, order_violation);
                    }

                    // 联合代价: 人行横道、同簇排序和横向偏好只在硬过滤后参与排序。
                    const double W_XWALK = 2000.0;
                    const double W_ORDER = 300.0;
                    const double W_ENDPOINT_SIDE = 500.0;
                    const double W_DEPTH = 50.0;
                    double cost = 0.0;
                    cost += W_CROSS * sib_x;
                    cost += W_OBST * obst_pen;
                    cost += W_BOUNDARY * boundary_pen;
                    cost += W_FENCE * fence_ovf;
                    cost += W_XWALK * xwalk_pen;
                    cost += W_ORDER * order_violation;
                    cost += W_ENDPOINT_SIDE * endpoint_side_violation;
                    cost += W_G1 * (1.0 - g1_min);
                    cost += W_MAXK * std::max(0.0, maxk - 1.0);
                    cost += W_ARCH * std::abs(arc_chord - 1.55);
                    // 软深度偏好: 惩罚"错误"深度方向
                    if (depth_preference > 0 && lead0_extra < 0.5) cost += W_DEPTH;
                    if (depth_preference < 0 && lead0_extra > 0.5) cost += W_DEPTH;
                    if (lead1_depth_preference > 0 && lead1_extra < 0.5) cost += W_DEPTH;
                    if (lead1_depth_preference < 0 && lead1_extra > 0.5) cost += W_DEPTH;
                    // 软横向偏好惩罚错误偏移方向。
                    if (lateral_preference > 0.5 && lat_bias < 0.5) cost += W_LATPREF;
                    if (lateral_preference < -0.5 && lat_bias > -0.5) cost += W_LATPREF;
                    // Crosswalk 场景降低深拱的紧凑度惩罚。
                    double lead0_penalty = input.crosswalks.empty() ? 0.3 : 0.05;
                    cost += W_COMPACT * (std::abs(lat_bias) + std::abs(handle_bias) + lead0_penalty * lead0_extra);
                    cost += W_LEN * std::abs(arc - ref_len);

                    if (better_uturn_candidate(
                            sib_x, obst_pen, boundary_pen, fence_ovf, cost, best)) {
                        best.curve = cand;
                        best.cost = cost;
                        best.sibling_crosses = sib_x;
                        best.obst_pen = obst_pen;
                        best.boundary_pen = boundary_pen;
                        best.fence_overflow = fence_ovf;
                        best.lateral_bias_used = lat_bias;
                        best.lead0_used = eff_min_lead0;
                    }
                }  // 把手偏置
                }  // 退出侧延长
            }  // 进入侧延长
        }  // 横向偏置
    }  // 尺度

    if (debug) {
        if (std::isfinite(best.cost)) {
            fprintf(stderr, "[UTURN-SOLVER] BEST: cost=%.2f sib_x=%d lat_bias=%.1f lead0=%.1f\n",
                    best.cost, best.sibling_crosses,
                    best.lateral_bias_used, best.lead0_used);
        } else {
            fprintf(stderr, "[UTURN-SOLVER] NO VALID CANDIDATE\n");
        }
    }

    // 有深度偏好但仍相交时，补回零延长候选；物理约束仍保持不变。
    if (std::isfinite(best.cost) && best.sibling_crosses > 0 &&
        (depth_preference != 0 || lead1_depth_preference != 0)) {
        UTurnSolverResult relaxed;
        Vec2d ut_axis_r = (T0 - T1).normalized();
        if (ut_axis_r.dot(T0) < 0) ut_axis_r = -ut_axis_r;
        std::vector<double> rl_lead0 = lead0_extras;
        if (std::find(rl_lead0.begin(), rl_lead0.end(), 0.0) == rl_lead0.end())
            rl_lead0.insert(rl_lead0.begin(), 0.0);
        std::vector<double> rl_lead1 = lead1_extras;
        if (std::find(rl_lead1.begin(), rl_lead1.end(), 0.0) == rl_lead1.end())
            rl_lead1.insert(rl_lead1.begin(), 0.0);
        for (double scale : scales) {
            for (double lat_bias : lat_biases) {
                for (double l0 : rl_lead0) {
                    for (double l1 : rl_lead1) {
                        for (double hb : handle_biases) {
                            double em0 = std::max(min_lead0_floor, min_lead0_floor + l0);
                            double em1 = std::max(min_lead1_floor, min_lead1_floor + l1);
                            BezierCurve cand;
                            cand.segs.push_back(makeAlignedUTurnCubic(
                                p0, T0, p1, T1, scale, hb, lat_bias, em0, em1));
                            if (cand.empty()) continue;
                            if (curveSelfIntersectsBusiness(cand, 1.0)) continue;
                            if (UTurnEnvelopeConstraint().exceeds(
                                    cand, reference, p0, p1)) continue;
                            if (UTurnEnvelopeConstraint().collapses(
                                    cand, reference, p0, p1)) continue;
                            if (escapes_endpoint_box(cand)) continue;
                            double mk = cand.maxCurvature(20);
                            if (mk > maxk_phys_bound) continue;
                            double arc_chord = chord_len > 1e-6
                                ? cand.arcLength() / chord_len : 1.0;
                            if (arc_chord <= min_arch_ratio || mk >= max_arch_curvature)
                                continue;

                            Vec2d st = cand.startTan().norm() > 1e-8 ? cand.startTan().normalized() : Vec2d(1, 0);
                            Vec2d et = cand.endTan().norm()   > 1e-8 ? cand.endTan().normalized()   : Vec2d(1, 0);
                            double g1_min = std::min(st.dot(T0), et.dot(T1));
                            if (g1_min < G1_HARD_MIN) continue;
                            SampledCurve cand_sample = backend.sample(cand);
                            int sx = backend.sibling_cross_count(cand);

                            // 人行横道穿透检查。
                            double xp = 0.0;
                            if (!input.crosswalks.empty()) {
                                auto pts_r = cand.sample(40);
                                Vec2d ap = p0; double mp = -1e18;
                                for (auto& pt : pts_r) {
                                    double pr = (pt - p0).dot(ut_axis_r);
                                    if (pr > mp) { mp = pr; ap = pt; }
                                }
                                for (auto& pt : pts_r) {
                                    if ((pt - ap).norm() > 3.0) continue;
                                    for (const auto& cw : input.crosswalks)
                                        if (polygonContains(cw.geometry, pt))
                                            xp += pointToPolygonDist(pt, cw.geometry);
                                }
                            }

                            double obst_pen = backend.obstacle_penalty(cand);

                            double boundary_pen = backend.boundary_penalty(cand);

                            double fence_ovf = fence_overflow(cand);

                            double endpoint_side_violation = backend.endpoint_side_violation(cand_sample);
                            double c = W_CROSS * sx
                                     + W_OBST * obst_pen
                                     + W_BOUNDARY * boundary_pen
                                     + W_FENCE * fence_ovf
                                     + 2000.0 * xp
                                     + 500.0 * endpoint_side_violation
                                     + W_MAXK * std::max(0.0, mk - 1.0)
                                     + 0.02 * cand.arcLength();
                            if (lateral_preference > 0.5 && lat_bias < 0.5) c += W_LATPREF;
                            if (lateral_preference < -0.5 && lat_bias > -0.5) c += W_LATPREF;
                            if (better_uturn_candidate(
                                    sx, obst_pen, boundary_pen, fence_ovf, c, relaxed)) {
                                relaxed.curve = cand;
                                relaxed.cost = c;
                                relaxed.sibling_crosses = sx;
                                relaxed.lead0_used = em0;
                                relaxed.lateral_bias_used = lat_bias;
                                relaxed.obst_pen = obst_pen;
                                relaxed.boundary_pen = boundary_pen;
                                relaxed.fence_overflow = fence_ovf;
                            }
                            if (sx == 0 && xp < 0.5 && boundary_pen <= 0.05) break;
                        }
                        if (relaxed.sibling_crosses == 0) break;
                    }
                    if (relaxed.sibling_crosses == 0) break;
                }
                if (relaxed.sibling_crosses == 0) break;
            }
            if (relaxed.sibling_crosses == 0) break;
        }
        if (std::isfinite(relaxed.cost) &&
            better_uturn_candidate(
                relaxed.sibling_crosses, relaxed.obst_pen, relaxed.boundary_pen,
                relaxed.fence_overflow, relaxed.cost, best)) {
            best = relaxed;
            if (debug) {
                fprintf(stderr, "[UTURN-SOLVER] RELAXED: sib_x=%d lead0=%.1f lat=%.1f\n",
                        best.sibling_crosses, best.lead0_used, best.lateral_bias_used);
            }
        }
    }

    return best;
}


}  // 命名空间 isg
