#include "penalty_cost.h"
#include "constraints/boundary_safety.h"
#include "constraints/fence_check.h"
#include "curve/curve_utils.h"
#include "utils.h"
#include <cmath>
#include <algorithm>
#include <limits>

namespace isg {

// ─── 缓存重建 ────────────────────────────────────────────────
/// 每外层迭代前重建缓存: 边界段、兄弟曲线采样点列、围栏非空标志
void PenaltyCostCache::rebuild(
        const std::vector<Boundary> &boundaries,
        const std::vector<SiblingCurve> &siblings, const Polygon2d &fence,
        const Vec2d& center) {
    bnd_segs.clear();
    for (auto& bnd : boundaries) {
        auto& pts = bnd.geometry.points;
        std::vector<BoundarySegment> boundary_segments;
        bool has_pos_center_side = false;
        bool has_neg_center_side = false;
        for (int i = 0; i + 1 < (int)pts.size(); ++i) {
            Vec2d d = pts[i + 1] - pts[i];
            double len = d.norm();
            if (len < 1e-8)
                continue;
            BoundarySegment seg;
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
            bnd_segs.push_back(seg);
        }
    }

    sib_polys.clear();
    for (auto& sib : siblings) {
        SiblingPoly sp;
        sp.exempt_a1 = sib.exempt_a1;
        sp.a2_radius = sib.exempt_a2_radius;
        sp.expected_side = sib.expected_side;
        sp.ref_perp = sib.ref_perp;
        sp.shared_endpoint = sib.shared_endpoint; // 仅标识共享连接点，不放行端点外重叠
        sp.fixed_shape = sib.fixed_shape;
        sp.pts = sib.curve.sampleByArcLength(32); // 32点采样(性能优化,原48)
        sib_polys.push_back(std::move(sp));
    }
    fence_empty = fence.outer.empty();
}

void PenaltyCost::buildCache() {
    cache_.rebuild(boundaries, siblings, fence, road_center);
}

static double roadEdgeOutsidePenaltyCached(
        const Vec2d& pt, const std::vector<BoundarySegment>& segments,
        double tol = 0.05, double max_influence_dist = 3.0) {
    const BoundarySegment* best = nullptr;
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

static double roadEdgeDistanceCached(
        const Vec2d& pt, const std::vector<BoundarySegment>& segments) {
    double best_dist = std::numeric_limits<double>::infinity();
    for (const auto& seg : segments) {
        if (!seg.road_edge)
            continue;
        Vec2d d = seg.b - seg.a;
        double len2 = d.squaredNorm();
        if (len2 < 1e-12)
            continue;
        double t = std::max(0.0, std::min(1.0, (pt - seg.a).dot(d) / len2));
        Vec2d closest = seg.a + t * d;
        best_dist = std::min(best_dist, (pt - closest).norm());
    }
    return best_dist;
}

// ─── 各项标量代价 ────────────────────────────────────────────
/// 平滑项: 曲率平方沿弧长积分 (能量法,促使曲线更直)
double PenaltyCost::evalSmooth(const BezierCurve &c) const {
    double cost = 0;
    constexpr int S = 8;  // 性能优化: 从16降至8
    for (auto& seg : c.segs) {
        double ds = seg.arcLength(S) / S;
        for (int i = 0; i <= S; ++i) {
            double k = seg.curvature((double)i / S);
            cost += k * k * ds;
        }
    }
    return cost;
}

/// 障碍物穿透惩罚: Σ max(0, clearance - sdf)²
double PenaltyCost::evalObstacle(const BezierCurve& c) const {
    if (!sdf) return 0;
    double cost = 0;
    constexpr int S = 12;  // 性能优化: 从24降至12
    for (auto& seg : c.segs)
        for (int i = 0; i <= S; ++i)
            cost += sdf->obstaclePenalty(seg.evaluate((double)i / S), obstacle_clearance);
    return cost;
}

/// 边界线穿越惩罚:
/// 1. 曲线非端点与任意Boundary相交为硬惩罚;
/// 2. RoadEdge按"路口中心所在侧"定义可行半平面, 曲线越到外侧时施加连续距离惩罚。
double PenaltyCost::evalBoundary(const BezierCurve& c) const {
    double cost = 0;
    constexpr int S = 24;
    std::vector<Vec2d> cp;
    cp.reserve(c.segs.size() * (S + 1));
    for (auto& seg : c.segs)
        for (int i = 0; i <= S; ++i)
            cp.push_back(seg.evaluate((double)i / S));

    for (auto& bs : cache_.bnd_segs)
        for (int ci = 0; ci + 1 < (int)cp.size(); ++ci)
            if (segmentsIntersect(cp[ci], cp[ci + 1], bs.a, bs.b)) {
                Vec2d isect;
                segmentsIntersect(cp[ci], cp[ci + 1], bs.a, bs.b, &isect);
                if ((isect - cp.front()).norm() <= 0.75 ||
                    (isect - cp.back()).norm() <= 0.75 ||
                    (bs.a_is_boundary_endpoint && (isect - bs.a).norm() <= 0.10) ||
                    (bs.b_is_boundary_endpoint && (isect - bs.b).norm() <= 0.10))
                    continue;
                cost += 25.0;
            }

    for (int i = 1; i + 1 < (int)cp.size(); ++i) {
        if ((cp[i] - cp.front()).norm() <= 0.75 ||
            (cp[i] - cp.back()).norm() <= 0.75)
            continue;
        double outside = roadEdgeOutsidePenaltyCached(cp[i], cache_.bnd_segs, 0.05);
        cost += outside * outside * 30.0;
    }

    if (road_edge_clearance > 0.0) {
        for (int i = 1; i + 1 < (int)cp.size(); ++i) {
            if ((cp[i] - cp.front()).norm() <= 0.75 ||
                (cp[i] - cp.back()).norm() <= 0.75)
                continue;
            double d_edge = roadEdgeDistanceCached(cp[i], cache_.bnd_segs);
            if (!std::isfinite(d_edge))
                continue;
            double deficit = road_edge_clearance - d_edge;
            if (deficit > 0.0)
                cost += deficit * deficit * 20.0;
        }
    }
    return cost;
}

/// 围栏越界惩罚: 越界点距围栏距离×2倍权重
double PenaltyCost::evalFence(const BezierCurve& c) const {
    if (cache_.fence_empty)
        return 0;
    double cost = 0;
    constexpr int S = 10;  // 性能优化: 从20降至10
    for (auto& seg : c.segs)
        for (int i = 0; i <= S; ++i) {
            Vec2d pt = seg.evaluate((double)i / S);
            if (!polygonContains(fence, pt))
                cost += pointToPolygonDist(pt, fence) * 2.0;
        }
    return cost;
}

// ─────────────────────────────────────────────────────────────────────────────
// evalCluster —— 同簇非端点不相交约束(保序拓扑约束)
//
// 针对同簇约束对(sibling)施加三重检测:
//   [A] 横向边距惩罚: 当前曲线相对兄弟曲线应保持 expected_side 指定的横向关系
//   [B] 符号变化惩罚: diffs[] 沿弧长方向不应出现符号翻转(穿越)
//   [C] 穷举段相交检测: 捕捉 nearest-point 漏检的不同弧分数交叉
// ─────────────────────────────────────────────────────────────────────────────
double PenaltyCost::evalCluster(const BezierCurve& c) const {
    if (cache_.sib_polys.empty())return 0;

    constexpr double MARGIN = 0.3;          // 米 — 横向最小间隔
    constexpr double SIGN_CHANGE_COST = 50.0; // 单次交叉基础代价
    constexpr double OBS_REDUCE = 0.4;      // 障碍物强制交叉时权重(非零,弱化)
    constexpr double MAX_DIST = 8.0;        // 米 — 最近点搜索半径(从12降至8)
    constexpr double MAX_DIST2 = MAX_DIST * MAX_DIST;
    constexpr double SKIP_FRAC = 0.02;      // 仅跳过连接点附近的数值噪声
    constexpr int N_SAMP = 32;              // 采样点数(从48降至32,性能优化)

    auto curve_pts = c.sampleByArcLength(N_SAMP);
    int N = (int) curve_pts.size();
    if (N < 4)return 0;

    // 预计算当前曲线包围盒用于兄弟曲线包围盒快速剔除
    BoundingBox2d cbox;
    for (auto& p : curve_pts) cbox.expand(p);
    cbox.min_pt -= Vec2d(MAX_DIST, MAX_DIST);
    cbox.max_pt += Vec2d(MAX_DIST, MAX_DIST);

    double cost = 0;

    for (auto& sp : cache_.sib_polys) {
        if (sp.exempt_a1)continue; // 结构性交叉或不同arm → 跳过

        double w = (sp.a2_radius > 0) ? OBS_REDUCE : 1.0; // 障碍物情形弱化但非零
        int M = (int) sp.pts.size();
        if (M < 2 || sp.ref_perp.norm() < 1e-9 || sp.expected_side == 0)continue;

        // 共端点固有形态(典型: 保留的掉头曲线)只约束非端点相交。
        // 固有掉头的局部法向不能拿来对左转/直行全程做横向排队,
        // 否则会把自然转弯弧压平或拉成S形。
        bool fixed_shared_endpoint = sp.fixed_shape && sp.shared_endpoint;
        double eff_skip = SKIP_FRAC;

        // 预计算兄弟曲线包围盒,与当前曲线包围盒比较,无重叠则跳过整个兄弟
        BoundingBox2d sbox;
        for (auto& p : sp.pts) sbox.expand(p);
        if (!cbox.intersects(sbox)) continue;

        // 预计算兄弟点列的 ref_perp 投影
        std::vector<double> sib_lat(M);
        for (int k = 0; k < M; ++k) sib_lat[k] = sp.pts[k].dot(sp.ref_perp);

        // 预计算每个采样点相对最近兄弟点的横向差
        std::vector<double> diffs(N, std::numeric_limits<double>::quiet_NaN());
        for (int i = 0; i < N; ++i) {
            double frac = (double) i / std::max(N - 1, 1);
            if (frac < eff_skip || frac > 1.0 - eff_skip)continue;
            Vec2d cp = curve_pts[i];
            double cur_lat = cp.dot(sp.ref_perp);
            // 在MAX_DIST范围内寻找最近兄弟采样点
            double best_d2 = MAX_DIST2 + 1;
            int best_k = -1;
            for (int k = 0; k < M; ++k) {
                double d2 = (sp.pts[k] - cp).squaredNorm();
                if (d2 < best_d2) {
                    best_d2 = d2;
                    best_k = k;
                }
            }
            if (best_k < 0)continue;
            // diff = latA - latB (正=当前在兄弟左侧)
            diffs[i] = cur_lat - sib_lat[best_k];

            if (!fixed_shared_endpoint) {
                // [A] 横向边距惩罚
                // expected_side=+1: 兄弟在左侧 → diff应为负 (sib>cur) → 违反=max(0, diff+MARGIN)
                // expected_side=-1: 兄弟在右侧 → diff应为正 (cur>sib) → 违反=max(0, -diff+MARGIN)
                double violation = 0;
                if (sp.expected_side == +1)
                    violation = std::max(0.0, diffs[i] + MARGIN); // diff应 < -MARGIN
                else
                    violation = std::max(0.0, -diffs[i] + MARGIN); // diff应 >  MARGIN
                cost += w * violation * violation;
            }
        }

        // [B] 符号变化(保序)惩罚: 沿弧长方向相邻采样点的横向差不应符号翻转
        if (!fixed_shared_endpoint) {
            double prev_diff = std::numeric_limits<double>::quiet_NaN();
            for (int i = 0; i < N; ++i) {
                if (std::isnan(diffs[i])) {
                    prev_diff = std::numeric_limits<double>::quiet_NaN();
                    continue;
                }
                if (!std::isnan(prev_diff) && std::abs(prev_diff) > 1e-6 && std::abs(diffs[i]) > 1e-6) {
                    if (prev_diff * diffs[i] < 0) {
                        // 符号翻转 → 检测到交叉
                        double depth = (std::abs(prev_diff) + std::abs(diffs[i])) * 0.5;
                        cost += w * SIGN_CHANGE_COST * depth;
                    }
                }
                prev_diff = diffs[i];
            }
        }

        // [C] 穷举段相交检测(跨弧分数交叉检测)
        // 捕捉曲线在截然不同弧分数位置相交的情形
        // O(N*M),N=M=32 → 1024次检查(原2304)
        constexpr double EP_TOL = 0.15; // 米 — 仅排除首尾连接点数值容差
        constexpr double EP_TOL2 = EP_TOL * EP_TOL;
        for (int ci = 0; ci + 1 < N; ++ci) {
            double frac_ci = (double) ci / std::max(N - 1, 1);
            if (frac_ci < eff_skip || frac_ci > 1.0 - eff_skip)
                continue;
            Vec2d A0 = curve_pts[ci], A1 = curve_pts[ci + 1];
            for (int si = 0; si + 1 < M; ++si) {
                Vec2d B0 = sp.pts[si], B1 = sp.pts[si + 1];
                // 包围盒距离快速剔除: 任一维度间距>MAX_DIST+2 跳过
                double dx = std::max({A0.x() - B1.x(), B0.x() - A1.x(), 0.0});
                double dy = std::max({A0.y() - B1.y(), B0.y() - A1.y(), 0.0});
                if (dx * dx + dy * dy > (MAX_DIST + 2.0) * (MAX_DIST + 2.0))
                    continue;
                Vec2d isect;
                if (!segmentsIntersect(A0, A1, B0, B1, &isect))
                    continue;
                // 跳过位于端点区域内的交点(允许端点连接)
                if ((isect - curve_pts.front()).squaredNorm() < EP_TOL2 ||
                    (isect - curve_pts.back()).squaredNorm() < EP_TOL2 ||
                    (isect - sp.pts.front()).squaredNorm() < EP_TOL2 ||
                    (isect - sp.pts.back()).squaredNorm() < EP_TOL2) continue;
                // 交叉深度 = 交叉前的横向分离量
                Vec2d near_B = B0;
                double bd2 = (B0 - A0).squaredNorm();
                double b1d2 = (B1 - A0).squaredNorm();
                if (b1d2 < bd2) {
                    bd2 = b1d2;
                    near_B = B1;
                }
                double lat_A0 = A0.dot(sp.ref_perp);
                double lat_NB = near_B.dot(sp.ref_perp);
                double depth = std::max(std::abs(lat_A0 - lat_NB), 0.3);
                cost += w * SIGN_CHANGE_COST * depth;
            }
        }
    }
    return cost;
}

/// 人行横道斜穿惩罚: 曲线在人行横道内偏离行人通行方向超过15°时惩罚
double PenaltyCost::evalCrosswalk(const BezierCurve& c) const {
    constexpr double DEG15 = 15.0 * M_PI / 180.0;
    double cost = 0;
    for (auto& cw : crosswalks) {
        for (auto& seg : c.segs)
            for (int i = 1; i < 20; ++i) {
                double t = (double)i / 20;
                Vec2d pt = seg.evaluate(t);
                if (!polygonContains(cw.geometry, pt))
                    continue;
                Vec2d tan = seg.tangent(t).normalized();
                double cosA = std::abs(tan.dot(cw.crossing_direction.normalized()));
                double exc = cosA - std::cos(M_PI / 2 - DEG15);
                if (exc > 0) cost += exc * exc;
            }
    }
    return cost;
}

// ─── G1端点强制 ──────────────────────────────────────────────
/// 将控制点投影到端点切向射线上,保证首末段G1连续
/// - 首段: ctrl[1] 必须位于射线 ctrl[0] + λ·start_dir (λ>0)
/// - 末段: ctrl[2] 必须位于射线 ctrl[3] - μ·end_dir   (μ>0)
static void enforceG1(
    BezierCurve& c, const Vec2d& start_dir, const Vec2d& end_dir,
    bool constrain_single_cubic_axes) {
    if (c.segs.empty())
        return;

    if (c.numSegments() == 1 && constrain_single_cubic_axes) {
        constrainOrdinarySingleCubicControls(
            c, c.startPt(), start_dir, c.endPt(), end_dir, true);
        return;
    }

    // 首段: ctrl[1] 投影到首点切向射线
    if (start_dir.norm() > 1e-10) {
        Vec2d& p0 = c.segs.front().ctrl[0];
        Vec2d& c1 = c.segs.front().ctrl[1];
        Vec2d sd = start_dir.normalized();
        // 当前ctrl[1]投影到射线上得到缩放λ
        double lam = (c1 - p0).dot(sd);
        lam = std::max(lam, 0.05); // 最小张力: 避免退化段
        c1 = p0 + lam * sd;
    }

    // 末段: ctrl[2] 投影到尾点切向射线
    if (end_dir.norm() > 1e-10) {
        Vec2d& p1 = c.segs.back().ctrl[3];
        Vec2d& c2 = c.segs.back().ctrl[2];
        Vec2d ed = end_dir.normalized();
        double mu = (p1 - c2).dot(ed);
        mu = std::max(mu, 0.05);
        c2 = p1 - mu * ed;
    }
}

// ─── 障碍物项解析梯度 ────────────────────────────────────────
void PenaltyCost::addObstacleGrad(
    const BezierCurve& c, const VecXd& params, double w, VecXd& grad) const {
    if (!sdf || w < 1e-15)
        return;
    constexpr int S = 12;  // 性能优化: 从24降至12
    for (int k = 0; k < c.segs.size(); ++k) {
        auto& seg = c.segs[k];
        for (int i = 0; i <= S; ++i) {
            double t = (double)i / S;
            Vec2d pt = seg.evaluate(t);
            auto kv = sdf->queryWithGrad(pt);
            double& d = kv.first;
            Vec2d& gd = kv.second;
            double slack = d - obstacle_clearance;
            if (slack >= 0)
                continue; // 无违反 → 零梯度

            double coeff = 2.0 * w * slack; // 负值(违反)

            // t处ctrl[1]与ctrl[2]的基函数导数
            double B1 = 3 * (1 - t) * (1 - t) * t; // ∂c/∂ctrl[1]
            double B2 = 3 * (1 - t) * t * t; // ∂c/∂ctrl[2]

            int base = 4 * k;
            // p_{4k}   = ctrl[1].x  →  ∂c.x/∂p = B1, ∂c.y/∂p = 0
            grad[base + 0] += coeff * (gd[0] * B1);
            // p_{4k+1} = ctrl[1].y  →  ∂c.x/∂p = 0, ∂c.y/∂p = B1
            grad[base + 1] += coeff * (gd[1] * B1);
            // p_{4k+2} = ctrl[2].x
            grad[base + 2] += coeff * (gd[0] * B2);
            // p_{4k+3} = ctrl[2].y
            grad[base + 3] += coeff * (gd[1] * B2);
        }
    }
}

//  平滑项解析梯度(中心差分)
//  性能优化: 直接修改控制点避免curveFromParams重建, S从12降至6
//  参数布局: 每段仅ctrl[1]与ctrl[2]为优化变量, 共4个值 [x1,y1,x2,y2]
void PenaltyCost::addSmoothGrad(
    const BezierCurve& c, const VecXd& params, double w, VecXd& grad) const {
    if (w < 1e-15)
        return;

    if (c.segs.size() == 1) return;
    constexpr int S = 6;
    const double h = 1e-4;
    const double inv_2h = 1.0 / (2 * h);
    (void)params;  // 直接基于曲线c操作,不重建
    for (int k = 0; k < c.segs.size(); ++k) {
        int base = 4 * k;
        double ds = c.segs[k].arcLength(S) / S;
        // 仅对ctrl[1](j=1)与ctrl[2](j=2)计算梯度,对应params[base+0..3]
        for (int jj = 0; jj < 2; ++jj) {
            int j = jj + 1;  // ctrl index: 1 or 2
            // x方向梯度
            BezierCurve cp_x = c, cm_x = c;
            cp_x.segs[k].ctrl[j].x() += h;
            cm_x.segs[k].ctrl[j].x() -= h;
            if (c.segs.size() == 1) {
                enforceG1(cp_x, start_tan_dir, end_tan_dir, constrain_single_cubic_axes);
                enforceG1(cm_x, start_tan_dir, end_tan_dir, constrain_single_cubic_axes);
            }
            double fp_x = 0, fm_x = 0;
            for (int i = 0; i <= S; ++i) {
                double t = (double)i / S;
                double kp = cp_x.segs[k].curvature(t);
                double km = cm_x.segs[k].curvature(t);
                fp_x += kp * kp * ds;
                fm_x += km * km * ds;
            }
            grad[base + jj * 2 + 0] += w * (fp_x - fm_x) * inv_2h;
            // y方向梯度
            BezierCurve cp_y = c, cm_y = c;
            cp_y.segs[k].ctrl[j].y() += h;
            cm_y.segs[k].ctrl[j].y() -= h;
            if (c.segs.size() == 1) {
                enforceG1(cp_y, start_tan_dir, end_tan_dir, constrain_single_cubic_axes);
                enforceG1(cm_y, start_tan_dir, end_tan_dir, constrain_single_cubic_axes);
            }
            double fp_y = 0, fm_y = 0;
            for (int i = 0; i <= S; ++i) {
                double t = (double)i / S;
                double kp = cp_y.segs[k].curvature(t);
                double km = cm_y.segs[k].curvature(t);
                fp_y += kp * kp * ds;
                fm_y += km * km * ds;
            }
            grad[base + jj * 2 + 1] += w * (fp_y - fm_y) * inv_2h;
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
//  次要项(边界/围栏/簇/人行横道)的数值梯度
//  仅当这些项代价非零时运行 → 收敛后的干净迭代中跳过。
//  性能优化: 直接修改控制点而非curveFromParams重建
// ─────────────────────────────────────────────────────────────────────────────
void PenaltyCost::addNumericGrad(
        const VecXd &params, double wb, double wf, double wc, double wx, VecXd &grad) {
    BezierCurve c0 = full_param_mode ? curveFromParamsFull(params, proto) : curveFromParams(params, proto);
    enforceG1(c0, start_tan_dir, end_tan_dir, constrain_single_cubic_axes);
    bool nb = (wb > 1e-6 && evalBoundary(c0) > 1e-9);
    bool nf = (wf > 1e-6 && evalFence(c0) > 1e-9);
    bool nc = (wc > 1e-6 && evalCluster(c0) > 1e-9);
    bool nx = (wx > 1e-6 && evalCrosswalk(c0) > 1e-9);
    if (!nb && !nf && !nc && !nx)return;
    const double h = 1e-4;
    const double inv_2h = 1.0 / (2 * h);

    // 性能优化: 直接对c0的控制点扰动,避免每次curveFromParams重建
    auto evalAt = [&](const BezierCurve& cv) -> double {
        double f = 0;
        if (nb)f += wb * evalBoundary(cv);
        if (nf)f += wf * evalFence(cv);
        if (nc)f += wc * evalCluster(cv);
        if (nx)f += wx * evalCrosswalk(cv);
        return f;
    };

    // 非full_param_mode: 仅ctrl[1]与ctrl[2]是优化变量, layout [x1,y1,x2,y2]
    if (!full_param_mode) {
        for (int k = 0; k < (int)c0.segs.size(); ++k) {
            int base = 4 * k;
            for (int jj = 0; jj < 2; ++jj) {
                int j = jj + 1;  // ctrl index
                // x方向
                BezierCurve cp = c0, cm = c0;
                cp.segs[k].ctrl[j].x() += h;
                cm.segs[k].ctrl[j].x() -= h;
                if (c0.segs.size() == 1) {
                    enforceG1(cp, start_tan_dir, end_tan_dir, constrain_single_cubic_axes);
                    enforceG1(cm, start_tan_dir, end_tan_dir, constrain_single_cubic_axes);
                }
                grad[base + jj * 2 + 0] += (evalAt(cp) - evalAt(cm)) * inv_2h;
                // y方向
                BezierCurve cp_y = c0, cm_y = c0;
                cp_y.segs[k].ctrl[j].y() += h;
                cm_y.segs[k].ctrl[j].y() -= h;
                if (c0.segs.size() == 1) {
                    enforceG1(cp_y, start_tan_dir, end_tan_dir, constrain_single_cubic_axes);
                    enforceG1(cm_y, start_tan_dir, end_tan_dir, constrain_single_cubic_axes);
                }
                grad[base + jj * 2 + 1] += (evalAt(cp_y) - evalAt(cm_y)) * inv_2h;
            }
        }
    } else {
        // full_param_mode: 保留原逻辑以正确处理多段
        const int n = (int)params.size();
        VecXd p = params;
        for (int i = 0; i < n; ++i) {
            double orig = p[i];
            p[i] = orig + h;
            auto cp = curveFromParamsFull(p, proto);
            enforceG1(cp, start_tan_dir, end_tan_dir, constrain_single_cubic_axes);
            double fp = evalAt(cp);
            p[i] = orig - h;
            auto cm = curveFromParamsFull(p, proto);
            enforceG1(cm, start_tan_dir, end_tan_dir, constrain_single_cubic_axes);
            double fm = evalAt(cm);
            grad[i] += (fp - fm) * inv_2h;
            p[i] = orig;
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
//  主operator(): 一次前向计算 → f + grad
//  总代价: O(n_segs × S × 解析项) 而非 O(2n × full_eval)
// ─────────────────────────────────────────────────────────────────────────────
double PenaltyCost::operator()(const VecXd& params, VecXd& grad) {
    const int n = (int)params.size();
    grad.resize(n);
    grad.setZero();
    BezierCurve c = full_param_mode ? curveFromParamsFull(params, proto) : curveFromParams(params, proto);

    enforceG1(c, start_tan_dir, end_tan_dir, constrain_single_cubic_axes);

    // ── 代价计算 ──────────────────────────────────────────────────
    double f_smooth = evalSmooth(c);
    double f_obstacle = evalObstacle(c);
    double f_boundary = evalBoundary(c);
    double f_fence = evalFence(c);
    double f_cluster = evalCluster(c);
    double f_xwalk = evalCrosswalk(c);

    double f = weights.smooth * f_smooth
        + weights.obstacle * f_obstacle
        + weights.boundary * f_boundary
        + weights.fence * f_fence
        + weights.cluster * f_cluster
        + weights.crosswalk * f_xwalk;

    // smooth — 稀疏解析(每段4参数,S采样点)
    addSmoothGrad(c, params, weights.smooth, grad);

    // obstacle — 基于SDF梯度的解析项(主导且低开销)
    addObstacleGrad(c, params, weights.obstacle, grad);

    // 数值差分,当代价为0时跳过(收敛后常见)
    addNumericGrad(params, weights.boundary, weights.fence, weights.cluster, weights.crosswalk, grad);

    return f;
}

// ─────────────────────────────────────────────────────────────────────────────
//  optimiseCurve — 外层自适应权重循环
// ─────────────────────────────────────────────────────────────────────────────
BezierCurve optimiseCurve(
    PenaltyCost& cost, LBFGSSolver& solver, const BezierCurve& initial, int outer_iters) {
    cost.proto = initial;
    cost.buildCache();
    VecXd params = cost.full_param_mode ? curveToParamsFull(initial) : curveToParams(initial);
    // 性能优化: 限制outer_iters最大为2
    outer_iters = std::min(outer_iters, 2);
    for (int outer = 0; outer < outer_iters; ++outer) {
        BezierCurve current = cost.full_param_mode ?
                              curveFromParamsFull(params, cost.proto) : curveFromParams(params, cost.proto);
        enforceG1(current, cost.start_tan_dir, cost.end_tan_dir,
                  cost.constrain_single_cubic_axes);
        double op0 = cost.evalObstacle(current);
        double bp0 = cost.evalBoundary(current);
        double fp0 = cost.evalFence(current);
        double cp0 = cost.evalCluster(current);
        if (op0 + bp0 + fp0 + cp0 < 1e-6)
            break;

        auto res = solver.solve([&](const VecXd& p, VecXd& g) { return cost(p, g); }, params);
        params = res.x;
        BezierCurve c = cost.full_param_mode ?
                        curveFromParamsFull(params, cost.proto) : curveFromParams(params, cost.proto);
        enforceG1(c, cost.start_tan_dir, cost.end_tan_dir,
                  cost.constrain_single_cubic_axes);
        double op = cost.evalObstacle(c);
        double bp = cost.evalBoundary(c);
        double fp = cost.evalFence(c);
        double cp = cost.evalCluster(c);
        cost.weights.update(op, bp, fp, cp);

        if (op + bp + fp + cp < 1e-6)
            break; // 所有约束满足
    }
    BezierCurve result = cost.full_param_mode ?
        curveFromParamsFull(params, cost.proto) : curveFromParams(params, cost.proto);
    enforceG1(result, cost.start_tan_dir, cost.end_tan_dir,
              cost.constrain_single_cubic_axes);
    return result;
}

}
