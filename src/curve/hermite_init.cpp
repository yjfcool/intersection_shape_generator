#include "hermite_init.h"
#include "optimizer/sdf_field.h"
#include "curve/curve_utils.h"
#include "utils.h"
#include <cmath>
#include <algorithm>

namespace isg {

// ─────────────────────────────────────────────────────────────────────────────
//  共享辅助函数
// ─────────────────────────────────────────────────────────────────────────────
static bool segmentsIntersect_internal(
    const Vec2d& a, const Vec2d& b, const Vec2d& c, const Vec2d& d, Vec2d* out = nullptr) {
    Vec2d r = b - a, s = d - c;
    double den = cross2d(r, s);
    if (std::abs(den) < 1e-12) return false;
    Vec2d ac = c - a;
    double t = cross2d(ac, s) / den;
    double u = cross2d(ac, r) / den;
    if (t >= 0 && t <= 1 && u >= 0 && u <= 1) {
        if (out) *out = a + t * r;
        return true;
    }
    return false;
}

static constexpr double MAX_TAN_DEV = 60.0 * M_PI / 180.0;

static Vec2d clampTangent(const Vec2d& tan, const Vec2d& ref, double max_a) {
    Vec2d t = tan.norm() > 1e-10 ? tan.normalized() : ref;
    if (t.dot(ref) >= std::cos(max_a)) return t;
    double sg = cross2d(ref, t) >= 0 ? 1.0 : -1.0;
    double ca = std::cos(max_a), sa = std::sin(max_a);
    return Vec2d(ref[0] * ca - sg * ref[1] * sa,
                 ref[0] * sa * sg + ref[1] * ca);
}

// ─────────────────────────────────────────────────────────────────────────────
//  Level-1: 几何直接构造
// ─────────────────────────────────────────────────────────────────────────────
// SDF场中所有障碍物缓冲几何的AABB并集。
// 沿直线路径探测SDF以发现障碍物包围盒,
// 无需障碍物多边形数据(SDF已可用)。
struct ObstacleAABB {
    double x_min, x_max, y_min, y_max;
    bool valid = false;
};

static ObstacleAABB probeObstacleAABB(
    const SDFField& sdf, const Vec2d& p0, const Vec2d& p1, double clearance = 0.0) {
    // 在直线路径附近采样稠密网格
    ObstacleAABB box;
    box.x_min = box.y_min = 1e18;
    box.x_max = box.y_max = -1e18;

    int nx = 60, ny = 30;
    Vec2d along = (p1 - p0);
    double len = along.norm();
    if (len < 1e-6) return box;
    along = along * (1.0 / len);
    Vec2d perp{-along[1], along[0]};
    // 扫描宽度: 固定窄带(每侧一个车道宽=3.5m)
    // 避免捕获远偏路径轴的障碍物。大扫描范围
    // 在多障碍物时使AABB跨越整个路口,
    
    double sweep = 4.0; // 路径中心线两侧各4米。

    for (int ix = 0; ix <= nx; ++ix) {
        for (int iy = -ny / 2; iy <= ny / 2; ++iy) {
            double s = (double)ix / nx * len;
            double t = (double)iy / (ny / 2) * sweep;
            Vec2d pt = p0 + s * along + t * perp;
            std::pair<double, Vec2d> _q = sdf.queryWithGrad(pt);
            if (_q.first < clearance) {
                box.x_min = std::min(box.x_min, pt.x());
                box.x_max = std::max(box.x_max, pt.x());
                box.y_min = std::min(box.y_min, pt.y());
                box.y_max = std::max(box.y_max, pt.y());
                box.valid = true;
            }
        }
    }
    return box;
}

// 计算指定侧的绕行apex点(+1=路径左侧, -1=右侧)。
// apex位于障碍物最大横向投影+间隙处,
// 沿路径方向的障碍物纵向中点。
static Vec2d computeApex(
    const ObstacleAABB& box, const Vec2d& p0, const Vec2d& p1, int side, double clearance) {
    Vec2d along = (p1 - p0);
    double len = along.norm();
    if (len < 1e-6) return 0.5 * (p0 + p1);
    along = along * (1.0 / len);
    Vec2d perp{-along[1], along[0]}; // 路径左法向单位向量。

    // 障碍物包围盒纵向中点(世界坐标 → 投影到路径轴)
    // 包围盒存储采样点的世界坐标,需将
    
    // 轴对齐世界包围盒的角点:
    double corners_x[4] = {box.x_min, box.x_max, box.x_max, box.x_min};
    double corners_y[4] = {box.y_min, box.y_min, box.y_max, box.y_max};

    double lon_min = 1e18, lon_max = -1e18;
    double lat_min = 1e18, lat_max = -1e18;
    for (int i = 0; i < 4; ++i) {
        Vec2d c(corners_x[i], corners_y[i]);
        double lon = c.dot(along);
        double lat = c.dot(perp);
        lon_min = std::min(lon_min, lon);
        lon_max = std::max(lon_max, lon);
        lat_min = std::min(lat_min, lat);
        lat_max = std::max(lat_max, lat);
    }

    // 障碍物纵向中点, 钳制到路径[20%, 80%]区间
    double path_lon_0 = p0.dot(along);
    double path_lon_1 = p1.dot(along);
    double lon_mid_obs = 0.5 * (lon_min + lon_max);
    double frac = (lon_mid_obs - path_lon_0) / (path_lon_1 - path_lon_0 + 1e-12);
    frac = std::max(0.2, std::min(0.8, frac));

    // 障碍物在选定侧的横向边缘 + 间隙
    double lat_edge = (side > 0) ? (lat_max + clearance) : (lat_min - clearance);

    // apex: 从路径沿法向移动到所需横向位置
    Vec2d base = p0 + frac * (p1 - p0);
    double base_lat = base.dot(perp);
    Vec2d apex = base + (lat_edge - base_lat) * perp;
    return apex;
}

// 选择路径哪一侧(左/右)可绕行。
//
// 坐标约定:
//   along = (p1-p0).normalized()
//   perp  = (-along.y, along.x)，即路径左法向。
//   side = +1 表示沿perp向左绕行，side = -1 表示向右绕行。
//
// 右侧通行规则下优先选择更靠近路口中心的内侧短弧；但如果障碍位于路径与
// 路口中心之间，则必须改走外侧。对左右转，内侧方向可由
// sign(cross2d(t0, p1-p0)) 推断: 左转为+1，右转为-1。直行或缺少切向时
// 回退到SDF梯度投票选择间隙更大的一侧。
//
// t0_hint为可选进入切向，用于判断转向类型；零向量表示只使用梯度投票。
static int chooseSide(
    const SDFField& sdf, const ObstacleAABB& box, const Vec2d& p0, const Vec2d& p1, double clearance,
    const std::vector<std::vector<Vec2d>>& sibling_polys, const Vec2d& t0_hint = Vec2d(0, 0)) {
    Vec2d along = (p1 - p0).normalized();
    Vec2d perp{-along[1], along[0]}; // 路径左法向。

    // ── 绕行侧选择：优先内侧，同时感知障碍相对中心的位置 ───────────────
    // 1. 计算路口中心相对路径的横向位置。
    // 2. 计算障碍AABB中心相对路径的横向位置。
    // 3. 中心与障碍分处路径两侧时，障碍在外侧，向中心侧绕行。
    // 4. 中心与障碍同侧时，障碍阻塞内侧，改向外侧绕行。
    // 若障碍或中心方向不明确，则回退到SDF间隙投票。
    if (t0_hint.norm() > 1e-8 || box.valid) {
        Vec2d junction_centre(0.0, 0.0);
        Vec2d mid_path = 0.5 * (p0 + p1);
        Vec2d to_centre = junction_centre - mid_path;
        double centre_lat = to_centre.dot(perp); // 正值表示中心在路径左侧。

        // 障碍AABB中心。
        double obs_cx = 0.5 * (box.x_min + box.x_max);
        double obs_cy = 0.5 * (box.y_min + box.y_max);
        Vec2d obs_centre_world(obs_cx, obs_cy);
        double obs_lat = (obs_centre_world - mid_path).dot(perp);

        const double LAT_THRESH = 0.5;
        bool centre_clear = std::abs(centre_lat) > LAT_THRESH;
        bool obs_clear = std::abs(obs_lat) > LAT_THRESH;

        if (centre_clear && obs_clear) {
            // 中心和障碍都能明确判定在路径哪一侧。
            bool same_side = (centre_lat > 0) == (obs_lat > 0);
            if (!same_side) {
                // 障碍在外侧，朝中心侧绕行。
                return (centre_lat > 0) ? +1 : -1;
            } else {
                // 障碍在内侧/中心侧，强制改向外侧绕行。
                return (centre_lat > 0) ? -1 : +1;
            }
        }
        if (centre_clear && !obs_clear) {
            // 障碍近似位于路径正前方且未明确阻塞内侧时，优先中心侧绕行。
            return (centre_lat > 0) ? +1 : -1;
        }
        if (centre_clear) {
            // 只有中心方向明确时默认内侧绕行。
            return (centre_lat > 0) ? +1 : -1;
        }
    }
    // SDF梯度指向远离障碍的自由空间；累计其横向分量，选择间隙更大的绕行侧。
    // 常规clearance为0，因此使用固定探测阈值而不是按clearance放大。
    double weighted_lat = 0.0;
    double total_weight = 0.0;
    constexpr double PROBE = 0.15;
    constexpr double SAMPLE_THRESH = 3.0; // 只探测障碍3米范围内的点。
    constexpr int N = 30;
    for (int i = 1; i < N; ++i) {
        double t = (double)i / N;
        Vec2d pt = p0 + t * (p1 - p0);
        std::pair<double, Vec2d> _q = sdf.queryWithGrad(pt);
        double d = _q.first;
        if (d > SAMPLE_THRESH) continue;
        std::pair<double, Vec2d> _ql = sdf.queryWithGrad(pt + PROBE * perp);
        std::pair<double, Vec2d> _qm = sdf.queryWithGrad(pt - PROBE * perp);
        double d_lp = _ql.first;
        double d_lm = _qm.first;
        double lat_grad = (d_lp - d_lm) / (2 * PROBE);
        double w = 1.0 / (d + 0.05);
        weighted_lat += w * lat_grad;
        total_weight += w;
    }

    // 选择SDF间隙更大的一侧。
    int primary_side = (total_weight > 1e-10 && weighted_lat > 0) ? +1 : -1;

    // 兄弟曲线占用侧惩罚：若首选侧已有明显更多兄弟点，切换到另一侧。
    Vec2d mid = 0.5 * (p0 + p1);
    int sib_left = 0, sib_right = 0;
    for (auto& poly : sibling_polys) {
        for (auto& pt : poly) {
            double lat = (pt - mid).dot(perp);
            if (lat > 0.2) sib_left++;
            if (lat < -0.2) sib_right++;
        }
    }
    int occ_primary = (primary_side > 0) ? sib_left : sib_right;
    int occ_opposite = (primary_side > 0) ? sib_right : sib_left;
    if (occ_primary > occ_opposite + 5) primary_side = -primary_side;

    return primary_side;
}

// 构造平滑两段拱形: p0(t0) → apex(apex_tan) → p1(t1)。
// 拱顶切向由两侧连线平均得到，必要时限制为接近垂直于总方向。
static BezierCurve buildArch(
    const Vec2d& p0, const Vec2d& t0, const Vec2d& apex,
    const Vec2d& p1, const Vec2d& t1, double alpha = 0.4) {
    // 拱顶切向取进入/退出两侧方向平均，保证两段拼接平顺。
    Vec2d leg0 = (apex - p0).norm() > 1e-8 ? (apex - p0).normalized() : t0.normalized();
    Vec2d leg1 = (p1 - apex).norm() > 1e-8 ? (p1 - apex).normalized() : t1.normalized();
    Vec2d apex_tan = (leg0 + leg1);
    if (apex_tan.norm() < 1e-8) apex_tan = leg0;
    apex_tan.normalize();

    BezierSegment s0 = makeCubicG1(p0, t0.normalized(), apex, apex_tan, alpha);
    BezierSegment s1 = makeCubicG1(apex, apex_tan, p1, t1.normalized(), alpha);

    BezierCurve c;
    c.segs.push_back(s0);
    c.segs.push_back(s1);
    return c;
}

// Level-1入口：基于障碍AABB直接构造绕行拱形。
static BezierCurve geometricBypass(
    const Vec2d& p0, const Vec2d& t0, const Vec2d& p1, const Vec2d& t1,
    const SDFField& sdf, const Polygon2d& fence, double clearance = 0.0,
    const std::vector<std::vector<Vec2d>>& sibling_polys = {}) {
    // 使用同一clearance探测，默认只检测真实穿透。
    auto box = probeObstacleAABB(sdf, p0, p1, clearance);
    if (!box.valid) {
        // 未检测到障碍穿透时直接返回单段自然弧。
        BezierCurve c;
        c.segs.push_back(makeCubicG1(p0, t0.normalized(), p1, t1.normalized(), 0.4));
        return c;
    }

    // 传入进入切向以便chooseSide按转向类型选择内/外侧。
    int side = chooseSide(sdf, box, p0, p1, clearance, sibling_polys, t0);
    // 拱顶放在障碍边缘之外，并增加少量缓冲。
    double apex_clearance = clearance + 0.3;
    Vec2d apex = computeApex(box, p0, p1, side, apex_clearance);

    // 若拱顶仍落在障碍内，改试另一侧。
    std::pair<double, Vec2d> _qa = sdf.queryWithGrad(apex);
    if (_qa.first < -0.1) {
        apex = computeApex(box, p0, p1, -side, apex_clearance);
    }

    // 采样校验拱形是否避开障碍。
    BezierCurve arch = buildArch(p0, t0, apex, p1, t1);
    bool arch_clear = true;
    for (auto& seg : arch.segs) {
        for (int i = 1; i < 20; ++i) {
            std::pair<double, Vec2d> _qs = sdf.queryWithGrad(seg.evaluate((double)i / 20));
            if (_qs.first < clearance * 0.5) {
                arch_clear = false;
                break;
            }
        }
        if (!arch_clear) break;
    }

    if (!arch_clear) {
        // 首选侧不通时尝试另一侧。
        BezierCurve arch2 = buildArch(p0, t0, computeApex(box, p0, p1, -side, apex_clearance), p1, t1);
        bool arch2_clear = true;
        for (auto& seg : arch2.segs) {
            for (int i = 1; i < 20; ++i) {
                std::pair<double, Vec2d> _qs2 = sdf.queryWithGrad(seg.evaluate((double)i / 20));
                if (_qs2.first < clearance * 0.5) {
                    arch2_clear = false;
                    break;
                }
            }
            if (!arch2_clear) break;
        }
        if (!arch2_clear) return {}; // 通知调用方回退到Level-2。
        arch = arch2;
    }

    // 校验是否与已生成兄弟曲线发生内部交叉。只排除拱弧真实首尾连接点
    // 附近的数值噪声，不允许按弧长裁掉共享前缀/后缀。
    if (!sibling_polys.empty()) {
        auto arch_pts = arch.sampleByArcLength(20);
        constexpr double endpoint_skip_dist = 0.15;
        for (auto& sp : sibling_polys) {
            for (int ai = 0; ai + 1 < (int)arch_pts.size(); ++ai) {
                for (int si = 0; si + 1 < (int)sp.size(); ++si) {
                    Vec2d isect;
                    if (!segmentsIntersect_internal(
                        arch_pts[ai], arch_pts[ai + 1],
                        sp[si], sp[si + 1], &isect))
                        continue;
                    // 跳过靠近拱弧任一端点的交叉。
                    double de = std::min(
                        (isect - arch_pts.front()).norm(),
                        (isect - arch_pts.back()).norm());
                    if (de < endpoint_skip_dist) continue;
                    // 真正内部交叉时交给Level-2初始化兜底。
                    return {};
                }
            }
        }
    }

    return arch;
}

// ─────────────────────────────────────────────────────────────────────────────
//  Level-2: 全控制点几何初始化
//
//  不走网格寻路，直接构造可优化的合理拱形初值。
//  所有控制点(包括拼接点)都交给后续优化；拱顶由切向延长线解析给出。
// ─────────────────────────────────────────────────────────────────────────────
static BezierCurve geometricInitLevel2(
    const Vec2d& p0, const Vec2d& t0, const Vec2d& p1, const Vec2d& t1, const SDFField& sdf) {
    // 步骤1：选择安全的横向绕行方向，尽量远离障碍物云。
    Vec2d along = (p1 - p0);
    double len = along.norm();
    if (len < 1e-6) {
        BezierCurve c;
        c.segs.push_back(makeCubicG1(p0, t0.normalized(), p1, t1.normalized(), 0.4));
        return c;
    }
    along = along * (1.0 / len);
    Vec2d perp{-along[1], along[0]};

    // 右侧通行场景优先朝路口中心绕行；若障碍本身在中心侧，则改向外侧。
    // 通过障碍AABB估计其相对路径的横向位置。
    double side = 0.0;
    {
        Vec2d junction_centre(0.0, 0.0);
        Vec2d mid_path = 0.5 * (p0 + p1);
        double centre_lat = (junction_centre - mid_path).dot(perp);

        // 复用Level-1辅助函数探测障碍AABB。
        auto box = probeObstacleAABB(sdf, p0, p1, 0.0);

        if (box.valid && std::abs(centre_lat) > 0.5) {
            // 计算障碍AABB中心相对路径的横向位置。
            Vec2d box_world_centre((box.x_min + box.x_max) * 0.5,
                                   (box.y_min + box.y_max) * 0.5);
            double obs_lat = (box_world_centre - mid_path).dot(perp);

            // 绕行侧选择:
            //   障碍在外侧时朝中心绕行；
            //   障碍在中心侧且分离明显时强制向外；
            //   障碍近似在路径正前方时优先选择中心侧，因为内侧弧更短。
            const double ON_PATH_THRESH = 1.5;
            if (std::abs(obs_lat) < ON_PATH_THRESH) {
                // 正前方障碍优先走内侧(中心侧)。
                side = (centre_lat > 0) ? 1.0 : -1.0;
            } else {
                bool same_side = (centre_lat > 0) == (obs_lat > 0);
                side = (!same_side)
                           ? ((centre_lat > 0) ? 1.0 : -1.0) // 外侧障碍 → 内侧绕行
                           : ((centre_lat > 0) ? -1.0 : 1.0); // 内侧障碍 → 外侧绕行
            }
        } else if (std::abs(centre_lat) > 0.5) {
            side = (centre_lat > 0) ? 1.0 : -1.0; // 默认内侧绕行。
        } else {
            // 近直行时选择SDF间隙更大的一侧。
            double d_left = 0, d_right = 0;
            int n_probe = 10;
            for (int i = 1; i <= n_probe; ++i) {
                double s = (double)i / (n_probe + 1) * len;
                Vec2d mid_pt = p0 + s * along;
                std::pair<double, Vec2d> _ql2 = sdf.queryWithGrad(mid_pt + 1.5 * perp);
                std::pair<double, Vec2d> _qr2 = sdf.queryWithGrad(mid_pt - 1.5 * perp);
                d_left += _ql2.first;
                d_right += _qr2.first;
            }
            side = (d_left >= d_right) ? 1.0 : -1.0;
        }
    }
    // 安全检查：若偏好侧位于障碍内部，则切换到另一侧。
    {
        Vec2d mid_pt = p0 + 0.5 * (p1 - p0);
        std::pair<double, Vec2d> _qc = sdf.queryWithGrad(mid_pt + side * 1.5 * perp);
        std::pair<double, Vec2d> _qo = sdf.queryWithGrad(mid_pt - side * 1.5 * perp);
        if (_qc.first < 0.0 && _qo.first > _qc.first) {
            side = -side; // 偏好侧在障碍内，改走另一侧。
        }
    }

    // 步骤2：放置三点拱形，拱顶位于50%纵向位置，横向偏移按直线最小SDF估计。
    double min_d = 1e18;
    for (int i = 1; i < 20; ++i) {
        double t = (double)i / 20;
        std::pair<double, Vec2d> _qd = sdf.queryWithGrad((1 - t) * p0 + t * p1);
        min_d = std::min(min_d, _qd.first);
    }
    double lateral_needed = std::max(1.0, -min_d + 1.5); // 所需绕行距离。

    Vec2d apex = p0 + 0.5 * (p1 - p0) + side * lateral_needed * perp;

    // 步骤3：构造经过 p0 → q1 → apex → q2 → p1 的四段拱形。
    // q1/q2位于25%/75%纵向位置，并向拱顶方向抬升一部分。
    Vec2d q1 = p0 + 0.25 * (p1 - p0) + side * (lateral_needed * 0.6) * perp;
    Vec2d q2 = p0 + 0.75 * (p1 - p0) + side * (lateral_needed * 0.6) * perp;

    // 中间切向由Catmull-Rom点列估计。
    std::vector<Vec2d> pts = {p0, q1, apex, q2, p1};
    std::vector<Vec2d> tans = {t0.normalized(), {}, {}, {}, t1.normalized()};
    for (int i = 1; i <= 3; ++i) {
        Vec2d d = 0.5 * (pts[i + 1] - pts[i - 1]);
        tans[i] = d.norm() > 1e-10 ? d.normalized() : (pts[i + 1] - pts[i]).normalized();
        // 限制相对总方向的偏转，避免生成回环。
        tans[i] = clampTangent(tans[i], along, MAX_TAN_DEV);
    }

    return makeCurveFromKnots(pts, tans, 0.35);
}

// ─────────────────────────────────────────────────────────────────────────────
//  对外接口
// ─────────────────────────────────────────────────────────────────────────────
BezierCurve buildInitialCurve(
    const Vec2d& p0, const Vec2d& t0, const Vec2d& p1, const Vec2d& t1,
    const SDFField& sdf, const Polygon2d& fence, const std::vector<std::vector<Vec2d>>& sibling_polys) {
    // U型调头使用专用构造。
    if (angleBetween(t0, t1) > M_PI * 0.85)
        return buildTwoSegmentUTurn(p0, t0, p1, t1, sdf, fence);

    // 初始曲线只在穿入当前SDF避让域(SDF < 0)时触发绕行。
    // SDF避让域由输入mode决定: mode=2包含1m绕障缓冲,其他mode仅为原始障碍。
    constexpr double INIT_CLEARANCE = 0.0;

    // ── 快路径：自然单段弧已无障碍时直接保留 ───────────────────────
    // 直接检查Bezier弧，而不是检查端点弦线；左右转弦线可能穿障碍，
    // 但真实转向弧本身是安全的。
    {
        BezierSegment trial = makeCubicG1(p0, t0.normalized(), p1, t1.normalized(), 0.4);
        bool bezier_clear = true;
        int samples = std::max(40, (int)std::ceil(trial.arcLength(24) / 0.15));
        samples = std::min(samples, 320);
        for (int i = 1; i < samples; ++i) {
            std::pair<double, Vec2d> _qt = sdf.queryWithGrad(trial.evaluate((double)i / samples));
            if (_qt.first < INIT_CLEARANCE) {
                bezier_clear = false;
                break;
            }
        }
        if (bezier_clear) {
            BezierCurve c;
            c.segs.push_back(trial);
            return c;
        }
    }

    // ── Level-1：几何直接绕行构造 ─────────────────────────────────
    {
        BezierCurve arch = geometricBypass(p0, t0, p1, t1, sdf, fence, INIT_CLEARANCE, sibling_polys);
        if (!arch.empty()) return arch;
    }

    // ── Level-2：全控制点几何初值 ─────────────────────────────────
    return geometricInitLevel2(p0, t0, p1, t1, sdf);
}

// ─────────────────────────────────────────────────────────────────────────────
//  U型调头：基于首尾射线对齐的单段三次Bezier。
//
//  旧实现会在p0和p1之间插入拱顶节点；近反向车道下控制多边形容易翻转，
//  形成可见折点或S形。当前实现只保留端点控制: 先沿进入/退出反向射线对齐
//  纵向延长，再用平齐横向间距的2/3作为近圆弧把手长度。
// ─────────────────────────────────────────────────────────────────────────────
BezierCurve buildTwoSegmentUTurn(
    const Vec2d& p0, const Vec2d& t0,
    const Vec2d& p1, const Vec2d& t1, const SDFField& sdf, const Polygon2d&) {
    (void)sdf;
    BezierCurve c;
    c.segs.push_back(makeAlignedUTurnCubic(p0, t0, p1, t1));
    return c;
}

}
