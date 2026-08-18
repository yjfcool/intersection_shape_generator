#include "cluster_order.h"
#include "optimizer/sdf_field.h"
#include "curve/curve_utils.h"
#include "utils.h"
#include <algorithm>
#include <cmath>
#include <numeric>

namespace isg {

// ─── 车道几何辅助 ────────────────────────────────────────────
static const Connectivity *findConn(const std::vector<Connectivity> &cs, const ConnId &id) {
    for (auto &c:cs) if (c.id == id)return &c;
    return nullptr;
}

static const Lane *findLane(const std::vector<Lane> &ls, const LaneId &id) {
    for (auto &l:ls) if (l.id == id)return &l;
    return nullptr;
}

/// 规范化连通关系的组ID：从LaneGroup与Lane的groupId字段反推enterGroupId/exitGroupId
static void normalizeConnectivityGroups(std::vector<Connectivity> &conns,
                                        const std::vector<Lane> &lanes,
                                        const std::vector<LaneGroup> &laneGroups) {
    std::unordered_map<LaneId, LaneGroupId> entry_gid;
    std::unordered_map<LaneId, LaneGroupId> exit_gid;
    for (auto &g : laneGroups) {
        for (auto &lid : g.lanes) {
            if (g.role == GroupRole::Entry)
                entry_gid[lid] = g.id;
            else
                exit_gid[lid] = g.id;
        }
    }
    for (auto &l : lanes) {
        if (!l.groupId.empty()) {
            if (!entry_gid.count(l.id))
                entry_gid[l.id] = l.groupId;
            if (!exit_gid.count(l.id))
                exit_gid[l.id] = l.groupId;
        }
    }
    for (auto &c : conns) {
        if (c.enterGroupId.empty()) {
            auto it = entry_gid.find(c.entry_lane_id);
            if (it != entry_gid.end())
                c.enterGroupId = it->second;
        }
        if (c.exitGroupId.empty()) {
            auto it = exit_gid.find(c.exit_lane_id);
            if (it != exit_gid.end())
                c.exitGroupId = it->second;
        }
    }
}

/// 进入车道端点：geometry末点
static Vec2d entryEndpoint(const Lane *l) {
    if (!l || l->geometry.points.empty())return Vec2d(0, 0);
    return l->geometry.points.back();
}

/// 退出车道端点：geometry首点
static Vec2d exitEndpoint(const Lane *l) {
    if (!l || l->geometry.points.empty())return Vec2d(0, 0);
    return l->geometry.points.front();
}

/// 进入车道末点切向（指向路口内）
static Vec2d entryTangent(const Lane *l) {
    if (!l || l->geometry.points.size() < 2)return Vec2d(0, 1);
    auto &p = l->geometry.points;
    Vec2d d = p.back() - p[p.size() - 2];
    return d.norm() > 1e-9 ? d.normalized() : Vec2d(0, 1);
}

/// 退出车道首点切向（指向路口外）
static Vec2d exitTangent(const Lane *l) {
    if (!l || l->geometry.points.size() < 2)return Vec2d(0, 1);
    auto &p = l->geometry.points;
    Vec2d d = p[1] - p[0];
    return d.norm() > 1e-9 ? d.normalized() : Vec2d(0, 1);
}

static bool nearPoint(const Vec2d& a, const Vec2d& b, double tol = 0.35) {
    return (a - b).squaredNorm() <= tol * tol;
}

static bool fixedGeometryEndpointDirection(
    const Connectivity& c, const Vec2d& endpoint, bool as_entry_dir, Vec2d& out_dir) {
    if (!c.fixed_shape)
        return false;
    const auto& pts = c.geometry.points;
    if (pts.size() < 2)
        return false;

    Vec2d dir(0, 0);
    if (nearPoint(endpoint, pts.front())) {
        dir = as_entry_dir ? (pts[1] - pts.front()) : (pts.front() - pts[1]);
    } else if (nearPoint(endpoint, pts.back())) {
        dir = as_entry_dir ? (pts[pts.size() - 2] - pts.back())
                           : (pts.back() - pts[pts.size() - 2]);
    } else {
        return false;
    }

    if (dir.norm() < 1e-9)
        return false;
    out_dir = dir.normalized();
    return true;
}

static Vec2d connectionDirection(
    const Connectivity& conn, const std::vector<Lane>& lanes, bool use_entry) {
    const Lane* lane = use_entry ? findLane(lanes, conn.entry_lane_id)
                                 : findLane(lanes, conn.exit_lane_id);
    Vec2d endpoint = use_entry ? entryEndpoint(lane) : exitEndpoint(lane);
    Vec2d fixed_dir;
    if (fixedGeometryEndpointDirection(conn, endpoint, use_entry, fixed_dir))
        return fixed_dir;
    return use_entry ? entryTangent(lane) : exitTangent(lane);
}

// ─── 组轴辅助 ────────────────────────────────────────────────
/// 计算指定簇内所有conn的进入/退出端点几何中心
static Vec2d meanPt(const std::vector<ConnId> &cids,
                    const std::vector<Connectivity> &conns,
                    const std::vector<Lane> &lanes, bool use_entry) {
    Vec2d s(0, 0);
    int n = 0;
    for (auto &cid:cids) {
        auto *c = findConn(conns, cid);
        if (!c) continue;
        const Lane *l = use_entry ? findLane(lanes, c->entry_lane_id) : findLane(lanes, c->exit_lane_id);
        s += use_entry ? entryEndpoint(l) : exitEndpoint(l);
        ++n;
    }
    return n > 0 ? (s / n) : Vec2d(0, 0);
}

/// 计算指定簇的平均切向左法向，作为横向参考轴
static Vec2d armLeftNormal(const std::vector<ConnId> &cids,
                           const std::vector<Connectivity> &conns,
                           const std::vector<Lane> &lanes, bool use_entry) {
    Vec2d fixed_mt(0, 0);
    int fixed_n = 0;
    Vec2d mt(0, 0);
    int n = 0;
    for (auto &cid:cids) {
        auto *c = findConn(conns, cid);
        if (!c) continue;
        const Lane *l = use_entry ? findLane(lanes, c->entry_lane_id) : findLane(lanes, c->exit_lane_id);
        Vec2d endpoint = use_entry ? entryEndpoint(l) : exitEndpoint(l);
        Vec2d fixed_dir;
        if (fixedGeometryEndpointDirection(*c, endpoint, use_entry, fixed_dir)) {
            fixed_mt += fixed_dir;
            ++fixed_n;
        }
        mt += use_entry ? entryTangent(l) : exitTangent(l);
        ++n;
    }
    if (fixed_n > 0 && fixed_mt.norm() > 1e-9) {
        fixed_mt.normalize();
        return Vec2d(-fixed_mt[1], fixed_mt[0]);
    }
    if (n == 0)return Vec2d(0, 1);
    mt /= n;
    if (mt.norm() < 1e-9)mt = Vec2d(0, 1);
    mt.normalize();
    return Vec2d(-mt[1], mt[0]);
}

// ─── 配对特定横向参考轴 ──────────────────────────────────────
/// 取两条曲线进入端点中点与退出端点中点的连线方向，返回其左法向
/// 作为该配对专用的横向参考轴，使排序约束在曲线对自身坐标系下度量
static Vec2d pairRefPerp(const Connectivity &ca, const Connectivity &cb,
                         const std::vector<Lane> &lanes) {
    const Lane* ea = findLane(lanes, ca.entry_lane_id);
    const Lane* xa = findLane(lanes, ca.exit_lane_id);

    Vec2d dir_sum(0, 0);
    int dir_count = 0;
    auto add_fixed_dir = [&](const Connectivity& conn, const Vec2d& endpoint, bool use_entry) {
        Vec2d d;
        if (fixedGeometryEndpointDirection(conn, endpoint, use_entry, d)) {
            dir_sum += d;
            ++dir_count;
        }
    };

    if (!ca.entry_lane_id.empty() && ca.entry_lane_id == cb.entry_lane_id) {
        Vec2d ep = entryEndpoint(ea);
        add_fixed_dir(ca, ep, true);
        add_fixed_dir(cb, ep, true);
    } else if (!ca.exit_lane_id.empty() && ca.exit_lane_id == cb.exit_lane_id) {
        Vec2d xp = exitEndpoint(xa);
        add_fixed_dir(ca, xp, false);
        add_fixed_dir(cb, xp, false);
    }
    if (dir_count > 0 && dir_sum.norm() > 1e-9) {
        dir_sum.normalize();
        return Vec2d(-dir_sum[1], dir_sum[0]);
    }

    Vec2d me = 0.5 * (entryEndpoint(findLane(lanes, ca.entry_lane_id))
                      + entryEndpoint(findLane(lanes, cb.entry_lane_id)));
    Vec2d mx = 0.5 * (exitEndpoint(findLane(lanes, ca.exit_lane_id))
                      + exitEndpoint(findLane(lanes, cb.exit_lane_id)));
    Vec2d al = mx - me;
    double n = al.norm();
    if (n < 1e-9)return Vec2d(0, 1);
    al.normalize();
    return Vec2d(-al[1], al[0]);
}

// ─── 排序键计算 ──────────────────────────────────────────────
/// 进入簇排序键：退出端点相对进入组中心点的方位角
/// 降序排列 → CCW在前 = 左侧在前
static double entryAngle(const Connectivity &c, const Vec2d &me, const std::vector<Lane> &lanes) {
    Vec2d diff = exitEndpoint(findLane(lanes, c.exit_lane_id)) - me;
    return std::atan2(diff[1], diff[0]);
}

/// 退出簇排序键：进入端点在退出arm左法向上的投影
/// 降序排列 → 左侧在前
static double exitEntryLat(const Connectivity &c, const Vec2d &ealn, const std::vector<Lane> &lanes) {
    return entryEndpoint(findLane(lanes, c.entry_lane_id)).dot(ealn);
}

// ─── 配对管理 ────────────────────────────────────────────────
bool ClusterOrderSolver::hasPair(
    const std::vector<CurvePair>& pairs, const ConnId& a, const ConnId& b) {
    for (auto& p : pairs)
        if ((p.id_a == a && p.id_b == b) || (p.id_a == b && p.id_b == a))
            return true;
    return false;
}

// ── O(1)查找的哈希索引 ───────────────────────────────────────
// key为配对中字典序较小的ConnId，value为该ConnId作为首key的配对索引列表
// 查询时仅需扫描该桶中的1~3个配对，避免热循环O(pairs_)线性扫描
void ClusterOrderSolver::buildPairIndex() {
    pair_index_.clear();
    for (int i = 0; i < (int)pairs_.size(); ++i) {
        const ConnId key = (pairs_[i].id_a < pairs_[i].id_b)
                           ? pairs_[i].id_a : pairs_[i].id_b;
        pair_index_[key].push_back(i);
    }
}

const CurvePair* ClusterOrderSolver::findPair(const ConnId& a, const ConnId& b) const {
    const ConnId key = (a < b) ? a : b;
    auto it = pair_index_.find(key);
    if (it == pair_index_.end())
        return nullptr;
    for (int idx : it->second) {
        const auto& p = pairs_[idx];
        if ((p.id_a == a && p.id_b == b) || (p.id_a == b && p.id_b == a))
            return &p;
    }
    return nullptr;
}

/// 从排序后的簇构造配对：每两个conn形成CurvePair
/// 排序在前(id_a)的预期位于左侧 → expected_side=+1
void ClusterOrderSolver::addPairsFromSortedCluster(
        const std::vector<ConnId> &sorted_ids,
        const std::vector<Connectivity> &conns,
        const std::vector<Lane> &lanes) {
    for (int i = 0; i < (int) sorted_ids.size(); ++i) {
        for (int j = i + 1; j < (int) sorted_ids.size(); ++j) {
            if (hasPair(pairs_, sorted_ids[i], sorted_ids[j]))continue;
            auto *ca = findConn(conns, sorted_ids[i]);
            auto *cb = findConn(conns, sorted_ids[j]);
            if (!ca || !cb)continue;
            CurvePair p;
            p.id_a = sorted_ids[i];
            p.id_b = sorted_ids[j];
            p.expected_side = +1;           // sorted[i] 在 sorted[j] 左侧
            p.ref_perp = pairRefPerp(*ca, *cb, lanes);

            // 共享退出车道：两条曲线在退出端点必然汇聚，
            // 仅允许在该端点重合，中段必须保持非交
            if (!ca->exit_lane_id.empty() &&
                ca->exit_lane_id == cb->exit_lane_id) {
                p.exempt = CrossExemption::None;
                p.shared_endpoint = true;   // 仅标识共享连接点，端点外仍受非交约束
                pairs_.push_back(p);
                continue;
            }

            // 共享进入车道：两条曲线从同一进入端点扇出，
            // 端点附近切向相同，约束梯度近零；仅标识合法连接点
            if (!ca->entry_lane_id.empty() &&
                ca->entry_lane_id == cb->entry_lane_id) {
                p.exempt = CrossExemption::None;
                p.shared_endpoint = true;
            } else {
                p.exempt = CrossExemption::None;
            }
            pairs_.push_back(p);
        }
    }
}

// ─── 拓扑倒置检测 ────────────────────────────────────────────
// 检测同簇配对中理论上"必须交叉"的拓扑倒置情形，标记为StructuralCross豁免
//
// 基本逻辑:
//   同进入组对(A,B):
//     rank_diff = rank_A - rank_B  (排序: 角度降序 → rank 0 = 最左)
//     lat_diff  = exit_lat_A - exit_lat_B  (进入arm法向坐标系下)
//     一致: rank_diff * lat_diff < 0  → 无需结构性交叉
//     倒置: rank_diff * lat_diff > 0  → StructuralCross
//
//   同退出组对(A,B):
//     rank_diff = rank_A - rank_B  (排序: entry_lat 降序 → rank 0 = 最左)
//     lat_diff  = exit_lat_in_exit_group_A - exit_lat_in_exit_group_B
//                 (退出端点在退出arm左法向投影，相对mean_exit点——
//                  与排序键独立，否则永远<0导致检测失效)
//     倒置: rank_diff * lat_diff > 0  → StructuralCross
void ClusterOrderSolver::detectTopologicalInversions(
        const std::vector<Connectivity> &conns, const std::vector<Lane> &lanes,
        const std::vector<Crosswalk> &crosswalks) {
    constexpr double LAT_EPS = 0.3;        // 米 — 忽略近相等横向距离
    constexpr double OPP_SIDE_THRESH = 3.0; // 米 — 对侧检测最小横向幅度
    (void)crosswalks;

    // 辅助函数: 从几何方位计算转向类型
    //   t0 = 进入车道切向, t1 = 退出车道切向, path = (p1-p0).normalized()
    //   cross2d(t0, path) > 0.5 → 左转
    //   cross2d(t0, path) < -0.5 → 右转
    //   |cross2d(t0, path)| < 0.5 且 t0·t1 > 0 → 直行
    //   t0·t1 < -0.5 → U型调头
    // 返回: 'U'=调头, 'L'=左转, 'R'=右转, 'S'=直行, '?'=未知
    auto geometricTurnType = [&](const Connectivity* c) -> char {
        const Lane* el = findLane(lanes, c->entry_lane_id);
        const Lane* xl = findLane(lanes, c->exit_lane_id);
        if (!el || !xl) return '?';
        Vec2d t0 = connectionDirection(*c, lanes, true);
        Vec2d t1 = connectionDirection(*c, lanes, false);
        Vec2d p0 = entryEndpoint(el);
        Vec2d p1 = exitEndpoint(xl);
        Vec2d path = p1 - p0;
        if (path.norm() < 1e-6) return '?';
        path.normalize();
        double dot = t0.dot(t1);
        if (dot < -0.5) return 'U';
        double cross = cross2d(t0, path);
        if (cross > 0.5) return 'L';
        if (cross < -0.5) return 'R';
        return 'S';
    };

    for (auto &pair:pairs_) {
        if (pair.exempt != CrossExemption::None)continue;
        auto *ca = findConn(conns, pair.id_a);
        auto *cb = findConn(conns, pair.id_b);
        if (!ca || !cb)continue;

        bool same_entry = (ca->enterGroupId == cb->enterGroupId);
        bool same_exit = (ca->exitGroupId == cb->exitGroupId);

        bool a_uturn = (is_uturn_.count(pair.id_a) && is_uturn_.at(pair.id_a));
        bool b_uturn = (is_uturn_.count(pair.id_b) && is_uturn_.at(pair.id_b));

        // 是否共享实际车道端点（同进入车道或同退出车道）
        bool share_entry_lane =
            (!ca->entry_lane_id.empty() && ca->entry_lane_id == cb->entry_lane_id);
        bool share_exit_lane =
            (!ca->exit_lane_id.empty() && ca->exit_lane_id == cb->exit_lane_id);
        bool share_lane_endpoint = share_entry_lane || share_exit_lane;

        char tA = geometricTurnType(ca);
        char tB = geometricTurnType(cb);
        bool a_is_uturn = a_uturn || tA == 'U';
        bool b_is_uturn = b_uturn || tB == 'U';
        bool a_is_left_right =
            ca->turn_type == ConnTurnType::TurnLeft ||
            ca->turn_type == ConnTurnType::TurnRight ||
            tA == 'L' || tA == 'R';
        bool b_is_left_right =
            cb->turn_type == ConnTurnType::TurnLeft ||
            cb->turn_type == ConnTurnType::TurnRight ||
            tB == 'L' || tB == 'R';
        bool uturn_turn_combo =
            (a_is_uturn && b_is_left_right) ||
            (b_is_uturn && a_is_left_right);

        // 业务定义唯一允许的同簇非端点相交组合: 掉头 + 左/右转。
        // 其它同连接点组合(直行/普通左右转/掉头之间)仍必须满足非端点不相交。
        if (uturn_turn_combo) {
            pair.exempt = CrossExemption::StructuralCross;
            continue;
        }

        // 共端点的: 仅端点相交,在共享进入处扇出或在共享退出处汇聚,
        // 不构成拓扑必须交叉,保持约束配对。U型调头与左/右转组合
        // 已在上方按业务定义作为唯一可非端点相交组合豁免。
        if (share_lane_endpoint) {
            // ── 共享端点的两个U型调头：可按"内层浅拱、外层深拱"有序排列，不属于结构性交叉 ──
            // 此处保持配对为受约束(shared_endpoint=true, 在addPairsFromSortedCluster中已设置)，
            // 端点附近只按连接点数值容差排除，端点外贴合/重叠/相交仍禁止；
            // U型调头多约束求解器主动搜索lead0/lateral_bias/handle_bias候选,保持两条拱弧的横向有序。
            continue;
        }

        // 非共享实际端点且入口行驶方向相反的 pair 属于对向交通流。
        // 这类曲线即便被同一进入/退出组归到同簇，也不属于同向线束的
        // “分散/合流后不得相交”约束；例如 100000643 的 78-97/99。
        // 同向的直行-掉头、掉头-掉头、同类型组合仍保留约束。
        Vec2d entry_dir_a = connectionDirection(*ca, lanes, true);
        Vec2d entry_dir_b = connectionDirection(*cb, lanes, true);
        bool opposite_entry_flow =
            entry_dir_a.norm() > 1e-8 && entry_dir_b.norm() > 1e-8 &&
            entry_dir_a.normalized().dot(entry_dir_b.normalized()) < -0.5;
        if (opposite_entry_flow) {
            pair.exempt = CrossExemption::StructuralCross;
            continue;
        }

        // ── 弦交叉检测 (连接关系逻辑交叉豁免) ──────────────────────────
        // 修复 BUG 1 衍生问题: 同退出组不同入口arm的右转(如 conn 21↔26, 23↔26
        // in intersection_jd.json) 从不同方向进入同一退出组, 路径必然交叉。
        // 现有 same_exit rank-lat 检测使用退出arm左法向投影, 无法捕捉跨arm
        // 几何交叉。弦交叉检测直接判断两条连通的直线弦是否相交, 若相交则
        // 任何连续曲线也必然交叉 (Jordan曲线定理推论), 属于"连接关系存在
        // 逻辑交叉"豁免场景。
        //
        // 仅对非共端点对应用, 弦交点需距所有端点 > 1.5m (排除端点附近相切)。
        if (!share_lane_endpoint) {
            const Lane* ea = findLane(lanes, ca->entry_lane_id);
            const Lane* xa = findLane(lanes, ca->exit_lane_id);
            const Lane* eb = findLane(lanes, cb->entry_lane_id);
            const Lane* xb = findLane(lanes, cb->exit_lane_id);
            if (ea && xa && eb && xb) {
                Vec2d p0a = entryEndpoint(ea), p1a = exitEndpoint(xa);
                Vec2d p0b = entryEndpoint(eb), p1b = exitEndpoint(xb);
                Vec2d isect;
                if (segmentsIntersect(p0a, p1a, p0b, p1b, &isect)) {
                    double min_ep_dist = std::min({
                        (isect - p0a).norm(), (isect - p1a).norm(),
                        (isect - p0b).norm(), (isect - p1b).norm()
                    });
                    if (min_ep_dist > 1.5) {
                        pair.exempt = CrossExemption::StructuralCross;
                        continue;
                    }
                }
            }
        }

        if (same_entry) {
            auto ira = entry_cluster_rank_.find(pair.id_a);
            auto irb = entry_cluster_rank_.find(pair.id_b);
            if (ira == entry_cluster_rank_.end() || irb == entry_cluster_rank_.end()) continue;
            double rank_diff = (double)(ira->second - irb->second);

            // ── 对侧退出 = 必然结构性交叉 ──
            // 如左转与右转永远交叉
            {
                auto ixa = exit_lat_sign_.find(pair.id_a);
                auto ixb = exit_lat_sign_.find(pair.id_b);
                if (ixa != exit_lat_sign_.end() && ixb != exit_lat_sign_.end()) {
                    double la = ixa->second, lb = ixb->second;
                    if (la * lb < 0.0 &&
                        std::abs(la) > OPP_SIDE_THRESH &&
                        std::abs(lb) > OPP_SIDE_THRESH) {
                        pair.exempt = CrossExemption::StructuralCross;
                        continue;
                    }
                }
            }

            // ── 进入-退出横向倒置检测 ──
            // 捕捉进入侧车道顺序与退出侧车道顺序倒置的情况
            // (如外右转来自内车道 必须与 内右转来自外车道 交叉;
            //  同west arm的直行东向 vs 右转南向)
            // 判据: (entry_lat_A - entry_lat_B) × (exit_lat_A - exit_lat_B) < 0
            {
                auto iea = entry_lat_in_entry_ref_.find(pair.id_a);
                auto ieb = entry_lat_in_entry_ref_.find(pair.id_b);
                bool share_exit = (!ca->exitGroupId.empty() &&
                                   ca->exitGroupId == cb->exitGroupId);
                if (iea != entry_lat_in_entry_ref_.end() &&
                    ieb != entry_lat_in_entry_ref_.end()) {
                    double entry_diff = iea->second - ieb->second;
                    if (std::abs(entry_diff) > LAT_EPS) {
                        double exit_diff = 0.0;
                        bool got_exit = false;
                        if (share_exit) {
                            // 同退出组: 用退出组内横向位置
                            auto ixa = exit_lat_in_exit_group_ref_.find(pair.id_a);
                            auto ixb = exit_lat_in_exit_group_ref_.find(pair.id_b);
                            if (ixa != exit_lat_in_exit_group_ref_.end() &&
                                ixb != exit_lat_in_exit_group_ref_.end()) {
                                exit_diff = ixa->second - ixb->second;
                                got_exit = true;
                            }
                        } else {
                            // 不同退出组: 用退出横向符号
                            auto ixa = exit_lat_sign_.find(pair.id_a);
                            auto ixb = exit_lat_sign_.find(pair.id_b);
                            if (ixa != exit_lat_sign_.end() &&
                                ixb != exit_lat_sign_.end()) {
                                exit_diff = ixa->second - ixb->second;
                                got_exit = true;
                            }
                        }
                        if (got_exit && std::abs(exit_diff) > LAT_EPS) {
                            if (entry_diff * exit_diff < 0.0) {
                                pair.exempt = CrossExemption::StructuralCross;
                                continue;
                            }
                        }
                    }
                }
            }

            // ── rank-lat 倒置检测 ──
            // 同退出组: 用退出组内横向位置(与排序键独立)
            // 不同退出组: 用退出横向(进入arm坐标系下)
            bool share_exit = (!ca->exitGroupId.empty() && ca->exitGroupId == cb->exitGroupId);
            double lat_diff;
            if (share_exit) {
                auto ixa = exit_lat_in_exit_group_ref_.find(pair.id_a);
                auto ixb = exit_lat_in_exit_group_ref_.find(pair.id_b);
                if (ixa == exit_lat_in_exit_group_ref_.end() ||
                    ixb == exit_lat_in_exit_group_ref_.end()) continue;
                lat_diff = ixa->second - ixb->second;
            } else {
                auto ixa = exit_lat_in_entry_ref_.find(pair.id_a);
                auto ixb = exit_lat_in_entry_ref_.find(pair.id_b);
                if (ixa == exit_lat_in_entry_ref_.end() ||
                    ixb == exit_lat_in_entry_ref_.end()) continue;
                lat_diff = ixa->second - ixb->second;
            }
            if (std::abs(lat_diff) < LAT_EPS) continue;
            if (rank_diff * lat_diff > 0)
                pair.exempt = CrossExemption::StructuralCross;

        } else if (same_exit) {
            // 同退出组: 用退出组内横向位置(与排序键独立)
            // 若用entry_lat_in_exit_ref_ (即排序键本身) → rank_diff×lat_diff 恒<0 → 永不触发
            auto ira = exit_cluster_rank_.find(pair.id_a);
            auto irb = exit_cluster_rank_.find(pair.id_b);
            auto ixa = exit_lat_in_exit_group_ref_.find(pair.id_a);
            auto ixb = exit_lat_in_exit_group_ref_.find(pair.id_b);
            if (ira == exit_cluster_rank_.end() || irb == exit_cluster_rank_.end()) continue;
            if (ixa == exit_lat_in_exit_group_ref_.end() ||
                ixb == exit_lat_in_exit_group_ref_.end()) continue;
            double rank_diff = (double)(ira->second - irb->second);
            double lat_diff  = ixa->second - ixb->second;
            if (std::abs(lat_diff) < LAT_EPS) continue;
            if (rank_diff * lat_diff > 0)
                pair.exempt = CrossExemption::StructuralCross;
        }
    }
}

// ─── 主构建入口 ──────────────────────────────────────────────
void ClusterOrderSolver::build(
    const std::vector<Connectivity>& conns,
    const std::vector<Lane>& lanes, const std::vector<LaneGroup>& laneGroups,
    const std::vector<Crosswalk>& crosswalks) {
    std::vector<Connectivity> norm_conns = conns;
    normalizeConnectivityGroups(norm_conns, lanes, laneGroups);

    entry_group_order_.clear();
    exit_group_order_.clear();
    pairs_.clear();
    is_uturn_.clear();
    exit_lat_sign_.clear();
    entry_lat_in_entry_ref_.clear();
    exit_lat_in_entry_ref_.clear();
    entry_lat_in_exit_ref_.clear();
    exit_lat_in_exit_group_ref_.clear();
    entry_cluster_rank_.clear();
    exit_cluster_rank_.clear();

    // 步骤1: 按enterGroupId/exitGroupId分组。
    for (auto& c : norm_conns) {
        if (!c.enterGroupId.empty())entry_group_order_[c.enterGroupId].push_back(c.id);
        if (!c.exitGroupId.empty())exit_group_order_[c.exitGroupId].push_back(c.id);
    }

    // 步骤2: 预计算每个退出组内退出车道的横向位置，
    // 用于进入簇排序时同向多车道的稳定并列打破。
    std::unordered_map<LaneId, double> exit_lane_lat_in_group;
    for (auto &kv:exit_group_order_) {
        auto &gid = kv.first;
        auto &cids = kv.second;
        if (cids.size() < 2)continue;
        Vec2d perp = armLeftNormal(cids, norm_conns, lanes, false);
        Vec2d ref = meanPt(cids, norm_conns, lanes, false);
        for (auto &cid:cids) {
            auto *c = findConn(norm_conns, cid);
            if (!c)continue;
            const Lane *xl = findLane(lanes, c->exit_lane_id);
            exit_lane_lat_in_group[c->exit_lane_id] = (exitEndpoint(xl) - ref).dot(perp);
        }
    }

    // 步骤3: 预计算每个进入组内进入车道的横向位置。
    std::unordered_map<LaneId, double> entry_lane_lat_in_group;
    for (auto &kv:entry_group_order_) {
        auto &gid = kv.first;
        auto &cids = kv.second;
        if (cids.size() < 2)continue;
        Vec2d perp = armLeftNormal(cids, norm_conns, lanes, true);
        Vec2d ref = meanPt(cids, norm_conns, lanes, true);
        for (auto &cid:cids) {
            auto *c = findConn(norm_conns, cid);
            if (!c)continue;
            const Lane *el = findLane(lanes, c->entry_lane_id);
            entry_lane_lat_in_group[c->entry_lane_id] = (entryEndpoint(el) - ref).dot(perp);
        }
    }

    // 步骤4: 进入簇排序。
    //   主键:   (exit_pt - mean_entry_pt) 角度降序 — 全局左→右排序
    //   并列1:  退出车道在退出组内横向位置降序 — 同向多车道稳定排序
    //   并列2:  进入车道在进入组内横向位置降序
    for (auto &kv:entry_group_order_) {
        auto &gid = kv.first;
        auto &cids = kv.second;
        if (cids.size() < 2)continue;
        Vec2d me = meanPt(cids, norm_conns, lanes, true);
        std::stable_sort(cids.begin(), cids.end(), [&](const ConnId &a, const ConnId &b) {
            auto *ca = findConn(norm_conns, a);
            auto *cb = findConn(norm_conns, b);
            if (!ca || !cb)return false;
            double ang_a = entryAngle(*ca, me, lanes);
            double ang_b = entryAngle(*cb, me, lanes);
            if (std::abs(ang_a - ang_b) > 0.035)return ang_a > ang_b; // >2° 角度决定
            // 同向多车道: 用退出组横向位置
            double xl_a = exit_lane_lat_in_group.count(ca->exit_lane_id) ?
                          exit_lane_lat_in_group.at(ca->exit_lane_id) : 0;
            double xl_b = exit_lane_lat_in_group.count(cb->exit_lane_id) ?
                          exit_lane_lat_in_group.at(cb->exit_lane_id) : 0;
            if (std::abs(xl_a - xl_b) > 0.1)return xl_a > xl_b;
            // 进入车道横向并列打破
            double el_a = entry_lane_lat_in_group.count(ca->entry_lane_id) ?
                          entry_lane_lat_in_group.at(ca->entry_lane_id) : 0;
            double el_b = entry_lane_lat_in_group.count(cb->entry_lane_id) ?
                          entry_lane_lat_in_group.at(cb->entry_lane_id) : 0;
            return el_a > el_b;
        });
    }

    // 步骤5: 退出簇排序。
    //   主键:  entry_pt 在退出arm左法向上投影降序
    //   并列:  退出车道在退出组内横向位置降序
    for (auto &kv:exit_group_order_) {
        auto &gid = kv.first;
        auto &cids = kv.second;
        if (cids.size() < 2)continue;
        Vec2d ealn = armLeftNormal(cids, norm_conns, lanes, false);
        std::stable_sort(cids.begin(), cids.end(), [&](const ConnId &a, const ConnId &b) {
            auto *ca = findConn(norm_conns, a);
            auto *cb = findConn(norm_conns, b);
            if (!ca || !cb)return false;
            double la = exitEntryLat(*ca, ealn, lanes);
            double lb = exitEntryLat(*cb, ealn, lanes);
            if (std::abs(la - lb) > 0.5)return la > lb;
            double xl_a = exit_lane_lat_in_group.count(ca->exit_lane_id) ?
                          exit_lane_lat_in_group.at(ca->exit_lane_id) : 0;
            double xl_b = exit_lane_lat_in_group.count(cb->exit_lane_id) ?
                          exit_lane_lat_in_group.at(cb->exit_lane_id) : 0;
            return xl_a > xl_b;
        });
    }

    // 步骤6: 存储rank与横向位置，供拓扑倒置检测使用。
    // 同时通过 t0·t1 点积预计算U型调头标志
    // (不依赖enterGroupId==exitGroupId,因为物理U型调头可能不同组ID)。
    for (auto &kv:entry_group_order_) {
        auto &gid = kv.first;
        auto &cids = kv.second;
        Vec2d me = meanPt(cids, norm_conns, lanes, true);
        Vec2d perp = armLeftNormal(cids, norm_conns, lanes, true);
        for (int i = 0; i < (int) cids.size(); ++i) {
            entry_cluster_rank_[cids[i]] = i;
            auto *c = findConn(norm_conns, cids[i]);
            if (!c)continue;
            const Lane *el = findLane(lanes, c->entry_lane_id);
            const Lane *xl = findLane(lanes, c->exit_lane_id);
            // 退出横向在进入arm坐标系下(值越大=越左)
            double elat = (exitEndpoint(xl) - me).dot(perp);
            exit_lat_in_entry_ref_[cids[i]] = elat;
            exit_lat_sign_[cids[i]] = elat; // 原值用于符号比较

            // 进入端点在自身进入组坐标系下的横向位置
            const Lane *entl = findLane(lanes, c->entry_lane_id);
            entry_lat_in_entry_ref_[cids[i]] = (entryEndpoint(entl) - me).dot(perp);
            // 用实际车道切向点积判断U型调头
            Vec2d t0 = connectionDirection(*c, lanes, true);
            Vec2d t1 = connectionDirection(*c, lanes, false);
            is_uturn_[cids[i]] = (t0.dot(t1) < -0.5);
        }
    }

    for (auto &kv:exit_group_order_) {
        auto &gid = kv.first;
        auto &cids = kv.second;
        Vec2d ealn = armLeftNormal(cids, norm_conns, lanes, false);
        // 退出组自身参考系下的退出横向位置:
        // entry_lat_in_exit_ref_ 与排序键相同 → rank_diff×lat_diff 恒<0 → 永不触发
        // exit_lat_in_exit_group_ref_ 用退出端点在退出arm坐标系下的投影,
        // 与排序键独立, 能揭示真正的几何倒置
        Vec2d exit_ref_pt = meanPt(cids, norm_conns, lanes, false);
        for (int i = 0; i < (int) cids.size(); ++i) {
            exit_cluster_rank_[cids[i]] = i;
            auto *c = findConn(norm_conns, cids[i]);
            if (!c)continue;
            const Lane *el = findLane(lanes, c->entry_lane_id);
            const Lane *xl = findLane(lanes, c->exit_lane_id);
            // 进入横向在退出arm坐标系下(兼容性保留)
            entry_lat_in_exit_ref_[cids[i]] = entryEndpoint(el).dot(ealn);
            // 退出横向在退出组自身坐标系下: (exit_endpoint - mean_exit) · exit_arm_left_normal
            exit_lat_in_exit_group_ref_[cids[i]] = (exitEndpoint(xl) - exit_ref_pt).dot(ealn);
            // 同时为未出现在任何进入组的conn填充is_uturn_
            if (!is_uturn_.count(cids[i])) {
                Vec2d t0 = connectionDirection(*c, lanes, true);
                Vec2d t1 = connectionDirection(*c, lanes, false);
                is_uturn_[cids[i]] = (t0.dot(t1) < -0.5);
            }
        }
    }

    // 步骤7: 从进入簇构造约束配对。
    for (auto &kv:entry_group_order_) {
        auto &gid = kv.first;
        auto &cids = kv.second;
        if (cids.size() > 1)addPairsFromSortedCluster(cids, norm_conns, lanes);
    }

    // 步骤8: 从退出簇构造约束配对(仅新增)。
    for (auto &kv:exit_group_order_) {
        auto &gid = kv.first;
        auto &cids = kv.second;
        if (cids.size() > 1) {
            addPairsFromSortedCluster(cids, norm_conns, lanes);
        }
    }

    // 步骤9: 检测并标记真正的拓扑倒置为StructuralCross。
    detectTopologicalInversions(norm_conns, lanes, crosswalks);

    // 步骤10: 构建哈希索引，使热循环中配对查找为O(1)。
    buildPairIndex();
}

// ─── 查询辅助 ────────────────────────────────────────────────
bool ClusterOrderSolver::pairExists(const ConnId& a, const ConnId& b) const {
    return findPair(a, b) != nullptr;
}

CrossExemption ClusterOrderSolver::exemptionOf(const ConnId& a, const ConnId& b) const {
    const CurvePair* p = findPair(a, b);
    return p ? p->exempt : CrossExemption::None;
}

int ClusterOrderSolver::expectedSideOf(const ConnId &a, const ConnId &b) const {
    const CurvePair* p = findPair(a, b);
    if (!p) return 0;
    // expected_side 以 id_a 视角存储; 查询顺序反向时取反
    return (p->id_a == a) ? p->expected_side : -p->expected_side;
}

Vec2d ClusterOrderSolver::refPerpOf(const ConnId &a, const ConnId &b) const {
    const CurvePair* p = findPair(a, b);
    return p ? p->ref_perp : Vec2d(0, 0);
}

/// 障碍物豁免标记: 当点pt的SDF值<半径r时,设置ObstacleCross
void ClusterOrderSolver::markObstacleExempt(CurvePair &pair, const Vec2d &pt,
                                            const SDFField &sdf, double r) {
    auto kv = sdf.queryWithGrad(pt);
    auto d = kv.first;
    auto _g = kv.second;
    (void)_g;
    if (d < r) {
        pair.exempt = CrossExemption::ObstacleCross;
        pair.exempt_zone_radius = r;
    }
}

/// 第二轮A2标记: 对实际相交且非豁免的配对,
/// 检查相交点是否在障碍物附近,若是则标记为ObstacleCross(弱化约束)
void ClusterOrderSolver::checkAndMarkA2(
    const std::unordered_map<ConnId, BezierCurve>& curves, const SDFField& sdf, double r) {
    if (!sdf.valid())
        return;
    for (auto& pair : pairs_) {
        if (pair.exempt != CrossExemption::None)
            continue;
        auto ia = curves.find(pair.id_a), ib = curves.find(pair.id_b);
        if (ia == curves.end() || ib == curves.end())
            continue;
        if (!curvesIntersectBusiness(ia->second, ib->second, 0.15))
            continue;
        for (auto& pt : curveCrossings(ia->second, ib->second, 0.3))
            markObstacleExempt(pair, pt, sdf, r);
    }
}

}
