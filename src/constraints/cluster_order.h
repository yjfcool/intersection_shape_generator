#pragma once

#include "types.h"
#include <unordered_map>

namespace isg {

/// 交叉豁免类型
enum class CrossExemption {
    None,              ///< 无豁免：正常受同簇非交约束
    StructuralCross,   ///< 结构性必须交叉（拓扑倒置、不同arm交叉交通等）
    ObstacleCross      ///< 障碍物强制交叉（穿障避障时允许）
};

/// 同簇曲线对：记录两条曲线间的约束关系
struct CurvePair {
    ConnId id_a, id_b;
    CrossExemption exempt = CrossExemption::None;  ///< 交叉豁免
    double exempt_zone_radius = 0.0;               ///< 豁免区半径
    int expected_side = 0;   ///< 预期方位：+1= id_a在id_b左侧; -1=id_a在id_b右侧; 0=未知
    Vec2d ref_perp{0, 0};    ///< 配对特定的横向参考轴（mean_entry→mean_exit的左法向）
    /// 共享端点：两条曲线共享同一退出车道或同一进入车道时为true，
    /// 仅用于识别合法连接点和端点侧向排序，连接点外仍禁止贴合/重叠/相交
    bool shared_endpoint = false;
};

/// 簇求解器的只读静态输出；运行时 A2 豁免状态不写回该快照。
struct ClusterTopology {
    std::vector<CurvePair> pairs;
    std::unordered_map<LaneGroupId, std::vector<ConnId> > entry_order;
    std::unordered_map<LaneGroupId, std::vector<ConnId> > exit_order;
};

class SDFField;

// ─────────────────────────────────────────────────────────────────────────────
// ClusterOrderSolver —— 拓扑优先的同簇排序求解器
//
// 排序策略:
//   进入簇: 按 (exit_pt - mean_entry_pt) 的角度降序排列 (CCW在前 = 左侧在前)
//           同向并列时按退出车道在退出组内的横向位置降序（稳定的同向多车道排序）
//   退出簇: 按 entry_pt 在退出arm左法向上的投影降序排列 (左侧在前)
//           并列时按退出车道在退出组内的横向位置降序
//
// 结构性交叉检测: 拓扑序 × 横向倒置
//   同进入组对(A,B):
//     entry_cluster_rank(A) < entry_cluster_rank(B) [A在进入簇中更左]
//     但 exit_lat(A) < exit_lat(B) [A退出更靠右]
//     → (rank_A-rank_B)*(exit_lat_A-exit_lat_B) > 0 → StructuralCross
//   同退出组对(A,B):
//     exit_cluster_rank(A) < exit_cluster_rank(B) [A在退出簇中更左]
//     但 entry_lat_A < entry_lat_B [A来自更右]
//     → (rank_A-rank_B)*(entry_lat_A-entry_lat_B) > 0 → StructuralCross
// ─────────────────────────────────────────────────────────────────────────────
class ClusterOrderSolver {
public:
    /// 构造簇排序与配对约束
    void build(const std::vector<Connectivity>&, const std::vector<Lane>&, const std::vector<LaneGroup>&,
               const std::vector<Crosswalk>& crosswalks = {});

    /// 标记障碍物豁免：当点pt位于障碍物SDF<半径r时设置ObstacleCross
    void markObstacleExempt(CurvePair&, const Vec2d&, const SDFField&, double r = 1.5);

    const std::vector<CurvePair>& pairs() const { return pairs_; }
    ClusterTopology topology() const {
        ClusterTopology result;
        result.pairs = pairs_;
        result.entry_order = entry_group_order_;
        result.exit_order = exit_group_order_;
        return result;
    }

    /// 当 (a,b) 共享同一退出/进入车道端点时返回true。
    bool isSharedEndpoint(const ConnId& a, const ConnId& b) const {
        const CurvePair* p = findPair(a, b);
        return p ? p->shared_endpoint : false;
    }

    CrossExemption exemptionOf(const ConnId&, const ConnId&) const;
    bool pairExists(const ConnId& a, const ConnId& b) const;
    int expectedSideOf(const ConnId &, const ConnId &) const;
    Vec2d refPerpOf(const ConnId &, const ConnId &) const;

    /// 第二轮A2标记：对实际相交且非豁免的配对，检查相交点是否在障碍物附近，
    /// 若是则标记为ObstacleCross（弱化约束）
    void checkAndMarkA2(const std::unordered_map<ConnId, BezierCurve>&, const SDFField&, double r = 1.5);

    /// 返回每个进入组按空间排序后的conn_id列表
    const std::unordered_map<LaneGroupId, std::vector<ConnId>>& entryGroupOrder() const {
        return entry_group_order_;
    }

    /// 返回每个退出组按空间排序后的conn_id列表
    const std::unordered_map<LaneGroupId, std::vector<ConnId>>& exitGroupOrder() const {
        return exit_group_order_;
    }

private:
    void addPairsFromSortedCluster(
            const std::vector<ConnId>& cids, const std::vector<Connectivity>& conns, const std::vector<Lane> &);

    static bool hasPair(const std::vector<CurvePair>& pairs, const ConnId& a, const ConnId& b);

    void detectTopologicalInversions(const std::vector<Connectivity> &, const std::vector<Lane> &,
                                     const std::vector<Crosswalk> &crosswalks);

    // ── 性能优化: 哈希索引，key为每个pair中字典序较小的ConnId
    // 在 build() 末尾构建，供 exemptionOf / expectedSideOf / refPerpOf /
    // pairExists / isSharedEndpoint 使用，避免热循环中O(pairs_)线性扫描
    void buildPairIndex();
    const CurvePair* findPair(const ConnId& a, const ConnId& b) const;

private:
    std::vector<CurvePair> pairs_;
    std::unordered_map<ConnId, std::vector<int>> pair_index_;
    std::unordered_map<LaneGroupId, std::vector<ConnId>> entry_group_order_;
    std::unordered_map<LaneGroupId, std::vector<ConnId>> exit_group_order_;

    // 通过 t0·t1 点积判断U型调头（不依赖enterGroupId==exitGroupId，
    // 因为物理U型调头可能存在不同的组ID）
    std::unordered_map<ConnId, bool> is_uturn_;

    // exit_lat_in_entry_ref 符号: +1=退出更左, -1=退出更右
    // 用于检测对侧结构性交叉（如左转与右转必然交叉）
    std::unordered_map<ConnId, double> exit_lat_sign_;

    // 进入车道端点在进入组坐标系下的横向位置
    // 用于进入-退出横向倒置检测:
    //   (entry_lat_A - entry_lat_B) × (exit_lat_A - exit_lat_B) < 0 → StructuralCross
    // 捕捉同侧退出但车道顺序倒置的情况
    std::unordered_map<ConnId, double> entry_lat_in_entry_ref_;

    // 拓扑倒置检测使用的横向位置
    // entry_cluster: exit_lat_in_entry_ref[cid] = exit_pt 投影到进入arm法向
    // exit_cluster:  entry_lat_in_exit_ref[cid]  = entry_pt 投影到退出arm左法向
    std::unordered_map<ConnId, double> exit_lat_in_entry_ref_;
    std::unordered_map<ConnId, double> entry_lat_in_exit_ref_;

    // 退出车道在自身退出组坐标系下的横向位置
    // = (exitEndpoint - mean_exit_pt) · exit_arm_left_normal
    // 与 entry_lat_in_exit_ref_ 不同（后者始终与sort key一致，
    // 使倒置检测对same_exit对永远失效）
    std::unordered_map<ConnId, double> exit_lat_in_exit_group_ref_;
    /// 簇内排序rank（索引越小=越左）
    std::unordered_map<ConnId, int> entry_cluster_rank_;
    std::unordered_map<ConnId, int> exit_cluster_rank_;
};

}
