#pragma once

#include "types.h"

#include <cstddef>
#include <string>
#include <unordered_set>
#include <vector>

namespace isg {

class ClusterOrderSolver;

struct UTurnFamilyRank {
    std::size_t family_size;
    std::size_t reverse_radius_rank;

    UTurnFamilyRank() : family_size(0), reverse_radius_rank(0) {}
};

struct UTurnLeadFloors {
    double lead0;
    double lead1;
    bool crosswalk0;
    bool crosswalk1;
    std::unordered_set<std::string> crosswalk_ids;

    UTurnLeadFloors()
        : lead0(0.0), lead1(0.0),
          crosswalk0(false), crosswalk1(false) {}
};

/// U-turn 家族的横向阶梯：家族统一步长 + 名义分档量 + 本成员两侧最终站位。
///
/// 入口侧与出口侧必须分别裁定。家族是一条交替由"共享入口端点"和"共享出口
/// 端点"连接起来的链（例：34—36 共享入口，36—35 共享出口，35—37 共享入口），
/// 而横向内缩在入口侧沿 +lateral、在出口侧沿 -lateral 施加。用同一个标量同时
/// 表达两侧，等于要求半径次序在两条不同的链上同时单调：把某成员压小以让开它
/// 共享入口的退化邻居，就会让它相对共享出口的邻居反号。实测 100000598 的
/// 58|57、60|59、64|63、15|14、17|16 与 100000699、100000643 的同型对全部由此
/// 产生（单标量口径下这些数据集违约共 +9）。
struct UTurnFamilyLadder {
    double step;              ///< 家族统一分档步长（米），相邻名次的名义站位差
    double nominal_stagger;   ///< 家族裁定前的名义分档量 = step * 逆序名次
    double stagger;           ///< 两侧站位的较大者，供只需单一标量的旧口径使用
    double entry_stagger;     ///< 入口侧最终横向站位（米，非负）
    double exit_stagger;      ///< 出口侧最终横向站位（米，非负）

    UTurnFamilyLadder()
        : step(0.0), nominal_stagger(0.0), stagger(0.0),
          entry_stagger(0.0), exit_stagger(0.0) {}
};

/// 单条 U-turn 使用的家族级只读几何快照。
struct UTurnFamilyInfo {
    bool geometric_uturn;
    bool boundary_alignment_required;
    Vec2d axis;
    double radius;
    double lead0;
    double lead1;
    double aligned_station;
    UTurnFamilyRank rank;
    std::vector<Crosswalk> clearance_crosswalks;

    UTurnFamilyInfo()
        : geometric_uturn(false), boundary_alignment_required(false),
          axis(1, 0), radius(0.0),
          lead0(0.0), lead1(0.0),
          aligned_station(0.0) {}
};

/// U-turn 家族共享的只读几何元数据。
class UTurnFamilyBuilder {
public:
    UTurnFamilyInfo build(
        const Connectivity& connectivity,
        const IntersectionInput& input,
        UTurnAlignmentScope scope,
        const ClusterOrderSolver* topology = nullptr,
        bool require_shared_endpoint_pair_for_rank = false) const;

    bool isGeometricUTurn(const Connectivity& connectivity, const IntersectionInput& input) const;

    double radiusKey(const Connectivity& connectivity, const IntersectionInput& input) const;

    Vec2d alignmentAxis(const Vec2d& entry_tangent, const Vec2d& exit_tangent) const;

    double requiredAlignedStation(const Vec2d& entry_point,
                                  const Vec2d& entry_tangent,
                                  const Vec2d& exit_point,
                                  const Vec2d& exit_tangent,
                                  double min_lead0,
                                  double min_lead1) const;

    double requiredAlignedStationOnAxis(const Vec2d& entry_point,
                                        const Vec2d& entry_tangent,
                                        const Vec2d& exit_point,
                                        const Vec2d& exit_tangent,
                                        const Vec2d& axis,
                                        double min_lead0,
                                        double min_lead1) const;

    std::vector<const Connectivity*> alignmentComponent(
        const Connectivity& connectivity,
        const IntersectionInput& input,
        UTurnAlignmentScope scope,
        const ClusterOrderSolver* topology = nullptr,
        bool require_shared_endpoint_pair = false) const;

    UTurnFamilyRank radiusRank(
        const Connectivity& connectivity,
        const IntersectionInput& input,
        UTurnAlignmentScope scope,
        const ClusterOrderSolver* topology = nullptr,
        bool require_shared_endpoint_pair = false) const;

    /// 家族统一横向阶梯：由家族全体成员的走廊宽度与共享端点拓扑共同裁定
    /// 本成员的横向分档站位，而不是让每个成员事后各自钳制。
    ///
    /// 名义阶梯为 step * reverse_radius_rank，名次越高（半径越小）越靠里。
    /// buildSegmented 会再按各成员自身走廊把分档量压到 0.25 * 走廊宽度以内；
    /// 一旦某个退化成员被压掉分档、而与它共享端点的外层邻居没有被压掉，阶梯
    /// 就会反号：外层成员的入口直段反而比内层成员更靠内，横扫过内层的窄走廊，
    /// 且是否相交还取决于两条曲线谁先生成，逐条做的同簇审计抓不住
    /// （100000412 的 34|36、16|18、8|10 全属此类，退化成员走廊只有 0.32~0.40m）。
    ///
    /// 这里把裁定提前到家族层面，分两步：
    ///   1. 自身走廊上限：cap = min(step * rank, 0.25 * radiusKey)；
    ///   2. 共享端点约束：若某成员与更内层成员共享真实端点坐标，则它在
    ///      "共享的那一侧"必须比该内层成员第 1 步的结果再小一个分离间距
    ///      sep = min(step, 0.5 * 内层 cap)。
    ///
    /// 第 2 步按侧独立裁定，返回 entry_stagger / exit_stagger 两个值。家族是一条
    /// 交替由共享入口、共享出口连接的链，内缩在两侧沿相反法向施加，单标量无法
    /// 同时让两条链单调——详见 UTurnFamilyLadder 的说明与实测清单。
    ///
    /// 第 2 步只看直接共享端点的邻居、不做链式传递：链式传递会把退化成员的
    /// 厘米级上限一路压到家族另一端，实测把 100000412 的 35|37 由不相交压成
    /// 相交（35 需要约 0.25m 的绝对站位才能避开退化为单段的 37）。
    /// 名义阶梯本身相邻名次恰好相差一个 step，因此在没有退化成员的家族里
    /// 第 2 步完全不生效，其他数据集的形态保持不变。
    ///
    /// 注：曾按需求原文尝试"家族全体同量同向平移"（把最紧成员的超出量整体
    /// 平移掉），实测不可行：100000412 家族 A 需要 0.915m 平移，而入口车道
    /// 两侧并没有这么多横向余量，靠外成员的直段被推到相邻直行/转向曲线上，
    /// 该数据集违约由 5 升到 20。因此改为在家族层面裁定阶梯：同样保证左右
    /// 次序与生成顺序无关，但不动用车道外侧不存在的横向空间。
    UTurnFamilyLadder familyLateralLadder(
        const Connectivity& connectivity,
        const IntersectionInput& input,
        UTurnAlignmentScope scope,
        double nominal_step,
        const ClusterOrderSolver* topology = nullptr,
        bool require_shared_endpoint_pair = false) const;

    UTurnLeadFloors leadFloors(const Vec2d& entry_point,
                               const Vec2d& entry_tangent,
                               const Vec2d& exit_point,
                               const Vec2d& exit_tangent,
                               const IntersectionInput& input,
                               double no_crosswalk_min_lead = 2.0) const;

    std::vector<Crosswalk> clearanceCrosswalks(
            const IntersectionInput& input, const UTurnLeadFloors& floors) const;

    double alignmentComponentStation(
        const Connectivity& connectivity,
        const IntersectionInput& input,
        UTurnAlignmentScope scope,
        const ClusterOrderSolver* topology = nullptr) const;

    void enforceAlignmentStation(const Vec2d& entry_point,
                                 const Vec2d& entry_tangent,
                                 const Vec2d& exit_point,
                                 const Vec2d& exit_tangent,
                                 double target_station,
                                 double& min_lead0,
                                 double& min_lead1) const;
};

}  // 命名空间 isg
