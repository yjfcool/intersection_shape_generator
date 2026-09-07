//
// RoadEdge 净距判定的统一口径：区分"曲线自己贴过去"与"连接点位置+切向已经锁定的亏欠"。
//
// 连接点位置由输入给定，端点切向又被 G1 连续性锁定在车道朝向上。若连接点自身就落在净距
// 要求之内（corpus 里 `110000741-u` 出口连接点距 RoadEdge 仅 0.9521m，而 mode=2 要求
// 1.0m），那么端点邻域的亏欠不是曲线选错了形状，而是输入强制的：任何满足 G1 的曲线离开
// 该连接点时都只能沿车道切向走，净距按 `d(s) ≈ d(0) + s·sinθ` 缓慢恢复。
//
// 无豁免时这笔还不掉的账会让生成侧把每一个候选都判成"物理风险"，从而放弃合法的单段三次
// 曲线、升级成两段绕行，再被形态约束判成非法拱形——一个不可满足的硬约束把三条本来只差
// 几毫米的曲线全部推成了真正的违规。
//
#ifndef ISG_ROAD_EDGE_CLEARANCE_H
#define ISG_ROAD_EDGE_CLEARANCE_H

#pragma once

#include "types.h"
#include <algorithm>
#include <limits>

namespace isg {

/// 净距判定的舍入容差（米）。与围栏连接点舍入同量级，吸收采样落位与折线离散误差。
static constexpr double kRoadEdgeClearanceRoundingTol = 0.0025;

/// 端点强制影响区的长度（米）：只有到较近连接点距离不超过它的采样点才按"沿切向直行"的
/// 可达下限判定，更远处一律按全额净距要求。
///
/// 这个上限不可省略。距离函数沿一条与路缘平行的射线几乎不增长，若不限定范围，任何与
/// RoadEdge 长距离平行且净距略欠的连通都会被整条豁免，硬约束等于作废。定标见文件末尾
/// 说明：corpus 里 `110000741-u` 出口的亏欠区延伸到离端点 0.58m，故取 1.0m——与切向合成
/// 点的插入距离同量级，仍属"端点邻域"，且区外只要求横向让出几厘米，任何曲线都做得到。
static constexpr double kRoadEdgeClearanceEndpointSpan = 1.0;

/// 曲线对 RoadEdge 的净距测量结果。距离均已排除首尾连接点自身的浮点邻域。
struct RoadEdgeClearanceMeasure {
    /// 全曲线最小净距。
    double minimum = std::numeric_limits<double>::infinity();
    /// 可归责于曲线的最大亏欠：逐采样点按"该点可达净距下限"判定后取最大值。
    /// 可达下限 = min(clearance, 沿端点切向直行到同等距离处的净距)，见文件头说明。
    double deficit = 0.0;
    /// 首/尾连接点自身到 RoadEdge 的距离。
    double start_distance = std::numeric_limits<double>::infinity();
    double end_distance = std::numeric_limits<double>::infinity();
    /// minimum 与 deficit 的取值位置，用于报告。
    Vec2d location{0, 0};
    Vec2d deficit_location{0, 0};
    /// 是否至少有一个有效采样点参与测量。
    bool valid = false;
};

/// 可归责于曲线的净距亏欠。测量阶段已按统一口径算好，此处只做取值，保证生成侧与
/// 审计侧共用同一判据。
inline double roadEdgeClearanceDeficit(
        const RoadEdgeClearanceMeasure& measure, double clearance) {
    if (clearance <= 0.0 || !measure.valid)
        return 0.0;
    return measure.deficit;
}

/// 是否属于"净距不足但完全由连接点位置与切向强制"的情形：此时曲线已做到端点允许的
/// 最好结果，审计应记为豁免而非违规，生成侧也不应据此判定物理风险并升级为绕行。
inline bool roadEdgeClearanceForcedByEndpoint(
        const RoadEdgeClearanceMeasure& measure, double clearance) {
    if (clearance <= 0.0 || !measure.valid)
        return false;
    return clearance - measure.minimum > kRoadEdgeClearanceRoundingTol &&
           measure.deficit <= kRoadEdgeClearanceRoundingTol;
}

/// 逐采样点的净距判定下限。`sample_distance_to_endpoint` 用采样点到较近端点的直线距离
/// （不超过弧长，方向偏保守），`ray_distance` 为沿该端点切向直行同等距离处的净距
/// （<0 表示调用方未计算，按全额判定）。
///
/// 判据只有一条：零曲率方案总是可达的。端点位置与切向都被输入和 G1 锁定，曲线离开连接点
/// 时只能先沿切向走，所以"沿切向直行到 s 处"的净距是这一点上任何合法曲线都逃不掉的量级；
/// 低于它才是曲线自身形态造成的。曲线额外弯离路缘最多再挣得 κs²/2（κ≤0.25、s≤0.5m 时约
/// 3cm），这部分刻意不追究——与围栏连接点舍入豁免同一取舍，理由是此处曲率预算已被转向
/// 半径本身占用，把理论弯曲余量算作"本该做到"会让唯一合法的单段三次曲线被判违规。
///
/// 注意不要反用 1-Lipschitz 界：`d(端点)+s ≥ clearance` 只说明全额净距不违反距离函数的
/// 增长上限，并不表示它可达——那需要曲线垂直于路缘离开，而切向通常近乎平行于路缘，实际
/// 只能按 `s·sinθ` 恢复。早先按该式提前返回全额 clearance，等于让判据比本注释描述的更严，
/// corpus 里 `110000741-u` 的合法单段曲线因此被判 0.0427m 亏欠并触发绕行升级。
inline double roadEdgeClearanceFloor(double clearance,
                                     double sample_distance_to_endpoint,
                                     double ray_distance) {
    if (sample_distance_to_endpoint > kRoadEdgeClearanceEndpointSpan ||
        ray_distance < 0.0)
        return clearance;
    return std::min(clearance, ray_distance);
}

/// 采样点是否落在端点强制影响区内：只有此时才需要额外计算切向直行参考点。
inline bool roadEdgeClearanceNeedsFloor(double sample_distance_to_endpoint) {
    return sample_distance_to_endpoint <= kRoadEdgeClearanceEndpointSpan;
}

}  // 命名空间 isg

#endif  // ISG_ROAD_EDGE_CLEARANCE_H
