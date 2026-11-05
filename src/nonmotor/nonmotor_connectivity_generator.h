#pragma once

#include "types.h"
#include "nonmotor_backbone_extractor.h"
#include "nonmotor_type.h"
#include <vector>

namespace isg {

struct NonmotorConifg {
    double road_edge_clearance      = 0.9;  ///< 一般非机动车道到 RoadEdge 最小净距（米）
    double motor_to_nm_clearance    = 1.0;  ///< 进入机动→退出非机动时的净距（米）
    double corner_radius            = 1.5;  ///< L型/U型折弯处圆角半径（米）
    double overlap_offset           = 0.1;  ///< 曲线重叠时的微偏移量（Q5，可配置）
    double island_max_dist          = 20.0; ///< 安全岛距端点的最大穿行距离（Q4）
    double l_turn_knee_max_offset   = 1.0;  ///< L型左转折弯膝点距两者中位的最大偏差（Q2）
};

/// 非机动车道挂接形态主生成器
/// 1.单条挂接路径生成
/// (1).按进入车道/退出车道截取有效范围骨干路径序列段
/// (2).对截取的骨干路径段做偏移处理：根据每段的引用次数N做相应N*0.1m的偏移处理
/// (3).生成挂接形态段序列：进入段 + 截取有效骨干段(逆时针序列) + 退出段
/// (4).平滑挂接形态段：连接段间存在大转弯角度时需做转角平滑处理
/// (5).构建完整挂接结构：ConnectivityCurve
/// 2.路口所有非机动挂接路径生成：generate
/// 依次遍历所有非机动挂接，并调用单条挂接路径生成，返回所有挂接路径集
class NonMotorConnectivityGenerator {
public:
    explicit NonMotorConnectivityGenerator(const NonmotorConifg& cfg = {});

    /// 生成所有按非机动逻辑处理的连通曲线
    /// @param input         经 ConnectivityMotorLogicFilter(false) 过滤的输入
    /// @param motor_curves  已生成的机动车道曲线（用于微偏移消重叠）
    std::vector<ConnectivityCurve> generate(
        const IntersectionInput& input,
        const std::vector<ConnectivityCurve>& motor_curves,
        NonMotorBackbone& backbone_output);

private:
    /// 生成单条完整非机动挂接线形态，并经过偏移、平滑处理
    LineString2d generateOne(
        const NonMotorBackbone& backbone,
        const Vec2d& entry_pt, const Vec2d& entry_dir,
        const Vec2d& exit_pt, const Vec2d& exit_dir,
        int path_index, double corner_radius = 1.5);

private:
    NonmotorConifg config_;

};

}  // namespace isg
