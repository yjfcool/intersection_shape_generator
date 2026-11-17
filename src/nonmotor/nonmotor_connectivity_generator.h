#pragma once

#include "types.h"
#include "nonmotor_backbone_extractor.h"
#include "nonmotor_type.h"
#include <vector>

namespace isg {

/// 非机动挂接生成配置。所有距离单位为米；进出骨干交点和内部点一起偏移，仅车道端点固定。
struct NonmotorConifg {
    double road_edge_clearance      = 0.9;  ///< 一般非机动车道到 RoadEdge 最小净距（米）
    double boundary_avoidance_clearance = 1.0; ///< 非机动绕行所有Boundary类型的最小净距（米）；RoadEdge取其与road_edge_clearance的较大值
    double motor_to_nm_clearance    = 1.0;  ///< 进入机动→退出非机动时的净距（米）
    double corner_radius            = 2.0;  ///< 折角两侧目标截距（米）；短段及相邻圆角相遇时缩短，0禁用圆角
    double overlap_offset           = 0.3;  ///< 共享骨干路径相邻横向轨位间距；0禁用分离但保留重叠验收
    double short_same_segment_max_length = 2.0;  ///< 非掉头同骨干同分段挂接的最大有效长度；范围内合并为一个公共交点（米）
    double source_surface_clearance = 0.3; ///< 来源面内中部目标边界净距；原骨干空间更窄时保持其已有净距（米）
    double island_max_dist          = 20.0; ///< 安全岛距端点的最大穿行距离（Q4）
    double l_turn_knee_max_offset   = 1.0;  ///< L型左转折弯膝点距两者中位的最大偏差（Q2）
};

/// 非机动车道挂接形态主生成器。
/// 1. 进入切向、退出反向切向的射线与主环/分支求最近交点，作为有效骨干路径的原始截取点。
/// 2. 普通直行/左转按主环逆时针截取，掉头取较短方向；同分支直接截取有效区间。
///    多点道路Gap的内部形态参与完整截取；分支按节点坐标映射环段，拓扑节点序与几何段序独立。
/// 3. 以实际几何段和截取站距构建占用区间，识别整段、部分及反向共享。
/// 4. 按每个局部共享区间的实际占用者分配居中等间距轨位；内部节点不强制回零。
///    overlap_offset 是相邻轨位的目标间距，进出骨干交点使用对应轨位，不锁定原坐标。
///    普通路径优先左右均衡；掉头优先占进出车道端点所在侧的外缘轨位，独立于骨干行驶方向。
///    横道/安全岛来源面限制偏移可用横向范围；首尾按偏移线与原面重新求交，面内中部不越界。
///    共享轨位整体平移或同比缩小间距适配可用范围，避免逐条截断造成重叠；优先靠近原骨干。
/// 5. 公共点沿前驱回退、沿后继前进corner_radius，按截取点和方向拟合贝塞尔，
///    以0.2米平面弧长采样（尾部余段补终点）；短段或相邻圆角相遇时缩短截距。
///    仅进入/退出车道端点固定，通过切向连续的进出段连接偏移后的骨干；fixed_shape不修改。
///    偏移段以支撑线交点拼接；进出骨干端段与车道射线重新求公共点，先做一次统一圆角，
///    再以受限小半径修正密集采样后仍残留的尖折；对端部短骨干段启用有界挂接窗口，
///    合并车道接入与邻近偏移转角并按解析曲率筛选独立控制柄；转角在骨干侧与对应偏移段相切，
///    中间直段保持轨位，不强制骨干交点采用车道切向。
///    尖折修复窗口按弧长精确截取（每侧最多min(0.8,corner_radius/2)米），窗口外点序和XYZ完整保留；
///    不能从入口或上一修复点直连到下一窗口，截取点使用保留路径切向拟合局部G1过渡。
///    公共点可能被圆角截去；仅固定车道端点。近平行或过长尖角采用有界回退。
///    生产拼接只接受沿路径前方的支撑线交点；真实折点前后的corner_radius保护带
///    不插入占用事件短段，但完整有效骨干范围始终保留，再在其末端做受限切向过渡，避免鼓肚和二次打折。
///    非掉头挂接若两端落在同一骨干同一几何分段且有效长度不超过
///    short_same_segment_max_length，则将两交点合并为两点间骨干中点，作为唯一公共交点。
///    圆角新增空间违规时在相同轨位上最多6次减半截距，保留0.2米采样规则。
/// 6. 统一验收最终几何的有效骨干重叠、完整路径自交及空间约束。
///    所有Boundary类型（含Other）均禁止穿越；进入、退出及骨干段的冲突在最终几何上局部绕行。
///    距边缘不足boundary_avoidance_clearance（默认1米）也触发绕行；外扩、可见边、平滑后及最终验收使用同一净距。
///    有限边界外扩后建立不超过512节点的局部可见图，每窗口最多尝试4个范围，每曲线最多修复12处；
///    平滑后重新验收边界、fence/洞、障碍、来源面、自交及原有空间约束。窗口外形态和车道端点保留。
///    进出段绕行直接接回各自偏移骨干，避免回到共同车道延长线新增重叠；无法避让的穿越为Infeasible。
///    可行性优先；同等约束候选再比较掉头侧别、两侧均衡与偏移幅度，空间受限时允许非对称避让。
///    非机动不验证不同曲线间的同簇非端点相交。
///    无法解决的固定形态或通道重叠冲突保留可审计几何，返回 Infeasible/违规原因；
///    障碍物与空间约束仍分别报告 WarnA2/Degraded，不将偏移执行过视为成功。
class NonMotorConnectivityGenerator {
public:
    explicit NonMotorConnectivityGenerator(const NonmotorConifg& cfg = {});

    /// 生成所有按非机动逻辑处理的连通曲线，输出顺序与 input.connectivities 一致。
    /// @param input         经 ConnectivityMotorLogicFilter(false) 过滤的输入
    /// @param motor_curves  已生成的机动车道曲线（参与最终重叠/净距验收）
    std::vector<ConnectivityCurve> generate(
        const IntersectionInput& input,
        const std::vector<ConnectivityCurve>& motor_curves,
        NonMotorBackbone& backbone_output);

private:
    NonmotorConifg config_;

};

}  // namespace isg
