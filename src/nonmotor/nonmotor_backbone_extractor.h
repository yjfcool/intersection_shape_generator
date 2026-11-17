#pragma once

#include "types.h"
#include "nonmotor_type.h"
#include "utils.h"
#include <string>
#include <vector>
#include <fstream>

namespace isg {

/// 非机动骨干路径提取器。
/// 从横道面、安全岛面和附近的直行车道推导骨干段，随后将骨干段组织成有序的图结构。
///  1.提取人行横道骨干线：搜索横道中心指定范围kLaneSearchRadius内的进入组/退出组直行车道，并沿其直行方向向横道发出射线，
///    先求有效直行方向的组平均垂向，再选择垂向投影最远的两条有效射线，以其相交弦的中点连线确定横道骨干方向；
///    两射线分别与人行横道相交于Ra1、Ra2和Rb1、Rb2，计算射线中点分别为Ra0=(Ra1+Ra2)/2与Rb0=(Rb1+Rb2)/2，
///    两个射线中点Ra0、Rb0连线并向两端延长与人行横道相交Cws、Cwe形成的有效线段CwSpine即是人行横道骨干线。
///    有效射线不足时回退到附近直行车道方向的垂向；附近没有有效直行车道时才使用横道外环最长边方向。
///    首尾优先按有效直行组方向统一；无有效直行基准时，使用附近可横穿横道、与骨干近似垂直（偏差不超过30度）的车道组，
///    不限制其转向类型，优先进入组并在同角色内选择最近端点，仅交换首尾而不改变骨干形状。
///    以进入组方向为基准，首点位于尾点左侧；以退出组方向为基准，首点位于尾点右侧；无有效组时保留几何首尾。
///  2.提取安全岛骨干线：
///    以安全岛面外扩指定范围kSearchRadius搜索与其相交的人行横道面，
///    多个候选横道面两两检查 1 米外扩相交关系；若相交则过滤其中距离安全岛面较远的横道面，避免误检横道参与后续骨干处理。
///    (1)搜到相交人行横道面：
///    仅一个邻接横道时，若其骨干线从端点向外延长与安全岛相交，连接该端点与延长方向的最远外环交点；否则沿用中心点连接方式。
///    恰有两个邻接横道且无法取得有效交点中心时，若两骨干线无向夹角小于30度且最近端点连线接触安全岛，则直接使用该连线作为安全岛骨干线。
///    过滤后仍有至少三个邻接横道时，若存在无向夹角小于10度的横道对，则取端点连线最短的一对；其中点必须在安全岛内且距边界至少0.5米，否则改用安全岛中心在该端点连线上的最近投影点，投影不满足安全条件时沿用原中心逻辑。
///    a.获取安全岛中心点：三条及以上邻接横道存在小于10度的骨干平行对时，取端点连线最短的一对，并优先以其中点为中心；若中点在岛外或距边界不足0.5米，则尝试将安全岛面中心投影到该线段，只有投影点位于安全岛内且距边界至少0.5米时才采用；否则计算所有邻接人行横道骨干线到路口面中心的距离，选取最近两条骨干线，将其双向延长交点作为中心；路口面缺失、交点退化或交点落在安全岛面外时回退到安全岛面中心；
///    b.分别连接安全岛中心点与每个相邻人行横道面骨干线的最近侧端点连线作为安全岛骨干线；
///    (2)未搜到相交人行横道面：在安全岛附近指定范围kLaneSearchRadius内搜索进入/退出车道组(可能有多个组)，
///    获取所有搜到车道组方向(以直行车道的方向为准)并将近似方向(方向差kDirectionTolerance不超过30度)进行聚类获得一个或多个方向组，并计算每组内所有方向的均值，
///    如果仅一个方向组时，以方向组均值方向的垂向且过安全岛面中心点做双向射线，射线与安全岛的交点构成的连线作为安全岛骨干线，
///    如果有多个方向组时，分别以每个方向组均值方向的垂向且过安全岛中心点做双向射线，射线与安全岛相交的最佳交点(
///      获取此方向组内与安全岛面中心点最近的车道组，与最近车道组较近的那个交点)与安全岛面中心点连线作为安全岛骨干线。
///    安全岛骨干线按对应方向组角色执行相同首尾规则；由横道端点或中心连线生成的骨干线按来源横道骨干线的原始向量方向排列。
///  3.生成非机动骨干路径图：
///    (1)空隙连接线补全：普通Gap在NonMotorGapOptions.search_radius（默认20米）范围内搜索；
///       有RoadEdge有效范围时使用独立的road_edge_search_radius（默认40米），端点分别靠近同一
///       RoadEdge网络即可进入长边缘候选，避免把两骨干端点直线距离误当作道路路径长度。
///       两骨干所属对象面各外扩surface_expansion（默认1米）并裁剪RoadEdge，在端点包围走廊内
///       沿真实相交/共端点构建的道路图以有界最短路截取连续有效范围；道路断缝或高程不匹配的交叉不得视为连通。
///       道路有效范围向路口中心一侧法向偏移inward_offset（默认0.9米），首尾支撑线延长与原骨干求公共端点；
///       普通道路候选端点调整不超过endpoint_adjustment（默认3米）；超出普通Gap范围的长RoadEdge
///       允许额外计入两份inward_offset作为偏移支撑线延长余量。两者均禁止反转骨干、过长尖角、自交、穿越Boundary/障碍及fence越界。
///       最近候选不可行时不抢占更远骨干的近邻关系，回退原逻辑：在3米范围内搜索其他骨干线端点；
///       若同一候选骨干线的两端均在范围内，选择距离孤立端点最近的一端。两端间距严格小于2米时，
///       优先沿两个端点的骨干方向双向延长求交，允许反向交点，并直接替换原有两个端点，避免生成过短Gap；
///       其他情况仅在两骨干线无向夹角大于45度且交点落在任一来源对象面内时合并端点，
///       否则检查孤立骨干线另一端是否能在孤立点3米搜索范围内，通过共享几何端点和一条或多条骨干线
///       连通到同一最近候选端点；若能则不生成局部Gap，保留孤立端点，
///       若不能则保留两端并连成Gap，连接方向按两端所属骨干线的首尾方向确定。
///    (2)有序路径及环平滑：保留骨干线端点逆时针排序形成的大环并补齐相邻Gap；孤立端点所属骨干不参与主环补Gap；
///       若相邻端点已由安全岛或其他骨干线连通，不重复生成大环Gap，此时输出骨干图不标记为连续大环；若过滤后不足两条主环骨干，
///       则回退全量排序以保留简单输入的大环闭合。3米搜索生成的局部连接另行加入；
///       若局部连接与大环间隔使用同一端点对，则复用该Gap以避免重复边；
///       按主环顺序和外围分枝顺序分段存储。骨干图中仅有两条骨干线共用且夹角小于限定角度的 terminal 端点可删除并合并；
///       合并后的骨干继续递归检查另一端的未处理骨干；两条骨干的远端 terminal 状态不影响该公共端点删除。
///       非terminal公共端点、局部Gap接入点、三线以上共点及夹角不小于限定角度的共点不调整。
///    (3)连通分量补全：对原规则留下的分离分量，按端点距离从近到远尝试有界道路Gap，失败时尝试原直线Gap。
///       每次仅合并不同分量，不重复加闭环边；接入目标骨干内部时拆边建立公共节点，完整保留主环两端及有效范围。
///       新连接及调整后的骨干需通过空间和拓扑相交检查；无安全候选时通过connected_components/diagnostics报告，
///       不以穿障、跨道路断缝或超出搜索范围的方式伪造连通。search_radius<=0禁用道路优先和跨分量补全。
///    输出主环时从最终骨干图剥离度为1的分支节点，仅将2-core中的闭合简单环写入ring_geometry；
///    ring_nodes只记录拓扑节点，ring_geometry按环遍历方向拼接每条边的完整折线，保留道路Gap所有内部形态点；
///    环外围的安全岛/横道及其他未进入主环的骨干边只写入分支图层，ring 与 branch 严格互斥。


    /// 安全岛与横道面的最大邻接距离，单位为米。
    const double kSearchRadius = 5.0; //距离过大会错误搜到转角处非真实邻接人行横道
    /// 邻接横道面两两误检过滤使用的外扩距离，单位为米。
    const double kCrosswalkMergeDistance = 1.0;
    /// 三条及以上邻接横道时识别近似平行骨干对的最大夹角，单位为弧度。
    const double kMultiCrosswalkParallelTolerance = 10.0 * DEG2RAD;
    /// 多横道对中点或投影点到安全岛边界的最小允许距离，单位为米。
    const double kParallelCenterBoundaryClearance = 0.5;
    /// 孤立骨干端点搜索其他骨干线端点的最大距离，单位为米。
    const double kIsolatedEndpointSearchRadius = 3.0;
    /// 短 Gap 直接改用两条骨干端点双向延长线交点的严格长度上限，单位为米。
    const double kShortGapLengthThreshold = 2.0;
    /// 孤立端点与近侧骨干线方向夹角超过此值时尝试延长求交，单位为弧度。
    const double kGapExtensionAngleThreshold = 45.0 * DEG2RAD;
    /// 骨干路径共端点简化的最大无向夹角，单位为弧度。
    const double kSpineSmoothAngleThreshold = 10.0 * DEG2RAD;
    /// 根据车道组方向推断骨干线时，允许采纳的车道端点搜索半径，单位为米。
    const double kLaneSearchRadius = 10.0;
    /// 安全岛车道方向无向聚类、横道定向回退时的垂向偏差限值，单位为弧度。
    const double kDirectionTolerance = 30.0 * DEG2RAD;
    /// 端点合并和零长度线判断使用的平面距离容差，单位为米。
    const double kPointTolerance = 1e-6;

/// 道路边缘优先补连的局部搜索配置；距离均为米，原3米Gap规则保持独立。
struct NonMotorGapOptions {
    double search_radius = 20.0;       ///< 孤立端点及最终跨分量补连的最大端点间距；<=0禁用新策略。
    double road_edge_search_radius = 40.0; ///< 道路边缘优先补连的端点搜索范围；道路连接允许超出普通Gap范围。
    double surface_expansion = 1.0;    ///< 横道/安全岛面与连续RoadEdge求交前的外扩量。
    double inward_offset = 0.9;        ///< RoadEdge有效范围向路口中心的法向偏移量。
    double endpoint_adjustment = 3.0;  ///< 公共交点相对原骨干端点的最大移动距离，禁止越过另一端。
};

class NonMotorBackboneExtractor {
public:
    /// 创建骨干提取器；可配置道路边缘补连搜索范围、外扩量、内移量和端点调整上限。
    NonMotorBackboneExtractor(std::string out_dir = "./",
                              NonMotorGapOptions gap_options = NonMotorGapOptions());
    ~NonMotorBackboneExtractor();

    /// 提取非机动骨干路径图
    /// 横道和安全岛的有效骨干段以带Z的LINESTRING WKT写入调试CSV，便于核查几何结果。
    /// @param input 路口输入数据，包含人行横道、安全岛、车道组、车道及连通关系。
    /// @return 包含节点、骨干边、Gap 连接边和主环索引的骨干图；没有可用骨干段时返回空图。
    NonMotorBackbone run(const IntersectionInput& input) const;

    /// 输出骨干路径图：使用shapefile格式输出
    /// 会写出逐边骨干图层、闭合主环图层及分支图层；输出目录由底层 Shapefile 写入器创建。
    /// @param backbone 已提取的骨干图。
    /// @param out_dir Shapefile 输出目录。
    /// @return 底层骨干 Shapefile 写入操作成功时返回 true。
    bool output(const NonMotorBackbone& backbone, std::string out_dir) const;

private:
    std::string out_dir_;
    NonMotorGapOptions gap_options_;
    mutable std::ofstream fs_wkt_;
    void logSpineWkt(const std::string& kind, const std::string& source_id, const Vec3d& first, const Vec3d& last) const;
    void logSpineWkt(const std::string& kind, const std::string& source_id, const std::vector<Vec3d>& pts) const;
};

}  // namespace isg
