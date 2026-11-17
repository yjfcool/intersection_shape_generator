#pragma once

namespace isg {

/// 非机动逻辑路径下的挂接类型（由分类器在 NonMotorConnectivityGenerator 内确定）
enum class NonMotorConnType {
    StraightThroughIsland,  ///< I.  直行穿安全岛（进出组均≤1，安全岛可穿，距端点≤20m）
    LeftTurnL,              ///< J.  左转L型（进入非机动，非提右车道）
    LeftTurnThroughIsland,  ///< K.  左转穿安全岛（进入非机动，提右车道，安全岛可穿）
    LeftTurnMotorToNM,      ///< L.  左转（进入机动，退出非机动，提右车道，安全岛可穿）
    UTurnU,                 ///< M.  掉头U型（进入非机动，或进入机动+安全岛可穿）
    FallbackMotor,          ///< 兜底：降级为机动逻辑（不应发生，但保留安全出口）
};

/// 路口内非机动骨干节点类型。
/// 骨干节点描述横道/安全岛穿行点、道路边缘转角以及少量外接车道端点。
enum class NonMotorBackboneNodeType {
    CrosswalkPt = 0,
    SafetyIslandPt = 1,
    RoadEdgeCorner = 2,
    LaneEndpoint = 3,
    Gap = 4,
};

/// 路口内非机动骨干边类型。
enum class NonMotorBackboneEdgeType {
    CrosswalkSpine = 0,
    SafetyIslandSpine = 1,
    Gap = 2,
    RoadEdge = 3,
    Branch = 4,
    Composite = 5,  ///< 删除环上共端点后由多段骨干/连接线合并形成的路径段。
};

/// 非机动骨干节点。
struct NonMotorBackboneNode {
    std::string id;
    NonMotorBackboneNodeType type = NonMotorBackboneNodeType::Gap;
    Vec3d point;
    std::string source_id;
};

/// 非机动骨干边。geometry 的首尾点与 from/to 节点一致。
struct NonMotorBackboneEdge {
    std::string id;
    std::string from;
    std::string to;
    NonMotorBackboneEdgeType type = NonMotorBackboneEdgeType::Gap;
    LineString2d geometry;
    /// 原始横道/安全岛面ID；合并、拆边时保留，Gap不冒充任何来源面。
    std::vector<std::string> source_surface_ids;
    bool ring = false;
    bool branch = false;
};

/// 整个路口的非机动骨干图。
/// ring_edges 表示闭合主环；branch_edges 表示所有未进入主环的骨干边（含局部 Gap）。
/// 两类边严格互斥，最终分类以 edge.ring 为准。
struct NonMotorBackbone {
    std::string id;
    std::vector<NonMotorBackboneNode> nodes;
    std::vector<NonMotorBackboneEdge> edges;
    /// 主环拓扑节点，不包括边内部的RoadEdge形态点；不可用其下标索引ring_geometry。
    std::vector<std::string> ring_nodes;
    std::vector<std::string> ring_edges;
    std::vector<std::string> branch_edges;
    LineString2d ring_geometry;
    bool closed = false;
    /// 最终无向边图的连通分量数；非空图期望为1，安全约束阻止补连时保留真实计数。
    size_t connected_components = 0;
    /// 无法在搜索范围内安全补连时的诊断，不通过越界/穿障伪造连通。
    std::vector<std::string> diagnostics;
};

}  // namespace isg
