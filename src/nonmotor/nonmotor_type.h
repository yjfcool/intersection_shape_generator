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
    bool ring = false;
    bool branch = false;
};

/// 整个路口的非机动骨干图。
/// ring_edges 表示主环，branch_edges 表示从环接出的少量分支。
struct NonMotorBackbone {
    std::string id;
    std::vector<NonMotorBackboneNode> nodes;
    std::vector<NonMotorBackboneEdge> edges;
    std::vector<std::string> ring_nodes;
    std::vector<std::string> ring_edges;
    std::vector<std::string> branch_edges;
    LineString2d ring_geometry;
    bool closed = false;
};

}  // namespace isg
