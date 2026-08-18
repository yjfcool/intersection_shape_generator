#pragma once

#include "domain/scene_view.h"

#include <unordered_map>

namespace isg {

class SDFField;
struct ClusterTopology;

/// 只读的转向元数据，供规划、初始化和审计阶段共享。
struct TurnMetadata {
    ConnTurnType turn = ConnTurnType::Unknown;
    Vec2d entry_point{0, 0};
    Vec2d exit_point{0, 0};
    Vec2d entry_tangent{1, 0};
    Vec2d exit_tangent{1, 0};
};

/// 生成阶段共享的只读场景对象。
struct SceneContext {
    SceneView view;
    const SDFField* sdf = nullptr;
    const SDFField* coarse_sdf = nullptr;
    const ClusterTopology* cluster_topology = nullptr;
    std::unordered_map<ConnId, TurnMetadata> turn_metadata;

    explicit SceneContext(const IntersectionInput& input) : view(input) {
        for (const auto& conn : input.connectivities) {
            TurnMetadata metadata;
            metadata.turn = conn.turn_type;
            const std::pair<Vec2d, Vec2d> entry = view.entryFrame(conn.entry_lane_id);
            const std::pair<Vec2d, Vec2d> exit = view.exitFrame(conn.exit_lane_id);
            metadata.entry_point = entry.first;
            metadata.entry_tangent = entry.second;
            metadata.exit_point = exit.first;
            metadata.exit_tangent = exit.second;
            turn_metadata[conn.id] = metadata;
        }
    }
    void bindSdf(const SDFField* fine, const SDFField* coarse = nullptr) {
        sdf = fine;
        coarse_sdf = coarse;
    }
    void bindClusterTopology(const ClusterTopology* topology) {
        cluster_topology = topology;
    }
};

}  // 命名空间 isg
