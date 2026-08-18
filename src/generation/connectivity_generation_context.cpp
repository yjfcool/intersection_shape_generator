#include "generation/connectivity_generation_context.h"

namespace isg {

CurveGenerationContext CurveGenerationContextBuilder::build(
    const SceneContext& scene, const Connectivity& connectivity,
    UTurnAlignmentScope scope, const ClusterOrderSolver* topology,
    const UTurnFamilyInfo* family_snapshot) const {
    CurveGenerationContext context;
    context.scene = &scene;
    context.connectivity = &connectivity;
    context.entry = scene.view.entryFrame(connectivity.entry_lane_id);
    context.exit = scene.view.exitFrame(connectivity.exit_lane_id);

    const std::unordered_map<ConnId, TurnMetadata>::const_iterator metadata =
        scene.turn_metadata.find(connectivity.id);
    context.turn = metadata == scene.turn_metadata.end()
        ? connectivity.turn_type : metadata->second.turn;
    context.uturn_family = family_snapshot
        ? *family_snapshot
        : UTurnFamilyBuilder().build(
              connectivity, scene.view.input(), scope, topology);
    return context;
}

}  // 命名空间 isg
