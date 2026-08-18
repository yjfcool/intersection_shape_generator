#pragma once

#include "domain/generation_types.h"
#include "domain/scene_context.h"
#include "preprocessing/uturn_family_builder.h"

#include <unordered_map>
#include <vector>

namespace isg {

struct CurveGenerationContext {
    const SceneContext* scene = nullptr;
    const Connectivity* connectivity = nullptr;
    std::pair<Vec2d, Vec2d> entry;
    std::pair<Vec2d, Vec2d> exit;
    ConnTurnType turn = ConnTurnType::Unknown;
    UTurnFamilyInfo uturn_family;
    ConstraintProfile profile;
};

/// 从只读场景集中组装单条曲线上下文，不持有跨 generate() 的状态。
class CurveGenerationContextBuilder {
public:
    CurveGenerationContext build(
        const SceneContext& scene,
        const Connectivity& connectivity,
        UTurnAlignmentScope scope = UTurnAlignmentScope::LaneEndpoint,
        const ClusterOrderSolver* topology = nullptr,
        const UTurnFamilyInfo* family_snapshot = nullptr) const;
};

struct GenerationState {
    std::vector<ConnectivityCurve> results;
    std::unordered_map<ConnId, BezierCurve> accepted_curves;
};

}  // 命名空间 isg
