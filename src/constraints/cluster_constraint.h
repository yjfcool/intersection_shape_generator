#pragma once

#include "domain/generation_types.h"
#include "generation/connectivity_generation_context.h"

namespace isg {

/// 使用静态簇拓扑和已接受曲线审计当前曲线的非端点交叉/贴合。
std::vector<ConstraintResult> evaluateClusterIntersections(
    const BezierCurve& curve,
    const CurveGenerationContext& context,
    const GenerationState& state);

}  // 命名空间 isg
