#pragma once

#include "domain/generation_types.h"
#include "generation/connectivity_generation_context.h"

namespace isg {

/// 检查普通直行/转向曲线的端点、控制轴和形态；严格避让回退允许两段
/// G1 路点曲线，普通自然候选仍保持单段 cubic。
ConstraintResult evaluateOrdinaryShape(const BezierCurve& curve,
                                       const CurveGenerationContext& context,
                                       const GenerationState& state);

/// 检查几何 U 型调头的三段式、平齐和首尾直行段。
ConstraintResult evaluateUTurnShape(const BezierCurve& curve,
                                    const CurveGenerationContext& context);

/// 检查 U 型首尾跨越相关 Crosswalk 且中弧不侵入。
ConstraintResult evaluateUTurnCrosswalk(const BezierCurve& curve,
                                        const CurveGenerationContext& context);

}  // 命名空间 isg
