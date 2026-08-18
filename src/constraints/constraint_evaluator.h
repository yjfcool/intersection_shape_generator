#pragma once

#include "domain/generation_types.h"
#include "domain/scene_context.h"
#include "generation/connectivity_generation_context.h"

namespace isg {

/// 候选和最终输出共用的轻量审计入口。
class ConstraintEvaluator {
public:
    ConstraintReport evaluate(const BezierCurve& curve,
                              const SceneContext& context,
                              const ConstraintProfile& profile) const;

    ConstraintReport evaluate(const BezierCurve& curve,
                              const CurveGenerationContext& context,
                              const GenerationState& state) const;

    ConstraintReport evaluate(const CurveCandidate& candidate,
                              const CurveGenerationContext& context,
                              const GenerationState& state) const;

    void audit(CurveCandidate& candidate,
               const CurveGenerationContext& context,
               const GenerationState& state) const;
};

}  // 命名空间 isg
