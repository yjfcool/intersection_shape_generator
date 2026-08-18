#pragma once

#include "domain/generation_types.h"
#include "generation/connectivity_generation_context.h"

#include <vector>

namespace isg {

struct UTurnFamilyInfo;

struct CurveInitializationOptions {
    std::vector<double> ordinary_alphas;
    bool allow_fixed_shape;
    double uturn_min_lead0;
    double uturn_min_lead1;
    double uturn_arc_alpha;
    double uturn_aligned_point_stagger;
    double uturn_lead0_extra_after_align;
    double uturn_lead1_extra_after_align;
    double uturn_aligned_entry_stagger;
    double uturn_aligned_exit_stagger;

    CurveInitializationOptions();

    /// 使用预处理快照填充 U-turn 初始化下限，不改变其他候选参数。
    void applyUTurnFamily(const UTurnFamilyInfo& family);
};

/// 按 U-turn、可保留 fixed shape、普通曲线的优先级构造未审计候选。
/// 场景相关参数由调用方提供，约束报告必须随后由 ConstraintEvaluator 写入。
class CurveInitializerRegistry {
public:
    std::vector<CurveCandidate> build(
        const CurveGenerationContext& context,
        const CurveInitializationOptions& options =
            CurveInitializationOptions()) const;
};

}  // 命名空间 isg
