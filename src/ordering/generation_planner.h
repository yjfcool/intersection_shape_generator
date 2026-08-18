#pragma once

#include "toolkits/toolkits.h"
#include "constraints/cluster_order.h"
#include "domain/scene_context.h"

#include <vector>

namespace isg {

struct GenerationBatch {
    int priority = 0;
    std::vector<ConnId> conn_ids;
};

struct GenerationPlan {
    std::vector<GenerationBatch> batches;
};

/// 稳定的几何优先生成顺序，与曲线构造方式解耦。
class GenerationPlanner {
public:
    GenerationPlan build(const IntersectionInput& input) const;
    GenerationPlan build(const IntersectionInput& input, const ClusterOrderSolver& topology) const;
    GenerationPlan build(const SceneContext& context, const ClusterOrderSolver& topology) const;
};

}  // 命名空间 isg
