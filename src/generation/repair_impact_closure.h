#pragma once

#include "constraints/cluster_order.h"

#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace isg {

struct RepairImpactClosure {
    std::unordered_set<ConnId> curve_ids;
    std::vector<CurvePair> pairs;
    std::unordered_map<ConnId, std::vector<ConnId> > neighbors;
};

/// 基于非结构性 Cluster 配对计算修复后的单跳影响闭包。
class RepairImpactClosureBuilder {
public:
    RepairImpactClosure build(
        const std::vector<CurvePair>& pairs,
        const std::unordered_set<ConnId>& changed_ids) const;

    std::unordered_map<ConnId, std::vector<ConnId> > neighbors(
        const std::vector<CurvePair>& pairs,
        const std::unordered_set<ConnId>& available_ids) const;
};

}  // 命名空间 isg
