#pragma once

#include "constraints/cluster_order.h"
#include "types.h"

#include <functional>
#include <vector>

namespace isg {

struct FinalPairFinding {
    ConnId id_a;
    ConnId id_b;
    CrossExemption exemption;
    Vec2d location;

    FinalPairFinding() : exemption(CrossExemption::None), location(0, 0) {}
};

typedef std::function<bool(
    const CurvePair&, const BezierCurve&, const BezierCurve&, Vec2d&)>
    FinalPairViolationDetector;

/// 全量遍历最终 Cluster 配对并将 finding 投影到兼容状态字段。
class FinalPairAuditor {
public:
    std::vector<FinalPairFinding> audit(
        const std::vector<CurvePair>& pairs,
        const std::vector<ConnectivityCurve>& results,
        const FinalPairViolationDetector& detector) const;

    void apply(const std::vector<FinalPairFinding>& findings,
               std::vector<ConnectivityCurve>& results) const;
};

}  // 命名空间 isg
