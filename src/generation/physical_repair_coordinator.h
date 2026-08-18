#pragma once

#include "types.h"

#include <functional>

namespace isg {

struct PhysicalRiskSnapshot {
    bool obstacle;
    bool boundary;
    bool fence;
    int sibling_crosses;

    PhysicalRiskSnapshot()
        : obstacle(false), boundary(false), fence(false),
          sibling_crosses(0) {}

    bool physical() const { return obstacle || boundary || fence; }
};

typedef std::function<PhysicalRiskSnapshot(const BezierCurve&)>
    PhysicalRiskAuditor;
typedef std::function<bool(const BezierCurve&, BezierCurve&)>
    PhysicalRepairAttempt;

struct PhysicalRepairResult {
    BezierCurve curve;
    PhysicalRiskSnapshot risk;
    bool boundary_repaired;

    PhysicalRepairResult() : boundary_repaired(false) {}
};

/// 固定顺序组织 U-turn Boundary 回退、普通 Boundary 修复和 Obstacle 修复。
class PhysicalRepairCoordinator {
public:
    PhysicalRepairResult repair(
        const BezierCurve& initial,
        const BezierCurve& current,
        bool geometric_uturn,
        const PhysicalRiskAuditor& auditor,
        const PhysicalRepairAttempt& boundary_repair,
        const PhysicalRepairAttempt& obstacle_repair) const;
};

}  // 命名空间 isg
