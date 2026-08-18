#pragma once

#include "types.h"

#include <functional>

namespace isg {

struct BoundedUTurnCandidateAudit {
    bool physical_violation;
    int sibling_crosses;

    BoundedUTurnCandidateAudit()
        : physical_violation(false), sibling_crosses(0) {}
};

typedef std::function<BoundedUTurnCandidateAudit(const BezierCurve&)>
    BoundedUTurnCandidateAuditor;

/// 有界单段 U-turn 候选枚举；约束语义由调用方审计回调提供。
class BoundedUTurnCandidateSearch {
public:
    bool search(const Vec2d& entry,
                const Vec2d& entry_tangent,
                const Vec2d& exit,
                const Vec2d& exit_tangent,
                const BezierCurve& reference,
                const BoundedUTurnCandidateAuditor& auditor,
                BezierCurve& result,
                double min_lead0 = 0.0,
                double min_lead1 = 0.0) const;
};

}  // 命名空间 isg
