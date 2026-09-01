#pragma once

#include "generation/curve_sampling.h"

#include <functional>
#include <limits>

namespace isg {

struct SegmentedUTurnAudit {
    bool physical_violation = false;
    int sibling_crosses = 0;
};

typedef std::function<SegmentedUTurnAudit(const BezierCurve&, bool)>
    SegmentedUTurnAuditor;

class SegmentedUTurnCandidateSearch {
public:
    bool search(
        const Vec2d& entry, const Vec2d& entry_tangent,
        const Vec2d& exit, const Vec2d& exit_tangent,
        const IntersectionInput& input,
        const std::vector<SampledSiblingCurve>& sampled_siblings,
        const SegmentedUTurnAuditor& audit, bool include_fence,
        BezierCurve& curve, double min_lead0 = 0.0, double min_lead1 = 0.0,
        const std::vector<Crosswalk>* crosswalks_for_clearance = nullptr,
        double aligned_point_stagger = 0.0,
        double base_lead0_extra_after_align = 0.0,
        double base_lead1_extra_after_align = 0.0,
        double aligned_family_station = std::numeric_limits<double>::quiet_NaN(),
        double aligned_entry_stagger = std::numeric_limits<double>::quiet_NaN(),
        double aligned_exit_stagger = std::numeric_limits<double>::quiet_NaN(),
        /// 家族分档步长(米)，透传给 UTurnCurveInitializer::buildSegmented，
        /// 用于把横向偏置钳制在本成员的家族槽位内。0 表示不做钳制。
        double family_stagger_step = 0.0) const;
};

}  // 命名空间 isg
