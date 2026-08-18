#include "constraints/cluster_constraint.h"

#include "constraints/cluster_order.h"
#include "curve/curve_utils.h"
#include "geometry/predicates.h"

namespace isg {

std::vector<ConstraintResult> evaluateClusterIntersections(
    const BezierCurve& curve,
    const CurveGenerationContext& context,
    const GenerationState& state) {
    std::vector<ConstraintResult> results;
    if (!context.connectivity || !context.scene || !context.scene->cluster_topology)
        return results;
    const ClusterTopology& topology = *context.scene->cluster_topology;
    for (const auto& pair : topology.pairs) {
        ConnId sibling_id;
        if (pair.id_a == context.connectivity->id) sibling_id = pair.id_b;
        else if (pair.id_b == context.connectivity->id) sibling_id = pair.id_a;
        else continue;
        std::unordered_map<ConnId, BezierCurve>::const_iterator sibling = state.accepted_curves.find(sibling_id);
        if (sibling == state.accepted_curves.end()) continue;

        const bool adherence = curvesHaveForbiddenAdherenceForAudit(
            curve, sibling->second, context.profile.cluster_endpoint_tolerance);
        const bool crossing = curvesIntersectBusiness(
            curve, sibling->second, context.profile.cluster_endpoint_tolerance);
        const bool structural_exempt = pair.exempt == CrossExemption::StructuralCross &&
                                       context.profile.allow_structural_cross;
        const bool obstacle_exempt = pair.exempt == CrossExemption::ObstacleCross &&
                                     context.profile.allow_obstacle_cross_exemption;
        if (!adherence && crossing && (structural_exempt || obstacle_exempt)) {
            ConstraintResult result;
            result.id = "cluster.intersection." + sibling_id;
            result.severity = ConstraintSeverity::Hard;
            result.state = ConstraintState::Exempt;
            result.reason = structural_exempt ? "structural crossing exemption" :
                                                "obstacle-adjacent crossing exemption";
            results.push_back(result);
            continue;
        }
        if (!adherence && !crossing)
            continue;

        ConstraintResult result;
        result.id = "cluster.intersection." + sibling_id;
        result.severity = ConstraintSeverity::Hard;
        result.state = ConstraintState::Violated;
        result.reason = adherence ? "same-cluster curves adhere or overlap away from endpoints" :
                                     "same-cluster curves cross away from connection endpoints";
        result.violation = adherence ? 2.0 : 1.0;
        if (crossing) {
            const std::vector<Vec2d> crossings = curveCrossings(curve, sibling->second, 0.05);
            if (!crossings.empty()) result.locations.push_back(crossings.front());
        }
        results.push_back(result);
    }
    return results;
}

}  // 命名空间 isg
