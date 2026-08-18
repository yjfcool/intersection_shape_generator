#pragma once

#include "constraints/cluster_order.h"
#include "generation/connectivity_generation_context.h"
#include "optimizer/lbfgs_solver.h"
#include "optimizer/penalty_cost.h"
#include "optimizer/sdf_field.h"
#include "types.h"

#include <unordered_set>

namespace isg {

/// 无状态连通曲线门面背后的单次调用可变会话。
class ConnectivityGenerationSession {
public:
    ConnectivityGenerationSession(const LBFGSConfig& config, const ConnectivityDirectionConfig& direction_config);

    std::vector<ConnectivityCurve> run(const IntersectionInput& input, SDFField& sdf, double* out_ms);

private:
    LBFGSSolver solver_;
    ClusterOrderSolver cluster_solver_;
    ConnectivityDirectionConfig direction_cfg_;

private:
    ConnectivityCurve generateOne(
        const CurveGenerationContext& context, const std::vector<SiblingCurve>& siblings, bool* physical_risk);

    std::vector<SiblingCurve> buildSiblings(
        const ConnId& id, const std::unordered_map<ConnId, BezierCurve>& completed,
        const ClusterOrderSolver& cluster_solver, const std::vector<Connectivity>& connectivities,
        bool constrained_only, const std::unordered_set<ConnId>* fixed_shape_ids) const;

    void validate(ConnectivityCurve& curve, const IntersectionInput& input, const SDFField& sdf) const;
};

}  // 命名空间 isg
