#pragma once

#include "types.h"

#include <string>
#include <cstddef>
#include <vector>

namespace isg {

enum class CandidateOrigin {
    FixedShape,
    NaturalSingleCubic,
    Hermite,
    SegmentedUTurn,
    BoundaryFit,
    ObstacleBypass,
    Optimized,
    Repaired,
    Fallback
};

enum class ConstraintSeverity { Hard, Soft, Advisory };
enum class ConstraintState { Satisfied, Violated, Exempt, NotApplicable };

struct ConstraintResult {
    std::string id;
    ConstraintSeverity severity = ConstraintSeverity::Advisory;
    ConstraintState state = ConstraintState::NotApplicable;
    double violation = 0.0;
    double penalty = 0.0;
    std::vector<Vec2d> locations;
    std::string reason;
};

struct ConstraintReport {
    std::vector<ConstraintResult> results;

    bool hasHardViolation() const {
        for (const auto& result : results)
            if (result.severity == ConstraintSeverity::Hard && result.state == ConstraintState::Violated)
                return true;
        return false;
    }
};

struct CurveCandidate {
    BezierCurve curve;
    CandidateOrigin origin = CandidateOrigin::Fallback;
    ConstraintReport report;
    int hard_violation_count = 0;
    double physical_violation = 0.0;
    int new_cluster_crosses = 0;
    int shape_hard_violations = 0;
    int fixed_shape_loss = 0;
    double soft_penalty = 0.0;
    double quality_score = 0.0;
    bool preserves_fixed_shape = false;
    std::size_t geometry_revision = 0;
    std::size_t audited_revision = static_cast<std::size_t>(-1);

    void replaceCurve(const BezierCurve& replacement) {
        curve = replacement;
        ++geometry_revision;
        audited_revision = static_cast<std::size_t>(-1);
        report.results.clear();
    }

    void acceptReport(const ConstraintReport& replacement) {
        report = replacement;
        audited_revision = geometry_revision;
    }

    bool hasCurrentReport() const {
        return audited_revision == geometry_revision;
    }
};

struct ConstraintProfile {
    bool require_single_segment = false;
    bool enforce_fence = true;
    bool check_self_intersection = true;
    bool check_obstacle = true;
    bool check_boundary = true;
    bool check_g1 = false;
    bool check_curvature = false;
    bool check_ordinary_shape = false;
    bool check_uturn_shape = false;
    bool check_crosswalk = false;
    bool check_cluster = false;
    bool require_fixed_shape_preservation = false;
    bool allow_structural_cross = false;
    bool allow_obstacle_cross_exemption = false;
    double obstacle_clearance = 0.0;
    double road_edge_clearance = 0.0;
    double max_curvature = 0.25;
    double g1_angle_deg = 15.0;
    double cluster_endpoint_tolerance = 1e-4;
    int samples = 25;
};

}  // 命名空间 isg
