#pragma once

#include "types.h"

#include <cstddef>
#include <string>
#include <unordered_set>
#include <vector>

namespace isg {

class ClusterOrderSolver;

struct UTurnFamilyRank {
    std::size_t family_size;
    std::size_t reverse_radius_rank;

    UTurnFamilyRank() : family_size(0), reverse_radius_rank(0) {}
};

struct UTurnLeadFloors {
    double lead0;
    double lead1;
    bool crosswalk0;
    bool crosswalk1;
    std::unordered_set<std::string> crosswalk_ids;

    UTurnLeadFloors()
        : lead0(0.0), lead1(0.0),
          crosswalk0(false), crosswalk1(false) {}
};

/// 单条 U-turn 使用的家族级只读几何快照。
struct UTurnFamilyInfo {
    bool geometric_uturn;
    Vec2d axis;
    double radius;
    double lead0;
    double lead1;
    double aligned_station;
    UTurnFamilyRank rank;
    std::vector<Crosswalk> clearance_crosswalks;

    UTurnFamilyInfo()
        : geometric_uturn(false), axis(1, 0), radius(0.0),
          lead0(0.0), lead1(0.0),
          aligned_station(0.0) {}
};

/// U-turn 家族共享的只读几何元数据。
class UTurnFamilyBuilder {
public:
    UTurnFamilyInfo build(
        const Connectivity& connectivity,
        const IntersectionInput& input,
        UTurnAlignmentScope scope,
        const ClusterOrderSolver* topology = nullptr,
        bool require_shared_endpoint_pair_for_rank = false) const;

    bool isGeometricUTurn(const Connectivity& connectivity, const IntersectionInput& input) const;

    double radiusKey(const Connectivity& connectivity, const IntersectionInput& input) const;

    Vec2d alignmentAxis(const Vec2d& entry_tangent, const Vec2d& exit_tangent) const;

    double requiredAlignedStation(const Vec2d& entry_point,
                                  const Vec2d& entry_tangent,
                                  const Vec2d& exit_point,
                                  const Vec2d& exit_tangent,
                                  double min_lead0,
                                  double min_lead1) const;

    double requiredAlignedStationOnAxis(const Vec2d& entry_point,
                                        const Vec2d& entry_tangent,
                                        const Vec2d& exit_point,
                                        const Vec2d& exit_tangent,
                                        const Vec2d& axis,
                                        double min_lead0,
                                        double min_lead1) const;

    std::vector<const Connectivity*> alignmentComponent(
        const Connectivity& connectivity,
        const IntersectionInput& input,
        UTurnAlignmentScope scope,
        const ClusterOrderSolver* topology = nullptr,
        bool require_shared_endpoint_pair = false) const;

    UTurnFamilyRank radiusRank(
        const Connectivity& connectivity,
        const IntersectionInput& input,
        UTurnAlignmentScope scope,
        const ClusterOrderSolver* topology = nullptr,
        bool require_shared_endpoint_pair = false) const;

    UTurnLeadFloors leadFloors(const Vec2d& entry_point,
                               const Vec2d& entry_tangent,
                               const Vec2d& exit_point,
                               const Vec2d& exit_tangent,
                               const IntersectionInput& input,
                               double no_crosswalk_min_lead = 2.0) const;

    std::vector<Crosswalk> clearanceCrosswalks(
            const IntersectionInput& input, const UTurnLeadFloors& floors) const;

    double alignmentComponentStation(
        const Connectivity& connectivity,
        const IntersectionInput& input,
        UTurnAlignmentScope scope,
        const ClusterOrderSolver* topology = nullptr) const;

    void enforceAlignmentStation(const Vec2d& entry_point,
                                 const Vec2d& entry_tangent,
                                 const Vec2d& exit_point,
                                 const Vec2d& exit_tangent,
                                 double target_station,
                                 double& min_lead0,
                                 double& min_lead1) const;
};

}  // 命名空间 isg
