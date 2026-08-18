#include "preprocessing/uturn_family_builder.h"

#include "constraints/cluster_order.h"
#include "preprocessing/crosswalk_clearance_calculator.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <unordered_set>

namespace isg {

UTurnFamilyInfo UTurnFamilyBuilder::build(
    const Connectivity& connectivity, const IntersectionInput& input,
    UTurnAlignmentScope scope, const ClusterOrderSolver* topology,
    bool require_shared_endpoint_pair_for_rank) const {
    UTurnFamilyInfo info;
    info.geometric_uturn = isGeometricUTurn(connectivity, input);
    if (!info.geometric_uturn)
        return info;

    const std::pair<Vec2d, Vec2d> entry =
        input.entryPtDir(connectivity.entry_lane_id);
    const std::pair<Vec2d, Vec2d> exit =
        input.exitPtDir(connectivity.exit_lane_id);
    info.axis = alignmentAxis(entry.second, exit.second);
    info.radius = radiusKey(connectivity, input);

    const UTurnLeadFloors floors = leadFloors(
        entry.first, entry.second, exit.first, exit.second, input);
    info.lead0 = floors.lead0;
    info.lead1 = floors.lead1;
    info.aligned_station = alignmentComponentStation(
        connectivity, input, scope, topology);
    enforceAlignmentStation(
        entry.first, entry.second, exit.first, exit.second,
        info.aligned_station, info.lead0, info.lead1);
    info.rank = radiusRank(
        connectivity, input, scope, topology,
        require_shared_endpoint_pair_for_rank);
    info.clearance_crosswalks = clearanceCrosswalks(input, floors);
    return info;
}

bool UTurnFamilyBuilder::isGeometricUTurn(
    const Connectivity& connectivity, const IntersectionInput& input) const {
    const std::pair<Vec2d, Vec2d> entry =
        input.entryPtDir(connectivity.entry_lane_id);
    const std::pair<Vec2d, Vec2d> exit =
        input.exitPtDir(connectivity.exit_lane_id);
    return entry.second.norm() > 1e-8 && exit.second.norm() > 1e-8 &&
           entry.second.normalized().dot(exit.second.normalized()) < -0.5;
}

Vec2d UTurnFamilyBuilder::alignmentAxis(
    const Vec2d& entry_tangent, const Vec2d& exit_tangent) const {
    Vec2d t0 = entry_tangent.norm() > 1e-8
        ? entry_tangent.normalized() : Vec2d(1, 0);
    Vec2d t1 = exit_tangent.norm() > 1e-8
        ? exit_tangent.normalized() : -t0;
    Vec2d axis = t0 - t1;
    if (axis.norm() < 1e-8)
        axis = t0;
    axis.normalize();
    if (axis.dot(t0) < 0.0)
        axis = -axis;
    return axis;
}

double UTurnFamilyBuilder::radiusKey(
    const Connectivity& connectivity, const IntersectionInput& input) const {
    const std::pair<Vec2d, Vec2d> entry =
        input.entryPtDir(connectivity.entry_lane_id);
    const std::pair<Vec2d, Vec2d> exit =
        input.exitPtDir(connectivity.exit_lane_id);
    const Vec2d axis = alignmentAxis(entry.second, exit.second);
    const Vec2d lateral(-axis.y(), axis.x());
    return std::max(0.0, std::abs((exit.first - entry.first).dot(lateral)));
}

double UTurnFamilyBuilder::requiredAlignedStation(
    const Vec2d& entry_point, const Vec2d& entry_tangent,
    const Vec2d& exit_point, const Vec2d& exit_tangent,
    double min_lead0, double min_lead1) const {
    return requiredAlignedStationOnAxis(
        entry_point, entry_tangent, exit_point, exit_tangent,
        alignmentAxis(entry_tangent, exit_tangent), min_lead0, min_lead1);
}

double UTurnFamilyBuilder::requiredAlignedStationOnAxis(
    const Vec2d& entry_point, const Vec2d& entry_tangent,
    const Vec2d& exit_point, const Vec2d& exit_tangent,
    const Vec2d& axis, double min_lead0, double min_lead1) const {
    if (axis.norm() < 1e-8)
        return -std::numeric_limits<double>::infinity();
    const Vec2d t0 = entry_tangent.norm() > 1e-8
        ? entry_tangent.normalized() : Vec2d(1, 0);
    const Vec2d t1 = exit_tangent.norm() > 1e-8
        ? exit_tangent.normalized() : -t0;
    const Vec2d exit_back = -t1;
    const double s0 = entry_point.dot(axis);
    const double s1 = exit_point.dot(axis);
    const double c0 = std::max(0.2, t0.dot(axis));
    const double c1 = std::max(0.2, exit_back.dot(axis));
    const double lead0 = std::max(min_lead0, std::max(0.0, (std::max(s0, s1) - s0) / c0));
    const double lead1 = std::max(min_lead1, std::max(0.0, (std::max(s0, s1) - s1) / c1));
    return std::max(s0 + lead0 * c0, s1 + lead1 * c1);
}

namespace {

LaneGroupId laneGroupIdForRole(const Connectivity& connectivity,
                               const IntersectionInput& input,
                               bool entry_side) {
    const LaneGroupId& explicit_id = entry_side
        ? connectivity.enterGroupId : connectivity.exitGroupId;
    if (!explicit_id.empty())
        return explicit_id;
    const GroupRole role = entry_side ? GroupRole::Entry : GroupRole::Exit;
    const LaneId& lane_id = entry_side
        ? connectivity.entry_lane_id : connectivity.exit_lane_id;
    for (const auto& group : input.lane_groups) {
        if (group.role == role &&
            std::find(group.lanes.begin(), group.lanes.end(), lane_id) !=
                group.lanes.end())
            return group.id;
    }
    const Lane* lane = input.findLane(lane_id);
    return lane ? lane->groupId : LaneGroupId();
}

bool sameAlignmentFamily(const Connectivity& a, const Connectivity& b,
                         const IntersectionInput& input, bool entry_side,
                         UTurnAlignmentScope scope) {
    if (scope == UTurnAlignmentScope::LaneGroup) {
        const LaneGroupId ga = laneGroupIdForRole(a, input, entry_side);
        const LaneGroupId gb = laneGroupIdForRole(b, input, entry_side);
        return !ga.empty() && ga == gb;
    }
    const LaneId& la = entry_side ? a.entry_lane_id : a.exit_lane_id;
    const LaneId& lb = entry_side ? b.entry_lane_id : b.exit_lane_id;
    if (la.empty() || lb.empty())
        return false;
    if (la == lb)
        return true;
    const std::pair<Vec2d, Vec2d> pa = entry_side
        ? input.entryPtDir(la) : input.exitPtDir(la);
    const std::pair<Vec2d, Vec2d> pb = entry_side
        ? input.entryPtDir(lb) : input.exitPtDir(lb);
    return pa.second.norm() > 1e-8 && pb.second.norm() > 1e-8 &&
           (pa.first - pb.first).norm() <= 0.15 &&
           pa.second.normalized().dot(pb.second.normalized()) > 0.99;
}

}  // namespace

std::vector<const Connectivity*> UTurnFamilyBuilder::alignmentComponent(
    const Connectivity& connectivity, const IntersectionInput& input,
    UTurnAlignmentScope scope, const ClusterOrderSolver* topology,
    bool require_shared_endpoint_pair) const {
    std::unordered_set<ConnId> ids;
    std::vector<const Connectivity*> component;
    ids.insert(connectivity.id);
    component.push_back(&connectivity);
    for (std::size_t i = 0; i < component.size(); ++i) {
        const Connectivity& current = *component[i];
        for (const auto& other : input.connectivities) {
            if (ids.count(other.id) || !isGeometricUTurn(other, input))
                continue;
            if (!sameAlignmentFamily(current, other, input, true, scope) &&
                !sameAlignmentFamily(current, other, input, false, scope))
                continue;
            if (require_shared_endpoint_pair && topology &&
                !topology->isSharedEndpoint(current.id, other.id))
                continue;
            ids.insert(other.id);
            component.push_back(&other);
        }
    }
    return component;
}

UTurnFamilyRank UTurnFamilyBuilder::radiusRank(
    const Connectivity& connectivity, const IntersectionInput& input,
    UTurnAlignmentScope scope, const ClusterOrderSolver* topology,
    bool require_shared_endpoint_pair) const {
    const std::vector<const Connectivity*> component = alignmentComponent(
        connectivity, input, scope, topology, require_shared_endpoint_pair);
    std::vector<std::pair<double, ConnId> > ordered;
    ordered.reserve(component.size());
    for (const Connectivity* member : component)
        ordered.push_back(std::make_pair(radiusKey(*member, input), member->id));
    std::sort(ordered.begin(), ordered.end(),
              [](const std::pair<double, ConnId>& a,
                 const std::pair<double, ConnId>& b) {
                  if (std::abs(a.first - b.first) > 1e-9)
                      return a.first < b.first;
                  return a.second < b.second;
              });
    UTurnFamilyRank result;
    result.family_size = ordered.size();
    for (std::size_t i = 0; i < ordered.size(); ++i) {
        if (ordered[i].second == connectivity.id) {
            result.reverse_radius_rank = ordered.size() - i;
            break;
        }
    }
    return result;
}

UTurnLeadFloors UTurnFamilyBuilder::leadFloors(
    const Vec2d& entry_point, const Vec2d& entry_tangent,
    const Vec2d& exit_point, const Vec2d& exit_tangent,
    const IntersectionInput& input, double no_crosswalk_min_lead) const {
    const CrosswalkClearanceCalculator calculator;
    const CrosswalkClearanceResult entry =
        calculator.ahead(entry_point, entry_tangent, input);
    const CrosswalkClearanceResult exit =
        calculator.behind(exit_point, exit_tangent, input);
    UTurnLeadFloors floors;
    floors.crosswalk0 = entry.found;
    floors.crosswalk1 = exit.found;
    if (entry.found)
        floors.crosswalk_ids.insert(entry.crosswalk_id);
    if (exit.found)
        floors.crosswalk_ids.insert(exit.crosswalk_id);
    if (floors.crosswalk0 || floors.crosswalk1) {
        floors.lead0 = floors.crosswalk0 ? entry.clearance : 0.0;
        floors.lead1 = floors.crosswalk1 ? exit.clearance : 0.0;
    } else {
        floors.lead0 = no_crosswalk_min_lead;
        floors.lead1 = no_crosswalk_min_lead;
    }
    return floors;
}

std::vector<Crosswalk> UTurnFamilyBuilder::clearanceCrosswalks(
    const IntersectionInput& input, const UTurnLeadFloors& floors) const {
    if (floors.crosswalk_ids.size() < 2)
        return input.crosswalks;
    std::vector<Crosswalk> selected;
    selected.reserve(floors.crosswalk_ids.size());
    for (const auto& crosswalk : input.crosswalks) {
        if (floors.crosswalk_ids.count(crosswalk.id))
            selected.push_back(crosswalk);
    }
    return selected.empty() ? input.crosswalks : selected;
}

double UTurnFamilyBuilder::alignmentComponentStation(
    const Connectivity& connectivity, const IntersectionInput& input,
    UTurnAlignmentScope scope, const ClusterOrderSolver* topology) const {
    const std::pair<Vec2d, Vec2d> entry =
        input.entryPtDir(connectivity.entry_lane_id);
    const std::pair<Vec2d, Vec2d> exit =
        input.exitPtDir(connectivity.exit_lane_id);
    const Vec2d axis = alignmentAxis(entry.second, exit.second);
    if (axis.norm() < 1e-8)
        return -std::numeric_limits<double>::infinity();

    const std::vector<const Connectivity*> component = alignmentComponent(
        connectivity, input, scope, topology, false);
    double target = -std::numeric_limits<double>::infinity();
    for (const Connectivity* member : component) {
        const std::pair<Vec2d, Vec2d> member_entry =
            input.entryPtDir(member->entry_lane_id);
        const std::pair<Vec2d, Vec2d> member_exit =
            input.exitPtDir(member->exit_lane_id);
        const UTurnLeadFloors floors = leadFloors(
            member_entry.first, member_entry.second,
            member_exit.first, member_exit.second, input);
        target = std::max(target, requiredAlignedStationOnAxis(
            member_entry.first, member_entry.second,
            member_exit.first, member_exit.second, axis,
            floors.lead0, floors.lead1));
    }
    return target;
}

void UTurnFamilyBuilder::enforceAlignmentStation(
    const Vec2d& entry_point, const Vec2d& entry_tangent,
    const Vec2d& exit_point, const Vec2d& exit_tangent,
    double target_station, double& min_lead0, double& min_lead1) const {
    if (!std::isfinite(target_station))
        return;
    const Vec2d t0 = entry_tangent.norm() > 1e-8
        ? entry_tangent.normalized() : Vec2d(1, 0);
    const Vec2d t1 = exit_tangent.norm() > 1e-8
        ? exit_tangent.normalized() : -t0;
    const Vec2d axis = alignmentAxis(t0, t1);
    const Vec2d exit_back = -t1;
    const double c0 = std::max(0.2, t0.dot(axis));
    const double c1 = std::max(0.2, exit_back.dot(axis));
    min_lead0 = std::max(min_lead0,
                         (target_station - entry_point.dot(axis)) / c0);
    min_lead1 = std::max(min_lead1,
                         (target_station - exit_point.dot(axis)) / c1);
    min_lead0 = std::max(0.0, min_lead0);
    min_lead1 = std::max(0.0, min_lead1);
}

}  // 命名空间 isg
