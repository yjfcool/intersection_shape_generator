#include "preprocessing/uturn_family_builder.h"

#include "constraints/cluster_order.h"
#include "preprocessing/crosswalk_clearance_calculator.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <unordered_set>

namespace isg {

namespace {

LaneGroupId laneGroupIdForRole(
        const Connectivity& connectivity, const IntersectionInput& input, bool entry_side);

bool sameAlignmentFamily(const Connectivity& a, const Connectivity& b,
        const IntersectionInput& input, bool entry_side, UTurnAlignmentScope scope);

}  // namespace

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

    // aligned_station 是沿当前连接自身 U 轴的标量。家族成员若共享
    // 同一入口/出口端点，而出口/入口切向略有差异，把这个标量分别投影
    // 回各自轴线会得到不同的 lead，首/尾平齐点随之产生可见台阶。
    // 对共享端点的成员先按各自快照求出所需 lead，再取家族最大值，
    // 让共享物理端点真正落在同一站位；不共享端点的另一侧仍保留各自站位。
    //
    // 判定"共享"必须用端点坐标重合，而不是 sameAlignmentFamily(scope)：
    // 后者在 LaneGroup 口径下把同车道组的不同车道也算作同族，而 lead 是
    // 从各自端点起算的长度，只有起点重合时才可比。对起点不同的成员取
    // lead 最大值，等价于把离路口更远的成员的 lead 套到更近的成员上，
    // 两者的轴向站位会被推开而不是对齐——LaneGroup 口径的"near/far 共站位"
    // 正是这样被破坏的（实测 near/far 由同站位变为相差 3m）。
    const double kSharedEndpointEps = 1e-3;
    const std::vector<const Connectivity*> component = alignmentComponent(
        connectivity, input, scope, topology, false);
    for (const Connectivity* member : component) {
        const std::pair<Vec2d, Vec2d> member_entry =
            input.entryPtDir(member->entry_lane_id);
        const std::pair<Vec2d, Vec2d> member_exit =
            input.exitPtDir(member->exit_lane_id);
        const bool shared_entry =
            (member_entry.first - entry.first).norm() < kSharedEndpointEps;
        const bool shared_exit =
            (member_exit.first - exit.first).norm() < kSharedEndpointEps;
        if (!shared_entry && !shared_exit)
            continue;
        const UTurnLeadFloors member_floors = leadFloors(
            member_entry.first, member_entry.second,
            member_exit.first, member_exit.second, input);
        double member_lead0 = member_floors.lead0;
        double member_lead1 = member_floors.lead1;
        const double member_station = alignmentComponentStation(
            *member, input, scope, topology);
        enforceAlignmentStation(
            member_entry.first, member_entry.second,
            member_exit.first, member_exit.second, member_station,
            member_lead0, member_lead1);
        if (shared_entry)
            info.lead0 = std::max(info.lead0, member_lead0);
        if (shared_exit)
            info.lead1 = std::max(info.lead1, member_lead1);
    }

    // 共享端点聚合会把 lead 抬到家族最大值：该最大值来自另一成员的
    // 端点与 U 轴快照，套用到本连接后，实际首/尾平齐点会越过上面记录
    // 的 aligned_station（cross 13 越过 0.054m，仅比站位门禁的 0.05m
    // 容差多出 4mm）。station 必须与最终 lead 同源重算，否则
    // SegmentedUTurnCandidateSearch 的站位门禁会否掉全部合规三段式候选，
    // U-turn 退化成没有首尾直行段的单段曲线。
    // 重算只会让站位沿轴前移（lead 只增不减），所以人行横道净距和
    // 2m 保底仍是下界，家族平齐由共享端点的同一 lead 保证。
    info.aligned_station = requiredAlignedStationOnAxis(
        entry.first, entry.second, exit.first, exit.second, info.axis,
        info.lead0, info.lead1);
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

namespace {

// 家族成员按半径升序（半径相同按 id）排列的稳定次序。
// radiusRank 与 familyStaggerStep 必须共用同一次序，否则分档量与分档上限
// 会依据不同的名次计算，阶梯单调性无从保证。
std::vector<std::pair<double, ConnId> > orderFamilyByRadius(
    const std::vector<const Connectivity*>& component,
    const IntersectionInput& input, const UTurnFamilyBuilder& builder) {
    std::vector<std::pair<double, ConnId> > ordered;
    ordered.reserve(component.size());
    for (const Connectivity* member : component)
        ordered.push_back(std::make_pair(builder.radiusKey(*member, input), member->id));
    std::sort(ordered.begin(), ordered.end(),
              [](const std::pair<double, ConnId>& a,
                 const std::pair<double, ConnId>& b) {
                  if (std::abs(a.first - b.first) > 1e-9)
                      return a.first < b.first;
                  return a.second < b.second;
              });
    return ordered;
}

// 成员在指定一侧的物理端点坐标。共享端点组的判定必须用坐标重合，而不是
// 车道 id 相等：同一物理端点可能被不同 id 的车道引用。
Vec2d endpointOnSide(const Connectivity& connectivity,
                     const IntersectionInput& input, bool entry_side) {
    return entry_side ? input.entryPtDir(connectivity.entry_lane_id).first
                      : input.exitPtDir(connectivity.exit_lane_id).first;
}

}  // namespace

UTurnFamilyRank UTurnFamilyBuilder::radiusRank(
    const Connectivity& connectivity, const IntersectionInput& input,
    UTurnAlignmentScope scope, const ClusterOrderSolver* topology,
    bool require_shared_endpoint_pair) const {
    const std::vector<const Connectivity*> component = alignmentComponent(
        connectivity, input, scope, topology, require_shared_endpoint_pair);
    const std::vector<std::pair<double, ConnId> > ordered =
        orderFamilyByRadius(component, input, *this);
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

UTurnFamilyLadder UTurnFamilyBuilder::familyLateralLadder(
    const Connectivity& connectivity, const IntersectionInput& input,
    UTurnAlignmentScope scope, double nominal_step,
    const ClusterOrderSolver* topology,
    bool require_shared_endpoint_pair) const {
    UTurnFamilyLadder ladder;
    ladder.step = std::max(0.0, nominal_step);
    if (ladder.step <= 0.0)
        return ladder;
    const std::vector<const Connectivity*> component = alignmentComponent(
        connectivity, input, scope, topology, require_shared_endpoint_pair);
    // 与 radiusRank 完全相同的次序：半径升序，半径相同按 id。
    std::vector<std::pair<double, const Connectivity*> > ordered;
    ordered.reserve(component.size());
    for (const Connectivity* member : component)
        ordered.push_back(std::make_pair(radiusKey(*member, input), member));
    std::sort(ordered.begin(), ordered.end(),
              [](const std::pair<double, const Connectivity*>& a,
                 const std::pair<double, const Connectivity*>& b) {
                  if (std::abs(a.first - b.first) > 1e-9)
                      return a.first < b.first;
                  return a.second->id < b.second->id;
              });
    const std::size_t n = ordered.size();
    // 第 1 步：名义阶梯 + 各成员自身走廊上限。0.25 与 buildSegmented 的
    // 0.25 * lateral_gap 同源：相向内缩两侧各占走廊四分之一，中弧保留一半。
    std::vector<double> capped(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        const double nominal =
            ladder.step * static_cast<double>(n - i);
        capped[i] = std::min(nominal, 0.25 * std::max(0.0, ordered[i].first));
        if (ordered[i].second->id == connectivity.id)
            ladder.nominal_stagger = nominal;
    }
    // 第 2 步：按"共享端点组"分侧裁定单调阶梯。
    //
    // 入口侧与出口侧必须分别裁定：家族是一条交替由共享入口、共享出口连接的链，
    // 内缩在两侧沿相反法向施加，用同一个标量同时表达两侧会要求半径次序在两条
    // 不同的链上同时单调，必然在其中一条上反号。
    //
    // 组内以最内层成员（半径最小）的自身走廊上限为预算：同组成员的首/尾直段都
    // 从同一个物理端点出发，只有全部嵌套在最内层成员的走廊里才不会横扫过它。
    // 组内步进量取 min(step, 预算 / 组员数)，因此三条以上共享同一端点时也能分到
    // 互不重合的档位（100000598 的 56/58/60 共享同一入口，退化成员 56 只有
    // 0.09m 走廊，若只按"内层 cap - sep"一次性裁剪，58 与 60 会拿到同一档位）。
    // 递推只在组内进行、不跨组传递：跨组传递会把退化成员的厘米级上限一路压到
    // 家族另一端，实测把 100000412 的 35|37 由不相交压成相交（35 需要约 0.25m
    // 的绝对站位才能避开退化为单段的 37）。
    // 组员数为 1 或组内无退化成员时该步完全不生效，其他数据集形态保持不变。
    std::vector<double> entry_resolved = capped;
    std::vector<double> exit_resolved = capped;
    for (int side = 0; side < 2; ++side) {
        const bool entry_side = side == 0;
        std::vector<double>& resolved = entry_side ? entry_resolved : exit_resolved;
        std::vector<bool> grouped(n, false);
        for (std::size_t i = 0; i < n; ++i) {
            if (grouped[i])
                continue;
            const Vec2d anchor = endpointOnSide(*ordered[i].second, input, entry_side);
            std::vector<std::size_t> group;
            group.push_back(i);
            grouped[i] = true;
            for (std::size_t j = i + 1; j < n; ++j) {
                if (grouped[j])
                    continue;
                const Vec2d other =
                    endpointOnSide(*ordered[j].second, input, entry_side);
                if ((anchor - other).norm() >= 1e-3)
                    continue;
                group.push_back(j);
                grouped[j] = true;
            }
            if (group.size() < 2)
                continue;
            const double budget = std::max(0.0, capped[group.front()]);
            const double sep = std::min(
                ladder.step, budget / static_cast<double>(group.size()));
            double previous = budget;
            resolved[group.front()] = previous;
            for (std::size_t m = 1; m < group.size(); ++m) {
                previous = std::min(capped[group[m]], previous - sep);
                resolved[group[m]] = previous;
            }
        }
    }
    for (std::size_t i = 0; i < n; ++i) {
        if (ordered[i].second->id == connectivity.id) {
            ladder.entry_stagger = std::max(0.0, entry_resolved[i]);
            ladder.exit_stagger = std::max(0.0, exit_resolved[i]);
            ladder.stagger =
                std::max(ladder.entry_stagger, ladder.exit_stagger);
            break;
        }
    }
    return ladder;
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
