//
// Created by yuanjinfa on 2026/8/14.
//

#include "toolkits.h"
#include "utils.h"
#include <unordered_set>

namespace isg {

    void ElevationInterpolator(std::vector<ConnectivityCurve> &curves, const IntersectionInput &input) {
        for (auto &curve : curves) {
            if (curve.fixed_shape || curve.geometry.points.empty()) continue;
            const Lane *entry = input.findLane(curve.entry_lane_id);
            const Lane *exit = input.findLane(curve.exit_lane_id);
            if (!entry || !exit || entry->geometry.points.empty() || exit->geometry.points.empty()) continue;
            const double z0 = entry->geometry.points.back().z();
            const double z1 = exit->geometry.points.front().z();
            std::vector<Vec2d> xy = toVec2dArray(curve.geometry.points);
            double total = 0.0;
            for (size_t i = 1; i < xy.size(); ++i) total += dist(xy[i - 1], xy[i]);
            double travelled = 0.0;
            for (size_t i = 0; i < curve.geometry.points.size(); ++i) {
                if (i > 0) travelled += dist(xy[i - 1], xy[i]);
                const double t = total > EPS ? travelled / total :
                                 (xy.size() > 1 ? static_cast<double>(i) / (xy.size() - 1) : 0.0);
                curve.geometry.points[i].z() = z0 * (1.0 - t) + z1 * t;
            }
        }
    }

    ConnTurnType TurnPreprocessor(
            const Connectivity &connectivity, const IntersectionInput &input) {
        const std::pair<Vec2d, Vec2d> entry = input.entryPtDir(connectivity.entry_lane_id);
        const std::pair<Vec2d, Vec2d> exit = input.exitPtDir(connectivity.exit_lane_id);
        if (entry.second.norm() < 1e-9 || exit.second.norm() < 1e-9)
            return connectivity.turn_type;
        const Vec2d entry_tangent = entry.second.normalized();
        const Vec2d exit_tangent = exit.second.normalized();
        if (entry_tangent.dot(exit_tangent) < -0.5)
            return cross2d(entry_tangent, exit.first - entry.first) >= 0.0 ?
                   ConnTurnType::UTurnLeft : ConnTurnType::UTurnRight;
        Vec2d chord = exit.first - entry.first;
        if (chord.norm() < 1e-9) return ConnTurnType::Straight;
        chord.normalize();
        const double cross = cross2d(entry_tangent, chord);
        if (cross > 0.25) return ConnTurnType::TurnLeft;
        if (cross < -0.25) return ConnTurnType::TurnRight;
        return ConnTurnType::Straight;
    }

    IntersectionInput InputNormalizer(const IntersectionInput& source) {
        IntersectionInput input = source;
        std::unordered_map<LaneId, LaneGroupId> entry_groups;
        std::unordered_map<LaneId, LaneGroupId> exit_groups;
        for (const auto& group : input.lane_groups) {
            for (const auto& lane_id : group.lanes) {
                (group.role == GroupRole::Entry ? entry_groups : exit_groups)[lane_id] = group.id;
            }
        }
        for (auto& lane : input.lanes) {
            if (lane.groupId.empty()) {
                std::unordered_map<LaneId, LaneGroupId>::const_iterator e = entry_groups.find(lane.id);
                std::unordered_map<LaneId, LaneGroupId>::const_iterator x = exit_groups.find(lane.id);
                if (e != entry_groups.end()) lane.groupId = e->second;
                else if (x != exit_groups.end()) lane.groupId = x->second;
            }
            if (!lane.groupId.empty()) {
                if (!entry_groups.count(lane.id)) entry_groups[lane.id] = lane.groupId;
                if (!exit_groups.count(lane.id)) exit_groups[lane.id] = lane.groupId;
            }
        }
        for (auto& conn : input.connectivities) {
            if (conn.enterGroupId.empty()) {
                std::unordered_map<LaneId, LaneGroupId>::const_iterator it = entry_groups.find(conn.entry_lane_id);
                if (it != entry_groups.end()) conn.enterGroupId = it->second;
            }
            if (conn.exitGroupId.empty()) {
                std::unordered_map<LaneId, LaneGroupId>::const_iterator it = exit_groups.find(conn.exit_lane_id);
                if (it != exit_groups.end()) conn.exitGroupId = it->second;
            }
        }
        return input;
    }


    using LaneIndex = std::unordered_map<LaneId, Lane*>;
    struct LaneDirectionSample {
        LaneId lane_id;
        Vec2d direction{1, 0};
        Vec2d endpoint{0, 0};
        int lane_order = 0;
        int group_index = 0;
    };
    void ConnectivityDirectionNormalizer(
            IntersectionInput& input, const ConnectivityDirectionConfig& config) {

        auto hasSingleConnectivity = [&](const Lane& lane) -> bool {
            const auto it = lane.attrs.find("CON_NUM");
            if (it == lane.attrs.end())
                return false;
            try {
                return std::abs(std::stod(it->second) - 1.0) < 1e-9;
            } catch (...) {
                return it->second == "1";
            }
        };

        auto directionForRole = [&](const Lane& lane, GroupRole role) -> Vec2d {
            Vec2d direction = role == GroupRole::Entry
                              ? entryLineTangent(lane.geometry.points) : exitLineTangent(lane.geometry.points);
            return direction.norm() > 1e-10 ? direction.normalized() : Vec2d(1, 0);
        };

        auto findLane = [&](const LaneIndex& lanes, const LaneId& id) -> Lane* {
            const auto it = lanes.find(id);
            return it == lanes.end() ? nullptr : it->second;
        };

        auto setEntryDirection = [&](Lane& lane, const Vec2d& direction) -> void {
            if (direction.norm() < 1e-10 || lane.geometry.points.empty())
                return;
            const Vec2d normalized = direction.normalized();
            auto& points = lane.geometry.points;
            const Vec2d endpoint = points.back();
            const double z = points.back().z();
            if (points.size() == 1)
                points.insert(points.begin(), Vec3d(endpoint - normalized, z));
            else
                points.insert(points.end() - 1, Vec3d(endpoint - normalized, z));
        };

        auto setExitDirection = [&](Lane& lane, const Vec2d& direction) -> void {
            if (direction.norm() < 1e-10 || lane.geometry.points.empty())
                return;
            const Vec2d normalized = direction.normalized();
            auto& points = lane.geometry.points;
            const Vec2d endpoint = points.front();
            const double z = points.front().z();
            if (points.size() == 1)
                points.push_back(Vec3d(endpoint + normalized, z));
            else
                points.insert(points.begin() + 1, Vec3d(endpoint + normalized, z));
        };

        auto unifiedDirection = [&](
                const LaneGroup& group, const LaneIndex& lanes,
                const std::vector<Connectivity>& connectivities) -> Vec2d {

            auto directionAngleDiff = [&](const Vec2d& a, const Vec2d& b) -> double {
                if (a.norm() < 1e-10 || b.norm() < 1e-10)
                    return M_PI;
                double cosine = a.normalized().dot(b.normalized());
                cosine = std::max(-1.0, std::min(1.0, cosine));
                return std::acos(cosine);
            };

            auto collectDirections = [&](
                    const LaneGroup& group, const LaneIndex& lanes) -> std::vector<LaneDirectionSample> {

                auto endpointForRole = [&](const Lane& lane, GroupRole role) -> Vec2d {
                    if (lane.geometry.points.empty())
                        return Vec2d(0, 0);
                    return role == GroupRole::Entry ?
                           entryLinePoint(lane.geometry.points) : exitLinePoint(lane.geometry.points);
                };

                std::vector<LaneDirectionSample> samples;
                samples.reserve(group.lanes.size());
                std::unordered_set<LaneId> seen;
                for (int i = 0; i < static_cast<int>(group.lanes.size()); ++i) {
                    if (!seen.insert(group.lanes[i]).second)
                        continue;
                    const Lane* lane = findLane(lanes, group.lanes[i]);
                    if (!lane || lane->geometry.points.empty())
                        continue;
                    LaneDirectionSample sample;
                    sample.lane_id = lane->id;
                    sample.direction = directionForRole(*lane, group.role);
                    sample.endpoint = endpointForRole(*lane, group.role);
                    sample.lane_order = lane->laneOrder;
                    sample.group_index = i;
                    samples.push_back(sample);
                }
                return samples;
            };

            auto leftmostDirection = [&](
                    const std::vector<LaneDirectionSample>& samples) -> const LaneDirectionSample*{
                if (samples.empty())
                    return nullptr;
                Vec2d sum(0, 0);
                Vec2d reference(0, 0);
                for (const auto& sample : samples) {
                    if (sample.direction.norm() > 1e-10)
                        sum += sample.direction.normalized();
                    reference += sample.endpoint;
                }
                const LaneDirectionSample* fallback = (samples.empty()) ? nullptr
                        : &*std::min_element(samples.begin(), samples.end(),
                                             [](const LaneDirectionSample& a, const LaneDirectionSample& b) {
                                                 if (a.lane_order != b.lane_order)
                                                     return a.lane_order < b.lane_order;
                                                 return a.group_index < b.group_index;
                                             });
                Vec2d direction = sum.norm() > 1e-10
                                  ? sum.normalized()
                                  : (fallback ? fallback->direction : Vec2d(1, 0));
                const Vec2d left(-direction.y(), direction.x());
                reference /= static_cast<double>(samples.size());
                return &*std::max_element(
                        samples.begin(), samples.end(),
                        [&](const LaneDirectionSample& a, const LaneDirectionSample& b) {
                            const double lateral_a = (a.endpoint - reference).dot(left);
                            const double lateral_b = (b.endpoint - reference).dot(left);
                            if (std::abs(lateral_a - lateral_b) > 1e-6)
                                return lateral_a < lateral_b;
                            if (a.lane_order != b.lane_order)
                                return a.lane_order > b.lane_order;
                            return a.group_index > b.group_index;
                        });
            };

            const auto samples = collectDirections(group, lanes);
            const auto* leftmost = leftmostDirection(samples);
            const Vec2d fallback = leftmost ? leftmost->direction : Vec2d(1, 0);

            std::unordered_map<LaneId, int> sample_by_lane;
            for (int i = 0; i < static_cast<int>(samples.size()); ++i)
                sample_by_lane[samples[i].lane_id] = i;

            int best_index = -1;
            double best_angle = std::numeric_limits<double>::infinity();
            double best_lateral = -std::numeric_limits<double>::infinity();
            Vec2d left(0, 1);
            Vec2d reference(0, 0);
            if (!samples.empty()) {
                Vec2d base_direction(0, 0);
                for (const auto& sample : samples) {
                    base_direction += sample.direction;
                    reference += sample.endpoint;
                }
                if (base_direction.norm() < 1e-10)
                    base_direction = fallback;
                base_direction.normalize();
                left = Vec2d(-base_direction.y(), base_direction.x());
                reference /= static_cast<double>(samples.size());
            }

            for (const auto& connectivity : connectivities) {
                if (connectivity.turn_type != ConnTurnType::Straight)
                    continue;
                const LaneId& lane_id = group.role == GroupRole::Entry
                                        ? connectivity.entry_lane_id : connectivity.exit_lane_id;
                const auto sample_it = sample_by_lane.find(lane_id);
                if (sample_it == sample_by_lane.end())
                    continue;
                const Lane* entry_lane = findLane(lanes, connectivity.entry_lane_id);
                const Lane* exit_lane = findLane(lanes, connectivity.exit_lane_id);
                if (!entry_lane || !exit_lane)
                    continue;
                const double angle = directionAngleDiff(
                        directionForRole(*entry_lane, GroupRole::Entry),
                        directionForRole(*exit_lane, GroupRole::Exit));
                const LaneDirectionSample& sample = samples[sample_it->second];
                const double lateral = (sample.endpoint - reference).dot(left);
                if (angle < best_angle - 1e-9 ||
                    (std::abs(angle - best_angle) <= 1e-9 &&
                     (lateral > best_lateral + 1e-6 ||
                      (std::abs(lateral - best_lateral) <= 1e-6 &&
                       (best_index < 0 ||
                        sample.group_index < samples[best_index].group_index))))) {
                    best_index = sample_it->second;
                    best_angle = angle;
                    best_lateral = lateral;
                }
            }
            return best_index >= 0 ? samples[best_index].direction : fallback;
        };

        if (config.mode != ConnectivityDirectionMode::GroupUnified)
            return;

        LaneIndex lanes;
        for (auto& lane : input.lanes)
            lanes[lane.id] = &lane;
        std::unordered_map<LaneGroupId, Vec2d> group_directions;
        for (const auto& group : input.lane_groups)
            group_directions[group.id] = unifiedDirection(group, lanes, input.connectivities);

        for (const auto& group : input.lane_groups) {
            const auto direction = group_directions.find(group.id);
            if (direction == group_directions.end())
                continue;
            for (const auto& lane_id : group.lanes) {
                Lane* lane = findLane(lanes, lane_id);
                if (!lane || hasSingleConnectivity(*lane))
                    continue;
                if (group.role == GroupRole::Entry)
                    setEntryDirection(*lane, direction->second);
                else
                    setExitDirection(*lane, direction->second);
            }
        }
    }
};
