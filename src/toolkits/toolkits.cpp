//
// Created by yuanjinfa on 2026/8/14.
//

#include "toolkits.h"
#include "preprocessing/uturn_family_builder.h"
#include "utils.h"
#include <algorithm>
#include <limits>
#include <unordered_set>

#include "nonmotor/nonmotor_utils.h"

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
        if (isGeometricUTurnByTurnType(
                connectivity.turn_type, entry_tangent, exit_tangent))
            {
            return cross2d(entry_tangent, exit.first - entry.first) >= 0.0 ?
                   ConnTurnType::UTurnLeft : ConnTurnType::UTurnRight;
            }
        if (connectivity.turn_type == ConnTurnType::TurnLeft ||
            connectivity.turn_type == ConnTurnType::TurnRight ||
            connectivity.turn_type == ConnTurnType::UTurnLeft ||
            connectivity.turn_type == ConnTurnType::UTurnRight)
            return connectivity.turn_type;
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

        // 抗抖动的"车道自身朝向"：自路口端点沿折线回溯至累计长度达到 baseline，取该点到端点
        // 的弦。端点两点式切向只看最后两点，末段短到几十厘米时（corpus 里最短 0.31m）方向被
        // 数字化噪声主导，不足以据此判定"这条车道属于另一条臂"。baseline <= 0 或折线长度不足
        // 时退回整条车道的弦。注意本函数只在插入合成点之前调用（normalizer 每次运行只处理
        // 每条车道一次，见下方 applied 去重）。
        auto robustDirectionForRole = [&](const Lane& lane, GroupRole role,
                                          double baseline) -> Vec2d {
            const auto& points = lane.geometry.points;
            if (baseline <= 0.0 || points.size() < 2)
                return Vec2d(0, 0);
            const bool from_back = role == GroupRole::Entry;
            const Vec2d endpoint = from_back ? Vec2d(points.back()) : Vec2d(points.front());
            Vec2d anchor = endpoint;
            double accumulated = 0.0;
            for (std::size_t step = 1; step < points.size(); ++step) {
                const std::size_t i = from_back ? points.size() - 1 - step : step;
                const Vec2d current(points[i]);
                const std::size_t prev = from_back ? i + 1 : i - 1;
                accumulated += (current - Vec2d(points[prev])).norm();
                anchor = current;
                if (accumulated >= baseline)
                    break;
            }
            const Vec2d chord = from_back ? endpoint - anchor : anchor - endpoint;
            return chord.norm() > 1e-10 ? chord.normalized() : Vec2d(0, 0);
        };

        auto directionAngleDiff = [&](const Vec2d& a, const Vec2d& b) -> double {
            if (a.norm() < 1e-10 || b.norm() < 1e-10)
                return M_PI;
            double cosine = a.normalized().dot(b.normalized());
            cosine = std::max(-1.0, std::min(1.0, cosine));
            return std::acos(cosine);
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
            std::unordered_set<LaneId> applied;
            for (const auto& lane_id : group.lanes) {
                if (!applied.insert(lane_id).second)
                    continue;  // group.lanes 存在重复 id，重复覆盖只会插入冗余点
                Lane* lane = findLane(lanes, lane_id);
                if (!lane || hasSingleConnectivity(*lane))
                    continue;
                // 跨臂误并保护：统一切向的用途是抹掉**同一路口臂内平行车道**的数字化抖动，
                // 不是把车道重新指向。若某条车道自身朝向与组方向差过大，说明这个
                // groupId 里混进了物理方向不同的另一条臂（或本车道是转向专用道），此时
                // 覆盖会让所有以它为端点的曲线在错误的切向上做 G1——曲线从连接点起就偏出
                // 车道数十度，随后必然自交、越路缘、与同簇兄弟成片相交，而审计读的是同一份
                // 规范化输入，反而看不出端点 G1 异常。这类车道保留自身朝向。
                //
                // 判据取两点式切向与长基线弦**二者与组方向夹角的较小值**，只有两个估计都
                // 认定"方向不同"才放弃统一。原因是两者各有失效模式：末段仅 0.31m 的车道
                // （`110000703-u` 的 `43106568`）两点式切向被噪声主导，虚报 23.8° 偏差；
                // 而弯曲的转向专用道整条弦又会偏离端点切向。取较小值使保护只在"无歧义的
                // 跨臂"上触发，其余一律沿用统一，避免把真实抖动误判成跨臂。
                //
                // corpus 定标（`diag_group_dir`，575 条车道）：被强制转过的角度
                // `<0.5° 420 条、<1° 55、<2° 59、<5° 19`，即 96.2% 在 5° 以内，属真实抖动；
                // 5°~20° 共 14 条；20° 以上 8 条，全部经查为"一个 groupId 装了两条臂"或
                // 转向专用道（`110004764` 的 `43107602` 把朝 -11.5° 的五条车道与朝
                // -83°/-88° 的两条并到一起，后者被整整扭转 71.8°/76.6°）。阈值取
                // `kGroupDirectionForceLimitDeg = 20°`：既是真实抖动上限的 4 倍，也低于
                // 实测最小的跨臂偏差；换算成几何量，3.5m 长的车道末段偏 20° 才产生 1.2m
                // 横向漂移，不可能来自数字化误差。
                //
                // 刻意不做"把离群车道再按方向聚类后组内统一"：`110000741-u` 的
                // `43100884` 只有两条车道且相差 34.7°，任何再统一都会把两条各扭 17°，
                // 比保留自身朝向更糟。
                const Vec2d own = directionForRole(*lane, group.role);
                const Vec2d robust = robustDirectionForRole(
                        *lane, group.role, config.group_robust_baseline_m);
                double deviation = directionAngleDiff(own, direction->second);
                if (robust.norm() > 1e-10)
                    deviation = std::min(
                            deviation, directionAngleDiff(robust, direction->second));
                if (deviation > config.group_force_limit_deg * M_PI / 180.0)
                    continue;
                if (group.role == GroupRole::Entry)
                    setEntryDirection(*lane, direction->second);
                else
                    setExitDirection(*lane, direction->second);
            }
        }
    };
    // ─────────────────────────────────────────────────────────────────────────
    // ConnectivityMotorLogicFilter
    //
    // 将输入的连通关系按"机动挂接逻辑"和"非机动挂接逻辑"精确分流。
    //
    // is_motorlogic=true  → 保留按机动逻辑处理的连通关系（机动车道 + 符合条件的非机动）
    // is_motorlogic=false → 保留按非机动逻辑处理的连通关系（其余非机动）
    //
    // 分流判据（GenNonmotorClassStategys）：
    //   机动逻辑：
    //     A. 所有机动车道（lane_type==MotorLane）
    //     B. 右转非机动车道（追加 0.9m 净距约束在生成时处理）
    //     C. 直行非机动：
    //        C1. 进或出只有一个是非机动车道
    //        C2. 右侧进入组≥2 OR 退出组≥2（右侧 = 进出组最右侧连线右侧，Q1）
    //        C3. 无安全岛
    //        C4. 安全岛距进入端点 > island_max_dist（Q4，默认20m）
    //        C5. 安全岛不可穿（穿行通道内点到RoadEdge < 0.9m）
    //     D. 左转非机动：
    //        D1. 进入组仅2车道（1机+1非机）且路口机动车道往返总数≤2（特殊情况F）
    //        D2. 进入机动，退出非机动，非提右车道（G）
    //        D3. 其他左转（进入非机动但安全岛不可穿时）
    //     E. 掉头非机动：路口无安全岛 且 进入机动退出非机动（H）
    //
    //   非机动逻辑（其余非机动连通，交由 NonMotorConnectivityGenerator 处理）
    // ─────────────────────────────────────────────────────────────────────────
    IntersectionInput ConnectivityMotorLogicFilter(const IntersectionInput& raw_input, bool is_motorlogic) {
        IntersectionInput input = raw_input;

        // 安全岛相关工具（用于直行/左转/掉头分流判断）
        static constexpr double kIslandMaxDist = 20.0;
        // TODO: 如何判断安全岛是否可穿（进入方向）
        auto islandPassable = [&](
            const SafetyIsland* si, const Vec2d& entry_dir, double clearance) -> bool {
            if (!si) return false;
            return true;
        };

        // 是否按机动车道处理逻辑
        auto isMotorProcess = [&](const Connectivity& conn) -> bool {
            // A. 所有机动车道
            if (conn.lane_type == MotorLaneType::MotorLane) return true;

            // 以下均为 NonMotorLane
            bool entry_is_nonmotor = isNonMotorLane(conn.entry_lane_id, input);
            bool exit_is_nonmotor  = isNonMotorLane(conn.exit_lane_id,  input);
            Vec2d entry_pt  = input.entryPtDir(conn.entry_lane_id).first;
            Vec2d entry_dir = input.entryPtDir(conn.entry_lane_id).second;

            // B. 右转非机动：全部按机动逻辑处理（追加 0.9m 净距约束在生成时注入）
            if (conn.turn_type == ConnTurnType::TurnRight) return true;

            // ── 直行非机动 ──────────────────────────────────────────────────
            if (conn.turn_type == ConnTurnType::Straight) {
                // C1. 进或出只有一个是非机动车道
                if (!entry_is_nonmotor || !exit_is_nonmotor) return true;
                // C2. 右侧进入组≥2 OR 退出组≥2
                if (countRightSideGroups(conn, GroupRole::Entry, input) >= 2) return true;
                if (countRightSideGroups(conn, GroupRole::Exit,  input) >= 2) return true;
                // C3. 无安全岛
                const SafetyIsland* si = findNearestSafetyIsland(entry_pt, input.safety_islands);
                if (!si) return true;
                // C4. 安全岛距进入端点 > 20m
                if (distIslandBoundaryToPoint(*si, entry_pt) > kIslandMaxDist) return true;
                // C5. 安全岛不可穿
                if (!islandPassable(si, entry_dir, 0.9)) return true;
                // → 安全岛可穿，非机动逻辑
                return false;
            }

            // ── 左转非机动 ──────────────────────────────────────────────────
            if (conn.turn_type == ConnTurnType::TurnLeft) {
                bool yield_right = isYieldRightLane(conn, input);
                const LaneGroup* eg = input.findGroup(conn.enterGroupId);

                // D1. 进入组仅2车道（1机+1非机）且路口机动车往返总数≤2 → 机动（特殊情况F）
                if (eg && eg->lanes.size() == 2 &&
                    countNonMotorLanesInGroup(conn.enterGroupId, input) == 1 &&
                    countMotorLanesTotalIntersection(input) <= 2)
                    return true;

                // D2. 进入机动，退出非机动，非提右车道 → 机动（G）
                if (!entry_is_nonmotor && exit_is_nonmotor && !yield_right)
                    return true;

                // 进入非机动，非提右车道 → J. 非机动L型
                if (entry_is_nonmotor && !yield_right) return false;

                // 进入非机动，提右车道 → K（安全岛可穿）or 机动（不可穿）
                if (entry_is_nonmotor && yield_right) {
                    const SafetyIsland* si = findNearestSafetyIsland(entry_pt, input.safety_islands);
                    if (si && distIslandBoundaryToPoint(*si, entry_pt) <= kIslandMaxDist
                           && islandPassable(si, entry_dir, 0.9))
                        return false;  // K. 非机动穿岛
                    return true;       // 无岛/不可穿 → 机动
                }

                // 进入机动，退出非机动，提右车道 → L（安全岛可穿）or 机动（不可穿）
                if (!entry_is_nonmotor && exit_is_nonmotor && yield_right) {
                    const SafetyIsland* si = findNearestSafetyIsland(entry_pt, input.safety_islands);
                    if (si && distIslandBoundaryToPoint(*si, entry_pt) <= kIslandMaxDist
                           && islandPassable(si, entry_dir, 1.0))
                        return false;  // L. 非机动处理
                    return true;       // 机动
                }

                return true;  // 其他左转 → 机动
            }

            // ── 掉头非机动 ──────────────────────────────────────────────────
            if (conn.turn_type == ConnTurnType::UTurnLeft ||
                conn.turn_type == ConnTurnType::UTurnRight) {
                // H. 进入机动，退出非机动，路口无安全岛 → 机动
                if (!entry_is_nonmotor && exit_is_nonmotor && input.safety_islands.empty())
                    return true;

                // 进入机动，退出非机动，有安全岛：可穿 → 非机动（M），不可穿 → 机动
                if (!entry_is_nonmotor && exit_is_nonmotor) {
                    const SafetyIsland* si = findNearestSafetyIsland(entry_pt, input.safety_islands);
                    if (si && distIslandBoundaryToPoint(*si, entry_pt) <= kIslandMaxDist
                           && islandPassable(si, entry_dir, 1.0))
                        return false;  // M. 非机动
                    return true;       // 机动
                }

                // 进入非机动 → M. 非机动U型
                if (entry_is_nonmotor) return false;

                return true;
            }

            // 未知转向类型 → 默认机动
            return true;
        };

        // 过滤 connectivities
        auto& conns_ref = input.connectivities;
        conns_ref.erase(
            std::remove_if(conns_ref.begin(), conns_ref.end(),
                [&](const Connectivity& c) {
                    bool motor = isMotorProcess(c);
                    return is_motorlogic ? !motor : motor;
                }),
            conns_ref.end());

        return input;
    };
};
