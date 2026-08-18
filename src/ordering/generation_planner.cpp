#include "ordering/generation_planner.h"
#include "preprocessing/uturn_family_builder.h"
#include "types.h"
#include "toolkits/toolkits.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <unordered_map>

namespace isg {
namespace {

int priorityOf(ConnTurnType turn) {
    return turn == ConnTurnType::Straight ? 0 :
        ((turn == ConnTurnType::UTurnLeft || turn == ConnTurnType::UTurnRight) ? 2 : 1);
}

}  // namespace

GenerationPlan GenerationPlanner::build(const IntersectionInput& input) const {
    ClusterOrderSolver topology;
    topology.build(input.connectivities, input.lanes, input.lane_groups, input.crosswalks);
    return build(input, topology);
}

GenerationPlan GenerationPlanner::build(
        const IntersectionInput& input, const ClusterOrderSolver& topology) const {
    std::unordered_map<ConnId, int> spatial_index;
    std::unordered_map<ConnId, int> group_size;
    for (const auto& item : topology.entryGroupOrder()) {
        const int size = static_cast<int>(item.second.size());
        for (int i = 0; i < size; ++i) {
            spatial_index[item.second[i]] = i;
            group_size[item.second[i]] = size;
        }
    }

    std::map<int, std::vector<const Connectivity*> > grouped;
    for (const auto& conn : input.connectivities)
        grouped[priorityOf(TurnPreprocessor(conn, input))].push_back(&conn);

    GenerationPlan plan;
    for (auto& item : grouped) {
        const int priority = item.first;
        std::vector<const Connectivity*>& connections = item.second;
        std::stable_sort(connections.begin(), connections.end(),
            [&](const Connectivity* a, const Connectivity* b) {
                const int ia = spatial_index.count(a->id) ? spatial_index[a->id] : 0;
                const int ib = spatial_index.count(b->id) ? spatial_index[b->id] : 0;
                const int sa = group_size.count(a->id) ? group_size[a->id] : 1;
                const int sb = group_size.count(b->id) ? group_size[b->id] : 1;
                const double da = std::abs(ia - (sa - 1) / 2.0);
                const double db = std::abs(ib - (sb - 1) / 2.0);
                if (std::abs(da - db) > 0.5) return da < db;
                if (a->enterGroupId != b->enterGroupId)
                    return a->enterGroupId < b->enterGroupId;
                return ia < ib;
            });

        if (priority == 2) {
            std::unordered_map<LaneId, std::vector<size_t> > family_slots;
            for (size_t i = 0; i < connections.size(); ++i)
                family_slots[connections[i]->entry_lane_id].push_back(i);
            for (const auto& family : family_slots) {
                if (family.second.size() < 2) continue;
                std::vector<const Connectivity*> ordered;
                for (size_t slot : family.second) ordered.push_back(connections[slot]);
                std::stable_sort(ordered.begin(), ordered.end(),
                    [&](const Connectivity* a, const Connectivity* b) {
                        const UTurnFamilyBuilder family_builder;
                        const double ra = family_builder.radiusKey(*a, input);
                        const double rb = family_builder.radiusKey(*b, input);
                        if (std::abs(ra - rb) > 1e-9) return ra < rb;
                        return a->id < b->id;
                    });
                for (size_t i = 0; i < family.second.size(); ++i)
                    connections[family.second[i]] = ordered[i];
            }
        }

        GenerationBatch batch;
        batch.priority = priority;
        for (const Connectivity* conn : connections) batch.conn_ids.push_back(conn->id);
        plan.batches.push_back(batch);
    }
    return plan;
}

GenerationPlan GenerationPlanner::build(const SceneContext& context,
                                        const ClusterOrderSolver& topology) const {
    return build(context.view.input(), topology);
}

}  // 命名空间 isg
