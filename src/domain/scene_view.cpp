#include "domain/scene_view.h"
#include "utils.h"

#include <algorithm>
#include <iostream>

namespace isg {

const bool IntersectionInput::IsEntryLane(const LaneId& id) const {
    for (const auto& group : lane_groups)
        if (std::find(group.lanes.begin(), group.lanes.end(), id) != group.lanes.end())
            return group.role == GroupRole::Entry;
    return false;
}

const bool IntersectionInput::IsEntryLaneEdge(const LaneEdgeId& id) const {
    for (const auto& group : lane_groups)
        if (std::find(group.boundaries.begin(), group.boundaries.end(), id) != group.boundaries.end())
            return group.role == GroupRole::Entry;
    for (const auto& group : lane_groups) {
        for (const auto& lane_id : group.lanes) {
            const Lane* item = findLane(lane_id);
            if (item && (item->left_edge_id == id || item->right_edge_id == id))
                return group.role == GroupRole::Entry;
        }
    }
    return false;
}

const Lane* IntersectionInput::findLane(const LaneId& id) const {
    for (const auto& item : lanes) if (item.id == id) return &item;
    return nullptr;
}

const LaneGroup* IntersectionInput::findGroup(const LaneGroupId& id) const {
    for (const auto& item : lane_groups) if (item.id == id) return &item;
    return nullptr;
}

const LaneEdge* IntersectionInput::findEdge(const LaneEdgeId& id) const {
    for (const auto& item : lane_edges) if (item.id == id) return &item;
    return nullptr;
}

bool IntersectionInput::laneGroupExists(const LaneGroupId& id) const {
    return findGroup(id) != nullptr;
}

std::pair<Vec2d, Vec2d> IntersectionInput::entryPtDir(const LaneId& id) const {
    const Lane* item = findLane(id);
    if (!item || item->geometry.points.empty()) {
        std::cout << "[WARN] entrylane:" << id << " no geometry!\n";
        return std::make_pair(Vec2d(0, 0), Vec2d(1, 0));
    }
    return std::make_pair(entryLinePoint(item->geometry.points),
                          entryLineTangent(item->geometry.points));
}

std::pair<Vec2d, Vec2d> IntersectionInput::exitPtDir(const LaneId& id) const {
    const Lane* item = findLane(id);
    if (!item || item->geometry.points.empty()) {
        std::cout << "[WARN] exitlane:" << id << " no geometry!\n";
        return std::make_pair(Vec2d(10, 0), Vec2d(1, 0));
    }
    return std::make_pair(exitLinePoint(item->geometry.points),
                          exitLineTangent(item->geometry.points));
}

SceneView::SceneView(const IntersectionInput& input) : input_(&input) {
    for (const auto& item : input.lanes)
        if (!lanes_.count(item.id)) lanes_[item.id] = &item;
    for (const auto& item : input.lane_groups) groups_[item.id] = &item;
    for (const auto& item : input.lane_edges) edges_[item.id] = &item;
    for (const auto& item : input.connectivities) connectivities_[item.id] = &item;
}

const Lane* SceneView::lane(const LaneId& id) const {
    std::unordered_map<LaneId, const Lane*>::const_iterator it = lanes_.find(id);
    return it == lanes_.end() ? nullptr : it->second;
}

const LaneGroup* SceneView::group(const LaneGroupId& id) const {
    std::unordered_map<LaneGroupId, const LaneGroup*>::const_iterator it = groups_.find(id);
    return it == groups_.end() ? nullptr : it->second;
}

const LaneEdge* SceneView::edge(const LaneEdgeId& id) const {
    std::unordered_map<LaneEdgeId, const LaneEdge*>::const_iterator it = edges_.find(id);
    return it == edges_.end() ? nullptr : it->second;
}

const Connectivity* SceneView::connectivity(const ConnId& id) const {
    std::unordered_map<ConnId, const Connectivity*>::const_iterator it = connectivities_.find(id);
    return it == connectivities_.end() ? nullptr : it->second;
}

std::pair<Vec2d, Vec2d> SceneView::entryFrame(const LaneId& id) const {
    const Lane* item = lane(id);
    if (!item || item->geometry.points.empty()) {
        std::cout << "[WARN] entrylane:" << id << " no geometry!\n";
        return std::make_pair(Vec2d(0, 0), Vec2d(1, 0));
    }
    return std::make_pair(entryLinePoint(item->geometry.points),
                          entryLineTangent(item->geometry.points));
}

std::pair<Vec2d, Vec2d> SceneView::exitFrame(const LaneId& id) const {
    const Lane* item = lane(id);
    if (!item || item->geometry.points.empty()) {
        std::cout << "[WARN] exitlane:" << id << " no geometry!\n";
        return std::make_pair(Vec2d(10, 0), Vec2d(1, 0));
    }
    return std::make_pair(exitLinePoint(item->geometry.points),
                          exitLineTangent(item->geometry.points));
}

}  // 命名空间 isg
