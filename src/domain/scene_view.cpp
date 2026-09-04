#include "domain/scene_view.h"
#include "utils.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <unordered_map>

namespace isg {

namespace {

// findLane 是全局最热的按 id 查询：entryPtDir/exitPtDir 经由
// UTurnFamilyBuilder::alignmentComponent、uturnSharedEndpointStagger 等路径
// 在候选评估里被调用数百万次，而原实现每次都要线性扫描 lanes 并逐个比较
// 字符串 id（100000643 实测 entryPtDir+exitPtDir 占 21% 自身耗时）。
//
// 这里缓存一张 id→Lane* 索引。索引只在「输入对象、lanes 缓冲区地址、lanes
// 数量」三者一致时复用，命中后还要校验 Lane::id 与查询 id 相等，任一条件
// 不成立就整表重建，因此返回值与线性扫描逐位一致；重复 id 只登记首次出现，
// 保持原来的「返回第一个匹配」语义。几何数据通过指针读取，即使输入被就地
// 改写也不会返回过期坐标。
struct LaneIndexCache {
    // 二级索引：id 内容 → Lane*。
    const void* owner = nullptr;
    const Lane* data = nullptr;
    std::size_t count = 0;
    std::unordered_map<LaneId, const Lane*> index;

    // 一级缓存：按 id 对象地址直查上次结果。热路径查询的都是 input 内
    // Connectivity/Lane 持有的稳定字符串对象，同一地址会被反复查询，
    // 这样连一次字符串哈希都不用做（实测哈希本身占 100000385-u 约 11%）。
    struct Slot {
        const LaneId* id_ptr = nullptr;
        const Lane* lane = nullptr;
    };
    static const std::size_t kSlots = 64;  // 2 的幂，直接映射
    Slot slots[kSlots];

    Slot& slot(const LaneId* id_ptr) {
        const std::uintptr_t bits = reinterpret_cast<std::uintptr_t>(id_ptr);
        return slots[(bits >> 4) & (kSlots - 1)];
    }

    void rebuild(const void* input, const std::vector<Lane>& lanes) {
        owner = input;
        data = lanes.data();
        count = lanes.size();
        index.clear();
        index.reserve(lanes.size() * 2);
        for (const Lane& item : lanes)
            index.emplace(item.id, &item);
        for (std::size_t i = 0; i < kSlots; ++i)
            slots[i] = Slot();
    }
};

}  // namespace

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
    static thread_local LaneIndexCache cache;
    if (cache.owner != this || cache.data != lanes.data() ||
        cache.count != lanes.size())
        cache.rebuild(this, lanes);
    // 地址命中后仍比较字符串内容：临时字符串复用同一栈地址时校验会失败，
    // 自动退回二级索引，因此结果与线性扫描一致。
    LaneIndexCache::Slot& slot = cache.slot(&id);
    if (slot.id_ptr == &id && slot.lane && slot.lane->id == id)
        return slot.lane;
    std::unordered_map<LaneId, const Lane*>::const_iterator it =
        cache.index.find(id);
    if (it != cache.index.end()) {
        if (it->second->id == id) {
            slot.id_ptr = &id;
            slot.lane = it->second;
            return it->second;
        }
        // 索引指向的 Lane 已换成别的 id：缓冲区被原地改写且规模未变。
        cache.rebuild(this, lanes);
        it = cache.index.find(id);
        return it == cache.index.end() ? nullptr : it->second;
    }
    // 未命中可能是真的不存在，也可能是索引过期。线性回退一次兜底：
    // 该分支只在缺失 id（原实现会打 WARN）时触发，代价与原实现相同。
    for (const auto& item : lanes) {
        if (item.id == id) {
            cache.rebuild(this, lanes);
            return &item;
        }
    }
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
