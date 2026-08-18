#pragma once

#include "types.h"

#include <unordered_map>

namespace isg {

/// 一个规范化路口输入的不可变索引视图。
class SceneView {
public:
    explicit SceneView(const IntersectionInput& input);

    const IntersectionInput& input() const { return *input_; }
    const Lane* lane(const LaneId& id) const;
    const LaneGroup* group(const LaneGroupId& id) const;
    const LaneEdge* edge(const LaneEdgeId& id) const;
    const Connectivity* connectivity(const ConnId& id) const;

    std::pair<Vec2d, Vec2d> entryFrame(const LaneId& id) const;
    std::pair<Vec2d, Vec2d> exitFrame(const LaneId& id) const;

private:
    const IntersectionInput* input_;
    std::unordered_map<LaneId, const Lane*> lanes_;
    std::unordered_map<LaneGroupId, const LaneGroup*> groups_;
    std::unordered_map<LaneEdgeId, const LaneEdge*> edges_;
    std::unordered_map<ConnId, const Connectivity*> connectivities_;
};

}  // 命名空间 isg
