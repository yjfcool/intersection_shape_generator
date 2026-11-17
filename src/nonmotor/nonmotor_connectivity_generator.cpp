#include "nonmotor/nonmotor_connectivity_generator.h"
#include "nonmotor/nonmotor_utils.h"
#include "nonmotor/nonmotor_corner_rounding.h"
#include "nonmotor/nonmotor_surface_constraints.h"
#include "utils.h"
#include "curve/curve_utils.h"
#include "constraints/fence_check.h"
#include "utils/clipper.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <map>
#include <queue>
#include <set>
#include <tuple>

namespace isg {
namespace {

struct Attachment {
    Vec3d point;
    // ring_segment is the synthetic segment index in ring_geometry for ring hits.
    size_t ring_segment = std::numeric_limits<size_t>::max();
    // Branch hits retain their source edge and local geometry segment so the
    // entry/exit connector can follow the actual branch back to the ring.
    size_t backbone_edge = std::numeric_limits<size_t>::max();
    size_t geometry_segment = std::numeric_limits<size_t>::max();
    double t = 0.0;
    double distance = std::numeric_limits<double>::infinity();

    bool onRing() const {
        return ring_segment != std::numeric_limits<size_t>::max();
    }
};

struct SegmentRef {
    Vec3d a, b;
    size_t edge;
    SegmentRef(const Vec3d& first, const Vec3d& last, size_t index)
        : a(first), b(last), edge(index) {}
};

struct BranchSegment {
    size_t edge_index;
    size_t geometry_segment;
    Vec3d a;
    Vec3d b;
};

const size_t kInvalidIndex = std::numeric_limits<size_t>::max();

double clamp01(double value) { return std::max(0.0, std::min(1.0, value)); }

Vec3d lerp3(const Vec3d& a, const Vec3d& b, double t) { return a + (b - a) * t; }

void appendUnique(std::vector<Vec3d>& points, const Vec3d& point);

double laneEndpointZ(const IntersectionInput& input, const LaneId& lane_id, bool entry) {
    const Lane* lane = input.findLane(lane_id);
    if (!lane || lane->geometry.points.empty()) return 0.0;
    return entry ? lane->geometry.points.back().z() : lane->geometry.points.front().z();
}

// 由车道端点切向与骨干边求挂接点；平行时退化为线段最近投影点。
Vec3d lineSegmentAttachment(const Vec2d& origin, const Vec2d& direction,
                            const Vec3d& a, const Vec3d& b, double* t_out,
                            double* line_station_out = nullptr) {
    const Vec2d d = direction.norm() > EPS ? direction.normalized() : Vec2d(1, 0);
    const Vec2d s = b.xy() - a.xy();
    const double denominator = cross2d(d, s);
    double u = 0.0;
    if (std::abs(denominator) > 1e-10) {
        u = cross2d(a.xy() - origin, d) / denominator;
    } else {
        u = (origin - a.xy()).dot(s) / std::max(1e-12, s.squaredNorm());
    }
    u = clamp01(u);
    if (t_out) *t_out = u;
    const Vec3d result = lerp3(a, b, u);
    if (line_station_out) *line_station_out = (result.xy() - origin).dot(d);
    return result;
}

std::vector<SegmentRef> ringSegments(const NonMotorBackbone& backbone) {
    std::vector<SegmentRef> result;
    const auto& points = backbone.ring_geometry.points;
    if (points.size() < 2) return result;
    for (size_t i = 0; i + 1 < points.size(); ++i) {
        if ((points[i].xy() - points[i + 1].xy()).norm() > 1e-8)
            result.push_back(SegmentRef{points[i], points[i + 1], i});
    }
    return result;
}

std::vector<BranchSegment> branchSegments(const NonMotorBackbone& backbone) {
    std::vector<BranchSegment> result;
    for (size_t edge_index = 0; edge_index < backbone.edges.size(); ++edge_index) {
        const NonMotorBackboneEdge& edge = backbone.edges[edge_index];
        if (edge.ring || edge.geometry.points.size() < 2) continue;
        for (size_t i = 1; i < edge.geometry.points.size(); ++i) {
            if ((edge.geometry.points[i - 1].xy() - edge.geometry.points[i].xy()).norm() < 1e-8)
                continue;
            result.push_back(BranchSegment{edge_index, i - 1,
                                           edge.geometry.points[i - 1], edge.geometry.points[i]});
        }
    }
    return result;
}

bool lineSegmentIntersection(const Vec2d& origin, const Vec2d& direction,
                             const Vec3d& a, const Vec3d& b,
                             Vec3d* point, double* t_out,
                             double* line_station_out = nullptr) {
    const Vec2d d = direction.norm() > EPS ? direction.normalized() : Vec2d(1, 0);
    const Vec2d s = b.xy() - a.xy();
    const double denominator = cross2d(d, s);
    if (std::abs(denominator) <= 1e-10) return false;
    const double u = cross2d(a.xy() - origin, d) / denominator;
    if (u < -1e-9 || u > 1.0 + 1e-9) return false;
    const double clamped = clamp01(u);
    const Vec3d result = lerp3(a, b, clamped);
    if (point) *point = result;
    if (t_out) *t_out = clamped;
    if (line_station_out) *line_station_out = (result.xy() - origin).dot(d);
    return true;
}

// 沿传入方向的延长射线寻找原始最近骨干交点；退出车道传入反向切向。
// 此方向约束用于截取原始骨干路径，交点在后续布局中参与偏移，最终仅车道端点固定。
Attachment attachToBackbone(const NonMotorBackbone& backbone, const Vec2d& point,
                            const Vec2d& direction, bool require_forward) {
    Attachment best;
    const std::vector<SegmentRef> ring = ringSegments(backbone);
    const std::vector<BranchSegment> branches = branchSegments(backbone);

    // First use real intersections of the forward lane extension ray with every
    // backbone segment. This is the required nearest-backbone attachment.
    for (const BranchSegment& segment : branches) {
        double t = 0.0;
        double station = 0.0;
        Vec3d candidate;
        if (!lineSegmentIntersection(point, direction, segment.a, segment.b,
                                     &candidate, &t, &station) ||
            (require_forward && station < -1e-9))
            continue;
        const double distance = (candidate.xy() - point).norm();
        if (distance < best.distance) {
            best.point = candidate;
            best.backbone_edge = segment.edge_index;
            best.geometry_segment = segment.geometry_segment;
            best.t = t;
            best.distance = distance;
        }
    }
    for (size_t i = 0; i < ring.size(); ++i) {
        double t = 0.0;
        double station = 0.0;
        Vec3d candidate;
        if (!lineSegmentIntersection(point, direction, ring[i].a, ring[i].b,
                                     &candidate, &t, &station) ||
            (require_forward && station < -1e-9))
            continue;
        const double distance = (candidate.xy() - point).norm();
        if (distance < best.distance) {
            best.point = candidate;
            best.ring_segment = i;
            best.t = t;
            best.distance = distance;
        }
    }

    // If no exact intersection exists, retain a bounded projection fallback so
    // malformed or parallel inputs still produce a deterministic path.
    if (!std::isfinite(best.distance)) {
        for (const BranchSegment& segment : branches) {
            double t = 0.0;
            double station = 0.0;
            const Vec3d candidate = lineSegmentAttachment(point, direction, segment.a, segment.b,
                                                          &t, &station);
            if (require_forward && station < -1e-9) continue;
            const double distance = (candidate.xy() - point).norm();
            if (distance < best.distance) {
                best.point = candidate;
                best.backbone_edge = segment.edge_index;
                best.geometry_segment = segment.geometry_segment;
                best.t = t;
                best.distance = distance;
            }
        }
        for (size_t i = 0; i < ring.size(); ++i) {
            double t = 0.0;
            double station = 0.0;
            const Vec3d candidate = lineSegmentAttachment(point, direction, ring[i].a, ring[i].b,
                                                          &t, &station);
            if (require_forward && station < -1e-9) continue;
            const double distance = (candidate.xy() - point).norm();
            if (distance < best.distance) {
                best.point = candidate;
                best.ring_segment = i;
                best.t = t;
                best.distance = distance;
            }
        }
    }
    return best;
}

struct RingPathResult {
    std::vector<Vec3d> points;
};

double pathLength(const std::vector<Vec3d>& points) {
    double length = 0.0;
    for (size_t i = 1; i < points.size(); ++i)
        length += (points[i].xy() - points[i - 1].xy()).norm();
    return length;
}

RingPathResult ringPathDirected(const NonMotorBackbone& backbone, const Attachment& from,
                                const Attachment& to, int direction) {
    const std::vector<SegmentRef> ring = ringSegments(backbone);
    if (ring.empty()) return RingPathResult();
    const size_t n = ring.size();
    RingPathResult result;
    result.points.push_back(from.point);
    if (!from.onRing() || !to.onRing()) return RingPathResult();
    size_t edge = from.ring_segment;
    double t = from.t;
    const int step_direction = direction >= 0 ? 1 : -1;
    // 正向沿提取器给出的逆时针环序；反向沿相邻骨干边反向行走。
    for (size_t step = 0; step <= n; ++step) {
        const bool reaches_target = edge == to.ring_segment &&
            (step_direction > 0 ? to.t >= t - 1e-9 : to.t <= t + 1e-9);
        if (reaches_target) {
            result.points.push_back(to.point);
            return result;
        }
        if (step_direction > 0) {
            result.points.push_back(ring[edge].b);
            edge = (edge + 1) % n;
            t = 0.0;
        } else {
            result.points.push_back(ring[edge].a);
            edge = (edge + n - 1) % n;
            t = 1.0;
        }
    }
    return RingPathResult();
}

RingPathResult ringPath(const NonMotorBackbone& backbone, const Attachment& from,
                        const Attachment& to) {
    // 普通直行/左转保持需求规定的逆时针截取。
    return ringPathDirected(backbone, from, to, 1);
}

RingPathResult shortestRingPath(const NonMotorBackbone& backbone, const Attachment& from,
                                const Attachment& to) {
    const RingPathResult counter_clockwise = ringPathDirected(backbone, from, to, 1);
    const RingPathResult clockwise = ringPathDirected(backbone, from, to, -1);
    if (counter_clockwise.points.empty()) return clockwise;
    if (clockwise.points.empty()) return counter_clockwise;
    // 长度相同保留逆时针结果，保证普通场景和历史输出稳定。
    return pathLength(clockwise.points) + 1e-9 < pathLength(counter_clockwise.points)
               ? clockwise : counter_clockwise;
}

struct BranchRoute {
    std::vector<Vec3d> points;
    std::string ring_node;
};

struct BranchParent {
    std::string previous;
    size_t edge;
    bool forward;
    BranchParent() : edge(kInvalidIndex), forward(true) {}
    BranchParent(const std::string& previous_node, size_t edge_index, bool is_forward)
        : previous(previous_node), edge(edge_index), forward(is_forward) {}
};

double backboneEdgeLength(const NonMotorBackboneEdge& edge) {
    return pathLength(edge.geometry.points);
}

std::vector<std::string> ringNodeIds(const NonMotorBackbone& backbone) {
    if (backbone.ring_nodes.size() >= 3) return backbone.ring_nodes;
    std::vector<std::string> result;
    const size_t count = backbone.ring_geometry.points.size() > 1
                             ? backbone.ring_geometry.points.size() - 1 : 0;
    for (size_t i = 0; i < count; ++i) {
        std::string id = "ring_" + std::to_string(i);
        double best = std::numeric_limits<double>::infinity();
        for (const auto& node : backbone.nodes) {
            const double distance = (node.point.xy() - backbone.ring_geometry.points[i].xy()).norm();
            if (distance < best && distance <= 1e-6) {
                best = distance;
                id = node.id;
            }
        }
        result.push_back(id);
    }
    return result;
}

size_t ringNodeIndex(const std::vector<std::string>& ids, const std::string& id) {
    for (size_t i = 0; i < ids.size(); ++i)
        if (ids[i] == id) return i;
    return kInvalidIndex;
}

void appendOrientedEdge(std::vector<Vec3d>& points, const NonMotorBackboneEdge& edge,
                        bool forward) {
    if (edge.geometry.points.empty()) return;
    if (forward) {
        for (size_t i = 1; i < edge.geometry.points.size(); ++i)
            appendUnique(points, edge.geometry.points[i]);
    } else {
        for (size_t i = edge.geometry.points.size() - 1; i > 0; --i)
            appendUnique(points, edge.geometry.points[i - 1]);
    }
}

bool shortestBranchToRing(const NonMotorBackbone& backbone, const std::string& start,
                          const std::set<std::string>& ring_nodes,
                          std::vector<Vec3d>& points,
                          std::string* ring_node, double* route_length) {
    std::map<std::string, std::vector<std::pair<std::string, BranchParent>>> adjacency;
    for (size_t i = 0; i < backbone.edges.size(); ++i) {
        const NonMotorBackboneEdge& edge = backbone.edges[i];
        if (edge.ring || edge.from.empty() || edge.to.empty()) continue;
        adjacency[edge.from].push_back(std::make_pair(
            edge.to, BranchParent{edge.from, i, true}));
        adjacency[edge.to].push_back(std::make_pair(
            edge.from, BranchParent{edge.to, i, false}));
    }

    std::map<std::string, double> distances;
    std::map<std::string, BranchParent> parents;
    typedef std::pair<double, std::string> QueueEntry;
    std::priority_queue<QueueEntry, std::vector<QueueEntry>, std::greater<QueueEntry>> queue;
    distances[start] = 0.0;
    queue.push(std::make_pair(0.0, start));
    std::string target;
    while (!queue.empty()) {
        const QueueEntry current = queue.top();
        queue.pop();
        if (current.first > distances[current.second] + 1e-9) continue;
        if (ring_nodes.count(current.second) != 0) {
            target = current.second;
            break;
        }
        const auto neighbors = adjacency.find(current.second);
        if (neighbors == adjacency.end()) continue;
        for (const auto& item : neighbors->second) {
            const NonMotorBackboneEdge& edge = backbone.edges[item.second.edge];
            const double candidate = current.first + backboneEdgeLength(edge);
            auto known = distances.find(item.first);
            if (known == distances.end() || candidate + 1e-9 < known->second) {
                distances[item.first] = candidate;
                parents[item.first] = item.second;
                queue.push(std::make_pair(candidate, item.first));
            }
        }
    }
    if (target.empty()) return false;

    std::vector<std::pair<size_t, bool>> route_edges;
    std::string node = target;
    while (node != start) {
        const auto parent = parents.find(node);
        if (parent == parents.end()) return false;
        route_edges.push_back(std::make_pair(parent->second.edge, parent->second.forward));
        node = parent->second.previous;
    }
    std::reverse(route_edges.begin(), route_edges.end());
    for (const auto& route_edge : route_edges) {
        appendOrientedEdge(points, backbone.edges[route_edge.first], route_edge.second);
    }
    if (ring_node) *ring_node = target;
    if (route_length) *route_length = distances[target];
    return true;
}

BranchRoute branchToRing(const NonMotorBackbone& backbone, const Attachment& attachment,
                         const std::set<std::string>& ring_nodes) {
    BranchRoute best;
    if (attachment.onRing() || attachment.backbone_edge >= backbone.edges.size()) return best;
    const NonMotorBackboneEdge& edge = backbone.edges[attachment.backbone_edge];
    if (edge.geometry.points.size() < 2) return best;

    for (int side = 0; side < 2; ++side) {
        const bool toward_from = side == 0;
        const std::string start_node = toward_from ? edge.from : edge.to;
        if (start_node.empty()) continue;
        BranchRoute candidate;
        candidate.points.push_back(attachment.point);
        if (toward_from) {
            for (size_t i = attachment.geometry_segment + 1; i > 0; --i)
                appendUnique(candidate.points, edge.geometry.points[i - 1]);
        } else {
            for (size_t i = attachment.geometry_segment + 1; i < edge.geometry.points.size(); ++i)
                appendUnique(candidate.points, edge.geometry.points[i]);
        }
        if (ring_nodes.count(start_node) != 0) {
            candidate.ring_node = start_node;
        } else if (!shortestBranchToRing(backbone, start_node, ring_nodes,
                                         candidate.points,
                                         &candidate.ring_node, nullptr)) {
            continue;
        }
        const double length = pathLength(candidate.points);
        if (best.points.empty() || length + 1e-9 < pathLength(best.points)) best = candidate;
    }
    return best;
}

Attachment ringNodeAnchor(const NonMotorBackbone& backbone, size_t node_index, int direction) {
    Attachment anchor;
    const std::vector<SegmentRef> ring = ringSegments(backbone);
    const auto ids = ringNodeIds(backbone);
    if (ring.empty() || node_index >= ids.size()) return anchor;
    // ring_nodes是拓扑节点序，RoadEdge Gap可含多个形态点，不能当作几何段序。
    const auto node = std::find_if(backbone.nodes.begin(), backbone.nodes.end(),
        [&](const NonMotorBackboneNode& value) { return value.id == ids[node_index]; });
    if (node != backbone.nodes.end()) {
        size_t segment = 0;
        while (segment < ring.size() && (ring[segment].a.xy() - node->point.xy()).norm() > 1e-6) ++segment;
        if (segment == ring.size()) return anchor;
        node_index = segment;
    } else if (node_index >= ring.size()) return anchor;
    anchor.point = ring[node_index].a;
    if (direction >= 0) {
        anchor.ring_segment = node_index;
        anchor.t = 0.0;
    } else {
        anchor.ring_segment = (node_index + ring.size() - 1) % ring.size();
        anchor.point = ring[anchor.ring_segment].b;
        anchor.t = 1.0;
    }
    anchor.distance = 0.0;
    return anchor;
}

void appendReversed(std::vector<Vec3d>& target, const std::vector<Vec3d>& source) {
    if (source.empty()) return;
    for (size_t i = source.size(); i > 1; --i) appendUnique(target, source[i - 2]);
}

RingPathResult backbonePath(const NonMotorBackbone& backbone, const Attachment& from,
                            const Attachment& to, bool shortest) {
    // 同一分段短距离合并后两端代表同一个公共交点，不再输出零长度骨干段。
    if ((from.point.xy() - to.point.xy()).norm() < 1e-8) {
        RingPathResult result;
        result.points.push_back(from.point);
        return result;
    }
    if (from.onRing() && to.onRing())
        return shortest ? shortestRingPath(backbone, from, to) : ringPath(backbone, from, to);

    // 同一分支上的两个挂接点直接截取边内区间，不能先到主环再原路返回。
    if (!from.onRing() && !to.onRing() && from.backbone_edge == to.backbone_edge &&
        from.backbone_edge < backbone.edges.size()) {
        const auto& geometry = backbone.edges[from.backbone_edge].geometry.points;
        const bool forward = from.geometry_segment < to.geometry_segment ||
            (from.geometry_segment == to.geometry_segment && from.t <= to.t);
        RingPathResult result; result.points.push_back(from.point);
        if (forward) {
            for (size_t i = from.geometry_segment + 1; i <= to.geometry_segment; ++i)
                appendUnique(result.points, geometry[i]);
        } else {
            for (size_t i = from.geometry_segment; i > to.geometry_segment; --i)
                appendUnique(result.points, geometry[i]);
        }
        appendUnique(result.points, to.point);
        return result;
    }

    const std::vector<std::string> node_ids = ringNodeIds(backbone);
    const std::set<std::string> ring_nodes(node_ids.begin(), node_ids.end());
    const BranchRoute start_branch = from.onRing() ? BranchRoute() : branchToRing(backbone, from, ring_nodes);
    const BranchRoute end_branch = to.onRing() ? BranchRoute() : branchToRing(backbone, to, ring_nodes);
    if ((!from.onRing() && start_branch.points.empty()) || (!to.onRing() && end_branch.points.empty()))
        return RingPathResult();

    const size_t start_index = from.onRing() ? kInvalidIndex : ringNodeIndex(node_ids, start_branch.ring_node);
    const size_t end_index = to.onRing() ? kInvalidIndex : ringNodeIndex(node_ids, end_branch.ring_node);
    if ((!from.onRing() && start_index == kInvalidIndex) || (!to.onRing() && end_index == kInvalidIndex))
        return RingPathResult();

    auto make_candidate = [&](int direction) {
        RingPathResult result;
        const Attachment start_anchor = from.onRing() ? from : ringNodeAnchor(backbone, start_index, direction);
        const Attachment end_anchor = to.onRing() ? to : ringNodeAnchor(backbone, end_index, direction);
        const RingPathResult middle = ringPathDirected(backbone, start_anchor, end_anchor, direction);
        if (middle.points.empty()) return result;
        if (from.onRing()) result.points.push_back(from.point);
        else result.points = start_branch.points;
        for (size_t i = 1; i < middle.points.size(); ++i) appendUnique(result.points, middle.points[i]);
        if (!to.onRing()) appendReversed(result.points, end_branch.points);
        return result;
    };

    const RingPathResult counter_clockwise = make_candidate(1);
    if (!shortest) return counter_clockwise;
    const RingPathResult clockwise = make_candidate(-1);
    if (counter_clockwise.points.empty()) return clockwise;
    if (clockwise.points.empty()) return counter_clockwise;
    return pathLength(clockwise.points) + 1e-9 < pathLength(counter_clockwise.points)
               ? clockwise : counter_clockwise;
}

void appendUnique(std::vector<Vec3d>& points, const Vec3d& point) {
    if (points.empty() || (points.back().xy() - point.xy()).norm() > 1e-7) points.push_back(point);
}

bool isUTurn(ConnTurnType turn_type) {
    return turn_type == ConnTurnType::UTurnLeft || turn_type == ConnTurnType::UTurnRight;
}

/// 非掉头两端若命中同一条骨干的同一几何分段，且两交点间有效长度很短，
/// 将两个交点收敛到其间骨干中点，避免短段偏移和端部圆角叠加成异常形态。
/// @param backbone 骨干图，用于读取命中分段长度。
/// @param from 进入车道延长线的骨干交点，成功时改写为公共交点。
/// @param to 退出车道延长线的骨干交点，成功时改写为公共交点。
/// @param max_length 允许合并的最大有效骨干长度，非正值禁用该规则。
/// @return 是否完成了公共交点合并。
bool collapseShortSameSegmentAttachments(const NonMotorBackbone& backbone,
                                         Attachment& from, Attachment& to,
                                         double max_length) {
    if (max_length <= 0.0 || from.onRing() != to.onRing()) return false;

    double segment_length = 0.0;
    bool same_segment = false;
    if (from.onRing()) {
        const auto segments = ringSegments(backbone);
        if (from.ring_segment < segments.size() &&
            from.ring_segment == to.ring_segment) {
            segment_length = (segments[from.ring_segment].b.xy() -
                              segments[from.ring_segment].a.xy()).norm();
            same_segment = true;
        }
    } else if (from.backbone_edge == to.backbone_edge &&
               from.backbone_edge < backbone.edges.size() &&
               from.geometry_segment == to.geometry_segment) {
        const auto& geometry = backbone.edges[from.backbone_edge].geometry.points;
        if (from.geometry_segment + 1 < geometry.size()) {
            segment_length = (geometry[from.geometry_segment + 1].xy() -
                              geometry[from.geometry_segment].xy()).norm();
            same_segment = true;
        }
    }
    if (!same_segment || segment_length < 1e-9) return false;

    const double effective_length = std::abs(from.t - to.t) * segment_length;
    if (!std::isfinite(effective_length) || effective_length > max_length + 1e-9)
        return false;

    const Vec3d midpoint = (from.point + to.point) * 0.5;
    const double middle_t = 0.5 * (from.t + to.t);
    from.point = midpoint;
    to.point = midpoint;
    from.t = middle_t;
    to.t = middle_t;
    return true;
}

// 某些数据的骨干提取结果可能没有填充 ring_geometry，但所有复合边按节点仍构成
// 一个简单闭环。若不恢复该拓扑，两个不同分支挂接点无法回溯到环，会退化为
// 进入段直接连接退出段。这里只在所有边端点无向度均为2且图连通时恢复环序，
// 不把含真实分支的图误判为主环。
bool recoverClosedRingFromEdgeGraph(NonMotorBackbone& backbone) {
    if (!backbone.ring_geometry.points.empty() || backbone.edges.size() < 3) return false;
    std::map<std::string, std::vector<size_t>> adjacency;
    for (size_t i = 0; i < backbone.edges.size(); ++i) {
        const auto& edge = backbone.edges[i];
        if (edge.from.empty() || edge.to.empty() || edge.geometry.points.size() < 2 || edge.from == edge.to)
            return false;
        adjacency[edge.from].push_back(i);
        adjacency[edge.to].push_back(i);
    }
    for (const auto& item : adjacency)
        if (item.second.size() != 2) return false;

    std::vector<bool> visited(backbone.edges.size(), false);
    std::vector<size_t> edge_order;
    std::vector<std::string> node_order;
    std::string start_node = backbone.edges.front().from;
    std::string current_node = start_node;
    size_t current_edge = 0;
    for (size_t step = 0; step < backbone.edges.size(); ++step) {
        if (visited[current_edge]) return false;
        visited[current_edge] = true;
        const auto& edge = backbone.edges[current_edge];
        const bool forward = edge.from == current_node;
        if (!forward && edge.to != current_node) return false;
        node_order.push_back(current_node);
        edge_order.push_back(current_edge);
        const std::string next_node = forward ? edge.to : edge.from;
        const auto& neighbors = adjacency[next_node];
        current_edge = neighbors[0] == current_edge ? neighbors[1] : neighbors[0];
        current_node = next_node;
    }
    if (current_node != start_node || std::any_of(visited.begin(), visited.end(), [](bool value) { return !value; }))
        return false;

    std::vector<Vec3d> ring_points;
    for (size_t i = 0; i < edge_order.size(); ++i) {
        const auto& edge = backbone.edges[edge_order[i]];
        const bool forward = edge.from == node_order[i];
        if (forward) {
            for (const auto& point : edge.geometry.points) appendUnique(ring_points, point);
        } else {
            for (size_t j = edge.geometry.points.size(); j > 0; --j)
                appendUnique(ring_points, edge.geometry.points[j - 1]);
        }
    }
    if (ring_points.size() < 3) return false;
    appendUnique(ring_points, ring_points.front());
    if (ring_points.size() < 4) return false;

    // 主环方向必须与普通左转的逆时针截取约定一致。边图首条边的输入方向
    // 不稳定，因此恢复后显式检查有向面积，必要时反转几何和节点序列。
    double signed_area_twice = 0.0;
    for (size_t i = 1; i < ring_points.size(); ++i)
        signed_area_twice += cross2d(ring_points[i - 1].xy(), ring_points[i].xy());
    if (signed_area_twice < -1e-9) {
        std::vector<Vec3d> counter_clockwise;
        counter_clockwise.reserve(ring_points.size());
        counter_clockwise.push_back(ring_points.front());
        for (size_t i = ring_points.size() - 1; i > 1; --i)
            counter_clockwise.push_back(ring_points[i - 1]);
        counter_clockwise.push_back(counter_clockwise.front());
        ring_points.swap(counter_clockwise);

        std::vector<std::string> counter_clockwise_nodes;
        counter_clockwise_nodes.reserve(node_order.size());
        counter_clockwise_nodes.push_back(node_order.front());
        for (size_t i = node_order.size(); i > 1; --i)
            counter_clockwise_nodes.push_back(node_order[i - 1]);
        node_order.swap(counter_clockwise_nodes);
    }
    backbone.ring_nodes = node_order;
    backbone.ring_geometry.points = ring_points;
    backbone.closed = true;
    backbone.ring_edges.clear();
    backbone.branch_edges.clear();
    for (auto& edge : backbone.edges) {
        edge.ring = true;
        edge.branch = false;
        backbone.ring_edges.push_back(edge.id);
    }
    return true;
}

bool segmentHitsPolygon(const Vec2d& a, const Vec2d& b, const Polygon2d& polygon) {
    if (polygonContains(polygon, a) || polygonContains(polygon, b)) return true;
    if (polygon.outer.size() < 2) return false;
    for (size_t i = 0; i < polygon.outer.size(); ++i) {
        const Vec2d c = polygon.outer[i].xy();
        const Vec2d d = polygon.outer[(i + 1) % polygon.outer.size()].xy();
        if (segmentsIntersect(a, b, c, d, nullptr)) return true;
    }
    return false;
}

bool pathHitsObstacle(const LineString2d& path, const std::vector<Obstacle>& obstacles) {
    for (size_t i = 1; i < path.points.size(); ++i)
        for (const auto& obstacle : obstacles)
            if (segmentHitsPolygon(path.points[i - 1].xy(), path.points[i].xy(), obstacle.geometry))
                return true;
    return false;
}

// 线段净距同时考虑内部相交和双方端点；只检查顶点会漏掉跨越短边的长段。
double segmentDistance(const Vec3d& a, const Vec3d& b, const Vec3d& c, const Vec3d& d) {
    if (segmentsIntersect(a.xy(), b.xy(), c.xy(), d.xy(), nullptr)) return 0;
    return std::min(std::min(pointToSegment(a.xy(), c, d).first, pointToSegment(b.xy(), c, d).first),
                    std::min(pointToSegment(c.xy(), a, b).first, pointToSegment(d.xy(), a, b).first));
}

bool pathOutsideFence(const LineString2d& path, const Polygon2d& fence) {
    if (fence.outer.empty()) return false;
    for (size_t i = 0; i < path.points.size(); ++i) {
        if (!polygonContains(fence, path.points[i].xy())) {
            const bool endpoint = i == 0 || i + 1 == path.points.size();
            if (!endpoint || pointToPolygonDist(path.points[i].xy(), fence) > 0.05) return true;
        }
        if (i == 0) continue;
        const Vec2d a = path.points[i - 1].xy();
        const Vec2d b = path.points[i].xy();
        for (int sample = 1; sample < 4; ++sample) {
            const Vec2d point = a + (b - a) * (sample / 4.0);
            if (polygonContains(fence, point)) continue;
            const double endpoint_distance = std::min((point - path.points.front().xy()).norm(),
                                                      (point - path.points.back().xy()).norm());
            if (endpoint_distance > 0.20 || pointToPolygonDist(point, fence) > 0.0025) return true;
        }
    }
    return false;
}

bool pathHitsBoundary(const LineString2d& path, const std::vector<Boundary>& boundaries) {
    if (path.points.size() < 2) return false;
    for (const auto& boundary : boundaries) {
        // Other、隔离带、绿化带同样是不可穿越边界；type不能成为绕行/验收的豁免。
        if (boundary.geometry.points.size() < 2) continue;
        for (size_t i = 1; i < path.points.size(); ++i)
            for (size_t j = 1; j < boundary.geometry.points.size(); ++j) {
                Vec2d intersection;
                if (!segmentsIntersect(path.points[i - 1].xy(), path.points[i].xy(),
                                       boundary.geometry.points[j - 1].xy(),
                                       boundary.geometry.points[j].xy(), &intersection))
                    continue;
                const bool endpoint = (intersection - path.points.front().xy()).norm() < 1e-4 ||
                                      (intersection - path.points.back().xy()).norm() < 1e-4;
                if (!endpoint) return true;
            }
    }
    return false;
}

double minimumPolylineDistance(const LineString2d& first, const LineString2d& second) {
    double result = std::numeric_limits<double>::infinity();
    if (first.points.empty() || second.points.size() < 2) return result;
    for (size_t i = 1; i < first.points.size(); ++i)
        for (size_t j = 1; j < second.points.size(); ++j)
            result = std::min(result, segmentDistance(first.points[i - 1], first.points[i], second.points[j - 1], second.points[j]));
    return result;
}

double minimumRoadEdgeDistance(const LineString2d& path, const std::vector<Boundary>& boundaries) {
    double result = std::numeric_limits<double>::infinity();
    for (size_t i = 1; i < path.points.size(); ++i)
        for (const auto& boundary : boundaries) {
            if (boundary.type != Boundary::Type::RoadEdge || boundary.geometry.points.size() < 2) continue;
            for (size_t j = 1; j < boundary.geometry.points.size(); ++j)
                result = std::min(result, segmentDistance(path.points[i - 1], path.points[i],
                                                         boundary.geometry.points[j - 1], boundary.geometry.points[j]));
        }
    return result;
}

// 后续布局以实际有效路径段为占用记录，不能再通过点下标猜测边引用。
#include "nonmotor/nonmotor_path_layout.inc"

}  // namespace

NonMotorConnectivityGenerator::NonMotorConnectivityGenerator(const NonmotorConifg& cfg)
    : config_(cfg) {}

std::vector<ConnectivityCurve> NonMotorConnectivityGenerator::generate(
    const IntersectionInput& input, const std::vector<ConnectivityCurve>& motor_curves,
    NonMotorBackbone& backbone) {
    return layoutNonMotorPaths(input, motor_curves, backbone, config_);
}

}  // namespace isg
