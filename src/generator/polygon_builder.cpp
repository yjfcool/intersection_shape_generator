#include "generator/polygon_builder.h"
#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <utility>
#include "utils/clipper.hpp"

#include "utils/shapefile.hpp"

namespace isg {

/**
 * 精细路口面构建
 * 处理：1.路口面 = 路口外轮廓有效范围的道路边缘线 + 进入组尾端(或退出组首端)的回缩点(或外扩点) / (暂时不用停止线);
 *      2.通过进入或退出组的缩扩线截取有效范围段的外轮廓道路边缘线;
 *      3.缩扩根据mode决定：mode=1|3进入组向路口内回缩、退出组向路口外扩, mode=2进入和退出组都向路口外扩;
 *      4.详细描述:按当前组中所有lane和laneedge各自去重后的几何线端点(进入组的线尾点或退出组线首点)的垂线按缩扩方向移动,
 *        lane、laneedge、boundary属于不同构面源类别,重复过滤只能在各自类别内部进行;
 *        laneedge同时被进入组和退出组共用时,只保留其矢量化方向与当前组方向同向的缩扩点;
 *        垂线与共端点所有连接线(进入|退出线和生成的路口内车道线)的交点(如果与连接线有多个交点时取距离端点最近的点)作为缩扩端点;
 *        缩扩线向两端延长后截取合并后的有效道路边缘, 同一道路边缘只保留最外侧/最内侧命中的连续有效段, 中间内部边缘直接过滤;
 * 输出：封闭可凹陷多边形（路口面）
 */
class IntersectionAreaBuilderImpl {
    // 进入退出组端侧边界向路口回缩或外扩的距离(米)
    double groupCutOffset = 0.05;
    // 端侧边界折线两端用于对齐道路边缘线的外延长度(米)
    double cutLineExtension = 4.0;
    // 路口模式: 1=bxn,2=jd,3=ds
    int mode = 1;
    // 路口面顶点吸附容差(米),距离小于此值的顶点会被合并。取值范围：0.1~2.0
    double snapTolerance = 0.5;
    // 路口面多边形绕向 (逆时针:"counterclockwise" | 顺时针:"clockwise")
    std::string winding = "clockwise";

public:
    explicit IntersectionAreaBuilderImpl(
        double _groupCutOffset,
        double _cutLineExtension = 7.0,
        double _snapTolerance = 0.5,
        std::string _winding = "clockwise")
        : snapTolerance(_snapTolerance),
          winding(_winding),
          groupCutOffset(_groupCutOffset),
          cutLineExtension(std::max(cutLineExtension, _cutLineExtension)) {}

    static std::vector<BoundaryLine> MergeBoundaryLines(
        std::vector<BoundaryLine> bnds) {
        auto sameEndpoint = [](const Vec3d& a, const Vec3d& b) {
            return dist(a, b) <= 1e-9;
        };
        auto sourceSegmentDir = [](const BoundaryLine& src, int segment_index) {
            if (segment_index >= 0 && segment_index < (int)src.segment_dirs.size()) {
                Vec2d d = src.segment_dirs[segment_index];
                if (d.norm() > 1e-8)
                    return d.normalized();
            }
            if (segment_index >= 0 && segment_index + 1 < (int)src.pts.size()) {
                Vec2d d = xyOf(src.pts[segment_index + 1]) - xyOf(src.pts[segment_index]);
                if (d.norm() > 1e-8)
                    return d.normalized();
            }
            return Vec2d(1, 0);
        };
        auto appendWithoutDuplicate = [&](BoundaryLine& dst,
                                          const BoundaryLine& src,
                                          int begin, int end, int step) {
            for (int i = begin; i != end; i += step) {
                if (i < 0 || i >= (int)src.pts.size())
                    continue;
                if (dst.pts.empty()) {
                    dst.pts.push_back(src.pts[i]);
                    continue;
                }
                if (dist(dst.pts.back(), src.pts[i]) <= 1e-9)
                    continue;
                int segment_index = step > 0 ? i - 1 : i;
                dst.segment_dirs.push_back(sourceSegmentDir(src, segment_index));
                dst.pts.push_back(src.pts[i]);
            }
        };
        auto addSourceIds = [](std::vector<std::string>& dst,
                               const std::vector<std::string>& src) {
            for (const auto& id : src) {
                if (!id.empty() && std::find(dst.begin(), dst.end(), id) == dst.end())
                    dst.push_back(id);
            }
        };
        auto addSourceEndpoints = [](
                std::vector<std::pair<Vec3d, Vec3d>>& dst,
                const std::vector<std::pair<Vec3d, Vec3d>>& src) {
            for (const auto& endpoints : src) {
                bool exists = false;
                for (const auto& existing : dst) {
                    if ((dist(existing.first, endpoints.first) <= 1e-9 &&
                         dist(existing.second, endpoints.second) <= 1e-9) ||
                        (dist(existing.first, endpoints.second) <= 1e-9 &&
                         dist(existing.second, endpoints.first) <= 1e-9)) {
                        exists = true;
                        break;
                    }
                }
                if (!exists)
                    dst.push_back(endpoints);
            }
        };
        auto mergedPair = [&](const BoundaryLine& a,
                              const BoundaryLine& b,
                              BoundaryLine& out) {
            if (a.pts.empty() || b.pts.empty())
                return false;
            const Vec3d& a0 = a.pts.front();
            const Vec3d& a1 = a.pts.back();
            const Vec3d& b0 = b.pts.front();
            const Vec3d& b1 = b.pts.back();

            out = {};
            if (sameEndpoint(a1, b0)) {
                out.pts = a.pts;
                out.segment_dirs = a.segment_dirs;
                appendWithoutDuplicate(out, b, 1, (int)b.pts.size(), 1);
            } else if (sameEndpoint(a1, b1)) {
                out.pts = a.pts;
                out.segment_dirs = a.segment_dirs;
                appendWithoutDuplicate(out, b, (int)b.pts.size() - 2, -1, -1);
            } else if (sameEndpoint(a0, b1)) {
                out.pts = b.pts;
                out.segment_dirs = b.segment_dirs;
                appendWithoutDuplicate(out, a, 1, (int)a.pts.size(), 1);
            } else if (sameEndpoint(a0, b0)) {
                appendWithoutDuplicate(out, b, (int)b.pts.size() - 1, -1, -1);
                appendWithoutDuplicate(out, a, 1, (int)a.pts.size(), 1);
            } else {
                return false;
            }
            addSourceIds(out.source_ids, a.source_ids);
            addSourceIds(out.source_ids, b.source_ids);
            addSourceEndpoints(out.source_endpoints, a.source_endpoints);
            addSourceEndpoints(out.source_endpoints, b.source_endpoints);
            return true;
        };

        bool changed = true;
        while (changed) {
            changed = false;
            for (int i = 0; i < (int)bnds.size() && !changed; ++i) {
                if (bnds[i].pts.size() < 2)
                    continue;
                for (int j = i + 1; j < (int)bnds.size(); ++j) {
                    if (bnds[j].pts.size() < 2)
                        continue;
                    BoundaryLine merged;
                    if (!mergedPair(bnds[i], bnds[j], merged))
                        continue;
                    bnds[i] = std::move(merged);
                    bnds.erase(bnds.begin() + j);
                    changed = true;
                    break;
                }
            }
        }

        std::vector<BoundaryLine> out;
        for (auto& bnd : bnds) {
            if (bnd.pts.size() >= 2)
                out.push_back(std::move(bnd));
        }
        return out;
    }

    IntersectionArea build(
        const IntersectionInput& inp,
        const std::vector<ConnectivityCurve>& centerlines,
        const std::vector<ConnectivityLaneEdge>& edgelines)
    {
        mode = inp.mode;

        IntersectionArea junctArea;
        junctArea.id = "int_polygon_" + inp.id;

        // 1. 确定路口中心和范围
        Vec2d center = computeCenter(inp, centerlines);
        double radius = computeRadius(inp, centerlines, center);

        std::vector<Vec3d> finePoly = buildFinePolygon(inp, centerlines, center, radius);
        if (finePoly.size() >= 4) {
            junctArea.geometry.outer = finePoly;
            junctArea.is_rough = false;
            return junctArea;
        }

        std::cout<<"!!!!! buildFinePolygon:"<<inp.id<<" failed !!!!!"<<std::endl;
        // 2. 裁剪道路边缘线，获取路口侧端点集
        std::vector<EdgeEndpoint> endpoints;
        for (int i = 0; i < inp.boundaries.size(); ++i) {
            const auto& re = inp.boundaries[i];
            // 获取有效路口道路边缘，其他类型视作路口内部非外轮廓边缘
            if (!isUsableRoadEdgeBoundary(inp, re))
                continue;

            // 找到路口侧的端点（靠近center的端点）
            double d0 = dist(re.geometry.points.front(), center);
            double d1 = dist(re.geometry.points.back(), center);

            Vec3d candidatePt;
            bool isStart;
            if (d0 < d1) {
                candidatePt = re.geometry.points.front();
                isStart = true;
            } else {
                candidatePt = re.geometry.points.back();
                isStart = false;
            }

            // 检查是否在路口范围内（在 rough_area 内或在 radius 内）
            bool inRange = false;
            if (!inp.area.geometry.empty() && !inp.area.geometry.outer.empty()) {
                inRange = pointInPolygon(candidatePt, inp.area.geometry.outer) ||
                    dist(candidatePt, center) < radius * 1.5;
            } else {
                inRange = dist(candidatePt, center) < radius * 1.5;
            }
            if (inRange) {
                // 吸附到最近的车道边线端点
                Vec3d snapped = snapToEdgePt(candidatePt, edgelines, inp, snapTolerance);
                endpoints.emplace_back(EdgeEndpoint{snapped, i, isStart});
            }
        }

        // 若道路边缘线端点不足，从连通关系连接点的外侧边线端点补充
        if (endpoints.size() < 3) {
            supplementEndpointsFromLaneEdges(endpoints, inp, edgelines, center, radius);
        }

        if (endpoints.size() < 3) {
            // 终极回退：用所有连接点的凸包
            junctArea.geometry.outer = buildConvexHullFallback(inp, centerlines);
            junctArea.is_rough = false;
            return junctArea;
        }

        // 3. 按极角排序
        std::sort(endpoints.begin(), endpoints.end(),
                  [&](const EdgeEndpoint& a, const EdgeEndpoint& b) {
                      double angA = std::atan2(a.pt[1] - center[1], a.pt[0] - center[0]);
                      double angB = std::atan2(b.pt[1] - center[1], b.pt[0] - center[0]);
                      if (winding == "clockwise")
                          return angA > angB;
                      else
                          return angA < angB;
                  });

        // 去重（距离过近的端点合并）
        deduplicateEndpoints(endpoints, 0.1);

        // 4. 连接端点构建多边形
        std::vector<Vec3d> rawPoly;
        for (int i = 0; i < (int)endpoints.size(); ++i) {
            const EdgeEndpoint& cur = endpoints[i];
            const EdgeEndpoint& next = endpoints[(i + 1) % endpoints.size()];

            rawPoly.push_back(cur.pt);

            // 若相邻端点来自同一条道路边缘线，沿边缘线几何连接
            if (cur.roadEdgeIdx >= 0 && cur.roadEdgeIdx == next.roadEdgeIdx) {
                std::vector<Vec3d> seg = extractEdgeSegment(
                    inp.boundaries[cur.roadEdgeIdx], cur.pt, next.pt, snapTolerance);
                for (auto& p : seg) rawPoly.push_back(p);
            }
            // 否则直连（或弧线连，按参数）
        }

        // 5. 闭合多边形
        if (!rawPoly.empty() && dist(rawPoly.front(), rawPoly.back()) > EPS)
            rawPoly.push_back(rawPoly.front());

        // 6. 多边形修复（移除共线点，确保逆/顺时针）
        junctArea.geometry.outer = repairPolygon(rawPoly, center);

        junctArea.is_rough = false;
        return junctArea;
    }

private:
    struct AreaChain {
        std::vector<Vec3d> pts;
    };

    struct LineIntersection {
        Vec3d pt{0, 0, 0};
        double station_a = 0.0;
        double station_b = 0.0;
        Vec2d dir_b{1, 0};
        Vec2d chain_dir_b{1, 0};
    };

    static std::vector<Vec3d> openRing(std::vector<Vec3d> ring) {
        while (ring.size() > 1 && dist(ring.front(), ring.back()) < 1e-6)
            ring.pop_back();
        return ring;
    }

    double polylineLength(const std::vector<Vec3d>& pts) const {
        double len = 0.0;
        for (int i = 0; i + 1 < (int)pts.size(); ++i)
            len += dist(pts[i], pts[i + 1]);
        return len;
    }

    Vec3d pointAtStation(const std::vector<Vec3d>& pts, double station) const {
        if (pts.empty())
            return Vec3d(0, 0, 0);
        if (pts.size() == 1)
            return pts.front();
        double total = polylineLength(pts);
        station = std::max(0.0, std::min(total, station));
        double acc = 0.0;
        for (int i = 0; i + 1 < (int)pts.size(); ++i) {
            double seg_len = dist(pts[i], pts[i + 1]);
            if (seg_len < 1e-9)
                continue;
            if (acc + seg_len >= station) {
                double t = (station - acc) / seg_len;
                return pts[i] * (1.0 - t) + pts[i + 1] * t;
            }
            acc += seg_len;
        }
        return pts.back();
    }

    std::vector<Vec3d> subPolylineByStation(
        const std::vector<Vec3d>& pts, double s0, double s1) const {
        std::vector<Vec3d> out;
        if (pts.size() < 2)
            return pts;
        double total = polylineLength(pts);
        s0 = std::max(0.0, std::min(total, s0));
        s1 = std::max(0.0, std::min(total, s1));
        if (s1 < s0)
            std::swap(s0, s1);
        pushUnique(out, pointAtStation(pts, s0));
        double acc = 0.0;
        for (int i = 0; i + 1 < (int)pts.size(); ++i) {
            double seg_len = dist(pts[i], pts[i + 1]);
            double next = acc + seg_len;
            if (seg_len > 1e-9 && next > s0 + 1e-9 && next < s1 - 1e-9)
                pushUnique(out, pts[i + 1]);
            acc = next;
        }
        pushUnique(out, pointAtStation(pts, s1));
        return out;
    }

    // 闭合 RoadEdge 链的两个缩扩命中点之间有两条拓扑弧。分别以各自
    // 首尾弦线为基准，计算折线段中点向路口中心侧的带权偏移；选择凹向
    // 路口内的一弧，不使用弧长作为业务判据。
    std::vector<Vec3d> subInwardClosedPolylineByStation(
        const std::vector<Vec3d>& pts, double s0, double s1,
        const Vec2d& center) const {
        const double total = polylineLength(pts);
        if (pts.size() < 3 || total <= 1e-9 || dist(pts.front(), pts.back()) > 0.03)
            return subPolylineByStation(pts, s0, s1);
        s0 = std::max(0.0, std::min(total, s0));
        s1 = std::max(0.0, std::min(total, s1));
        if (s1 < s0)
            std::swap(s0, s1);
        std::vector<Vec3d> direct = subPolylineByStation(pts, s0, s1);
        std::vector<Vec3d> wrapped = subPolylineByStation(pts, s1, total);
        std::vector<Vec3d> head = subPolylineByStation(pts, 0.0, s0);
        for (const auto& pt : head)
            pushUnique(wrapped, pt);

        auto inwardScore = [&](const std::vector<Vec3d>& arc) {
            if (arc.size() < 2)
                return -std::numeric_limits<double>::infinity();
            Vec2d chord = xyOf(arc.back()) - xyOf(arc.front());
            if (chord.norm() < 1e-9)
                return -std::numeric_limits<double>::infinity();
            chord.normalize();
            double center_side = cross2d(chord, center - xyOf(arc.front()));
            const bool center_on_chord = std::abs(center_side) < 1e-9;
            const double side = center_side >= 0.0 ? 1.0 : -1.0;
            double weighted_offset = 0.0;
            double length = 0.0;
            for (int i = 0; i + 1 < (int)arc.size(); ++i) {
                const double segment_length = dist(arc[i], arc[i + 1]);
                if (segment_length < 1e-9)
                    continue;
                const Vec2d midpoint = 0.5 * (xyOf(arc[i]) + xyOf(arc[i + 1]));
                weighted_offset += (center_on_chord
                    ? -dist(midpoint, center)
                    : side * cross2d(chord, midpoint - xyOf(arc.front()))) * segment_length;
                length += segment_length;
            }
            return length > 1e-9
                ? weighted_offset / length
                : -std::numeric_limits<double>::infinity();
        };
        return inwardScore(wrapped) > inwardScore(direct) ? wrapped : direct;
    }

    void pushUnique(std::vector<Vec3d>& out, const Vec3d& pt, double tol = 0.03) const {
        if (out.empty() || dist(out.back(), pt) > tol)
            out.push_back(pt);
    }

    double angleOf(const Vec3d& pt, const Vec2d& center) const {
        return std::atan2(pt[1] - center[1], pt[0] - center[0]);
    }

    double normalizedAngleDelta(double a0, double a1) const {
        double d = a1 - a0;
        while (d > M_PI)
            d -= 2.0 * M_PI;
        while (d < -M_PI)
            d += 2.0 * M_PI;
        return d;
    }

    AreaChain orientChain(AreaChain chain, const Vec2d& center) const {
        if (chain.pts.size() < 2)
            return chain;
        double a0 = angleOf(chain.pts.front(), center);
        double a1 = angleOf(chain.pts.back(), center);
        double delta = normalizedAngleDelta(a0, a1);
        bool should_reverse = (winding == "clockwise") ? (delta > 0.0) : (delta < 0.0);
        if (should_reverse)
            std::reverse(chain.pts.begin(), chain.pts.end());
        return chain;
    }

    Vec3d chainMidpoint(const AreaChain& chain) const {
        Vec3d mid(0, 0, 0);
        if (chain.pts.empty())
            return mid;
        for (const auto& pt : chain.pts)
            mid += pt;
        return mid / (double)chain.pts.size();
    }

    bool laneInGroup(const LaneGroup& group, const LaneId& lane_id) const {
        return std::find(group.lanes.begin(), group.lanes.end(), lane_id) != group.lanes.end();
    }

    bool isRoadEdgeBridgeBoundary(const IntersectionInput& inp, const Boundary& boundary) const {
        if (boundary.type != Boundary::Type::Other || boundary.geometry.points.size() < 2)
            return false;
        auto endpointTouchesRoadEdge = [&](const Vec3d& pt) {
            for (const auto& other : inp.boundaries) {
                if (other.id == boundary.id ||
                    other.type != Boundary::Type::RoadEdge ||
                    other.geometry.points.size() < 2)
                    continue;
                if (dist(pt, other.geometry.points.front()) <= 1e-9 ||
                    dist(pt, other.geometry.points.back()) <= 1e-9)
                    return true;
            }
            return false;
        };
        return endpointTouchesRoadEdge(boundary.geometry.points.front()) ||
               endpointTouchesRoadEdge(boundary.geometry.points.back());
    }

    bool isUsableRoadEdgeBoundary(const IntersectionInput& inp, const Boundary& boundary) const {
        if (boundary.geometry.points.size() < 2)
            return false;
        if (boundary.type != Boundary::Type::RoadEdge &&
            !isRoadEdgeBridgeBoundary(inp, boundary))
            return false;
        if (inp.findLane(boundary.id))
            return false;
        return true;
    }

    bool groupContainsLaneId(const LaneGroup& group, const LaneId& lane_id) const {
        return laneInGroup(group, lane_id);
    }

    std::vector<LaneId> uniqueGroupLaneIds(const LaneGroup& group) const {
        std::vector<LaneId> ids;
        for (const auto& lane_id : group.lanes) {
            if (std::find(ids.begin(), ids.end(), lane_id) == ids.end())
                ids.push_back(lane_id);
        }
        return ids;
    }

    bool containsLaneId(const std::vector<LaneId>& ids, const LaneId& id) const {
        return std::find(ids.begin(), ids.end(), id) != ids.end();
    }

    void addUniqueLaneId(std::vector<LaneId>& ids, const LaneId& id) const {
        if (!id.empty() && !containsLaneId(ids, id))
            ids.push_back(id);
    }

    struct GroupCutPoint {
        Vec3d pt{0, 0, 0};
        Vec2d shift_dir{1, 0};
        double lateral = 0.0;
        int order = 0;
        int sequence = 0;
        bool has_order = false;
        AreaSourceKind source_kind = AreaSourceKind::LaneEdge;
    };

    struct CutSourceLine {
        bool is_lane = false;
        AreaSourceKind source_kind = AreaSourceKind::LaneEdge;
        LaneId lane_id;
        LaneEdgeId edge_id;
        std::vector<Vec3d> pts;
        int order = 0;
        bool has_order = false;
    };

    struct CutLineBoundaryProbe {
        std::vector<Vec3d> line;
        std::vector<LaneEdgeId> boundary_ids;
        Vec3d prefer{0, 0, 0};
    };

    struct LaneGroupCutLine {
        LaneGroupId groupid;
        std::vector<Vec3d> line;
        std::vector<AreaSourceKind> point_sources;
        std::vector<CutLineBoundaryProbe> boundary_probes;
    };

    struct LineCutIntersection {
        bool found = false;
        Vec3d pt{0, 0, 0};
        double anchor_distance = 1e18;
        double base_distance = 1e18;
    };

    // 通过mode获取内缩或外扩方向
    Vec2d applyModeToGroupShiftDir(GroupRole role, Vec2d dir) const {
        if (dir.norm() < 1e-8) dir = Vec2d(1, 0);
        else dir.normalize();
        if (role == GroupRole::Entry && mode == 2)
            dir = -dir;
        return dir;
    }

    bool shouldFlipEndpointDirToCenter(Vec2d dir, Vec2d wanted) const {
        if (dir.norm() < 1e-8 || wanted.norm() < 1e-8)
            return false;
        dir.normalize();
        wanted.normalize();
        // 仅在端点切向与中心侧预期方向明显相反时纠正。斜向道路的端点
        // 切向常接近垂直于中心径向，不能因微小负点积把mode=2外扩翻成内缩。
        return dir.dot(wanted) < -0.35;
    }

    // 统一车道/边线在组端侧的连接方向：进入组指向路口内，退出组指向路口外。
    // 先取连接端切向，再仅在明显反向时用路口中心校正方向，避免斜向道路被误翻转。
    template <typename P>
    Vec2d directedConnDir(
        const std::vector<P>& pts, GroupRole role, const Vec2d& center) const {
        bool is_entry = (role == GroupRole::Entry);
        Vec2d p = getConnPoint(pts, is_entry);
        Vec2d d = getConnTangent(pts, is_entry);
        Vec2d wanted = is_entry ? (center - p) : (p - center);
        if (shouldFlipEndpointDirToCenter(d, wanted))
            d = -d;
        if (d.norm() < 1e-8)
            return Vec2d(1, 0);
        return d.normalized();
    }

    template <typename P>
    Vec2d centerSideEndpointDir(
        const std::vector<P>& pts, GroupRole role, const Vec2d& center) const {
        if (pts.size() < 2)
            return Vec2d(1, 0);
        bool use_front = dist(pts.front(), center) <= dist(pts.back(), center);
        Vec2d p = use_front ? xyOf(pts.front()) : xyOf(pts.back());
        Vec2d d = use_front
            ? (xyOf(pts[1]) - xyOf(pts.front()))
            : (xyOf(pts.back()) - xyOf(pts[pts.size() - 2]));
        Vec2d wanted = (role == GroupRole::Entry) ? (center - p) : (p - center);
        if (shouldFlipEndpointDirToCenter(d, wanted))
            d = -d;
        if (d.norm() < 1e-8)
            return Vec2d(1, 0);
        return d.normalized();
    }

    Vec2d laneEdgeVectorDir(const std::vector<Vec3d>& pts) const {
        if (pts.size() < 2)
            return Vec2d(1, 0);
        Vec2d d = xyOf(pts.back()) - xyOf(pts.front());
        if (d.norm() < 1e-8)
            return Vec2d(1, 0);
        return d.normalized();
    }

    bool laneEdgeOpposesGroupDirection(
        const IntersectionInput& inp, const LaneGroup& group,
        const LaneEdgeId& edge_id, const std::vector<Vec3d>& edge_pts,
        const Vec2d& center) const {
        Vec2d edge_dir = laneEdgeVectorDir(edge_pts);
        Vec2d group_dir = laneGroupDirectionForEdge(inp, group, edge_id, center);
        if (edge_dir.norm() < 1e-8 || group_dir.norm() < 1e-8)
            return false;
        return edge_dir.normalized().dot(group_dir.normalized()) < -0.35;
    }

    template <typename P>
    Vec3d laneEdgeCutAnchor(
        const std::vector<P>& pts, GroupRole role, bool edge_opposes_group) const {
        if (pts.empty())
            return Vec3d(0, 0, 0);
        if (pts.size() == 1)
            return xyzOf(pts.front());
        bool use_entry_endpoint = (role == GroupRole::Entry);
        if (edge_opposes_group)
            use_entry_endpoint = !use_entry_endpoint;
        return getConnPoint3d(pts, use_entry_endpoint);
    }

    template <typename P>
    Vec2d laneEdgeCutDir(
        const std::vector<P>& pts, GroupRole role, const Vec2d& center,
        bool edge_opposes_group, const Vec2d& group_dir) const {
        if (pts.size() < 2)
            return Vec2d(1, 0);
        if (edge_opposes_group) {
            bool use_entry_endpoint = (role == GroupRole::Entry);
            use_entry_endpoint = !use_entry_endpoint;
            Vec2d d = -getConnTangent(pts, use_entry_endpoint);
            if (group_dir.norm() > 1e-8 && d.dot(group_dir) < 0.0)
                d = -d;
            if (d.norm() < 1e-8)
                return group_dir.norm() > 1e-8 ? group_dir.normalized() : Vec2d(1, 0);
            return d.normalized();
        }
        if (mode == 2 && role == GroupRole::Entry) {
            Vec2d d = getConnTangent(pts, true);
            if (d.norm() < 1e-8)
                return Vec2d(1, 0);
            return d.normalized();
        }
        return directedConnDir(pts, role, center);
    }

    template <typename P>
    Vec2d laneEdgeCutDir(
        const std::vector<P>& pts, GroupRole role, const Vec2d& center) const {
        if (pts.size() < 2)
            return Vec2d(1, 0);
        if (mode == 2 && role == GroupRole::Entry) {
            Vec2d d = getConnTangent(pts, true);
            if (d.norm() < 1e-8)
                return Vec2d(1, 0);
            return d.normalized();
        }
        return directedConnDir(pts, role, center);
    }

    bool laneEdgeReferencedByGroup(
        const IntersectionInput& inp, const LaneGroup& group, const LaneEdgeId& edge_id) const {
        if (std::find(group.boundaries.begin(), group.boundaries.end(), edge_id) != group.boundaries.end())
            return true;
        for (const auto& lane_id : uniqueGroupLaneIds(group)) {
            const Lane* lane = inp.findLane(lane_id);
            if (lane && (lane->left_edge_id == edge_id || lane->right_edge_id == edge_id))
                return true;
        }
        return false;
    }

    bool laneEdgeReferencedByRole(
        const IntersectionInput& inp, const LaneEdgeId& edge_id, GroupRole role) const {
        for (const auto& group : inp.lane_groups) {
            if (group.role == role && laneEdgeReferencedByGroup(inp, group, edge_id))
                return true;
        }
        return false;
    }

    bool laneEdgeSharedByEntryAndExit(
        const IntersectionInput& inp, const LaneEdgeId& edge_id) const {
        return laneEdgeReferencedByRole(inp, edge_id, GroupRole::Entry) &&
               laneEdgeReferencedByRole(inp, edge_id, GroupRole::Exit);
    }

    template <typename P>
    Vec2d rawLaneEdgeDirForRole(const std::vector<P>& pts, GroupRole role) const {
        Vec2d d = getConnTangent(pts, role == GroupRole::Entry);
        if (d.norm() < 1e-8)
            return Vec2d(1, 0);
        return d.normalized();
    }

    Vec2d laneGroupDirectionForEdge(
        const IntersectionInput& inp, const LaneGroup& group,
        const LaneEdgeId& edge_id, const Vec2d& center) const {
        Vec2d dir(0, 0);
        int count = 0;
        for (const auto& lane_id : uniqueGroupLaneIds(group)) {
            const Lane* lane = inp.findLane(lane_id);
            if (!lane || (lane->left_edge_id != edge_id && lane->right_edge_id != edge_id))
                continue;
            dir += directedConnDir(lane->geometry.points, group.role, center);
            ++count;
        }
        if (count == 0) {
            for (const auto& lane_id : uniqueGroupLaneIds(group)) {
                const Lane* lane = inp.findLane(lane_id);
                if (!lane)
                    continue;
                dir += directedConnDir(lane->geometry.points, group.role, center);
                ++count;
            }
        }
        if (dir.norm() < 1e-8)
            return Vec2d(1, 0);
        return dir.normalized();
    }

    bool sharedLaneEdgeMatchesGroupDirection(
        const IntersectionInput& inp, const LaneGroup& group,
        const LaneEdgeId& edge_id, const std::vector<Vec3d>& edge_pts,
        const Vec2d& center) const {
        if (!laneEdgeSharedByEntryAndExit(inp, edge_id))
            return true;
        Vec2d edge_dir = rawLaneEdgeDirForRole(edge_pts, group.role);
        Vec2d group_dir = laneGroupDirectionForEdge(inp, group, edge_id, center);
        return edge_dir.dot(group_dir) > -0.15;
    }

    Vec2d laneGroupTrafficDir(
        const IntersectionInput& inp, const LaneGroup& group, const Vec2d& center) const {
        Vec2d dir(0, 0);
        int count = 0;
        for (const auto& lane_id : uniqueGroupLaneIds(group)) {
            const Lane* lane = inp.findLane(lane_id);
            if (!lane || lane->geometry.points.size() < 2)
                continue;
            dir += directedConnDir(lane->geometry.points, group.role, center);
            ++count;
        }
        if (count == 0) {
            for (const auto& edge_id : laneGroupEdgeIds(inp, group)) {
                const LaneEdge* edge = inp.findEdge(edge_id);
                if (!edge || edge->geometry.points.size() < 2)
                    continue;
                dir += laneEdgeCutDir(edge->geometry.points, group.role, center);
                ++count;
            }
        }
        if (dir.norm() < 1e-8)
            return Vec2d(1, 0);
        return dir.normalized();
    }

    // 收集车道边线 id 时去重，避免 group.boundaries 与车道左右边线重复。
    void addUniqueLaneEdgeId(std::vector<LaneEdgeId>& ids, const LaneEdgeId& id) const {
        if (!id.empty() && std::find(ids.begin(), ids.end(), id) == ids.end())
            ids.push_back(id);
    }

    std::vector<LaneGroup> effectiveLaneGroups(
        const IntersectionInput& inp, const std::vector<ConnectivityCurve>& centerlines) const {
        std::vector<LaneGroup> groups = inp.lane_groups;
        std::map<LaneGroupId, size_t> index;
        for (size_t i = 0; i < groups.size(); ++i)
            index[groups[i].id] = i;

        auto ensure_group = [&](const LaneGroupId& id, GroupRole role) -> LaneGroup& {
            auto it = index.find(id);
            if (it != index.end())
                return groups[it->second];
            LaneGroup group;
            group.id = id;
            group.role = role;
            groups.push_back(group);
            index[id] = groups.size() - 1;
            return groups.back();
        };

        auto add_lane_to_group = [&](const LaneGroupId& id, GroupRole role, const LaneId& lane_id) {
            if (id.empty() || lane_id.empty())
                return;
            LaneGroup& group = ensure_group(id, role);
            addUniqueLaneId(group.lanes, lane_id);
            const Lane* lane = inp.findLane(lane_id);
            if (!lane)
                return;
            addUniqueLaneEdgeId(group.boundaries, lane->left_edge_id);
            addUniqueLaneEdgeId(group.boundaries, lane->right_edge_id);
        };

        for (const auto& conn : inp.connectivities) {
            add_lane_to_group(conn.enterGroupId, GroupRole::Entry, conn.entry_lane_id);
            add_lane_to_group(conn.exitGroupId, GroupRole::Exit, conn.exit_lane_id);
        }

        for (const auto& cc : centerlines) {
            const Connectivity* conn = nullptr;
            for (const auto& c : inp.connectivities) {
                if (c.id == cc.id) {
                    conn = &c;
                    break;
                }
            }
            if (conn) {
                add_lane_to_group(conn->enterGroupId, GroupRole::Entry, cc.entry_lane_id);
                add_lane_to_group(conn->exitGroupId, GroupRole::Exit, cc.exit_lane_id);
                continue;
            }
            const Lane* entry_lane = inp.findLane(cc.entry_lane_id);
            const Lane* exit_lane = inp.findLane(cc.exit_lane_id);
            if (entry_lane)
                add_lane_to_group(entry_lane->groupId, GroupRole::Entry, cc.entry_lane_id);
            if (exit_lane)
                add_lane_to_group(exit_lane->groupId, GroupRole::Exit, cc.exit_lane_id);
        }

        return groups;
    }

    // 提取组内参与端侧边界的车道边线：优先组边界，再补充组内车道左右边线。
    std::vector<LaneEdgeId> laneGroupEdgeIds(
        const IntersectionInput& inp, const LaneGroup& group) const {
        std::vector<LaneEdgeId> ids;
        for (const auto& edge_id : group.boundaries)
            addUniqueLaneEdgeId(ids, edge_id);
        for (const auto& lane_id : uniqueGroupLaneIds(group)) {
            const Lane* lane = inp.findLane(lane_id);
            if (!lane) continue;
            addUniqueLaneEdgeId(ids, lane->left_edge_id);
            addUniqueLaneEdgeId(ids, lane->right_edge_id);
        }
        return ids;
    }

    // 将组端侧形点按容差去重后加入候选集；只在同一数据类别内过滤重复，
    // lane、laneedge、boundary 属于不同构面源，不能互相吞掉。
    void addGroupCutPoint(std::vector<GroupCutPoint>& pts, GroupCutPoint p, double tol = 0.03) const {
        for (const auto& existing : pts) {
            if (existing.source_kind != p.source_kind)
                continue;
            if (dist(existing.pt, p.pt) <= tol)
                return;
        }
        p.sequence = (int)pts.size();
        pts.push_back(p);
    }

    bool sameCutSourceGeometry(
        const std::vector<Vec3d>& a,
        const std::vector<Vec3d>& b,
        double tol = 0.03) const {
        if (a.size() < 2 || b.size() < 2)
            return false;
        if (std::abs(polylineLength(a) - polylineLength(b)) > tol)
            return false;
        for (const auto& pt : a) {
            if (pointToPolyline(pt, b) > tol)
                return false;
        }
        for (const auto& pt : b) {
            if (pointToPolyline(pt, a) > tol)
                return false;
        }
        return true;
    }

    void addUniqueCutSourceLine(
        std::vector<CutSourceLine>& sources,
        CutSourceLine source) const {
        if (source.pts.size() < 2)
            return;
        for (auto& existing : sources) {
            if (existing.source_kind != source.source_kind)
                continue;
            if (!sameCutSourceGeometry(existing.pts, source.pts))
                continue;
            return;
        }
        sources.push_back(std::move(source));
    }

    std::vector<CutSourceLine> laneGroupCutSourceLines(
        const IntersectionInput& inp,
        const LaneGroup& group,
        const Vec2d& center) const {
        std::vector<CutSourceLine> sources;
        for (const auto& edge_id : laneGroupEdgeIds(inp, group)) {
            const LaneEdge* edge = inp.findEdge(edge_id);
            if (!edge || edge->geometry.points.size() < 2)
                continue;
            if (!sharedLaneEdgeMatchesGroupDirection(
                    inp, group, edge_id, edge->geometry.points, center))
                continue;
            CutSourceLine source;
            source.is_lane = false;
            source.source_kind = AreaSourceKind::LaneEdge;
            source.edge_id = edge_id;
            source.pts = edge->geometry.points;
            source.order = edge->lineOrder;
            source.has_order = true;
            addUniqueCutSourceLine(sources, std::move(source));
        }
        for (const auto& lane_id : uniqueGroupLaneIds(group)) {
            const Lane* lane = inp.findLane(lane_id);
            if (!lane || lane->geometry.points.size() < 2)
                continue;
            CutSourceLine source;
            source.is_lane = true;
            source.source_kind = AreaSourceKind::Lane;
            source.lane_id = lane_id;
            source.pts = lane->geometry.points;
            source.order = lane->laneOrder * 2 + 1;
            source.has_order = true;
            addUniqueCutSourceLine(sources, std::move(source));
        }
        return sources;
    }

    // 计算组端侧折线的缩扩方向：汇总边线、中心线和生成线簇端侧切向，
    // 失败时用组端点到中心兜底。mode=2时进入组取反为路口外扩，其他模式进入组向路口内缩。
    Vec2d laneGroupShiftDir(
        const IntersectionInput& inp,
        const LaneGroup& group,
        const Vec2d& center,
        const std::vector<ConnectivityCurve>* centerlines = nullptr) const {
        Vec2d dir(0, 0);
        int count = 0;
        if (centerlines) {
            for (const auto& cc : *centerlines) {
                if (!cc.curve)
                    continue;
                if (group.role == GroupRole::Entry) {
                    const Lane* lane = inp.findLane(cc.entry_lane_id);
                    if (!lane || (lane->groupId != group.id && !groupContainsLaneId(group, cc.entry_lane_id)))
                        continue;
                    Vec2d t = cc.curve->startTan();
                    if (t.norm() > 1e-8) {
                        dir += t.normalized();
                        ++count;
                    }
                } else {
                    const Lane* lane = inp.findLane(cc.exit_lane_id);
                    if (!lane || (lane->groupId != group.id && !groupContainsLaneId(group, cc.exit_lane_id)))
                        continue;
                    Vec2d t = cc.curve->endTan();
                    if (t.norm() > 1e-8) {
                        dir += t.normalized();
                        ++count;
                    }
                }
            }
            if (count > 0 && dir.norm() > 1e-8) {
                dir.normalize();
                Vec2d mid(0, 0);
                int mid_count = 0;
                for (const auto& cc : *centerlines) {
                    if (!cc.curve)
                        continue;
                    if (group.role == GroupRole::Entry) {
                        const Lane* lane = inp.findLane(cc.entry_lane_id);
                        if (!lane || (lane->groupId != group.id && !groupContainsLaneId(group, cc.entry_lane_id)))
                            continue;
                        mid += cc.curve->startPt();
                    } else {
                        const Lane* lane = inp.findLane(cc.exit_lane_id);
                        if (!lane || (lane->groupId != group.id && !groupContainsLaneId(group, cc.exit_lane_id)))
                            continue;
                        mid += cc.curve->endPt();
                    }
                    ++mid_count;
                }
                if (mid_count > 0) {
                    mid /= (double)mid_count;
                    Vec2d wanted = (group.role == GroupRole::Entry) ? (center - mid) : (mid - center);
                    if (shouldFlipEndpointDirToCenter(dir, wanted))
                        dir = -dir;
                }
                return applyModeToGroupShiftDir(group.role, dir);
            }
        }

        for (const auto& edge_id : laneGroupEdgeIds(inp, group)) {
            const LaneEdge* edge = inp.findEdge(edge_id);
            if (!edge || edge->geometry.points.size() < 2)
                continue;
            dir += directedConnDir(edge->geometry.points, group.role, center);
            ++count;
        }
        for (const auto& lane_id : uniqueGroupLaneIds(group)) {
            const Lane* lane = inp.findLane(lane_id);
            if (!lane || lane->geometry.points.size() < 2)
                continue;
            dir += directedConnDir(lane->geometry.points, group.role, center);
            ++count;
        }
        if (count > 0 && dir.norm() > 1e-8)
            return applyModeToGroupShiftDir(group.role, dir);

        Vec2d mid(0, 0);
        count = 0;
        for (const auto& lane_id : uniqueGroupLaneIds(group)) {
            const Lane* lane = inp.findLane(lane_id);
            if (!lane || lane->geometry.points.empty())
                continue;
            mid += getConnPoint(lane->geometry.points, group.role == GroupRole::Entry);
            ++count;
        }
        if (count > 0) {
            mid /= (double)count;
            dir = (group.role == GroupRole::Entry) ? (center - mid) : (mid - center);
            if (dir.norm() > 1e-8)
                return applyModeToGroupShiftDir(group.role, dir);
        }
        return Vec2d(1, 0);
    }

    Vec3d shiftedGroupCutPoint(const GroupCutPoint& p) const {
        Vec2d dir = p.shift_dir.norm() > 1e-8 ? p.shift_dir.normalized() : Vec2d(1, 0);
        return offsetXY(p.pt, groupCutOffset * dir);
    }

    CutLineBoundaryProbe endpointCutBoundaryProbe(
        const Vec3d& anchor,
        const Vec2d& base_dir,
        GroupRole role,
        std::vector<LaneEdgeId> boundary_ids) const {
        CutLineBoundaryProbe probe;
        probe.boundary_ids = std::move(boundary_ids);
        Vec2d tangent = base_dir.norm() > 1e-8 ? base_dir.normalized() : Vec2d(1, 0);
        Vec3d cut_base = offsetXY(anchor, groupCutOffset * applyModeToGroupShiftDir(role, tangent));
        Vec2d cut_dir = rotLeft(tangent);
        if (cut_dir.norm() < 1e-8)
            return probe;
        cut_dir.normalize();
        double len = std::max(cutLineExtension, 3.5);
        probe.prefer = cut_base;
        probe.line = {offsetXY(cut_base, -len * cut_dir),
                      offsetXY(cut_base,  len * cut_dir)};
        return probe;
    }

    Vec3d connectivityEndpoint3d(const ConnectivityCurve& cc, bool entry) const {
        if (!cc.geometry.points.empty())
            return entry ? cc.geometry.points.front() : cc.geometry.points.back();
        if (cc.curve)
            return Vec3d(entry ? cc.curve->startPt() : cc.curve->endPt(), 0.0);
        return Vec3d(0, 0, 0);
    }

    // 由已生成的连通曲线端点构建端侧边界折线：
    // 每个缩扩点沿自身进入/退出线切向移动，保证点落在原线或其延长线上。
    std::vector<Vec3d> getLaneGroupCurveCutLine(
        const IntersectionInput& inp, const LaneGroup& group,
        const std::vector<ConnectivityCurve>& centerlines, const Vec2d& center) const {

        auto curveEndpointShiftDir = [&](
                const ConnectivityCurve& cc, const Lane* lane,
                const LaneGroup& group, const Vec2d& center) -> Vec2d {
            const bool is_entry_group = group.role == GroupRole::Entry;
            Vec2d dir = is_entry_group ? cc.curve->startTan() : cc.curve->endTan();
            Vec2d conn_dir = lane ? directedConnDir(lane->geometry.points, group.role, center)
                                  : ((group.role == GroupRole::Entry)
                                     ? (center - cc.curve->startPt()) : (cc.curve->endPt() - center));
            if (dir.norm() < 1e-8)
                dir = conn_dir;
            else {
                dir.normalize();
                if (conn_dir.norm() > 1e-8 && dir.dot(conn_dir) < 0.0)
                    dir = -dir;
            }
            return applyModeToGroupShiftDir(group.role, dir);
        };

        std::vector<GroupCutPoint> cut_pts;
        Vec2d shift_dir = laneGroupShiftDir(inp, group, center, &centerlines);

        for (const auto& cc : centerlines) {
            if (!cc.curve)
                continue;
            const bool is_entry_group = group.role == GroupRole::Entry;
            const Lane* lane = inp.findLane(is_entry_group ? cc.entry_lane_id : cc.exit_lane_id);
            const LaneId& lane_id = is_entry_group ? cc.entry_lane_id : cc.exit_lane_id;
            if (!lane || (lane->groupId != group.id && !groupContainsLaneId(group, lane_id)))
                continue;

            GroupCutPoint cp;
            cp.pt = connectivityEndpoint3d(cc, is_entry_group);
            cp.shift_dir = curveEndpointShiftDir(cc, lane, group, center);
            cp.order = lane->laneOrder * 2 + 1;
            cp.has_order = true;
            cp.source_kind = AreaSourceKind::Lane;
            addGroupCutPoint(cut_pts, cp);
        }

        if (cut_pts.size() < 2)
            return {};

        Vec2d lateral_dir = rotLeft(shift_dir);
        for (auto& p : cut_pts)
            p.lateral = p.pt.dot(lateral_dir);

        std::sort(cut_pts.begin(), cut_pts.end(), [](const GroupCutPoint& a, const GroupCutPoint& b) {
            if (a.has_order != b.has_order)
                return a.has_order > b.has_order;
            if (a.has_order && a.order != b.order)
                return a.order < b.order;
            if (std::abs(a.lateral - b.lateral) > 1e-9)
                return a.lateral < b.lateral;
            return a.sequence < b.sequence;
        });

        std::vector<Vec3d> out;
        for (const auto& p : cut_pts)
            pushUnique(out, shiftedGroupCutPoint(p), 0.03);
        return out;
    }

    // 构建进入/退出组的端侧边界折线：逐条收集组内车道边线端点和中心线端点，按线序/车道序连接成有序折线。
    // 每个中心线缩扩点由端点垂线平移后与共端点连接线求交得到，保证形点落在实际连接线上。
    LaneGroupCutLine getLaneGroupCutLine(
        const IntersectionInput& inp,
        const LaneGroup& group,
        const std::vector<ConnectivityCurve>& centerlines,
        const Vec2d& center) const {

        auto shiftedGroupCutPointByConnections = [&](
                const Vec3d &anchor, const Vec2d &base_dir, GroupRole role,
                const std::vector<std::vector<Vec3d>> &connection_lines) -> Vec3d {

            auto nearestIntersectionOnConnectionLine = [&](
                    const std::vector<Vec3d> &connection, const Vec3d &cut_base,
                    const Vec2d &cut_dir, const Vec3d &anchor) -> LineCutIntersection {

                auto lineSegmentIntersection = [&](
                        const Vec3d& line_pt, const Vec2d& line_dir,
                        const Vec3d& a, const Vec3d& b, const Vec3d& anchor, Vec3d* out) -> bool {
                    Vec2d seg = xyOf(b) - xyOf(a);
                    if (line_dir.norm() < 1e-9 || seg.norm() < 1e-9)
                        return false;
                    double den = cross2d(seg, line_dir);
                    if (std::abs(den) < 1e-12) {
                        if (std::abs(cross2d(xyOf(a) - xyOf(line_pt), line_dir)) > 1e-8)
                            return false;
                        *out = dist(a, anchor) <= dist(b, anchor) ? a : b;
                        return true;
                    }
                    double t = cross2d(xyOf(line_pt) - xyOf(a), line_dir) / den;
                    if (t < -1e-9 || t > 1.0 + 1e-9)
                        return false;
                    t = std::max(0.0, std::min(1.0, t));
                    *out = a + t * (b - a);
                    return true;
                };

                LineCutIntersection best;
                if (connection.size() < 2 || cut_dir.norm() < 1e-9)
                    return best;
                for (int i = 0; i + 1 < (int) connection.size(); ++i) {
                    Vec3d hit;
                    if (!lineSegmentIntersection(cut_base, cut_dir, connection[i], connection[i + 1], anchor, &hit))
                        continue;
                    double ad = dist(hit, anchor);
                    double bd = dist(hit, cut_base);
                    if (!best.found ||
                        ad < best.anchor_distance - 1e-8 ||
                        (std::abs(ad - best.anchor_distance) <= 1e-8 && bd < best.base_distance)) {
                        best.found = true;
                        best.pt = hit;
                        best.anchor_distance = ad;
                        best.base_distance = bd;
                    }
                }
                return best;
            };

            Vec2d tangent = base_dir.norm() > 1e-8 ? base_dir.normalized() : Vec2d(1, 0);
            Vec2d shift_dir = applyModeToGroupShiftDir(role, tangent);
            Vec3d cut_base = offsetXY(anchor, groupCutOffset * shift_dir);
            Vec2d cut_dir = rotLeft(tangent);
            if (cut_dir.norm() < 1e-8)
                return cut_base;
            cut_dir.normalize();

            LineCutIntersection best;
            for (const auto &line : connection_lines) {
                auto hit = nearestIntersectionOnConnectionLine(line, cut_base, cut_dir, anchor);
                if (!hit.found)
                    continue;
                if (!best.found ||
                    hit.anchor_distance < best.anchor_distance - 1e-8 ||
                    (std::abs(hit.anchor_distance - best.anchor_distance) <= 1e-8 &&
                     hit.base_distance < best.base_distance)) {
                    best = hit;
                }
            }
            return best.found ? best.pt : cut_base;
        };

        auto laneEndpointConnectionLines = [&](
                const LaneId& lane_id, GroupRole role, const std::vector<Vec3d>& lane_pts,
                const std::vector<ConnectivityCurve>& centerlines) -> std::vector<std::vector<Vec3d>> {
            std::vector<std::vector<Vec3d>> lines;
            if (lane_pts.size() >= 2)
                lines.push_back(lane_pts);
            for (const auto& cc : centerlines) {
                if ((role == GroupRole::Entry && cc.entry_lane_id != lane_id) ||
                    (role == GroupRole::Exit && cc.exit_lane_id != lane_id))
                    continue;
                std::vector<Vec3d> pts;
                if (!cc.geometry.points.empty())
                    pts = cc.geometry.points;
                else if (cc.curve)
                    pts = toVec3dArray(cc.curve->sampleByArcLength(160));
                if (pts.size() >= 2)
                    lines.push_back(pts);
            }
            return lines;
        };

        LaneGroupCutLine result;
        result.groupid = group.id;
        std::vector<GroupCutPoint> cut_pts;
        Vec2d shift_dir = laneGroupShiftDir(inp, group, center);

        auto laneHasRoleCenterline = [&](
                const LaneId& lane_id, GroupRole role,
                const std::vector<ConnectivityCurve>& centerlines) -> bool {
            for (const auto& cc : centerlines) {
                if ((role == GroupRole::Entry && cc.entry_lane_id == lane_id) ||
                    (role == GroupRole::Exit && cc.exit_lane_id == lane_id))
                    return true;
            }
            return false;
        };

        for (const auto& source : laneGroupCutSourceLines(inp, group, center)) {
            GroupCutPoint cp;
            Vec3d anchor;
            Vec2d base_dir;
            std::vector<std::vector<Vec3d>> connection_lines;
            if (source.is_lane) {
                anchor = getConnPoint3d(source.pts, group.role == GroupRole::Entry);
                base_dir = directedConnDir(source.pts, group.role, center);
                connection_lines = laneEndpointConnectionLines(
                    source.lane_id, group.role, source.pts, centerlines);
            } else {
                // LaneEdge 通常与 lane 一样按组角色取进入尾端/退出首端；若其
                // 矢量化方向与当前组交通方向相反，则翻转端点语义，避免把远端
                // 边线点拉入组端侧缩扩线。
                Vec2d group_dir = laneGroupDirectionForEdge(
                    inp, group, source.edge_id, center);
                bool edge_opposes_group = laneEdgeOpposesGroupDirection(
                    inp, group, source.edge_id, source.pts, center);
                anchor = laneEdgeCutAnchor(source.pts, group.role, edge_opposes_group);
                base_dir = laneEdgeCutDir(
                    source.pts, group.role, center, edge_opposes_group, group_dir);
                connection_lines = {source.pts};
            }
            cp.pt = shiftedGroupCutPointByConnections(
                anchor, base_dir, group.role, connection_lines);
            cp.shift_dir = applyModeToGroupShiftDir(group.role, base_dir);
            cp.order = source.order;
            cp.has_order = source.has_order;
            cp.source_kind = source.source_kind;
            addGroupCutPoint(cut_pts, cp);

            if (source.is_lane &&
                !centerlines.empty() &&
                group.role == GroupRole::Exit &&
                !laneHasRoleCenterline(source.lane_id, group.role, centerlines)) {
                const Lane* lane = inp.findLane(source.lane_id);
                if (!lane) continue;
                CutLineBoundaryProbe probe = endpointCutBoundaryProbe(
                    anchor, base_dir, group.role,
                    std::vector<LaneEdgeId>{lane->left_edge_id, lane->right_edge_id});
                if (probe.line.size() >= 2)
                    result.boundary_probes.push_back(std::move(probe));
            }
        }

        if (cut_pts.size() < 2)
            return result;
        Vec2d lateral_dir = rotLeft(shift_dir);
        for (auto& p : cut_pts)
            p.lateral = p.pt.dot(lateral_dir);

        std::sort(cut_pts.begin(), cut_pts.end(), [](const GroupCutPoint& a, const GroupCutPoint& b) {
            if (a.has_order != b.has_order)
                return a.has_order > b.has_order;
            if (a.has_order && a.order != b.order)
                return a.order < b.order;
            if (std::abs(a.lateral - b.lateral) > 1e-9)
                return a.lateral < b.lateral;
            return a.sequence < b.sequence;
        });

        for (const auto& p : cut_pts) {
            if (!result.line.empty() &&
                result.point_sources.back() == p.source_kind &&
                dist(result.line.back(), p.pt) <= 0.03)
                continue;
            result.line.push_back(p.pt);
            result.point_sources.push_back(p.source_kind);
        }
        return result;
    }

    struct BoundaryHit {
        bool found = false;
        int boundary_index = -1;
        Vec3d pt{0, 0, 0};
        double station = 0.0;
        double cut_station = 0.0;
        double distance = 1e18;
        Vec2d boundary_dir{1, 0};
        Vec2d chain_dir{1, 0};
        bool has_boundarynode = true; //附近是否有Boundary的端点
    };

    struct SelectedBoundaryHit {
        BoundaryHit hit;
        bool from_first_endpoint = true;
    };

    struct BoundaryHitStation {
        double station = 0.0;
        std::vector<int> cut_indices;
        std::vector<GroupRole> cut_roles;
        Vec3d pt{0, 0, 0};//额外记录位置
        Vec2d boundary_dir{1, 0};
        Vec2d chain_dir{1, 0};
    };

    struct AlignedAreaParts {
        std::vector<std::vector<Vec3d>> cut_lines;
        std::vector<std::vector<Vec3d>> boundary_lines;
    };

    struct EdgePart {
        std::vector<Vec3d> pts;
        int node0 = -1;
        int node1 = -1;
        bool is_artificial_closure = false;
    };

    struct EdgeNode {
        Vec3d pt{0, 0, 0};
        std::vector<std::pair<int, bool>> incident;
    };

    // 对齐端侧边界折线与道路边缘：端侧折线两端外延求交，并同步截取对应道路边缘。
    // 同一道路边缘只保留最外侧到最内侧命中的连续有效段；完整生成已有连通曲线时，
    // 退出组内未生成连通曲线的退出线还用首端缩扩探测线补充对应RoadEdge的station，
    // 避免单命中时把整条边缘纳入路口面。
    AlignedAreaParts alignCutLinesAndBoundaries(
        const std::vector<BoundaryLine>& boundary_lines,
        std::vector<std::vector<Vec3d>> cut_lines,
        const std::vector<GroupRole>& cut_roles,
        const std::vector<LaneGroupCutLine>& cut_meta,
        const std::vector<Vec2d>& cut_group_dirs,
        const std::vector<std::vector<LaneEdgeId>>& cut_group_edge_ids,
        const Vec2d& center,
        const std::vector<LaneGroupId> group_ids) const {

        // 将端侧边界折线两端按接口参数外延，用外延段与道路边缘求交来消除未对齐尖角。仅延长首尾端点，中间形点保持原折线形态。
        auto extendCutLineEnds = [&](const std::vector<Vec3d>& line, double extension_len) -> std::vector<Vec3d> {
            std::vector<Vec3d> out = line;
            if (out.size() < 2)
                return out;
            double len = std::max(cutLineExtension, extension_len);
            Vec2d first_dir = xyOf(out.front()) - xyOf(out[1]);
            if (first_dir.norm() > 1e-8)
                out.front() = offsetXY(out.front(), len * first_dir.normalized());
            Vec2d last_dir = xyOf(out.back()) - xyOf(out[out.size() - 2]);
            if (last_dir.norm() > 1e-8)
                out.back() = offsetXY(out.back(), len * last_dir.normalized());
            return out;
        };

        // 查找外延探测线与每条道路边缘线的最近交点：同一条道路边缘有多个交点时，
        // 取距离原始未延长缩扩线端点 prefer 最近的点。命中点是否靠近该道路边缘
        // 原始端点会被记录下来；后续只允许端点邻近命中写入道路边缘截取 station，
        // 避免其他组的 RoadEdge 被当前组的缩扩延长线误归属。道路边缘截取 station
        // 另行按局部 RoadEdge 方向与组方向过滤，避免掉头形边缘被反向交点打穿。
        auto boundaryHitsByBoundary = [&](LaneGroupId groupid, const std::vector<Vec3d>& probe,
                const std::vector<BoundaryLine>& boundaries, const Vec3d& prefer,
                const Vec2d& group_dir,
                bool require_direction_match = false,
                double cut_station_offset = 0.0) -> std::vector<BoundaryHit> {
            auto hasBoundarynodeNearHit = [&](const BoundaryLine& boundary, const Vec3d& hit) {
                if (boundary.source_endpoints.empty())
                    return true;
                bool ret = false;
                double endpoint_tol = groupCutOffset + 0.5;
                for (const auto& eps : boundary.source_endpoints) {
                    double sdist = dist(hit, eps.first);
                    double edist = dist(hit, eps.second);
                    if (sdist <= endpoint_tol || edist <= endpoint_tol) {
                        ret = true;
                        break;
                    }
                }
                return ret;
            };

            // 获取两折线交点
            auto polylineIntersections = [&](
                    const std::vector<Vec3d>& a, const BoundaryLine& b) -> std::vector<LineIntersection> {
                // 获取线段交点
                auto segmentIntersection = [&](
                        const Vec3d& a, const Vec3d& b, const Vec3d& c, const Vec3d& d,
                        Vec3d* out = nullptr, double* ta = nullptr, double* tb = nullptr) -> bool {
                    Vec2d r = xyOf(b) - xyOf(a);
                    Vec2d s = xyOf(d) - xyOf(c);
                    double den = cross2d(r, s);
                    if (std::abs(den) < 1e-12)
                        return false;
                    Vec2d ac = xyOf(c) - xyOf(a);
                    double t = cross2d(ac, s) / den;
                    double u = cross2d(ac, r) / den;
                    if (t < -1e-9 || t > 1.0 + 1e-9 || u < -1e-9 || u > 1.0 + 1e-9)
                        return false;
                    t = std::max(0.0, std::min(1.0, t));
                    u = std::max(0.0, std::min(1.0, u));
                    if (out) *out = a + t * (b - a);
                    if (ta) *ta = t;
                    if (tb) *tb = u;
                    return true;
                };
                std::vector<LineIntersection> out;
                double acc_a = 0.0;
                for (int i = 0; i + 1 < (int)a.size(); ++i) {
                    double seg_a = dist(a[i], a[i + 1]);
                    double acc_b = 0.0;
                    for (int j = 0; j + 1 < (int)b.pts.size(); ++j) {
                        double seg_b = dist(b.pts[j], b.pts[j + 1]);
                        Vec3d pt;
                        double ta = 0.0, tb = 0.0;
                        if (segmentIntersection(a[i], a[i + 1], b.pts[j], b.pts[j + 1], &pt, &ta, &tb)) {
                            LineIntersection hit;
                            hit.pt = pt;
                            hit.station_a = acc_a + ta * seg_a;
                            hit.station_b = acc_b + tb * seg_b;
                            Vec2d dir_b = j < (int)b.segment_dirs.size()
                                ? b.segment_dirs[j]
                                : xyOf(b.pts[j + 1]) - xyOf(b.pts[j]);
                            hit.dir_b = dir_b.norm() > 1e-8 ? dir_b.normalized() : Vec2d(1, 0);
                            Vec2d chain_dir_b = xyOf(b.pts[j + 1]) - xyOf(b.pts[j]);
                            hit.chain_dir_b = chain_dir_b.norm() > 1e-8
                                ? chain_dir_b.normalized() : hit.dir_b;
                            bool duplicate = false;
                            for (const auto& existing : out) {
                                if (dist(existing.pt, hit.pt) < 0.03) {
                                    duplicate = true;
                                    break;
                                }
                            }
                            if (!duplicate)
                                out.push_back(hit);
                        }
                        acc_b += seg_b;
                    }
                    acc_a += seg_a;
                }
                return out;
            };

            std::vector<BoundaryHit> hits;
            for (int i = 0; i < (int)boundaries.size(); ++i) {
                BoundaryHit best;
                for (const auto& hit : polylineIntersections(probe, boundaries[i])) {
                    Vec2d boundary_dir = hit.dir_b.norm() > 1e-8 ? hit.dir_b.normalized() : Vec2d(1, 0);
                    Vec2d normalized_group_dir = group_dir.norm() > 1e-8 ? group_dir.normalized() : Vec2d(1, 0);
                    if (require_direction_match && boundary_dir.dot(normalized_group_dir) <= 0.15)
                        continue;
                    double d = dist(hit.pt, prefer);
                    bool has_boundarynode = hasBoundarynodeNearHit(boundaries[i], hit.pt);
                    if (!has_boundarynode) //未找到边缘端点（不属于当前组的边缘）直接跳过
                        continue;
                    if (d < best.distance) {
                        best.found = true;
                        best.boundary_index = i;
                        best.pt = hit.pt;
                        best.station = hit.station_b;
                        best.cut_station = cut_station_offset + hit.station_a;
                        best.distance = d;
                        best.boundary_dir = boundary_dir;
                        best.chain_dir = hit.chain_dir_b;
                        best.has_boundarynode = has_boundarynode;
                    }
                }
                if (best.found)
                    hits.push_back(best);
            }
            return hits;
        };

        auto nearestOverallBoundaryHit = [&](const std::vector<BoundaryHit>& hits) -> BoundaryHit {
            BoundaryHit best;
            for (const auto& hit : hits) {
                if (hit.distance < best.distance)
                    best = hit;
            }
            return best;
        };

        auto boundaryMatchesSourceIds = [&](const BoundaryLine& boundary, const std::vector<LaneEdgeId>& ids) -> bool {
            if (ids.empty())
                return true;
            for (const auto& id : ids) {
                if (std::find(boundary.source_ids.begin(), boundary.source_ids.end(), id) != boundary.source_ids.end())
                    return true;
            }
            return false;
        };

        auto boundaryHitMatchesGroupDir = [&](const BoundaryHit& hit, const Vec2d& group_dir) -> bool {
            Vec2d boundary_dir = hit.boundary_dir.norm() > 1e-8
                ? hit.boundary_dir.normalized() : Vec2d(1, 0);
            Vec2d normalized_group_dir = group_dir.norm() > 1e-8
                ? group_dir.normalized() : Vec2d(1, 0);
            return boundary_dir.dot(normalized_group_dir) > 0.15;
        };

        // 记录道路边缘被端侧边界折线命中的里程位置，后续据此截取有效道路边缘段。
        auto addBoundaryHitStation = [&](
                std::vector<std::vector<BoundaryHitStation>>& hit_stations,
                const BoundaryHit& hit, int cut_index, GroupRole cut_role) -> void {
            if (!hit.found || hit.boundary_index < 0 ||
                hit.boundary_index >= (int)hit_stations.size())
                return;
            if (!hit.has_boundarynode)
                return;
            for (auto& existing : hit_stations[hit.boundary_index]) {
                if (std::abs(existing.station - hit.station) < 0.03) {
                    if (std::find(existing.cut_indices.begin(), existing.cut_indices.end(), cut_index) ==
                        existing.cut_indices.end()) {
                        existing.cut_indices.push_back(cut_index);
                        existing.cut_roles.push_back(cut_role);
                    }
                    return;
                }
            }
            BoundaryHitStation station;
            station.station = hit.station;
            station.cut_indices.push_back(cut_index);
            station.cut_roles.push_back(cut_role);
            station.pt = hit.pt;
            station.boundary_dir = hit.boundary_dir;
            station.chain_dir = hit.chain_dir;
            hit_stations[hit.boundary_index].push_back(station);
        };

        AlignedAreaParts result;
        std::vector<std::vector<BoundaryHitStation>> hit_stations(boundary_lines.size());

        auto uniqueHitStationCount = [](const std::vector<BoundaryHitStation>& stations) -> int {
            std::vector<double> unique;
            for (const auto& hit : stations) {
                bool exists = false;
                for (double station : unique) {
                    if (std::abs(station - hit.station) <= 0.03) {
                        exists = true;
                        break;
                    }
                }
                if (!exists)
                    unique.push_back(hit.station);
            }
            return (int)unique.size();
        };

        auto addProbeHits = [&](const LaneGroupCutLine& meta, int cut_index, GroupRole cut_role) {
            if (cut_role != GroupRole::Exit)
                return;
            Vec2d group_dir = cut_index < (int)cut_group_dirs.size()
                ? cut_group_dirs[cut_index] : Vec2d(1, 0);
            LaneGroupId groupid = group_ids[cut_index];
            for (const auto& probe : meta.boundary_probes) {
                if (probe.line.size() < 2)
                    continue;
                std::vector<BoundaryHit> hits = boundaryHitsByBoundary(groupid,
                    probe.line, boundary_lines, probe.prefer, group_dir,
                    true, polylineLength(probe.line) * 0.5);
                for (const auto& hit : hits) {
                    if (hit.boundary_index < 0 ||
                        hit.boundary_index >= (int)boundary_lines.size())
                        continue;
                    if (!boundaryMatchesSourceIds(boundary_lines[hit.boundary_index], probe.boundary_ids))
                        continue;
                    if (uniqueHitStationCount(hit_stations[hit.boundary_index]) >= 2)
                        continue;
                    addBoundaryHitStation(hit_stations, hit, cut_index, cut_role);
                }
            }
        };

        result.cut_lines.reserve(cut_lines.size());
        for (auto& cut : cut_lines) {
            if (cut.size() < 2)
                continue;
            int cut_index = (int)result.cut_lines.size();
            LaneGroupId groupid = group_ids[cut_index];
            // 组切割线两端延长
            std::vector<Vec3d> extended = extendCutLineEnds(cut, cutLineExtension);

            double front_original_station = dist(extended.front(), cut.front());
            double back_original_station = polylineLength(extended) - dist(cut.back(), extended.back());
            std::vector<Vec3d> first_probe = {extended.front(), extended[1]};//切割线首端
            std::vector<Vec3d> last_probe = {extended[extended.size() - 2], extended.back()};//切割线首端
            double last_probe_offset = polylineLength(extended) - dist(extended[extended.size() - 2], extended.back());
            GroupRole cut_role = cut_index < (int)cut_roles.size() ? cut_roles[cut_index] : GroupRole::Entry;
            Vec2d group_dir = cut_index < (int)cut_group_dirs.size()
                ? cut_group_dirs[cut_index] : Vec2d(1, 0);
            std::vector<BoundaryHit> first_hits = boundaryHitsByBoundary(groupid,
                first_probe, boundary_lines, cut.front(), group_dir);
            std::vector<BoundaryHit> last_hits = boundaryHitsByBoundary(groupid,
                last_probe, boundary_lines, cut.back(), group_dir, false, last_probe_offset);
            std::vector<SelectedBoundaryHit> selected_hits;
            auto hitMatchesCurrentGroupEdge = [&](const BoundaryHit& hit) {
                if (!hit.found || hit.boundary_index < 0 ||
                    hit.boundary_index >= (int)boundary_lines.size() ||
                    cut_index >= (int)cut_group_edge_ids.size())
                    return false;
                return boundaryMatchesSourceIds(
                    boundary_lines[hit.boundary_index], cut_group_edge_ids[cut_index]);
            };
            auto selectNearestHitForBoundary = [&](const BoundaryHit& hit, bool from_first_endpoint) {
                if (!hit.found)
                    return;
                for (auto& selected : selected_hits) {
                    if (selected.hit.boundary_index == hit.boundary_index) {
                        if (hit.distance < selected.hit.distance) {
                            selected.hit = hit;
                            selected.from_first_endpoint = from_first_endpoint;
                        }
                        return;
                    }
                }
                SelectedBoundaryHit selected;
                selected.hit = hit;
                selected.from_first_endpoint = from_first_endpoint;
                selected_hits.push_back(selected);
            };
            for (const auto& hit : first_hits)
                selectNearestHitForBoundary(hit, true);
            for (const auto& hit : last_hits)
                selectNearestHitForBoundary(hit, false);

            std::vector<BoundaryHit> selected_first_hits;
            std::vector<BoundaryHit> selected_last_hits;
            for (const auto& selected : selected_hits) {
                if (selected.from_first_endpoint)
                    selected_first_hits.push_back(selected.hit);
                else
                    selected_last_hits.push_back(selected.hit);
            }
            BoundaryHit first_hit = nearestOverallBoundaryHit(selected_first_hits);
            BoundaryHit last_hit = nearestOverallBoundaryHit(selected_last_hits);
            auto endpointSource = [&](bool first_endpoint) {
                if (cut_index >= (int)cut_meta.size() || cut_meta[cut_index].point_sources.empty())
                    return AreaSourceKind::Boundary;
                const auto& sources = cut_meta[cut_index].point_sources;
                return first_endpoint ? sources.front() : sources.back();
            };
            double first_trim_station = first_hit.cut_station;
            double last_trim_station = last_hit.cut_station;
            // 当前组自己的 LaneEdge 既是缩扩点来源，又可能作为 RoadEdge
            // 链参与对齐。命中这类同源边缘时，组切线已经在正确的端点
            // 处生成，不能再被内部边缘的交点截短，否则会丢掉相邻车道
            // 的连续缩扩折线（100000536 的 43112504 即此情形）。
            if (first_hit.found && first_hit.distance > 1.0 &&
                (endpointSource(true) == AreaSourceKind::LaneEdge ||
                 hitMatchesCurrentGroupEdge(first_hit)))
                first_trim_station = front_original_station;
            if (last_hit.found && last_hit.distance > 1.0 &&
                (endpointSource(false) == AreaSourceKind::LaneEdge ||
                 hitMatchesCurrentGroupEdge(last_hit)))
                last_trim_station = back_original_station;
            if (first_hit.found && last_hit.found) {
                double s0 = std::min(first_trim_station, last_trim_station);
                double s1 = std::max(first_trim_station, last_trim_station);
                cut = subPolylineByStation(extended, s0, s1);
            } else if (first_hit.found) {
                cut = subPolylineByStation(extended, first_trim_station, back_original_station);
            } else if (last_hit.found) {
                cut = subPolylineByStation(extended, front_original_station, last_trim_station);
            }
            for (const auto& hit : selected_hits) {
                // 首/尾探针包含“外延段 + 第一个原始 cut 线段”。当第二个
                // cut 形点本身来自当前组 LaneEdge 时，会在离端点较远处
                // 命中同源 RoadEdge。该命中只用于禁止 cut 被截短，不能再
                // 写成 RoadEdge 裁剪 station，否则同一几何会以 cut 和
                // RoadEdge 两种身份入环并形成折返副环。
                const bool internal_same_group_hit =
                    hit.hit.distance > 1.0 && hitMatchesCurrentGroupEdge(hit.hit);
                if (!internal_same_group_hit &&
                    boundaryHitMatchesGroupDir(hit.hit, group_dir)) {
                    addBoundaryHitStation(hit_stations, hit.hit, cut_index, cut_role);
                }
            }
            result.cut_lines.push_back(cut);
        }

        for (int i = 0; i < (int)cut_meta.size() && i < (int)result.cut_lines.size(); ++i)
            addProbeHits(cut_meta[i], i, i < (int)cut_roles.size() ? cut_roles[i] : GroupRole::Entry);

        for (int i = 0; i < (int)boundary_lines.size(); ++i) {
            if (boundary_lines[i].pts.size() < 2)
                continue;
            if (hit_stations[i].empty())
                continue;
            double total = polylineLength(boundary_lines[i].pts);
            std::sort(hit_stations[i].begin(), hit_stations[i].end(),
                      [](const BoundaryHitStation& a, const BoundaryHitStation& b) {
                          return a.station < b.station;
                      });

            auto hasRole = [](const BoundaryHitStation& hit, GroupRole role) {
                return std::find(hit.cut_roles.begin(), hit.cut_roles.end(), role) != hit.cut_roles.end();
            };
            std::vector<double> unique_stations;
            for (const auto& hit : hit_stations[i]) {
                if (unique_stations.empty() || std::abs(unique_stations.back() - hit.station) > 0.03)
                    unique_stations.push_back(hit.station);
            }
            if (unique_stations.size() >= 2) {
                double s0 = std::max(0.0, std::min(total, unique_stations.front()));
                double s1 = std::max(0.0, std::min(total, unique_stations.back()));
                if (s1 < s0)
                    std::swap(s0, s1);
                if (s1 - s0 > 0.05) {
                    std::vector<double> cuts;
                    cuts.push_back(s0);
                    for (double station : unique_stations) {
                        if (station > s0 + 0.03 && station < s1 - 0.03)
                            cuts.push_back(station);
                    }
                    cuts.push_back(s1);
                    std::sort(cuts.begin(), cuts.end());
                    std::vector<Vec3d> clipped;
                    if (cuts.size() == 2) {
                        clipped = subInwardClosedPolylineByStation(
                            boundary_lines[i].pts, cuts.front(), cuts.back(), center);
                    } else {
                        for (int k = 0; k + 1 < (int)cuts.size(); ++k) {
                            std::vector<Vec3d> part = subPolylineByStation(
                                boundary_lines[i].pts, cuts[k], cuts[k + 1]);
                            for (const auto& pt : part)
                                pushUnique(clipped, pt);
                        }
                    }
                    result.boundary_lines.push_back(std::move(clipped));
                }
                continue;
            }

            const auto& hit = hit_stations[i].front();
            if (hasRole(hit, GroupRole::Entry) && hasRole(hit, GroupRole::Exit)) {
                if (total > 0.05)
                    result.boundary_lines.push_back(boundary_lines[i].pts);
            } else {
                // 合并链在多个 RoadEdge 共享端点处可能发生回折。若单命中点
                // 紧邻这样的共享节点，链方向只能告诉我们交通方向上的前后，
                // 不能决定应该把共享节点另一侧的整条远端分支收入路口面。
                // 先识别命中点附近的共享端点；该局部弧是道路边缘进入节点后
                // 向路口内凹回的部分，应优先于 H->远端 的整条分支。
                auto localSharedEndpointBranch = [&]() {
                    std::vector<Vec3d> empty;
                    if (boundary_lines[i].source_ids.size() <= 1 ||
                        boundary_lines[i].source_endpoints.size() < 2)
                        return empty;

                    auto stationOfPoint = [&](const Vec3d& point) {
                        double best_station = 0.0;
                        double best_distance = 1e18;
                        double station = 0.0;
                        for (int k = 0; k + 1 < (int)boundary_lines[i].pts.size(); ++k) {
                            Vec2d a = xyOf(boundary_lines[i].pts[k]);
                            Vec2d b = xyOf(boundary_lines[i].pts[k + 1]);
                            Vec2d ab = b - a;
                            double length = ab.norm();
                            if (length < 1e-9)
                                continue;
                            double t = (xyOf(point) - a).dot(ab) / ab.squaredNorm();
                            t = std::max(0.0, std::min(1.0, t));
                            Vec2d projected = a + t * ab;
                            double current_distance = dist(projected, xyOf(point));
                            if (current_distance < best_distance) {
                                best_distance = current_distance;
                                best_station = station + t * length;
                            }
                            station += length;
                        }
                        return best_station;
                    };

                    std::vector<Vec3d> unique_endpoints;
                    std::vector<int> endpoint_counts;
                    for (const auto& source : boundary_lines[i].source_endpoints) {
                        for (int endpoint_index = 0; endpoint_index < 2; ++endpoint_index) {
                            const Vec3d& endpoint = endpoint_index == 0
                                ? source.first : source.second;
                            if (endpoint_index == 1 && dist(source.first, source.second) <= 0.03)
                                continue;
                            int existing = -1;
                            for (int k = 0; k < (int)unique_endpoints.size(); ++k) {
                                if (dist(unique_endpoints[k], endpoint) <= 0.03) {
                                    existing = k;
                                    break;
                                }
                            }
                            if (existing < 0) {
                                unique_endpoints.push_back(endpoint);
                                endpoint_counts.push_back(1);
                            } else {
                                ++endpoint_counts[existing];
                            }
                        }
                    }

                    const double hit_station = std::max(0.0, std::min(total, hit.station));
                    double best_distance = 1e18;
                    double best_station = 0.0;
                    for (int k = 0; k < (int)unique_endpoints.size(); ++k) {
                        if (endpoint_counts[k] < 2)
                            continue;
                        const double endpoint_station = stationOfPoint(unique_endpoints[k]);
                        const double distance = std::abs(endpoint_station - hit_station);
                        if (distance < best_distance) {
                            best_distance = distance;
                            best_station = endpoint_station;
                        }
                    }

                    if (best_distance >= 1e17)
                        return empty;
                    const double sample_distance = 0.10;
                    const Vec2d shared_point = xyOf(pointAtStation(
                        boundary_lines[i].pts, best_station));
                    const Vec2d before = xyOf(pointAtStation(
                        boundary_lines[i].pts,
                        std::max(0.0, best_station - sample_distance))) - shared_point;
                    const Vec2d after = xyOf(pointAtStation(
                        boundary_lines[i].pts,
                        std::min(total, best_station + sample_distance))) - shared_point;
                    // 近共享点的局部弧只有在共享点确实是一个道路边缘转角时
                    // 才能代表内凹返回段。近 180 度的连续/掉头形链虽然也有
                    // 共享端点，但应继续保留其原有切点方向（100000610）。
                    if (before.norm() <= 1e-8 || after.norm() <= 1e-8 ||
                        before.normalized().dot(after.normalized()) <
                            std::cos(5.0 * M_PI / 6.0))
                        return empty;

                    // 1m 是现有“道路边缘端点命中”判定使用的同一物理邻近范围；
                    // 超出该范围仍按整条合并链的既有方向规则处理。
                    const double branch_threshold = groupCutOffset + 0.5;
                    if (best_distance <= 0.05 || best_distance > branch_threshold)
                        return empty;
                    return subPolylineByStation(
                        boundary_lines[i].pts, hit_station, best_station);
                };

                std::vector<Vec3d> local_branch = localSharedEndpointBranch();
                if (local_branch.size() >= 2 && polylineLength(local_branch) > 0.05) {
                    result.boundary_lines.push_back(std::move(local_branch));
                    continue;
                }

                // 单命中时，保留命中点朝道路外侧的一段：Entry 组的道路外侧
                // 在链方向的反侧，Exit 组的道路外侧在链方向的同侧。不能只按
                // 角色固定取 prefix；MergeBoundaryLines 可能把相邻 RoadEdge
                // 合成带回折的链（100002465 的 485|703），固定取 prefix
                // 会把另一条边的远端突尖段带入路口面。
                bool keep_prefix = false;
                bool direction_resolved = false;
                for (size_t k = 0; k < hit.cut_indices.size() &&
                                   k < hit.cut_roles.size(); ++k) {
                    const int cut_index = hit.cut_indices[k];
                    if (cut_index < 0 || cut_index >= (int)cut_group_dirs.size())
                        continue;
                    if (boundary_lines[i].source_ids.size() <= 1)
                        continue;
                    Vec2d group_dir = cut_group_dirs[cut_index];
                    if (group_dir.norm() <= 1e-8 || hit.chain_dir.norm() <= 1e-8)
                        continue;
                    // 原始边段与合并链反向时，旧的 source-direction 规则已经
                    // 能可靠识别链首/链尾；仅对同向且存在回折风险的链启用组方向
                    // 判定，避免改变 100000610 等正常的反向合并链。
                    if (hit.boundary_dir.norm() <= 1e-8 ||
                        hit.boundary_dir.dot(hit.chain_dir) <= 0.15)
                        continue;
                    const double alignment = hit.chain_dir.dot(group_dir.normalized());
                    if (std::abs(alignment) <= 0.15)
                        continue;
                    const bool chain_toward_outside = alignment > 0.0;
                    keep_prefix = hit.cut_roles[k] == GroupRole::Entry
                        ? chain_toward_outside : !chain_toward_outside;
                    direction_resolved = true;
                    break;
                }
                if (!direction_resolved) {
                    // 兼容缺少组方向元数据的旧输入，仍使用原始 RoadEdge
                    // 方向和组角色推断保留侧。
                    keep_prefix = hasRole(hit, GroupRole::Exit);
                    if (hit.boundary_dir.dot(hit.chain_dir) < 0.0)
                        keep_prefix = !keep_prefix;
                }
                double station = std::max(0.0, std::min(total, hit.station));
                if (keep_prefix && station > 0.05) {
                    result.boundary_lines.push_back(
                        subPolylineByStation(boundary_lines[i].pts, 0.0, station));
                } else if (!keep_prefix && total - station > 0.05) {
                    result.boundary_lines.push_back(
                        subPolylineByStation(boundary_lines[i].pts, station, total));
                }
            }
        }
        return result;
    }


    Vec2d partOutgoingDir(const EdgePart& part, bool from_front) const {
        const auto& pts = part.pts;
        if (pts.size() < 2)
            return Vec2d(1, 0);
        if (from_front) {
            for (int i = 1; i < (int)pts.size(); ++i) {
                Vec2d d = xyOf(pts[i]) - xyOf(pts.front());
                if (d.norm() > 1e-8)
                    return d.normalized();
            }
        } else {
            for (int i = (int)pts.size() - 2; i >= 0; --i) {
                Vec2d d = xyOf(pts[i]) - xyOf(pts.back());
                if (d.norm() > 1e-8)
                    return d.normalized();
            }
        }
        return Vec2d(1, 0);
    }

    void appendPartPoints(std::vector<Vec3d>& raw, const EdgePart& part, bool forward) const {
        if (forward) {
            for (const auto& pt : part.pts)
                pushUnique(raw, pt);
        } else {
            for (auto it = part.pts.rbegin(); it != part.pts.rend(); ++it)
                pushUnique(raw, *it);
        }
    }

    bool isSimpleRing(const std::vector<Vec3d>& raw) const {
        std::vector<Vec3d> ring = openRing(raw);
        int n = (int)ring.size();
        if (n < 4)
            return n >= 3;
        auto onSegment = [](const Vec2d& a, const Vec2d& b, const Vec2d& p) {
            return std::abs(cross2d(b - a, p - a)) <= 1e-8 &&
                   p.x() >= std::min(a.x(), b.x()) - 1e-8 &&
                   p.x() <= std::max(a.x(), b.x()) + 1e-8 &&
                   p.y() >= std::min(a.y(), b.y()) - 1e-8 &&
                   p.y() <= std::max(a.y(), b.y()) + 1e-8;
        };
        auto intersectsInclusive = [&](const Vec3d& a3, const Vec3d& b3,
                                        const Vec3d& c3, const Vec3d& d3) {
            const Vec2d a = xyOf(a3), b = xyOf(b3), c = xyOf(c3), d = xyOf(d3);
            const Vec2d ab = b - a, cd = d - c;
            const double o1 = cross2d(ab, c - a);
            const double o2 = cross2d(ab, d - a);
            const double o3 = cross2d(cd, a - c);
            const double o4 = cross2d(cd, b - c);
            auto opposite = [](double x, double y) {
                return (x > 1e-8 && y < -1e-8) || (x < -1e-8 && y > 1e-8);
            };
            if (opposite(o1, o2) && opposite(o3, o4))
                return true;
            return (std::abs(o1) <= 1e-8 && onSegment(a, b, c)) ||
                   (std::abs(o2) <= 1e-8 && onSegment(a, b, d)) ||
                   (std::abs(o3) <= 1e-8 && onSegment(c, d, a)) ||
                   (std::abs(o4) <= 1e-8 && onSegment(c, d, b));
        };
        for (int i = 0; i < n; ++i) {
            for (int j = i + 2; j < n; ++j) {
                if (i == 0 && j == n - 1)
                    continue;
                if (intersectsInclusive(ring[i], ring[(i + 1) % n],
                                         ring[j], ring[(j + 1) % n]))
                    return false;
            }
        }
        return true;
    }

    // 将连接图产生的自接触环按偶奇规则拆成简单环。连接图中的人工闭合边
    // 可能把一条端点折返 RoadEdge 围成一个很小的副环；直接返回原环会把
    // 非相邻端点接触漏过。Clipper 只负责 XY 拆环，Z 从最近的原始点继承。
    std::vector<Vec3d> simplifyNonSimpleRing(
        const std::vector<Vec3d>& raw, const Vec2d& center) const {
        const std::vector<Vec3d> ring = openRing(raw);
        if (ring.size() < 3)
            return {};
        ClipperLib::Path path;
        path.reserve(ring.size());
        const double scale = 1e6;
        for (const auto& p : ring)
            path.emplace_back((ClipperLib::cInt)std::llround(p.x() * scale),
                              (ClipperLib::cInt)std::llround(p.y() * scale));
        ClipperLib::Paths split;
        ClipperLib::SimplifyPolygon(path, split, ClipperLib::pftEvenOdd);
        if (split.empty())
            return {};

        // 记录原始环中由非相邻边相碰形成的顶点。Clipper 拆环后这类
        // 连接锚点可能仍留在外环上；它们只是自接触副环的分叉点，
        // 不应继续作为路口面顶点输出。
        std::vector<Vec2d> touch_vertices;
        auto onSegment = [](const Vec2d& a, const Vec2d& b, const Vec2d& p) {
            return std::abs(cross2d(b - a, p - a)) <= 1e-8 &&
                   p.x() >= std::min(a.x(), b.x()) - 1e-8 &&
                   p.x() <= std::max(a.x(), b.x()) + 1e-8 &&
                   p.y() >= std::min(a.y(), b.y()) - 1e-8 &&
                   p.y() <= std::max(a.y(), b.y()) + 1e-8;
        };
        for (int i = 0; i < (int)ring.size(); ++i) {
            const Vec2d vertex = xyOf(ring[i]);
            for (int j = i + 2; j < (int)ring.size(); ++j) {
                if (i == 0 && j == (int)ring.size() - 1)
                    continue;
                const Vec2d other_vertex = xyOf(ring[j]);
                if (onSegment(xyOf(ring[j]), xyOf(ring[(j + 1) % ring.size()]), vertex) ||
                    onSegment(xyOf(ring[i]), xyOf(ring[(i + 1) % ring.size()]), other_vertex)) {
                    touch_vertices.push_back(vertex);
                    if (dist(other_vertex, xyOf(ring[i])) > 1e-8 &&
                        dist(other_vertex, xyOf(ring[(i + 1) % ring.size()])) > 1e-8)
                        touch_vertices.push_back(other_vertex);
                    break;
                }
            }
        }

        int selected = -1;
        double selected_area = -1.0;
        const ClipperLib::IntPoint center_i(
            (ClipperLib::cInt)std::llround(center.x() * scale),
            (ClipperLib::cInt)std::llround(center.y() * scale));
        for (int i = 0; i < (int)split.size(); ++i) {
            const double area = std::abs(ClipperLib::Area(split[i]));
            const bool contains = ClipperLib::PointInPolygon(center_i, split[i]) != 0;
            const bool selected_contains = selected >= 0 &&
                ClipperLib::PointInPolygon(center_i, split[selected]) != 0;
            if (selected < 0 ||
                (contains && !selected_contains) ||
                (contains == selected_contains && area > selected_area)) {
                selected = i;
                selected_area = area;
            }
        }
        if (selected < 0 || split[selected].size() < 3)
            return {};

        std::vector<Vec3d> result;
        result.reserve(split[selected].size() + 1);
        for (const auto& ip : split[selected]) {
            const Vec2d xy(ip.X / scale, ip.Y / scale);
            bool is_touch_vertex = false;
            for (const auto& touch : touch_vertices) {
                if (dist(touch, xy) <= 1e-6) {
                    is_touch_vertex = true;
                    break;
                }
            }
            if (is_touch_vertex)
                continue;
            int nearest = 0;
            double nearest_d = 1e18;
            for (int i = 0; i < (int)ring.size(); ++i) {
                const double d = dist(ring[i], xy);
                if (d < nearest_d) {
                    nearest_d = d;
                    nearest = i;
                }
            }
            result.emplace_back(xy, ring[nearest].z());
        }
        return repairPolygon(result, center);
    }

    struct RingCandidate {
        std::vector<Vec3d> pts;
        int used_count = 0;
        int real_used_count = 0;
        bool closed = false;
        bool simple = false;
        bool contains_center = false;
        double abs_area = 0.0;
    };

    RingCandidate traceEdgePartsRing(
        const std::vector<EdgePart>& parts,
        const std::vector<EdgeNode>& nodes,
        int start_part,
        bool start_forward,
        int turn_policy,
        const Vec2d& center) const {
        RingCandidate out;
        if (start_part < 0 || start_part >= (int)parts.size())
            return out;

        std::vector<bool> used(parts.size(), false);
        const EdgePart& first = parts[start_part];
        int start_node = start_forward ? first.node0 : first.node1;
        int current_node = start_forward ? first.node1 : first.node0;
        int current_part = start_part;
        used[start_part] = true;
        out.used_count = 1;
        out.real_used_count = first.is_artificial_closure ? 0 : 1;
        appendPartPoints(out.pts, first, start_forward);

        while (out.used_count < (int)parts.size()) {
            int best_part = -1;
            bool best_at_front = true;
            double best_score = (turn_policy == 2) ? -1e18 : 1e18;
            Vec2d incoming(1, 0);
            if (out.pts.size() >= 2) {
                incoming = xyOf(out.pts.back()) - xyOf(out.pts[out.pts.size() - 2]);
                if (incoming.norm() > 1e-8)
                    incoming.normalize();
            }
            double incoming_ang = std::atan2(incoming[1], incoming[0]);

            bool has_real_candidate = false;
            for (const auto& incident : nodes[current_node].incident) {
                int candidate_part = incident.first;
                if (candidate_part != current_part && !used[candidate_part] &&
                    !parts[candidate_part].is_artificial_closure) {
                    has_real_candidate = true;
                    break;
                }
            }
            if (current_node == start_node && out.used_count > 1 && !has_real_candidate)
                break;

            for (const auto& incident : nodes[current_node].incident) {
                int candidate_part = incident.first;
                bool at_front = incident.second;
                if (candidate_part == current_part || used[candidate_part])
                    continue;
                if (has_real_candidate && parts[candidate_part].is_artificial_closure)
                    continue;
                Vec2d outgoing = partOutgoingDir(parts[candidate_part], at_front);
                double outgoing_ang = std::atan2(outgoing[1], outgoing[0]);
                double turn = normalizedAngleDelta(incoming_ang, outgoing_ang);
                double score = std::abs(turn);
                if (turn_policy == 1)
                    score = turn;
                else if (turn_policy == 2)
                    score = turn;

                bool better = false;
                if (turn_policy == 2)
                    better = score > best_score;
                else
                    better = score < best_score;
                if (better) {
                    best_score = score;
                    best_part = candidate_part;
                    best_at_front = at_front;
                }
            }
            if (best_part < 0)
                break;

            bool forward = best_at_front;
            appendPartPoints(out.pts, parts[best_part], forward);
            used[best_part] = true;
            ++out.used_count;
            if (!parts[best_part].is_artificial_closure)
                ++out.real_used_count;
            current_part = best_part;
            current_node = forward ? parts[best_part].node1 : parts[best_part].node0;
        }

        out.closed = (current_node == start_node && out.used_count > 1);
        out.simple = out.closed && isSimpleRing(out.pts);
        out.abs_area = polygonArea(openRing(out.pts));
        out.contains_center = out.closed && pointInPolygon(center, openRing(out.pts));
        return out;
    }

    std::vector<Vec3d> buildPolygonByAngularChains(
        const std::vector<std::vector<Vec3d>>& boundary_lines,
        const std::vector<std::vector<Vec3d>>& cut_lines, const Vec2d& center) const {
        std::vector<AreaChain> chains;
        for (const auto& line : boundary_lines) {
            if (line.size() < 2)
                continue;
            AreaChain chain;
            chain.pts = line;
            chains.push_back(orientChain(chain, center));
        }
        for (const auto& line : cut_lines) {
            if (line.size() < 2)
                continue;
            AreaChain chain;
            chain.pts = line;
            chains.push_back(orientChain(chain, center));
        }
        if (chains.size() < 3)
            return {};

        std::sort(chains.begin(), chains.end(), [&](const AreaChain& a, const AreaChain& b) {
            double angA = angleOf(chainMidpoint(a), center);
            double angB = angleOf(chainMidpoint(b), center);
            return winding == "clockwise" ? angA > angB : angA < angB;
        });

        std::vector<Vec3d> raw;
        for (const auto& chain : chains)
            for (const auto& pt : chain.pts)
                pushUnique(raw, pt);
        if (raw.size() < 3)
            return {};
        return repairPolygon(raw, center);
    }

    std::vector<Vec3d> buildPolygonByConnectedParts(
        const std::vector<std::vector<Vec3d>>& boundary_lines,
        const std::vector<std::vector<Vec3d>>& cut_lines, const Vec2d& center,
        bool prefer_center_containment = false) const {

        std::vector<EdgePart> parts;
        // 添加边缘段
        auto add_part = [&](const std::vector<Vec3d>& line) {
            if (line.size() < 2 || dist(line.front(), line.back()) < 1e-6) return;
            EdgePart part;
            part.pts = line;
            parts.push_back(part);
        };

        // 获取端点索引号
        auto endpointNode = [&](std::vector<EdgeNode>& nodes, const Vec3d& pt, double tol) -> int {
            for (int i = 0; i < (int)nodes.size(); ++i) {
                if (dist(nodes[i].pt, pt) <= tol)
                    return i;
            }
            EdgeNode node;
            node.pt = pt;
            nodes.push_back(node);
            return (int)nodes.size() - 1;
        };

        //
        auto betterRingCandidate = [&](const RingCandidate& a, const RingCandidate& b) -> bool {
            if (a.simple != b.simple)
                return a.simple;
            if (a.closed != b.closed)
                return a.closed;
            if (prefer_center_containment && a.contains_center != b.contains_center)
                return a.contains_center;
            if (a.real_used_count != b.real_used_count)
                return a.real_used_count > b.real_used_count;
            if (a.used_count != b.used_count)
                return a.used_count > b.used_count;
            return a.abs_area > b.abs_area;
        };

        // 存储组端边缘线和道路边缘线
        for (const auto& line : boundary_lines)
            add_part(line);
        for (const auto& line : cut_lines)
            add_part(line);
        if (parts.size() < 2)
            return {};

        // 获取所有边缘端点
        std::vector<EdgeNode> nodes;
        double node_tol = std::max(0.05, std::min(0.25, snapTolerance));
        for (int i = 0; i < (int)parts.size(); ++i) {
            parts[i].node0 = endpointNode(nodes, parts[i].pts.front(), node_tol);
            parts[i].node1 = endpointNode(nodes, parts[i].pts.back(), node_tol);
            if (parts[i].node0 == parts[i].node1)
                continue;
            parts[i].pts.front() = nodes[parts[i].node0].pt;
            parts[i].pts.back() = nodes[parts[i].node1].pt;
            nodes[parts[i].node0].incident.push_back({i, true});
            nodes[parts[i].node1].incident.push_back({i, false});
        }

        std::vector<int> dangling_nodes;
        for (int i = 0; i < (int)nodes.size(); ++i) {
            if (nodes[i].incident.size() == 1)
                dangling_nodes.push_back(i);
        }
        if (dangling_nodes.size() >= 2) {
            std::sort(dangling_nodes.begin(), dangling_nodes.end(), [&](int a, int b) {
                double aa = angleOf(nodes[a].pt, center);
                double ab = angleOf(nodes[b].pt, center);
                return winding == "clockwise" ? aa > ab : aa < ab;
            });
            auto has_existing_real_part = [&](int n0, int n1) {
                for (const auto& part : parts) {
                    if (part.is_artificial_closure)
                        continue;
                    if ((part.node0 == n0 && part.node1 == n1) ||
                        (part.node0 == n1 && part.node1 == n0))
                        return true;
                }
                return false;
            };
            for (int i = 0; i < (int)dangling_nodes.size(); ++i) {
                int n0 = dangling_nodes[i];
                int n1 = dangling_nodes[(i + 1) % dangling_nodes.size()];
                if (n0 == n1 || dist(nodes[n0].pt, nodes[n1].pt) < 1e-6)
                    continue;
                if (has_existing_real_part(n0, n1))
                    continue;
                EdgePart part;
                part.pts = {nodes[n0].pt, nodes[n1].pt};
                part.node0 = n0;
                part.node1 = n1;
                part.is_artificial_closure = true;
                int idx = (int)parts.size();
                parts.push_back(part);
                nodes[n0].incident.push_back({idx, true});
                nodes[n1].incident.push_back({idx, false});
            }
        }

        RingCandidate best;
        for (int i = 0; i < (int)parts.size(); ++i) {
            if (parts[i].node0 < 0 || parts[i].node1 < 0 || parts[i].node0 == parts[i].node1)
                continue;
            for (int policy = 0; policy < 3; ++policy) {
                RingCandidate a = traceEdgePartsRing(parts, nodes, i, true, policy, center);
                if (betterRingCandidate(a, best))
                    best = a;
                RingCandidate b = traceEdgePartsRing(parts, nodes, i, false, policy, center);
                if (betterRingCandidate(b, best))
                    best = b;
            }
        }

        if (!best.closed || best.used_count < 3 || best.pts.size() < 3)
            return {};
        if (prefer_center_containment && !best.contains_center)
            return {};
        return repairPolygon(best.pts, center);
    }

    bool polygonUsesLineSegment(
        const std::vector<Vec3d>& polygon,
        const std::vector<Vec3d>& line,
        double tol = 0.08,
        double min_station_span = 0.20) const {
        std::vector<double> stations;
        double acc = 0.0;
        for (int i = 0; i + 1 < (int)line.size(); ++i) {
            Vec2d a = xyOf(line[i]);
            Vec2d b = xyOf(line[i + 1]);
            Vec2d ab = b - a;
            double seg_len = ab.norm();
            if (seg_len < 1e-9)
                continue;
            for (const auto& pt : polygon) {
                double t = (xyOf(pt) - a).dot(ab) / ab.squaredNorm();
                if (t < -1e-6 || t > 1.0 + 1e-6)
                    continue;
                t = std::max(0.0, std::min(1.0, t));
                Vec2d proj = a + t * ab;
                if (dist(xyOf(pt), proj) <= tol)
                    stations.push_back(acc + t * seg_len);
            }
            acc += seg_len;
        }
        if (stations.size() < 2)
            return false;
        auto mm = std::minmax_element(stations.begin(), stations.end());
        return *mm.second - *mm.first >= min_station_span;
    }

    int usedBoundarySegmentCount(
        const std::vector<Vec3d>& polygon,
        const std::vector<std::vector<Vec3d>>& boundary_lines) const {
        int count = 0;
        for (const auto& line : boundary_lines) {
            if (polygonUsesLineSegment(polygon, line))
                ++count;
        }
        return count;
    }

    // 组装精细路口面：道路边缘线与所有进入/退出组端侧边界折线共同排序成可凹多边形。
    // 端侧边界先按组回缩/外扩，再与道路边缘对齐，最后统一修复多边形绕向与闭合。
    std::vector<Vec3d> buildFromLaneGroupCuts(const IntersectionInput& inp,
        const std::vector<ConnectivityCurve>& centerlines, const Vec2d& center) const {

        std::vector<LaneGroupId> group_ids;
        std::vector<std::vector<Vec3d>> cut_lines; //进入组尾端(或退出组首端)缩扩端点折线
        std::vector<GroupRole> cut_roles;
        std::vector<LaneGroupCutLine> cut_meta;
        std::vector<Vec2d> cut_group_dirs;
        std::vector<std::vector<LaneEdgeId>> cut_group_edge_ids;
        std::vector<BoundaryLine> boundary_lines; //道路边缘线
        auto fillBoundarySegmentDirs = [](BoundaryLine& line) {
            line.segment_dirs.clear();
            for (int i = 0; i + 1 < (int)line.pts.size(); ++i) {
                Vec2d d = xyOf(line.pts[i + 1]) - xyOf(line.pts[i]);
                line.segment_dirs.push_back(d.norm() > 1e-8 ? d.normalized() : Vec2d(1, 0));
            }
        };

        // 获取构成路口外轮廓的进入/退出组的缩扩端点折线
        std::vector<LaneGroup> groups = effectiveLaneGroups(inp, centerlines);
        for (const auto& group : groups) {
            LaneGroupCutLine cut = getLaneGroupCutLine(inp, group, centerlines, center);
            std::vector<Vec3d> line = cut.line;
            if (line.size() < 2)
                line = getLaneGroupCurveCutLine(inp, group, centerlines, center);
            if (line.size() >= 2) {
                cut_lines.push_back(line);
                cut_roles.push_back(group.role);
                cut_group_dirs.push_back(laneGroupTrafficDir(inp, group, center));
                cut_group_edge_ids.push_back(laneGroupEdgeIds(inp, group));
                cut_meta.push_back(std::move(cut));
                group_ids.push_back(group.id);
            }
        }

        // 获取构成路口外轮廓的道路边缘线
        if (centerlines.empty()) {
            std::vector<std::vector<Vec3d>> raw_boundary_lines;
            for (const auto& bnd : inp.boundaries) {
                if (!isUsableRoadEdgeBoundary(inp, bnd))
                    continue;
                raw_boundary_lines.push_back(bnd.geometry.points);
            }
            raw_boundary_lines = IntersectionAreaBuilder::MergeBoundaries(
                std::move(raw_boundary_lines));
            for (auto& pts : raw_boundary_lines) {
                BoundaryLine line;
                line.pts = std::move(pts);
                line.source_kind = AreaSourceKind::Boundary;
                fillBoundarySegmentDirs(line);
                boundary_lines.push_back(std::move(line));
            }
        } else {
            for (const auto& bnd : inp.boundaries) {
                if (!isUsableRoadEdgeBoundary(inp, bnd))
                    continue;
                BoundaryLine line;
                line.pts = bnd.geometry.points;
                line.source_ids.push_back(bnd.id);
                line.source_kind = AreaSourceKind::Boundary;
                line.source_endpoints.push_back({bnd.geometry.points.front(), bnd.geometry.points.back()});
                fillBoundarySegmentDirs(line);
                boundary_lines.push_back(std::move(line));
            }
            boundary_lines = MergeBoundaryLines(std::move(boundary_lines));
        }

        // 打断处理后集合段：组端边缘线+道路边缘
        AlignedAreaParts parts = alignCutLinesAndBoundaries(
            boundary_lines, cut_lines, cut_roles, cut_meta, cut_group_dirs,
            cut_group_edge_ids, center, group_ids);

        bool prefer_center_containment = parts.boundary_lines.size() >= 3;
        std::vector<Vec3d> connected = buildPolygonByConnectedParts(
            parts.boundary_lines, parts.cut_lines, center, prefer_center_containment);
        std::vector<Vec3d> angular = buildPolygonByAngularChains(
            parts.boundary_lines, parts.cut_lines, center);
        // 先按包含非相邻端点接触的严格口径验环。历史连接图在人工闭合边
        // 参与时可能生成自接触环；仅在候选不简单时拆分，避免改变正常凹环。
        if (!isSimpleRing(connected))
            connected = simplifyNonSimpleRing(connected, center);
        if (!isSimpleRing(angular))
            angular = simplifyNonSimpleRing(angular, center);
        if (connected.size() >= 4 && angular.size() >= 4 && isSimpleRing(angular)) {
            int connected_boundary_count = usedBoundarySegmentCount(connected, parts.boundary_lines);
            int angular_boundary_count = usedBoundarySegmentCount(angular, parts.boundary_lines);
            if (angular_boundary_count > connected_boundary_count)
                return angular;
            return connected;
        }
        if (connected.size() >= 4)
            return connected;
        return angular;
    }

    // 精细路口面入口：优先使用进入/退出组端侧边界折线围成可凹多边形。
    // 道路边缘存在时参与端侧线对齐和边缘截断；不存在时仍保留组端侧线的回缩/外扩结果。
    std::vector<Vec3d> buildFinePolygon(const IntersectionInput& inp,
        const std::vector<ConnectivityCurve>& centerlines, const Vec2d& center, double radius) const {
        std::vector<Vec3d> by_group_cuts = buildFromLaneGroupCuts(inp, centerlines, center);
        if (by_group_cuts.size() >= 4)
            return by_group_cuts;
        return {};
    }

    // 计算路口中心点（连接点的质心）
    Vec2d computeCenter(const IntersectionInput& inp,
                        const std::vector<ConnectivityCurve>& cls) const {
        Vec2d c{0, 0};
        int cnt = 0;
        for (auto& gcl : cls) {
            if(gcl.curve) {
                c += gcl.curve->startPt(); ++cnt;
                c += gcl.curve->endPt();  ++cnt;
            }
        }
        if (cnt > 0) return c * (1.0 / cnt);

        // 回退：用中心线连接点
        for (auto& cl : inp.lanes) {
            c += getConnPoint(cl.geometry.points, inp.IsEntryLane(cl.id));
            ++cnt;
        }
        return cnt > 0 ? c * (1.0 / cnt) : Vec2d{0, 0};
    }

    // 估算路口半径（连接点到中心的最大距离）
    double computeRadius(const IntersectionInput& inp,
                         const std::vector<ConnectivityCurve>& cls,
                         const Vec2d& center) const {
        double maxD = 5.0;
        for (auto& gcl : cls) {
            if (gcl.curve) {
                maxD = std::max(maxD, dist(gcl.curve->startPt(), center));
                maxD = std::max(maxD, dist(gcl.curve->endPt(), center));
            }
        }
        for (auto& cl : inp.lanes) {
            auto connPt = getConnPoint(cl.geometry.points, inp.IsEntryLane(cl.id));
            maxD = std::max(maxD, dist(connPt, center));
        }
        return maxD;
    }

    // 将端点吸附到最近的车道边线端点（若在容差范围内）
    Vec3d snapToEdgePt(
        const Vec3d& pt, const std::vector<ConnectivityLaneEdge>& edgelines,
        const IntersectionInput& inp, double tol) const {
        double minD = tol;
        Vec3d best = pt;
        // 检查生成的路口内边线端点
        for (auto& el : edgelines) {
            if (el.geometry.points.empty()) continue;
            double d0 = dist(pt, el.geometry.points.front());
            double d1 = dist(pt, el.geometry.points.back());
            if (d0 < minD) {
                minD = d0;
                best = el.geometry.points.front();
            }
            if (d1 < minD) {
                minD = d1;
                best = el.geometry.points.back();
            }
        }
        // 检查路口外边线连接点
        for (auto& el : inp.lane_edges) {
            if (el.geometry.points.empty()) continue;
            auto connPt = getConnPoint3d(el.geometry.points, true);
            double d = dist(pt, connPt);
            if (d < minD) {
                minD = d;
                best = connPt;
            }
        }
        return best;
    }

    // 从车道边线端点补充路口面顶点
    void supplementEndpointsFromLaneEdges(
        std::vector<struct EdgeEndpoint>& endpoints,
        const IntersectionInput& inp,
        const std::vector<ConnectivityLaneEdge>& edgelines,
        const Vec2d& center, double radius) const {
        // 收集路口外边线的连接点（这些是路口面的边界顶点）
        for (auto& el : inp.lane_edges) {
            if (el.geometry.points.empty()) continue;
            Vec3d cp = getConnPoint3d(el.geometry.points, inp.IsEntryLaneEdge(el.id));
            if (dist(cp, center) > radius * 2.0) continue;

            // 检查是否已有相近的端点
            bool dup = false;
            for (auto& ep : endpoints) {
                if (dist(ep.pt, cp) < 0.5) {
                    dup = true;
                    break;
                }
            }
            if (!dup) endpoints.emplace_back(EdgeEndpoint{cp, -1, true});
        }

        // 补充中心线连接点的法线方向偏移点（近似边界点）
        for (auto& cl : inp.lanes) {
            if (cl.geometry.points.size() < 2) continue;
            Vec3d cp3 = getConnPoint3d(cl.geometry.points, inp.IsEntryLane(cl.id));
            Vec2d ct = getConnTangent(cl.geometry.points, inp.IsEntryLane(cl.id));
            Vec2d n = rotLeft(ct);
            // 左侧点
            Vec3d lp = offsetXY(cp3, n * snapTolerance * 2);
            Vec3d rp = offsetXY(cp3, -n * snapTolerance * 2);
            bool dupL = false, dupR = false;
            for (auto& ep : endpoints) {
                if (dist(ep.pt, lp) < 0.5) dupL = true;
                if (dist(ep.pt, rp) < 0.5) dupR = true;
            }
            if (!dupL) endpoints.emplace_back(EdgeEndpoint{lp, -1, true});
            if (!dupR) endpoints.emplace_back(EdgeEndpoint{rp, -1, true});
        }
    }

    // 从道路边缘线提取两点间的几何段
    std::vector<Vec3d> extractEdgeSegment(
        const Boundary& re, const Vec3d& a, const Vec3d& b, double tol) const {
        if (re.geometry.points.size() < 2) return {};

        // 找到 a 和 b 在 re.geometry.points 上最近的索引
        int idxA = -1, idxB = -1;
        double minDA = 1e18, minDB = 1e18;
        for (int i = 0; i < (int)re.geometry.points.size(); ++i) {
            double da = dist(re.geometry.points[i], a);
            double db = dist(re.geometry.points[i], b);
            if (da < minDA) {
                minDA = da, idxA = i;
            }
            if (db < minDB) {
                minDB = db, idxB = i;
            }
        }
        if (idxA < 0 || idxB < 0 || idxA == idxB) return {};
        std::vector<Vec3d> seg;
        int step = (idxA < idxB) ? 1 : -1;
        for (int i = idxA; i != idxB; i += step) {
            seg.push_back(re.geometry.points[i]);
        }
        // 不包含 b 本身（由下一个端点自己加）
        return seg;
    }

    // 去重端点（合并距离过近的点）
    void deduplicateEndpoints(std::vector<struct EdgeEndpoint>& eps, double tol) const {
        std::vector<bool> remove(eps.size(), false);
        for (size_t i = 0; i < eps.size(); ++i) {
            if (remove[i]) continue;
            for (size_t j = i + 1; j < eps.size(); ++j) {
                if (!remove[j] && dist(eps[i].pt, eps[j].pt) < tol) {
                    remove[j] = true;
                }
            }
        }
        std::vector<struct EdgeEndpoint> out;
        for (size_t i = 0; i < eps.size(); ++i)
            if (!remove[i]) out.push_back(eps[i]);
        eps = out;
    }

    // 多边形修复：移除共线点，确保方向正确
    std::vector<Vec3d> repairPolygon(const std::vector<Vec3d>& raw, const Vec2d& center) const {
        if (raw.size() < 3) return raw;
        // 移除重复点，保留参与路口面输出的三维坐标高程。
        std::vector<Vec3d> pts;
        for (auto& p : raw) {
            if (pts.empty() || dist(pts.back(), p) > EPS) pts.push_back(p);
        }
        // 移除首尾近重复点，避免闭合处保留厘米级回钩造成自交。
        while (pts.size() > 1 && dist(pts.front(), pts.back()) < 0.03)
            pts.pop_back();
        if (pts.size() < 3) return raw;
        bool removed_loop_anchor = true;
        while (removed_loop_anchor) {
            removed_loop_anchor = false;
            for (int i = 0; i < (int)pts.size() && !removed_loop_anchor; ++i) {
                for (int j = i + 2; j < (int)pts.size(); ++j) {
                    if (i == 0 && j == (int)pts.size() - 1)
                        continue;
                    if (dist(pts[i], pts[j]) < 0.03) {
                        pts.erase(pts.begin() + i);
                        removed_loop_anchor = true;
                        break;
                    }
                }
            }
        }
        if (pts.size() < 3) return raw;
        // 检查方向（有符号面积）
        double area = polygonSignedArea(pts);
        bool isCCW = area > 0;
        if (winding == "clockwise" && isCCW) {
            std::reverse(pts.begin(), pts.end());
        } else if (winding == "counter_clockwise" && !isCCW) {
            std::reverse(pts.begin(), pts.end());
        }
        // 闭合
        if (dist(pts.front(), pts.back()) > EPS)
            pts.push_back(pts.front());
        return pts;
    }

    // 最终回退：所有连接点的凸包
    std::vector<Vec3d> buildConvexHullFallback(
        const IntersectionInput& inp, const std::vector<ConnectivityCurve>& cls) const {
        std::vector<Vec3d> pts;
        for (auto& gcl : cls) {
            if (gcl.curve) {
                pts.push_back(connectivityEndpoint3d(gcl, true));
                pts.push_back(connectivityEndpoint3d(gcl, false));
            }
        }
        for (auto& cl : inp.lanes) {
            Vec3d cp = getConnPoint3d(cl.geometry.points, inp.IsEntryLane(cl.id));
            pts.push_back(cp);
        }
        if (pts.empty()) return {};
        return convexHull(pts);
    }

    // Andrew单调链凸包算法。排序与叉积使用 XY，输出保留每个凸包顶点的 Z。
    std::vector<Vec3d> convexHull(std::vector<Vec3d> pts) const {
        int n = pts.size();
        if (n < 3) return pts;
        std::sort(pts.begin(), pts.end(), [](const Vec3d& a, const Vec3d& b) {
            return a[0] < b[0] || (a[0] == b[0] && a[1] < b[1]);
        });

        std::vector<Vec3d> hull;
        // 下凸包
        for (int i = 0; i < n; ++i) {
            while (hull.size() >= 2) {
                Vec2d a = xyOf(hull[hull.size() - 2]), b = xyOf(hull.back());
                if (cross2d((b - a), xyOf(pts[i]) - a) <= 0) hull.pop_back();
                else break;
            }
            hull.push_back(pts[i]);
        }
        // 上凸包
        int lower = hull.size();
        for (int i = n - 2; i >= 0; --i) {
            while ((int)hull.size() > lower) {
                Vec2d a = xyOf(hull[hull.size() - 2]), b = xyOf(hull.back());
                if (cross2d((b - a), xyOf(pts[i]) - a) <= 0) hull.pop_back();
                else break;
            }
            hull.push_back(pts[i]);
        }
        hull.pop_back();
        if (!hull.empty()) hull.push_back(hull.front()); // 闭合
        return hull;
    }
};

std::vector<BoundaryLine> IntersectionAreaBuilder::MergeBoundaryLines(
    std::vector<BoundaryLine> boundaries) {
    return IntersectionAreaBuilderImpl::MergeBoundaryLines(std::move(boundaries));
}

IntersectionAreaBuilder::IntersectionAreaBuilder(
    double group_cut_offset,
    double cut_line_extension,
    double snap_tolerance,
    std::string winding)
    : group_cut_offset_(group_cut_offset),
      cut_line_extension_(cut_line_extension),
      snap_tolerance_(snap_tolerance),
      winding_(std::move(winding)) {}

IntersectionArea IntersectionAreaBuilder::build(
    const IntersectionInput& input,
    const std::vector<ConnectivityCurve>& centerlines,
    const std::vector<ConnectivityLaneEdge>& edge_lines) {
    IntersectionAreaBuilderImpl impl(
        group_cut_offset_, cut_line_extension_, snap_tolerance_, winding_);
    return impl.build(input, centerlines, edge_lines);
}

}
