#include "nonmotor_backbone_extractor.h"
#include "io/iodata_shapefile.h"
#include "utils.h"
#include "utils/clipper.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <queue>
#include <set>
#include <sstream>
#include <utility>
#include <vector>


namespace isg {

/// 无限直线与多边形边界的单个交点，以及该点沿射线方向的有符号站距。
struct LineHit {
    /// 交点的三维坐标，平面相交计算使用 XY。
    Vec3d point;
    /// 交点相对直线起点沿单位方向的有符号距离。
    double station = 0.0;

    /// 用交点坐标和站距构造一个边界相交记录。
    /// @param p 直线与多边形边界的交点。
    /// @param s 交点在直线方向上的有符号站距。
    LineHit(const Vec3d& p, double s) : point(p), station(s) {}
};

/// 尚未写入骨干图的线段记录，同时保存其来源以及两端是否属于主环端点。
struct Primitive {
    /// 输出骨干图中使用的边 ID。
    std::string id;
    /// 生成该线段的横道或安全岛 ID。
    std::string source_id;
    /// 线段类别，用于构建图边及区分横道/安全岛骨干。
    NonMotorBackboneEdgeType type = NonMotorBackboneEdgeType::Gap;
    /// 有向线段起点。
    Vec3d first;
    /// 有向线段终点。
    Vec3d last;
    /// 起点是否参与主环端点排序和 Gap 闭环。
    bool first_terminal = true;
    /// 终点是否参与主环端点排序和 Gap 闭环。
    bool last_terminal = true;
    /// 来源线具有明确车道组方向时的参考向量。
    Vec2d group_direction{1, 0};
    /// 来源线的进入/退出组角色。
    GroupRole group_role = GroupRole::Entry;
    /// 来源线是否具有车道组方向基准。
    bool has_group_direction = false;

    /// 使用边标识、来源、类别和端点属性构造待入图线段。
    /// @param edge_id 输出边 ID。
    /// @param source 来源横道或安全岛 ID。
    /// @param edge_type 输出边类型。
    /// @param a 线段起点。
    /// @param b 线段终点。
    /// @param first_is_terminal 起点是否是主环端点。
    /// @param last_is_terminal 终点是否是主环端点。
    /// @param direction 来源车道组的矢量方向。
    /// @param role 来源车道组角色。
    /// @param has_direction 是否存在有效车道组方向。
    Primitive(const std::string& edge_id, const std::string& source,
              NonMotorBackboneEdgeType edge_type, const Vec3d& a, const Vec3d& b,
              bool first_is_terminal, bool last_is_terminal,
              const Vec2d& direction = Vec2d(1, 0),
              GroupRole role = GroupRole::Entry, bool has_direction = false)
        : id(edge_id), source_id(source), type(edge_type), first(a), last(b),
          first_terminal(first_is_terminal), last_terminal(last_is_terminal),
          group_direction(direction), group_role(role), has_group_direction(has_direction) {}
};

/// 一条裁剪后的骨干线及其定向基准。
struct SpineGeometry {
    /// 统一方向后的骨干线首点。
    Vec3d first;
    /// 统一方向后的骨干线尾点。
    Vec3d last;
    /// 用于判定首尾左右关系的组矢量方向。
    Vec2d group_direction{1, 0};
    /// 组矢量方向所属的进入/退出角色。
    GroupRole group_role = GroupRole::Entry;
    /// 是否存在可用的车道组方向基准。
    bool has_group_direction = false;
};

/// 按组矢量方向和组角色统一骨干线首尾点。
/// 进入组要求首点位于尾点左侧，退出组要求首点位于尾点右侧。
bool orientByGroupDirection(Vec3d& first, Vec3d& last,
                            const Vec2d& group_direction, GroupRole role) {
    if (group_direction.norm() < 1e-10) return false;
    const Vec2d first_from_last = first.xy() - last.xy();
    const double side = group_direction.normalized().x() * first_from_last.y() -
                        group_direction.normalized().y() * first_from_last.x();
    if ((role == GroupRole::Entry && side < -1e-10) ||
        (role == GroupRole::Exit && side > 1e-10)) {
        std::swap(first, last);
        return true;
    }
    return false;
}

/// 按来源骨干线的原始向量方向决定端点/中点连接线的首尾。
bool orientBySourceSpine(Vec3d& first, Vec3d& last, const Primitive& source) {
    const Vec2d source_direction = source.last.xy() - source.first.xy();
    if (source_direction.norm() > 1e-10 &&
        (last.xy() - first.xy()).dot(source_direction) < 0.0) {
        std::swap(first, last);
        return true;
    }
    return false;
}

/// 计算多边形外环的平面面积质心，并以外环顶点高程的平均值作为 Z。
/// 空多边形返回原点；面积退化时改用外环顶点的算术平均，避免除以零。
/// @param polygon 输入多边形；洞不参与本函数的中心计算。
/// @return 多边形外环的三维中心点。
Vec3d polygonCenter(const Polygon2d& polygon) {
    if (polygon.outer.empty()) return Vec3d(0, 0, 0);
    // 对非自交面使用面积质心；退化面退回顶点平均值，避免产生 NaN。
    double twice_area = 0.0;
    Vec2d weighted = Vec2d::Zero();
    double z = 0.0;
    for (size_t i = 0; i < polygon.outer.size(); ++i) {
        const Vec3d& a = polygon.outer[i];
        const Vec3d& b = polygon.outer[(i + 1) % polygon.outer.size()];
        const double cross = a.x() * b.y() - b.x() * a.y();
        twice_area += cross;
        weighted += (a.xy() + b.xy()) * cross;
        z += a.z();
    }
    if (std::abs(twice_area) > 1e-12) {
        return Vec3d(weighted / (3.0 * twice_area), z / polygon.outer.size());
    }
    Vec3d center(0, 0, 0);
    for (const Vec3d& p : polygon.outer) center += p;
    return center / static_cast<double>(polygon.outer.size());
}

/// 从多边形外环选取最长边方向，作为缺少明确方向时的主方向估计。
/// @param polygon 输入多边形。
/// @param fallback 外环没有有效边时使用的回退方向。
/// @return 单位化的最长边方向；若其退化则返回单位化的回退方向或 +X。
Vec2d polygonPrincipalDirection(const Polygon2d& polygon, const Vec2d& fallback) {
    double best_len = 0.0;
    Vec2d direction = fallback;
    for (size_t i = 0; i < polygon.outer.size(); ++i) {
        const Vec2d d = polygon.outer[(i + 1) % polygon.outer.size()].xy() -
                        polygon.outer[i].xy();
        if (d.squaredNorm() > best_len) {
            best_len = d.squaredNorm();
            direction = d;
        }
    }
    if (direction.norm() < 1e-10) direction = Vec2d(1, 0);
    return direction.normalized();
}

/// 求一条无限直线与多边形外环各边的交点，并按沿线站距排序、合并重复顶点。
/// 交点高程按所在边的端点高程线性插值；洞不参与相交计算。
/// @param polygon 提供待相交外环的多边形。
/// @param origin 无限直线上的基准点。
/// @param direction 无限直线方向，函数内部会单位化。
/// @return 按站距递增排列且已去重的边界交点；输入退化时返回空数组。
std::vector<LineHit> linePolygonHits(
    const Polygon2d& polygon, const Vec2d& origin, const Vec2d& direction) {
    std::vector<LineHit> hits;
    if (polygon.outer.size() < 2 || direction.norm() < 1e-10) return hits;
    const Vec2d d = direction.normalized();
    for (size_t i = 0; i < polygon.outer.size(); ++i) {
        const Vec3d& a3 = polygon.outer[i];
        const Vec3d& b3 = polygon.outer[(i + 1) % polygon.outer.size()];
        const Vec2d a = a3.xy();
        const Vec2d edge = b3.xy() - a;
        const double denominator = cross2d(d, edge);
        if (std::abs(denominator) < 1e-10) {
            if (std::abs(cross2d(a - origin, d)) < 1e-8) {
                hits.push_back({a3, (a - origin).dot(d)});
                hits.push_back({b3, (b3.xy() - origin).dot(d)});
            }
            continue;
        }
        const double t = cross2d(a - origin, edge) / denominator;
        const double u = cross2d(a - origin, d) / denominator;
        if (u < -1e-9 || u > 1.0 + 1e-9) continue;
        const Vec2d p = origin + t * d;
        const double z = a3.z() + (b3.z() - a3.z()) * std::max(0.0, std::min(1.0, u));
        hits.push_back({Vec3d(p, z), t});
    }
    std::sort(hits.begin(), hits.end(), [](const LineHit& a, const LineHit& b) {
        return a.station < b.station;
    });
    std::vector<LineHit> unique;
    for (const LineHit& hit : hits) {
        if (unique.empty() || dist(unique.back().point, hit.point) > 1e-7)
            unique.push_back(hit);
    }
    return unique;
}

/// 两条骨干线双向延长后的唯一交点；平行或退化线段没有唯一交点。
bool extendedSpinesIntersection(const Primitive& first, const Primitive& second,
                                Vec3d& intersection) {
    const Vec2d p = first.first.xy();
    const Vec2d r = first.last.xy() - p;
    const Vec2d q = second.first.xy();
    const Vec2d s = second.last.xy() - q;
    const double denominator = cross2d(r, s);
    if (r.squaredNorm() < 1e-20 || s.squaredNorm() < 1e-20 ||
        std::abs(denominator) < 1e-10) return false;

    const double t = cross2d(q - p, s) / denominator;
    const double u = cross2d(q - p, r) / denominator;
    const Vec2d point = p + t * r;
    const double z1 = first.first.z() + t * (first.last.z() - first.first.z());
    const double z2 = second.first.z() + u * (second.last.z() - second.first.z());
    intersection = Vec3d(point, 0.5 * (z1 + z2));
    return true;
}

/// 判断点是否在多边形外环内或边界上。
bool pointInRingOrBoundary(const Vec2d& point, const std::vector<Vec3d>& ring) {
    if (ring.size() < 3) return false;
    if (pointInPolygon(point, ring)) return true;
    for (size_t i = 0; i < ring.size(); ++i) {
        if (pointToSegment(point, ring[i], ring[(i + 1) % ring.size()]).first <= 1e-7)
            return true;
    }
    return false;
}

/// 判断点是否落在多边形面内或外环边界上，洞内部除外。
bool pointInPolygonSurface(const Vec2d& point, const Polygon2d& polygon) {
    if (!pointInRingOrBoundary(point, polygon.outer)) return false;
    for (const std::vector<Vec3d>& hole : polygon.holes) {
        if (pointInRingOrBoundary(point, hole)) return false;
    }
    return true;
}

/// 求线段与多边形各边的相交边编号，外环和洞分别编号以避免编号冲突。
std::vector<std::pair<size_t, size_t>> segmentPolygonBoundaryEdges(
    const Vec3d& start, const Vec3d& end, const Polygon2d& polygon) {
    std::vector<std::pair<size_t, size_t>> result;
    const Vec2d p = start.xy();
    const Vec2d r = end.xy() - p;
    if (r.squaredNorm() < 1e-20) return result;

    const auto check_ring = [&](const std::vector<Vec3d>& ring, size_t ring_index) {
        for (size_t i = 0; i < ring.size(); ++i) {
            const Vec2d q = ring[i].xy();
            const Vec2d s = ring[(i + 1) % ring.size()].xy() - q;
            const double denominator = cross2d(r, s);
            bool intersects = false;
            if (std::abs(denominator) > 1e-10) {
                const double t = cross2d(q - p, s) / denominator;
                const double u = cross2d(q - p, r) / denominator;
                intersects = t >= -1e-8 && t <= 1.0 + 1e-8 &&
                             u >= -1e-8 && u <= 1.0 + 1e-8;
            } else if (std::abs(cross2d(q - p, r)) <= 1e-8 * r.norm()) {
                intersects = pointToSegment(start, ring[i], ring[(i + 1) % ring.size()]).first <= 1e-7 ||
                             pointToSegment(end, ring[i], ring[(i + 1) % ring.size()]).first <= 1e-7 ||
                             pointToSegment(ring[i], start, end).first <= 1e-7 ||
                             pointToSegment(ring[(i + 1) % ring.size()], start, end).first <= 1e-7;
            }
            if (intersects) result.push_back({ring_index, i});
        }
    };

    check_ring(polygon.outer, 0);
    for (size_t i = 0; i < polygon.holes.size(); ++i)
        check_ring(polygon.holes[i], i + 1);
    std::sort(result.begin(), result.end());
    result.erase(std::unique(result.begin(), result.end()), result.end());
    return result;
}

/// 判断点是否属于安全岛面，外环边界算在面内，洞不算。
bool pointInIslandSurface(const Vec2d& point, const Polygon2d& polygon) {
    return pointInPolygonSurface(point, polygon);
}

/// 一条骨干的主环端点及其来源 Primitive。
struct TerminalEndpointOwner {
    size_t primitive_index = 0;
    bool is_first = false;

    TerminalEndpointOwner(size_t index, bool first)
        : primitive_index(index), is_first(first) {}
};

/// 一条孤立端点补全产生的有向 Gap 连接。
struct GapConnection {
    size_t from_node = 0;
    size_t to_node = 0;
    std::string id;
    /// 道路边缘Gap保留完整折线；旧Gap为空，入图时按两个节点构造直线。
    std::vector<Vec3d> points;

    GapConnection(size_t from, size_t to, const std::string& edge_id)
        : from_node(from), to_node(to), id(edge_id) {}
};

/// 根据骨干来源查找其对应的横道面或安全岛面。
const Polygon2d* backboneObjectSurface(const Primitive& spine,
                                       const IntersectionInput& input) {
    if (spine.type == NonMotorBackboneEdgeType::CrosswalkSpine) {
        const auto found = std::find_if(input.crosswalks.begin(), input.crosswalks.end(),
            [&](const Crosswalk& crosswalk) { return crosswalk.id == spine.source_id; });
        return found == input.crosswalks.end() ? nullptr : &found->geometry;
    }
    if (spine.type == NonMotorBackboneEdgeType::SafetyIslandSpine) {
        const auto found = std::find_if(input.safety_islands.begin(), input.safety_islands.end(),
            [&](const SafetyIsland& island) { return island.id == spine.source_id; });
        return found == input.safety_islands.end() ? nullptr : &found->geometry;
    }
    return nullptr;
}

// 道路边缘图、精确外扩面裁剪和有界补连，独立于原骨干提取规则。
#include "nonmotor_road_gap.inc"

/// 从两个最近的主环端点沿各自骨干方向求延长交点；短 Gap 可双向延长。
/// 默认仅接受端点向外的射线交点，允许反向时接受两条无限直线交点；平行或退化时失败。
/// @param first 第一条骨干线及其端点角色。
/// @param first_is_first 第一条骨干线是否取首端点。
/// @param second 第二条骨干线及其端点角色。
/// @param second_is_first 第二条骨干线是否取首端点。
/// @param allow_reverse 是否允许任一端点反向延长。
/// @param intersection 返回两条延长线的交点，高程沿原骨干线外推后取均值。
/// @param extension_distance 返回两个端点到交点的平面距离之和，反向距离也取正值。
/// @return 存在唯一且满足方向要求的交点时返回 true。
bool extendedTerminalIntersection(const Primitive& first, bool first_is_first,
                                  const Primitive& second, bool second_is_first,
                                  bool allow_reverse,
                                  Vec3d& intersection, double& extension_distance) {
    const Vec3d& first_end = first_is_first ? first.first : first.last;
    const Vec3d& first_other = first_is_first ? first.last : first.first;
    const Vec3d& second_end = second_is_first ? second.first : second.last;
    const Vec3d& second_other = second_is_first ? second.last : second.first;
    const Vec2d p = first_end.xy();
    const Vec2d r = p - first_other.xy();
    const Vec2d q = second_end.xy();
    const Vec2d s = q - second_other.xy();
    const double denominator = cross2d(r, s);
    if (r.squaredNorm() < 1e-20 || s.squaredNorm() < 1e-20 ||
        std::abs(denominator) < 1e-10) return false;

    double t = cross2d(q - p, s) / denominator;
    double u = cross2d(q - p, r) / denominator;
    if (!allow_reverse) {
        if (t < -1e-9 || u < -1e-9) return false;
        t = std::max(0.0, t);
        u = std::max(0.0, u);
    }
    const Vec2d point = p + t * r;
    const double first_z = first_end.z() + t * (first_end.z() - first_other.z());
    const double second_z = second_end.z() + u * (second_end.z() - second_other.z());
    intersection = Vec3d(point, 0.5 * (first_z + second_z));
    extension_distance = std::abs(t) * r.norm() + std::abs(u) * s.norm();
    return true;
}

/// 计算点到安全岛外环及洞边界的最短距离。
double islandBoundaryDistance(const Vec2d& point, const Polygon2d& polygon) {
    double distance = std::numeric_limits<double>::infinity();
    const auto update_ring_distance = [&](const std::vector<Vec3d>& ring) {
        for (size_t i = 0; i < ring.size(); ++i)
            distance = std::min(distance,
                pointToSegment(point, ring[i], ring[(i + 1) % ring.size()]).first);
    };
    update_ring_distance(polygon.outer);
    for (const std::vector<Vec3d>& hole : polygon.holes) update_ring_distance(hole);
    return distance;
}

/// 检查安全岛中心到连线投影策略，并要求结果在岛内且离边界至少指定距离。
bool safeIslandCenterCandidate(const Vec3d& candidate, const Polygon2d& island) {
    return pointInIslandSurface(candidate.xy(), island) &&
           islandBoundaryDistance(candidate.xy(), island) + kPointTolerance >=
               kParallelCenterBoundaryClearance;
}

/// 以最近端点连线的中点为中心；中点不安全时尝试投影安全岛面中心。
/// 投影仍不在岛内或离边界过近时返回 false，由调用方执行原中心逻辑。
bool parallelPairIslandCenter(const Vec3d& start, const Vec3d& end,
                              const Polygon2d& island, Vec3d& center) {
    const Vec3d midpoint = 0.5 * (start + end);
    if (safeIslandCenterCandidate(midpoint, island)) {
        center = midpoint;
        return true;
    }

    const Vec2d segment = end.xy() - start.xy();
    if (segment.squaredNorm() < 1e-20) return false;
    const Vec2d island_center = polygonCenter(island).xy();
    const double t = std::max(0.0, std::min(1.0,
        (island_center - start.xy()).dot(segment) / segment.squaredNorm()));
    const double z = start.z() + t * (end.z() - start.z());
    const Vec3d projected(start.xy() + t * segment, z);
    if (!safeIslandCenterCandidate(projected, island)) return false;
    center = projected;
    return true;
}

/// 沿横道骨干两端向外延长，选择首先遇到安全岛的一端及该方向最远的外环交点。
/// @param spine 邻接横道骨干线。
/// @param island 安全岛面。
/// @param endpoint 输出选中的横道端点。
/// @param far_hit 输出延长方向上的最远边界交点。
/// @return 找到非零长度的延长线时返回 true。
bool extendedSpineThroughIsland(const Primitive& spine, const Polygon2d& island,
                                Vec3d& endpoint, Vec3d& far_hit) {
    const Vec3d ends[] = {spine.first, spine.last};
    double nearest_hit = std::numeric_limits<double>::infinity();
    bool found = false;
    for (int i = 0; i < 2; ++i) {
        const Vec2d direction = ends[i].xy() - ends[1 - i].xy();
        const std::vector<LineHit> hits = linePolygonHits(island, ends[i].xy(), direction);
        double first_station = std::numeric_limits<double>::infinity();
        const LineHit* last = nullptr;
        for (const LineHit& hit : hits) {
            if (hit.station < -1e-7) continue;
            first_station = std::min(first_station, std::max(0.0, hit.station));
            if (hit.station > kPointTolerance) last = &hit;
        }
        if (last && first_station < nearest_hit) {
            nearest_hit = first_station;
            endpoint = ends[i];
            far_hit = last->point;
            found = true;
        }
    }
    return found;
}

/// 判断线段是否接触安全岛面或其外环边界，完全位于洞内的线段不算相交。
/// @param first 线段起点。
/// @param last 线段终点。
/// @param island 安全岛面。
/// @return 非退化线段接触安全岛面时返回 true。
bool segmentTouchesIsland(const Vec3d& first, const Vec3d& last, const Polygon2d& island) {
    const Vec2d direction = last.xy() - first.xy();
    const double length = direction.norm();
    if (length <= kPointTolerance) return false;
    if (pointInIslandSurface(first.xy(), island) || pointInIslandSurface(last.xy(), island))
        return true;
    for (const LineHit& hit : linePolygonHits(island, first.xy(), direction)) {
        if (hit.station >= -1e-7 && hit.station <= length + 1e-7) return true;
    }
    return false;
}

/// 将经过给定点和方向的无限直线裁剪到多边形外环范围内。
/// 若边界交点不足两个，返回同一个回退点构成的零长度线段，供调用方识别失败。
/// @param polygon 待裁剪的多边形。
/// @param origin 无限直线上的基准点。
/// @param direction 无限直线方向。
/// @return 多边形边界上的两个端点，顺序沿方向从负站距到正站距排列。
std::pair<Vec3d, Vec3d> clippedLine(
    const Polygon2d& polygon, const Vec2d& origin, const Vec2d& direction) {
    const std::vector<LineHit> hits = linePolygonHits(polygon, origin, direction);
    if (hits.size() >= 2) return std::make_pair(hits.front().point, hits.back().point);
    const Vec3d fallback(origin.x(), origin.y(), polygonCenter(polygon).z());
    return std::make_pair(fallback, fallback);
}

/// 计算两个多边形外环之间的最短平面距离；相交、接触或互相包含时距离为零。
/// @param a 第一个多边形。
/// @param b 第二个多边形。
/// @return 外环之间的最短距离；任一外环为空时返回正无穷。
double polygonDistance(const Polygon2d& a, const Polygon2d& b) {
    if (a.outer.empty() || b.outer.empty()) return std::numeric_limits<double>::infinity();
    for (const Vec3d& p : a.outer)
        if (pointInPolygon(p, b.outer)) return 0.0;
    for (const Vec3d& p : b.outer)
        if (pointInPolygon(p, a.outer)) return 0.0;
    // 两个边界可能在边段内部相交，不能只依赖顶点包含关系。
    for (size_t i = 0; i < a.outer.size(); ++i) {
        const Vec3d& a0 = a.outer[i];
        const Vec3d& a1 = a.outer[(i + 1) % a.outer.size()];
        for (size_t j = 0; j < b.outer.size(); ++j) {
            const Vec3d& b0 = b.outer[j];
            const Vec3d& b1 = b.outer[(j + 1) % b.outer.size()];
            if (segmentsIntersectStrict(a0, a1, b0, b1) ||
                pointToSegment(a0, b0, b1).first < 1e-9 ||
                pointToSegment(b0, a0, a1).first < 1e-9) return 0.0;
        }
    }
    double best = std::numeric_limits<double>::infinity();
    for (size_t i = 0; i < a.outer.size(); ++i) {
        const Vec3d& p = a.outer[i];
        for (size_t j = 0; j < b.outer.size(); ++j) {
            best = std::min(best, pointToSegment(p, b.outer[j],
                b.outer[(j + 1) % b.outer.size()]).first);
        }
    }
    for (size_t i = 0; i < b.outer.size(); ++i) {
        const Vec3d& p = b.outer[i];
        for (size_t j = 0; j < a.outer.size(); ++j) {
            best = std::min(best, pointToSegment(p, a.outer[j],
                a.outer[(j + 1) % a.outer.size()]).first);
        }
    }
    return best;
}

/// 判断一个横道面外扩 1 米后是否与另一个横道面相交。
/// 多边形相交、接触或互相包含时同样返回 true。
bool crosswalkExpandedIntersects(const Crosswalk& first, const Crosswalk& second) {
    return polygonDistance(first.geometry, second.geometry) <=
           kCrosswalkMergeDistance + kPointTolerance;
}

/// 求两条横道骨干线端点间的最短连接线段。
/// @param first 第一条横道骨干线。
/// @param second 第二条横道骨干线。
/// @param start 输出第一条骨干线上的最近端点。
/// @param end 输出第二条骨干线上的最近端点。
/// @return 最近端点间距。
double nearestSpineEndpointSegment(const Primitive& first, const Primitive& second,
                                  Vec3d& start, Vec3d& end) {
    const Vec3d first_ends[] = {first.first, first.last};
    const Vec3d second_ends[] = {second.first, second.last};
    double best_distance = std::numeric_limits<double>::infinity();
    for (const Vec3d& first_end : first_ends) {
        for (const Vec3d& second_end : second_ends) {
            const double distance = dist(first_end, second_end);
            if (distance < best_distance) {
                best_distance = distance;
                start = first_end;
                end = second_end;
            }
        }
    }
    return best_distance;
}

/// 判断两条横道骨干线的无向夹角是否严格小于指定阈值。
/// @param first 第一条横道骨干线。
/// @param second 第二条横道骨干线。
/// @param tolerance 最大无向夹角，单位为弧度。
/// @return 两条线均有效且夹角严格小于阈值时返回 true。
bool spinesNearlyParallel(const Primitive& first, const Primitive& second, double tolerance) {
    const Vec2d a = first.last.xy() - first.first.xy();
    const Vec2d b = second.last.xy() - second.first.xy();
    return a.norm() > 1e-10 && b.norm() > 1e-10 &&
           std::abs(a.normalized().dot(b.normalized())) > std::cos(tolerance);
}

/// 过滤与其他候选横道外扩面相交且距离安全岛更远的误检横道。
std::vector<const Primitive*> filterAdjacentSpines(
    const std::vector<std::pair<const Crosswalk*, const Primitive*>>& candidates,
    const SafetyIsland& island) {
    std::vector<double> island_distances;
    island_distances.reserve(candidates.size());
    for (const auto& candidate : candidates)
        island_distances.push_back(polygonDistance(island.geometry, candidate.first->geometry));

    std::vector<bool> discarded(candidates.size(), false);
    for (size_t i = 0; i < candidates.size(); ++i) {
        for (size_t j = i + 1; j < candidates.size(); ++j) {
            if (!crosswalkExpandedIntersects(*candidates[i].first, *candidates[j].first)) continue;
            if (island_distances[i] > island_distances[j] + kPointTolerance)
                discarded[i] = true;
            else if (island_distances[j] > island_distances[i] + kPointTolerance)
                discarded[j] = true;
        }
    }

    std::vector<const Primitive*> result;
    for (size_t i = 0; i < candidates.size(); ++i)
        if (!discarded[i]) result.push_back(candidates[i].second);
    return result;
}

/// 判断两个方向是否可归入同一无向方向组，正反方向按等价方向处理。
/// @param a 第一个方向。
/// @param b 第二个方向。
/// @return 两方向夹角不超过 30 度时返回 true；任一方向退化时返回 false。
bool closeUndirected(const Vec2d& a, const Vec2d& b) {
    if (a.norm() < 1e-10 || b.norm() < 1e-10) return false;
    return std::acos(std::max(-1.0, std::min(1.0,
        std::abs(a.normalized().dot(b.normalized()))))) <= kDirectionTolerance;
}

/// 翻转方向符号，使其与参考方向同向，便于无向方向进行向量均值计算。
/// @param direction 待对齐方向。
/// @param reference 作为符号基准的参考方向。
/// @return 与参考方向同向的单位方向；输入方向退化时回退为 +X。
Vec2d alignedDirection(const Vec2d& direction, const Vec2d& reference) {
    Vec2d result = direction.norm() > 1e-10 ? direction.normalized() : Vec2d(1, 0);
    if (result.dot(reference) < 0.0) result = -result;
    return result;
}

/// 判断指定车道是否至少出现在一条标记为直行的连通关系中。
/// @param input 路口输入数据。
/// @param lane_id 待检查车道 ID。
/// @return 车道作为直行连通关系的进入或退出车道时返回 true。
bool isStraightLane(const IntersectionInput& input, const LaneId& lane_id) {
    for (const Connectivity& conn : input.connectivities) {
        if (conn.turn_type != ConnTurnType::Straight) continue;
        if (conn.entry_lane_id == lane_id || conn.exit_lane_id == lane_id) return true;
    }
    return false;
}

/// 车道在路口侧的端点、行驶方向及其所属组角色。
struct LaneDirectionSample {
    Vec2d point;
    Vec2d direction;
    GroupRole role;

    LaneDirectionSample(const Vec2d& p, const Vec2d& d, GroupRole r)
        : point(p), direction(d), role(r) {}
};

/// 收集搜索半径内的车道组方向；形状提取只收集直行车道，首尾定向回退允许转弯车道。
/// @param input 路口输入数据，包含车道组、车道和连通关系。
/// @param center 搜索中心点。
/// @param straight_only 是否只保留参与直行连通关系的车道。
/// @return 路口侧端点距离中心不超过 kLaneSearchRadius 的方向样本。
std::vector<LaneDirectionSample> nearbyLaneDirections(
    const IntersectionInput& input, const Vec2d& center, bool straight_only) {
    std::vector<LaneDirectionSample> result;
    for (const LaneGroup& group : input.lane_groups) {
        for (const LaneId& lane_id : group.lanes) {
            const Lane* lane = input.findLane(lane_id);
            if (!lane || lane->geometry.points.empty() ||
                (straight_only && !isStraightLane(input, lane_id))) continue;
            const std::pair<Vec2d, Vec2d> endpoint_direction =
                group.role == GroupRole::Entry ? input.entryPtDir(lane_id)
                                               : input.exitPtDir(lane_id);
            if (dist(endpoint_direction.first, center) <= kLaneSearchRadius)
                result.push_back({endpoint_direction.first, endpoint_direction.second, group.role});
        }
    }
    return result;
}

/// 根据附近直行车道方向射线与横道的交点提取横道骨干线。
/// 每条有效射线与横道相交形成弦，先对有效直行方向求平均并取其垂向，
/// 再选择垂向投影最远的两条有效弦，以两个弦中点的实际连线确定骨干方向，
/// 再延长并裁剪到横道边界。有效中点不足两个时，
/// 使用可用直行车道方向的垂向；没有有效车道方向时回退到横道外环最长边方向。
/// 首尾定向优先使用有效直行组；缺失时单独搜索附近横穿横道、且与骨干近似垂直的车道组，
/// 不限制其转向类型，优先进入组并在同角色内选择最近端点；此回退只交换首尾，不改变裁剪几何。
/// 进入组令首点在尾点左侧，退出组令首点在尾点右侧；完全没有有效组时保留几何首尾。
/// @param crosswalk 待提取骨干线的人行横道。
/// @param input 提供横道附近直行形状样本和车道组定向候选的路口输入数据。
/// @return 横道骨干几何及用于统一方向的组矢量和角色；无法裁剪时两点重合。
SpineGeometry crosswalkSpine(
    const Crosswalk& crosswalk, const IntersectionInput& input) {
    const Vec3d center3 = polygonCenter(crosswalk.geometry);
    const std::vector<LaneDirectionSample> lanes =
        nearbyLaneDirections(input, center3.xy(), true);

    std::vector<Vec2d> midpoints;
    std::vector<Vec2d> valid_lane_directions;
    std::vector<LaneDirectionSample> valid_lane_samples;
    for (const LaneDirectionSample& lane : lanes) {
        if (lane.direction.norm() < 1e-10) continue;
        const std::vector<LineHit> hits = linePolygonHits(
            crosswalk.geometry, lane.point, lane.direction);
        if (hits.size() < 2) continue;
        valid_lane_samples.push_back(lane);
        const Vec2d midpoint = 0.5 * (hits.front().point.xy() + hits.back().point.xy());
        bool duplicate = false;
        for (const Vec2d& existing : midpoints)
            duplicate = duplicate || dist(existing, midpoint) < 1e-5;
        if (!duplicate) {
            midpoints.push_back(midpoint);
            valid_lane_directions.push_back(lane.direction.normalized());
        }
    }

    Vec2d line_origin = center3.xy();
    Vec2d spine_direction(0, 0);
    Vec2d group_perpendicular(0, 0);
    SpineGeometry result;
    if (!valid_lane_samples.empty()) {
        const auto reference = std::find_if(valid_lane_samples.begin(), valid_lane_samples.end(),
            [](const LaneDirectionSample& sample) { return sample.role == GroupRole::Entry; });
        const LaneDirectionSample& basis = reference == valid_lane_samples.end()
                                                ? valid_lane_samples.front() : *reference;
        result.group_direction = basis.direction.normalized();
        result.group_role = basis.role;
        result.has_group_direction = true;
    }
    if (!valid_lane_directions.empty()) {
        Vec2d group_direction = Vec2d::Zero();
        const Vec2d reference = result.has_group_direction
                                    ? result.group_direction : valid_lane_directions.front();
        for (const Vec2d& direction : valid_lane_directions)
            group_direction += alignedDirection(direction, reference);
        if (group_direction.norm() > 1e-10)
            group_perpendicular = rotLeft(group_direction.normalized());
    }
    if (midpoints.size() >= 2) {
        Vec2d selection_axis = group_perpendicular;
        if (selection_axis.norm() < 1e-10)
            selection_axis = polygonPrincipalDirection(crosswalk.geometry, Vec2d(1, 0));
        auto first_midpoint = std::min_element(midpoints.begin(), midpoints.end(),
            [&](const Vec2d& a, const Vec2d& b) {
                return a.dot(selection_axis) < b.dot(selection_axis);
            });
        auto last_midpoint = std::max_element(midpoints.begin(), midpoints.end(),
            [&](const Vec2d& a, const Vec2d& b) {
                return a.dot(selection_axis) < b.dot(selection_axis);
            });
        if (first_midpoint != midpoints.end() && last_midpoint != midpoints.end() &&
            dist(*first_midpoint, *last_midpoint) > 1e-10) {
            line_origin = 0.5 * (*first_midpoint + *last_midpoint);
            spine_direction = *last_midpoint - *first_midpoint;
        }
    }

    if (spine_direction.norm() < 1e-10 && group_perpendicular.norm() > 1e-10)
        spine_direction = group_perpendicular;
    if (spine_direction.norm() < 1e-10)
        spine_direction = polygonPrincipalDirection(crosswalk.geometry, Vec2d(1, 0));

    const std::pair<Vec3d, Vec3d> clipped = clippedLine(
        crosswalk.geometry, line_origin, spine_direction);
    result.first = clipped.first;
    result.last = clipped.second;
    if (!result.has_group_direction && dist(result.first, result.last) > kPointTolerance) {
        // T形路口支路可能没有直行连通，仍需按实际车道组方向定向；
        // 只接纳延长线横穿横道且与已提取骨干的垂向偏差不超过30度的候选，避免邻路误定向。
        const std::vector<LaneDirectionSample> orientation_samples =
            nearbyLaneDirections(input, center3.xy(), false);
        const Vec2d spine_axis = (result.last.xy() - result.first.xy()).normalized();
        const LaneDirectionSample* basis = nullptr;
        double nearest_distance = std::numeric_limits<double>::infinity();
        for (const LaneDirectionSample& sample : orientation_samples) {
            if (sample.direction.norm() < 1e-10) continue;
            if (!closeUndirected(rotLeft(sample.direction), spine_axis)) continue;
            const std::vector<LineHit> hits =
                linePolygonHits(crosswalk.geometry, sample.point, sample.direction);
            if (hits.size() < 2 || dist(hits.front().point, hits.back().point) <= kPointTolerance)
                continue;
            if (basis && basis->role == GroupRole::Entry && sample.role != GroupRole::Entry)
                continue;
            const double sample_distance = dist(sample.point, center3.xy());
            if (!basis || (sample.role == GroupRole::Entry && basis->role != GroupRole::Entry) ||
                sample_distance < nearest_distance) {
                basis = &sample;
                nearest_distance = sample_distance;
            }
        }
        if (basis) {
            result.group_direction = basis->direction.normalized();
            result.group_role = basis->role;
            result.has_group_direction = true;
        }
    }
    if (result.has_group_direction)
        orientByGroupDirection(result.first, result.last, result.group_direction, result.group_role);
    return result;
}

/// 安全岛方向组：保存同向车道方向的均值和对应车道端点代表位置。
struct IslandDirectionGroup {
    /// 方向组的单位化平均行驶方向。
    Vec2d direction{1, 0};
    /// 方向组内车道端点的代表位置，用于确定安全岛骨干的近侧交点。
    Vec2d representative_point{0, 0};
    /// 产生该方向的进入/退出组角色。
    GroupRole role = GroupRole::Entry;
    /// 方向是否来自有效车道组；纯几何回退不应用组角色定向。
    bool has_group_direction = true;

    /// 以平均行驶方向及车道端点代表位置初始化方向组。
    /// @param d 方向组的平均行驶方向。
    /// @param p 方向组车道端点的代表位置。
    /// @param r 方向组角色。
    IslandDirectionGroup(const Vec2d& d, const Vec2d& p, GroupRole r,
                         bool has_direction = true)
        : direction(d), representative_point(p), role(r), has_group_direction(has_direction) {}
};

/// 按车道组收集安全岛附近的直行方向，并将相差不超过 30 度的无向方向合并。
/// 同向方向优先保留进入组角色；仅有退出组方向时保留退出组角色。
/// @param island 当前安全岛；其面积中心作为车道端点搜索中心。
/// @param input 路口输入数据，提供车道组、车道端点和直行连通关系。
/// @return 每个方向组的平均方向和代表端点；未找到有效车道方向时返回空数组。
std::vector<IslandDirectionGroup> islandDirectionGroups(
    const SafetyIsland& island, const IntersectionInput& input) {
    const Vec2d center = polygonCenter(island.geometry).xy();
    std::vector<IslandDirectionGroup> groups;
    for (GroupRole role : {GroupRole::Entry, GroupRole::Exit}) {
        for (const LaneGroup& lane_group : input.lane_groups) {
            if (lane_group.role != role) continue;
            std::vector<std::pair<Vec2d, Vec2d>> samples;
            for (const LaneId& lane_id : lane_group.lanes) {
                const Lane* lane = input.findLane(lane_id);
                if (!lane || lane->geometry.points.empty()) continue;
                if (!isStraightLane(input, lane_id)) continue;
                const std::pair<Vec2d, Vec2d> sample =
                    lane_group.role == GroupRole::Entry ? input.entryPtDir(lane_id)
                                                        : input.exitPtDir(lane_id);
                if (dist(sample.first, center) <= kLaneSearchRadius) samples.push_back(sample);
            }
            if (samples.empty()) continue;
            Vec2d direction = samples.front().second.normalized();
            Vec2d point = Vec2d::Zero();
            for (const auto& sample : samples) {
                direction += alignedDirection(sample.second, direction);
                point += sample.first;
            }
            if (direction.norm() < 1e-10) continue;
            direction.normalize();
            point /= static_cast<double>(samples.size());
            bool merged = false;
            for (IslandDirectionGroup& group : groups) {
                if (!closeUndirected(group.direction, direction)) continue;
                group.direction = alignedDirection(group.direction + direction, group.direction).normalized();
                group.representative_point = 0.5 * (group.representative_point + point);
                merged = true;
                break;
            }
            if (!merged) groups.push_back({direction, point, role});
        }
    }
    return groups;
}

/// 将横道或安全岛骨干线以带Z的WKT写入调试CSV，便于核查提取结果。
/// @param kind 骨干类别中文名称，例如“人行横道骨干线”或“安全岛骨干线”。
/// @param source_id 生成该线段的横道或安全岛 ID。
/// @param first 骨干线起点。
/// @param last 骨干线终点。
void NonMotorBackboneExtractor::logSpineWkt(
    const std::string& kind, const std::string& source_id, const Vec3d& first, const Vec3d& last) const {
    std::ostringstream wkt;
    wkt << std::setprecision(12) << "LINESTRING Z ("
        << first.x() << " " << first.y() << " " << first.z()
        << ", " << last.x() << " " << last.y() << " " << last.z()
        << ")";
    //std::cout << kind << " source_id=" << source_id << " WKT=" << wkt.str() << std::endl;
    fs_wkt_ << kind << ", " << source_id << ", \"" << wkt.str() << "\"" << std::endl;
}
void NonMotorBackboneExtractor::logSpineWkt(
        const std::string& kind, const std::string& source_id, const std::vector<Vec3d>& pts) const {
    std::ostringstream wkt;
    wkt << std::setprecision(12) << "LINESTRING Z (";
    for (int i = 0; i < pts.size(); ++i) {
        if (i != 0) wkt << ", ";
        wkt << pts[i].x() << ' ' << pts[i].y() << ' ' << pts[i].z();
    }
    wkt << ")";
    std::cout << kind << " source_id=" << source_id << " WKT=" << wkt.str() << std::endl;
    fs_wkt_ << kind << ", " << source_id << ", \"" << wkt.str() << "\"" << std::endl;
}

///////////////////////////////////////////////////////////////////
NonMotorBackboneExtractor::NonMotorBackboneExtractor(std::string out_dir,
                                                   NonMotorGapOptions gap_options) {
    out_dir_ = out_dir;
    gap_options_ = gap_options;
}
NonMotorBackboneExtractor::~NonMotorBackboneExtractor() {
    if (fs_wkt_.is_open()) fs_wkt_.close();
}
NonMotorBackbone NonMotorBackboneExtractor::run(const IntersectionInput& input) const {
    fs_wkt_.open(out_dir_ + "/nonmotor_backbone.csv");

    NonMotorBackbone bb;
    bb.id = input.id.empty() ? "nonmotor_backbone" : input.id + "_nonmotor";

    std::vector<Primitive> primitives;
    /// 获取人行横道骨干线
    for (const Crosswalk& crosswalk : input.crosswalks) {
        if (crosswalk.geometry.outer.size() < 3) continue;
        const SpineGeometry spine = crosswalkSpine(crosswalk, input);
        if (dist(spine.first, spine.last) <= kPointTolerance) continue;
        logSpineWkt("人行横道骨干线", crosswalk.id, spine.first, spine.last);
        primitives.push_back({"cw_" + crosswalk.id + "_spine", crosswalk.id,
                              NonMotorBackboneEdgeType::CrosswalkSpine,
                              spine.first, spine.last, true, true,
                              spine.group_direction, spine.group_role, spine.has_group_direction});
    }

    /// 获取安全岛骨干线
    for (const SafetyIsland& island : input.safety_islands) {
        if (island.geometry.outer.size() < 3) continue;
        // 获取安全岛相邻接的人行横道集合
        std::vector<std::pair<const Crosswalk*, const Primitive*>> adjacent_candidates;
        for (const Crosswalk& crosswalk : input.crosswalks) {
            if (polygonDistance(island.geometry, crosswalk.geometry) > kSearchRadius) continue;
            const auto edge_it = std::find_if(primitives.begin(), primitives.end(),
                [&](const Primitive& primitive) { return primitive.source_id == crosswalk.id; });
            if (edge_it == primitives.end()) continue;
            adjacent_candidates.push_back(std::make_pair(&crosswalk, &*edge_it));
        }
        // 邻接横道面两两检查 1 米外扩相交关系，过滤距离安全岛更远的误检面。
        const std::vector<const Primitive*> adjacent_spines =
            filterAdjacentSpines(adjacent_candidates, island);

        // 单横道优先沿其朝向安全岛的端点延伸至远侧边界。
        if (adjacent_spines.size() == 1) {
            Vec3d endpoint, far_hit;
            if (extendedSpineThroughIsland(*adjacent_spines.front(), island.geometry,
                                           endpoint, far_hit)) {
                const Primitive& source = *adjacent_spines.front();
                orientBySourceSpine(endpoint, far_hit, source);
                logSpineWkt("安全岛骨干线", island.id, endpoint, far_hit);
                primitives.push_back({"si_" + island.id + "_spine", island.id,
                                      NonMotorBackboneEdgeType::SafetyIslandSpine,
                                      endpoint, far_hit, true, true,
                                      source.group_direction, source.group_role,
                                      source.has_group_direction});
                continue;
            }
        }

        // 安全岛中心优先级：三条及以上邻接横道存在小于10度的骨干平行对时，
        // 取端点连线最短的一对；其中点离岛边界不足0.5米时，改用岛中心在线段上的投影，
        // 投影也不安全时沿用最近两条骨干线交点逻辑。
        // 原逻辑缺少有效路口面、交点退化或交点越出安全岛时使用安全岛面中心。
        Vec3d center = polygonCenter(island.geometry);
        bool has_intersection_center = false;
        bool has_parallel_pair_center = false;
        if (adjacent_spines.size() >= 3) {
            double shortest_connection = std::numeric_limits<double>::infinity();
            Vec3d nearest_start, nearest_end;
            for (size_t i = 0; i < adjacent_spines.size(); ++i) {
                for (size_t j = i + 1; j < adjacent_spines.size(); ++j) {
                    if (!spinesNearlyParallel(*adjacent_spines[i], *adjacent_spines[j],
                                              kMultiCrosswalkParallelTolerance)) continue;
                    Vec3d start, end;
                    const double distance = nearestSpineEndpointSegment(
                        *adjacent_spines[i], *adjacent_spines[j], start, end);
                    if (distance < shortest_connection) {
                        shortest_connection = distance;
                        nearest_start = start;
                        nearest_end = end;
                    }
                }
            }
            if (std::isfinite(shortest_connection))
                has_parallel_pair_center = parallelPairIslandCenter(
                    nearest_start, nearest_end, island.geometry, center);
        }
        if (!has_parallel_pair_center && adjacent_spines.size() >= 2 &&
            input.area.geometry.outer.size() >= 3) {
            const Vec3d intersection_center = polygonCenter(input.area.geometry);
            std::vector<const Primitive*> nearest_spines = adjacent_spines;
            std::sort(nearest_spines.begin(), nearest_spines.end(),
                [&](const Primitive* first, const Primitive* second) {
                    const double first_distance = pointToSegment(
                        intersection_center, first->first, first->last).first;
                    const double second_distance = pointToSegment(
                        intersection_center, second->first, second->last).first;
                    if (first_distance != second_distance)
                        return first_distance < second_distance;
                    return first->source_id < second->source_id;
                });
            Vec3d intersection;
            if (extendedSpinesIntersection(*nearest_spines[0], *nearest_spines[1], intersection) &&
                pointInIslandSurface(intersection.xy(), island.geometry)) {
                center = intersection;
                has_intersection_center = true;
            }
        }
        // 双横道无法取得有效交点中心时，近似平行且连线接触安全岛则直接连通。
        if (adjacent_spines.size() == 2 && !has_intersection_center) {
            const Primitive& first = *adjacent_spines[0];
            const Primitive& second = *adjacent_spines[1];
            if (spinesNearlyParallel(first, second, kDirectionTolerance)) {
                Vec3d start, end;
                nearestSpineEndpointSegment(first, second, start, end);
                if (segmentTouchesIsland(start, end, island.geometry)) {
                    orientBySourceSpine(start, end, first);
                    logSpineWkt("安全岛骨干线", island.id, start, end);
                    primitives.push_back({"si_" + island.id + "_spine", island.id,
                                          NonMotorBackboneEdgeType::SafetyIslandSpine,
                                          start, end, true, true,
                                          first.group_direction, first.group_role,
                                          first.has_group_direction});
                    continue;
                }
            }
        }

        // 获取安全岛骨干线：
        // 1.连接邻接的人行横道
        if (!adjacent_spines.empty()) {
            // 获取安全岛中心与邻接人行横道骨干线最近的端点
            struct AdjacentEndpoint {
                Vec3d point;
                const Primitive* spine;
            };
            std::vector<AdjacentEndpoint> adjacent;
            for (const Primitive* spine : adjacent_spines) {
                const Vec3d endpoint = dist(center, spine->first) <= dist(center, spine->last)
                                           ? spine->first : spine->last;
                adjacent.push_back({endpoint, spine});
            }
            // 连接岛中心与邻接人行横道的最近端点构成安全岛骨干线
            if (!adjacent.empty()) {
                for (size_t i = 0; i < adjacent.size(); ++i) {
                    if (dist(center, adjacent[i].point) <= kPointTolerance) continue;
                    Vec3d first = center;
                    Vec3d last = adjacent[i].point;
                    const bool reversed = orientBySourceSpine(first, last, *adjacent[i].spine);
                    // 只有沿用平行对中点/投影作为安全岛中心时固定中心端点，避免主环平滑
                    // 将已经通过岛内及边界净距检查的中心位置再次移出安全区域。
                    logSpineWkt("安全岛骨干线", island.id, first, last);
                    primitives.push_back({"si_" + island.id + "_spine_" + std::to_string(i), island.id,
                                          NonMotorBackboneEdgeType::SafetyIslandSpine,
                                          first, last, reversed, !reversed,
                                          adjacent[i].spine->group_direction,
                                          adjacent[i].spine->group_role,
                                          adjacent[i].spine->has_group_direction});
                }
                continue;
            }
        }
        // 2.邻近组射线法：先通过组方向射线与安全岛交点确认方向，方向聚类，构建岛骨干线
        else {
            std::vector<IslandDirectionGroup> direction_groups = islandDirectionGroups(island, input);
            if (direction_groups.empty()) {
                direction_groups.push_back({polygonPrincipalDirection(island.geometry, Vec2d(1, 0)),
                                            center.xy(), GroupRole::Entry, false});
            }
            const Vec2d center_xy = center.xy();
            if (direction_groups.size() == 1) {
                const std::pair<Vec3d, Vec3d> line = clippedLine(
                    island.geometry, center_xy, rotLeft(direction_groups.front().direction));
                Vec3d first = line.first;
                Vec3d last = line.second;
                if (direction_groups.front().has_group_direction)
                    orientByGroupDirection(first, last, direction_groups.front().direction,
                                           direction_groups.front().role);
                if (dist(first, last) > kPointTolerance) {
                    logSpineWkt("安全岛骨干线", island.id, first, last);
                    primitives.push_back({"si_" + island.id + "_spine", island.id,
                                          NonMotorBackboneEdgeType::SafetyIslandSpine,
                                          first, last, true, true,
                                          direction_groups.front().direction,
                                          direction_groups.front().role,
                                          direction_groups.front().has_group_direction});
                }
            } else {
                for (size_t i = 0; i < direction_groups.size(); ++i) {
                    const IslandDirectionGroup& group = direction_groups[i];
                    const std::pair<Vec3d, Vec3d> line = clippedLine(
                        island.geometry, center_xy, rotLeft(group.direction));
                    if (dist(line.first, line.second) <= kPointTolerance) continue;
                    Vec3d endpoint = dist(group.representative_point, line.first) <=
                                                   dist(group.representative_point, line.second)
                                               ? line.first : line.second;
                    Vec3d first = center;
                    Vec3d last = endpoint;
                    const bool reversed = group.has_group_direction &&
                        orientByGroupDirection(first, last, group.direction, group.role);
                    logSpineWkt("安全岛骨干线", island.id, first, last);
                    primitives.push_back({"spine_si_" + island.id + "_" + std::to_string(i),
                        island.id, NonMotorBackboneEdgeType::SafetyIslandSpine,
                        first, last, reversed, !reversed,
                        group.direction, group.role, group.has_group_direction});
                }
            }
        }
    }

    /// 构建骨干图
    std::vector<Vec3d> node_points;
    std::vector<NonMotorBackboneNodeType> node_types;
    std::vector<std::string> node_sources;
    if (primitives.empty()) return bb;

    // 按平面坐标容差查找或新建节点，合并各骨干边共享的端点。
    // @param point 待登记节点的三维坐标。
    // @param type 节点类别；合并到既有节点时会保留已有类别，Gap 类型可被具体类别替换。
    // @param source 该节点对应的横道或安全岛来源 ID。
    // @return 既有或新建节点的稳定 ID，格式为 n<序号>。
    auto nodeFor = [&](const Vec3d& point, NonMotorBackboneNodeType type,
                       const std::string& source) -> std::string {
        for (size_t i = 0; i < node_points.size(); ++i) {
            if (dist(node_points[i], point) <= kPointTolerance) {
                if (node_types[i] == NonMotorBackboneNodeType::Gap && type != NonMotorBackboneNodeType::Gap)
                    node_types[i] = type;
                if (node_sources[i].empty()) node_sources[i] = source;
                return "n" + std::to_string(i);
            }
        }
        node_points.push_back(point);
        node_types.push_back(type);
        node_sources.push_back(source);
        return "n" + std::to_string(node_points.size() - 1);
    };

    std::vector<std::string> primitive_from;
    std::vector<std::string> primitive_to;
    for (const Primitive& primitive : primitives) {
        const NonMotorBackboneNodeType type =
            primitive.type == NonMotorBackboneEdgeType::CrosswalkSpine
                ? NonMotorBackboneNodeType::CrosswalkPt : NonMotorBackboneNodeType::SafetyIslandPt;
        primitive_from.push_back(nodeFor(primitive.first, type, primitive.source_id));
        primitive_to.push_back(nodeFor(primitive.last, type, primitive.source_id));
    }

    std::vector<std::string> terminal_nodes;
    for (size_t i = 0; i < primitives.size(); ++i) {
        if (primitives[i].first_terminal) terminal_nodes.push_back(primitive_from[i]);
        if (primitives[i].last_terminal) terminal_nodes.push_back(primitive_to[i]);
    }
    const auto sortTerminalNodes = [&]() {
        std::sort(terminal_nodes.begin(), terminal_nodes.end());
        terminal_nodes.erase(std::unique(terminal_nodes.begin(), terminal_nodes.end()), terminal_nodes.end());
        Vec2d sort_center = Vec2d::Zero();
        for (const std::string& id : terminal_nodes) {
            const size_t index = static_cast<size_t>(std::stoul(id.substr(1)));
            sort_center += node_points[index].xy();
        }
        if (!terminal_nodes.empty()) sort_center /= static_cast<double>(terminal_nodes.size());
        std::sort(terminal_nodes.begin(), terminal_nodes.end(), [&](const std::string& a, const std::string& b) {
            const size_t ia = static_cast<size_t>(std::stoul(a.substr(1)));
            const size_t ib = static_cast<size_t>(std::stoul(b.substr(1)));
            return std::atan2(node_points[ia].y() - sort_center.y(), node_points[ia].x() - sort_center.x()) <
                   std::atan2(node_points[ib].y() - sort_center.y(), node_points[ib].x() - sort_center.x());
        });
    };
    sortTerminalNodes();

    // 记录每个骨干端点的来源及其在骨干线上的方向。
    std::vector<std::vector<TerminalEndpointOwner>> endpoint_owners(node_points.size());
    for (size_t i = 0; i < primitives.size(); ++i) {
        const size_t first_index = static_cast<size_t>(std::stoul(primitive_from[i].substr(1)));
        const size_t last_index = static_cast<size_t>(std::stoul(primitive_to[i].substr(1)));
        endpoint_owners[first_index].push_back({i, true});
        endpoint_owners[last_index].push_back({i, false});
    }

    // 仅统计横道/安全岛骨干线的端点关联；连接线不参与孤立端点发现。
    std::vector<size_t> incident_spines(node_points.size(), 0);
    for (size_t i = 0; i < endpoint_owners.size(); ++i) {
        std::vector<size_t> unique_owners;
        for (const TerminalEndpointOwner& owner : endpoint_owners[i])
            if (std::find(unique_owners.begin(), unique_owners.end(), owner.primitive_index) ==
                unique_owners.end()) unique_owners.push_back(owner.primitive_index);
        incident_spines[i] = unique_owners.size();
    }

    std::vector<size_t> parent(node_points.size());
    for (size_t i = 0; i < parent.size(); ++i) parent[i] = i;
    const auto rootOf = [&](size_t node) {
        while (parent[node] != node) node = parent[node];
        return node;
    };

    // 骨干端点构成无向连通图；仅通过共享几何端点连接，并限定在孤立点的搜索半径内。
    std::vector<std::vector<size_t>> spine_neighbors(node_points.size());
    for (size_t i = 0; i < primitives.size(); ++i) {
        const size_t first_node = static_cast<size_t>(std::stoul(primitive_from[i].substr(1)));
        const size_t last_node = static_cast<size_t>(std::stoul(primitive_to[i].substr(1)));
        if (first_node == last_node) continue;
        spine_neighbors[first_node].push_back(last_node);
        spine_neighbors[last_node].push_back(first_node);
    }
    const auto connectedWithinSearchRadius = [&](size_t start_node, size_t target_node,
                                                  const Vec3d& search_origin) {
        start_node = rootOf(start_node);
        target_node = rootOf(target_node);
        if (dist(node_points[start_node], search_origin) >
                kIsolatedEndpointSearchRadius + kPointTolerance ||
            dist(node_points[target_node], search_origin) >
                kIsolatedEndpointSearchRadius + kPointTolerance) return false;

        std::vector<bool> visited(node_points.size(), false);
        std::vector<size_t> pending(1, start_node);
        visited[start_node] = true;
        for (size_t cursor = 0; cursor < pending.size(); ++cursor) {
            const size_t node = rootOf(pending[cursor]);
            if (node == target_node) return true;
            for (const size_t neighbor_id : spine_neighbors[node]) {
                const size_t neighbor = rootOf(neighbor_id);
                if (visited[neighbor] ||
                    dist(node_points[neighbor], search_origin) >
                        kIsolatedEndpointSearchRadius + kPointTolerance) continue;
                visited[neighbor] = true;
                pending.push_back(neighbor);
            }
        }
        return false;
    };

    std::vector<GapConnection> isolated_gaps;
    // 道路边缘优先于原3米Gap；只移动断头端点，避免改变共享节点上其他骨干的方向。
    const RoadGapNetwork road_network(input);
    Vec2d gap_center = Vec2d::Zero();
    if (!input.area.geometry.outer.empty()) gap_center = polygonCenter(input.area.geometry).xy();
    else {
        for (const Vec3d& p : node_points) gap_center += p.xy();
        if (!node_points.empty()) gap_center /= static_cast<double>(node_points.size());
    }
    std::vector<std::vector<Polygon2d>> expanded_surfaces(primitives.size());
    if (!road_network.segments.empty() && gap_options_.search_radius > 0)
        for (size_t i = 0; i < primitives.size(); ++i) {
            const Polygon2d* surface = backboneObjectSurface(primitives[i], input);
            if (surface) expanded_surfaces[i] = roadGapExpandedSurface(*surface, gap_options_.surface_expansion);
        }
    // 道路边缘有效段可能连接相距超过普通端点Gap搜索范围的骨干；道路优先阶段使用独立上限，
    // 但仍受search_radius<=0总开关约束，避免把普通直线补连的范围意外放大。
    const double road_edge_radius = gap_options_.search_radius > 0
        ? std::max(gap_options_.search_radius, gap_options_.road_edge_search_radius)
        : 0.0;
    const bool road_extension_enabled = road_edge_radius > gap_options_.search_radius + kPointTolerance;
    // 道路Gap的搜索半径约束端点到RoadEdge网络的局部接近范围，不能继续误用
    // 两个骨干端点之间的直线距离；长道路边缘的有效路径可能远长于该局部半径。
    const auto road_edge_near = [&](const Vec3d& point) {
        double nearest = std::numeric_limits<double>::infinity();
        for (const auto& segment : road_network.segments)
            nearest = std::min(nearest, pointToSegment(point, segment.a, segment.b).first);
        return nearest;
    };
    std::vector<bool> road_joined(node_points.size(), false);
    for (size_t node = 0; node < node_points.size() && !road_network.segments.empty(); ++node) {
        if (incident_spines[node] != 1 || road_joined[node] || gap_options_.search_radius <= 0) continue;
        const auto owner = endpoint_owners[node].front();
        struct RoadCandidate { double distance; size_t node; bool road_supported; };
        std::vector<RoadCandidate> candidates;
        double nearest_endpoint = road_edge_radius;
        for (size_t other = 0; other < node_points.size(); ++other) {
            if (node == other ||
                owner.primitive_index == endpoint_owners[other].front().primitive_index) continue;
            const double distance = dist(node_points[node], node_points[other]);
            nearest_endpoint = std::min(nearest_endpoint, distance);
            const bool road_supported = road_extension_enabled &&
                distance > gap_options_.search_radius + kPointTolerance &&
                road_edge_near(node_points[node]) <= road_edge_radius + kPointTolerance &&
                road_edge_near(node_points[other]) <= road_edge_radius + kPointTolerance;
            if (distance <= road_edge_radius + kPointTolerance || road_supported)
                candidates.push_back({distance, other, road_supported});
        }
        std::sort(candidates.begin(), candidates.end(), [](const RoadCandidate& a, const RoadCandidate& b) {
            if (a.road_supported != b.road_supported) return a.road_supported > b.road_supported;
            if (a.distance != b.distance) return a.distance < b.distance;
            return a.node < b.node;
        });
        for (const auto& candidate : candidates) {
            const size_t other = candidate.node;
            // 最近骨干不可移动/不可沿RoadEdge连接时回退原规则，不越过它劫持更远端点。
            if (!candidate.road_supported && candidate.distance > nearest_endpoint + kPointTolerance) break;
            if (incident_spines[other] != 1 || road_joined[other]) continue;
            const auto other_owner = endpoint_owners[other].front();
            double other_nearest = candidate.distance;
            for (size_t nearby = 0; nearby < node_points.size(); ++nearby) {
                if (nearby == other || endpoint_owners[nearby].front().primitive_index == other_owner.primitive_index) continue;
                other_nearest = std::min(other_nearest, dist(node_points[nearby], node_points[other]));
            }
            // 两侧必须互为最近端点；否则先保留目标端与其更近骨干的原合并关系。
            if (!candidate.road_supported && other_nearest + kPointTolerance < candidate.distance) continue;
            const size_t opposite = static_cast<size_t>(std::stoul(
                (owner.is_first ? primitive_to[owner.primitive_index] : primitive_from[owner.primitive_index]).substr(1)));
            if (connectedWithinSearchRadius(opposite, other, node_points[node])) continue;
            Primitive& first = primitives[owner.primitive_index];
            Primitive& second = primitives[other_owner.primitive_index];
            auto points = road_network.path(node_points[node], node_points[other],
                expanded_surfaces[owner.primitive_index], expanded_surfaces[other_owner.primitive_index],
                road_edge_radius);
            // 偏移后道路线端点和骨干支撑线交点各承担一份法向偏移；公共点允许的
            // 延长量相应增加两份 inward_offset。只有超出普通Gap范围的长道路连接
            // 才启用这项道路延长余量，普通候选仍严格服从显式endpoint_adjustment。
            NonMotorGapOptions road_options = gap_options_;
            if (candidate.distance > gap_options_.search_radius + kPointTolerance)
                road_options.endpoint_adjustment += 2.0 * gap_options_.inward_offset;
            if (!shiftedRoadGap(points, gap_center, first, owner.is_first, second, other_owner.is_first, road_options)) continue;
            // 改端点后的原骨干与新增Gap联合验收，禁止穿障、越界或自交。
            std::vector<Vec3d> whole{owner.is_first ? first.last : first.first};
            whole.insert(whole.end(), points.begin(), points.end());
            whole.push_back(other_owner.is_first ? second.last : second.first);
            // 骨干端点与道路边缘的既有公共接触均合法；新增Gap只能禁止其它
            // 位置的Boundary穿越，不能因检查完整骨干而误拒绝端点接触。
            const std::vector<Vec3d> allowed_contacts{
                first.first, first.last, second.first, second.last};
            if (!safeRoadGapPath(whole, input, allowed_contacts)) continue;
            bool crosses_other = false;
            for (size_t i = 0; i < primitives.size(); ++i) {
                if (i == owner.primitive_index || i == other_owner.primitive_index) continue;
                if (roadGapCrossesOther(whole, {primitives[i].first, primitives[i].last})) crosses_other = true;
            }
            for (const auto& gap : isolated_gaps)
                if (roadGapCrossesOther(whole, gap.points)) crosses_other = true;
            if (crosses_other) continue;
            node_points[node] = points.front(); node_points[other] = points.back();
            (owner.is_first ? first.first : first.last) = points.front();
            (other_owner.is_first ? second.first : second.last) = points.back();
            GapConnection gap(node, other, "gap_road_" + std::to_string(isolated_gaps.size()));
            gap.points = points;
            if (owner.is_first) {
                std::swap(gap.from_node, gap.to_node);
                std::reverse(gap.points.begin(), gap.points.end());
            }
            isolated_gaps.push_back(gap);
            road_joined[node] = road_joined[other] = true;
            spine_neighbors[node].push_back(other); spine_neighbors[other].push_back(node);
            break;
        }
    }
    std::set<std::pair<size_t, size_t>> processed_pairs;
    std::vector<bool> processed_isolated = road_joined;
    // 未找到近侧候选，或因已有骨干连通而跳过局部 Gap 的端点继续保留为孤立端点；
    // 后续大环排序不会再强制连接这些端点。
    std::set<size_t> unresolved_isolated_nodes;
    for (size_t isolated_node = 0; isolated_node < node_points.size(); ++isolated_node) {
        if (incident_spines[isolated_node] != 1 || processed_isolated[isolated_node]) continue;
        const TerminalEndpointOwner isolated_owner = endpoint_owners[isolated_node].front();
        const Primitive& isolated_spine = primitives[isolated_owner.primitive_index];
        const Vec3d& isolated_point = isolated_owner.is_first
                                          ? isolated_spine.first : isolated_spine.last;

        size_t near_node = node_points.size();
        TerminalEndpointOwner near_owner(0, false);
        double nearest_distance = std::numeric_limits<double>::infinity();
        for (size_t candidate_spine_index = 0;
             candidate_spine_index < primitives.size(); ++candidate_spine_index) {
            if (candidate_spine_index == isolated_owner.primitive_index) continue;
            const size_t candidate_nodes[] = {
                static_cast<size_t>(std::stoul(primitive_from[candidate_spine_index].substr(1))),
                static_cast<size_t>(std::stoul(primitive_to[candidate_spine_index].substr(1)))};
            for (int endpoint = 0; endpoint < 2; ++endpoint) {
                const size_t candidate_node = candidate_nodes[endpoint];
                if (candidate_node == isolated_node || road_joined[candidate_node]) continue;
                const Primitive& candidate_spine = primitives[candidate_spine_index];
                const Vec3d& candidate_point = endpoint == 0
                                                   ? candidate_spine.first : candidate_spine.last;
                const double distance = dist(isolated_point, candidate_point);
                if (distance > kIsolatedEndpointSearchRadius + kPointTolerance) continue;
                if (distance < nearest_distance) {
                    nearest_distance = distance;
                    near_node = candidate_node;
                    near_owner = TerminalEndpointOwner(candidate_spine_index, endpoint == 0);
                }
            }
        }
        if (near_node == node_points.size()) {
            unresolved_isolated_nodes.insert(isolated_node);
            continue;
        }
        const std::pair<size_t, size_t> pair_key(std::min(isolated_node, near_node),
                                                  std::max(isolated_node, near_node));
        if (!processed_pairs.insert(pair_key).second) {
            unresolved_isolated_nodes.insert(isolated_node);
            continue;
        }
        processed_isolated[isolated_node] = true;
        if (incident_spines[near_node] == 1) processed_isolated[near_node] = true;

        const Primitive& near_spine = primitives[near_owner.primitive_index];
        const Vec2d isolated_direction = isolated_spine.last.xy() - isolated_spine.first.xy();
        const Vec2d near_direction = near_spine.last.xy() - near_spine.first.xy();
        double spine_angle = 0.0;
        if (isolated_direction.norm() > 1e-10 && near_direction.norm() > 1e-10) {
            const double cosine = std::abs(isolated_direction.normalized().dot(
                near_direction.normalized()));
            spine_angle = std::acos(std::max(-1.0, std::min(1.0, cosine)));
        }

        // 短于2米的端点连线优先使用双向延长交点，允许反向且无需位于来源面内，
        // 以避免输出过短的 Gap；达到阈值时仍沿用原有射线方向、夹角和来源面约束。
        const bool short_gap = nearest_distance < kShortGapLengthThreshold;
        bool found_common_intersection = false;
        Vec3d selected_intersection;
        double shortest_extension = std::numeric_limits<double>::infinity();
        if (short_gap || spine_angle > kGapExtensionAngleThreshold) {
            const Polygon2d* isolated_surface = backboneObjectSurface(isolated_spine, input);
            for (const TerminalEndpointOwner& candidate_owner : endpoint_owners[near_node]) {
                if (candidate_owner.primitive_index == isolated_owner.primitive_index) continue;
                const Primitive& candidate_spine = primitives[candidate_owner.primitive_index];
                const Polygon2d* candidate_surface = backboneObjectSurface(candidate_spine, input);
                if (!isolated_surface || !candidate_surface) continue;

                Vec3d intersection;
                double extension_distance = 0.0;
                if (!extendedTerminalIntersection(isolated_spine, isolated_owner.is_first,
                                                  candidate_spine, candidate_owner.is_first,
                                                  short_gap, intersection, extension_distance)) continue;
                const bool in_source_surface =
                    isolated_surface && candidate_surface &&
                    (pointInPolygonSurface(intersection.xy(), *isolated_surface) ||
                     pointInPolygonSurface(intersection.xy(), *candidate_surface));
                if (!short_gap && !in_source_surface) continue;
                if (extension_distance < shortest_extension) {
                    shortest_extension = extension_distance;
                    selected_intersection = intersection;
                    found_common_intersection = true;
                }
            }
        }

        if (found_common_intersection) {
            const size_t isolated_root = rootOf(isolated_node);
            const size_t near_root = rootOf(near_node);
            if (isolated_root == near_root) continue;
            // 端点并入同一公共点后，合并两侧邻接表，供后续局部连通搜索使用。
            for (const size_t neighbor : spine_neighbors[near_root])
                spine_neighbors[isolated_root].push_back(neighbor);
            parent[near_root] = isolated_root;
            node_points[isolated_root] = selected_intersection;
        } else {
            // 若骨干线另一端能在孤立点的搜索范围内，通过共享端点和一条或多条骨干线
            // 连通到同一最近候选端点，则不额外连Gap，并保留原孤立端点。
            const size_t opposite_isolated_node = static_cast<size_t>(std::stoul(
                (isolated_owner.is_first ? primitive_to[isolated_owner.primitive_index]
                                         : primitive_from[isolated_owner.primitive_index]).substr(1)));
            if (connectedWithinSearchRadius(opposite_isolated_node, near_node, isolated_point)) {
                unresolved_isolated_nodes.insert(isolated_node);
                continue;
            }

            // 一端为骨干线尾、另一端为骨干线首时可同时延续两条线的方向；
            // 若两端的首尾角色相同而无法同时满足，则优先延续孤立端点所属骨干线的方向。
            const bool isolated_prefers_from = !isolated_owner.is_first;
            const bool from_isolated = isolated_prefers_from;
            isolated_gaps.push_back(from_isolated
                ? GapConnection(isolated_node, near_node, "gap_near_" + std::to_string(isolated_gaps.size()))
                : GapConnection(near_node, isolated_node, "gap_near_" + std::to_string(isolated_gaps.size())));
        }
    }

    // 递归平滑骨干共点：每轮只在当前仍恰有两条活动骨干共端点、且夹角小于阈值时合并；
    // 合并后的 Composite 段会重新进入下一轮，因此可继续处理其另一端的未处理骨干。
    // 仅固定非 terminal 的公共端点和局部 Gap 接入点；两条骨干的远端 terminal 状态不再阻止合并。
    bool merged_spine = true;
    while (merged_spine) {
        merged_spine = false;
        std::vector<std::vector<TerminalEndpointOwner>> root_owners(node_points.size());
        for (size_t i = 0; i < primitives.size(); ++i) {
            const size_t first_root = rootOf(static_cast<size_t>(
                std::stoul(primitive_from[i].substr(1))));
            const size_t last_root = rootOf(static_cast<size_t>(
                std::stoul(primitive_to[i].substr(1))));
            root_owners[first_root].push_back({i, true});
            root_owners[last_root].push_back({i, false});
        }

        for (size_t shared = 0; shared < root_owners.size(); ++shared) {
            std::vector<TerminalEndpointOwner> owners;
            for (const TerminalEndpointOwner& owner : root_owners[shared]) {
                if (std::find_if(owners.begin(), owners.end(), [&](const TerminalEndpointOwner& item) {
                        return item.primitive_index == owner.primitive_index;
                    }) == owners.end()) owners.push_back(owner);
            }
            if (owners.size() != 2) continue;

            const TerminalEndpointOwner first_owner = owners[0];
            const TerminalEndpointOwner second_owner = owners[1];
            // 只保护当前公共端点；合并后的 Composite 可以保留另一端的 non-terminal 属性。
            const bool first_shared_terminal = first_owner.is_first
                ? primitives[first_owner.primitive_index].first_terminal
                : primitives[first_owner.primitive_index].last_terminal;
            const bool second_shared_terminal = second_owner.is_first
                ? primitives[second_owner.primitive_index].first_terminal
                : primitives[second_owner.primitive_index].last_terminal;
            if (!first_shared_terminal || !second_shared_terminal) continue;

            bool has_local_gap = false;
            for (const GapConnection& gap : isolated_gaps) {
                if (rootOf(gap.from_node) == shared || rootOf(gap.to_node) == shared) {
                    has_local_gap = true;
                    break;
                }
            }
            if (has_local_gap) continue;

            const Primitive& first = primitives[first_owner.primitive_index];
            const Primitive& second = primitives[second_owner.primitive_index];

            const size_t first_far_root = rootOf(static_cast<size_t>(std::stoul(
                (first_owner.is_first ? primitive_to[first_owner.primitive_index] :
                                        primitive_from[first_owner.primitive_index]).substr(1))));
            const size_t second_far_root = rootOf(static_cast<size_t>(std::stoul(
                (second_owner.is_first ? primitive_to[second_owner.primitive_index] :
                                         primitive_from[second_owner.primitive_index]).substr(1))));
            if (first_far_root == shared || second_far_root == shared ||
                first_far_root == second_far_root) continue;

            const Vec2d first_ray = node_points[first_far_root].xy() - node_points[shared].xy();
            const Vec2d second_ray = node_points[second_far_root].xy() - node_points[shared].xy();
            if (first_ray.norm() < 1e-10 || second_ray.norm() < 1e-10) continue;
            const double cosine = std::abs(first_ray.normalized().dot(second_ray.normalized()));
            const double angle = std::acos(std::max(-1.0, std::min(1.0, cosine)));
            if (angle >= kSpineSmoothAngleThreshold) continue;

            size_t from_root = first_far_root;
            size_t to_root = second_far_root;
            if (!first_owner.is_first && second_owner.is_first) {
                from_root = first_far_root;
                to_root = second_far_root;
            } else if (!second_owner.is_first && first_owner.is_first) {
                from_root = second_far_root;
                to_root = first_far_root;
            } else {
                const Vec2d source_direction = first.last.xy() - first.first.xy();
                if ((node_points[to_root].xy() - node_points[from_root].xy()).dot(source_direction) < 0.0)
                    std::swap(from_root, to_root);
            }
            const bool from_terminal = from_root == first_far_root
                ? (first_owner.is_first ? first.last_terminal : first.first_terminal)
                : (second_owner.is_first ? second.last_terminal : second.first_terminal);
            const bool to_terminal = to_root == first_far_root
                ? (first_owner.is_first ? first.last_terminal : first.first_terminal)
                : (second_owner.is_first ? second.last_terminal : second.first_terminal);
            const std::string merged_id = "smooth_" + first.id + "_" + second.id;
            Primitive merged(merged_id, first.source_id + "+" + second.source_id,
                             NonMotorBackboneEdgeType::Composite,
                             node_points[from_root], node_points[to_root],
                             from_terminal, to_terminal, Vec2d(1, 0),
                             GroupRole::Entry, false);

            const size_t first_index = first_owner.primitive_index;
            const size_t second_index = second_owner.primitive_index;
            std::vector<Primitive> merged_primitives;
            std::vector<std::string> merged_from;
            std::vector<std::string> merged_to;
            merged_primitives.reserve(primitives.size() - 1);
            merged_from.reserve(primitive_from.size() - 1);
            merged_to.reserve(primitive_to.size() - 1);
            for (size_t i = 0; i < primitives.size(); ++i) {
                if (i == first_index) {
                    merged_primitives.push_back(merged);
                    merged_from.push_back("n" + std::to_string(from_root));
                    merged_to.push_back("n" + std::to_string(to_root));
                }
                if (i == first_index || i == second_index) continue;
                merged_primitives.push_back(primitives[i]);
                merged_from.push_back(primitive_from[i]);
                merged_to.push_back(primitive_to[i]);
            }
            primitives.swap(merged_primitives);
            primitive_from.swap(merged_from);
            primitive_to.swap(merged_to);
            merged_spine = true;
            break;
        }
    }
    terminal_nodes.clear();
    for (size_t i = 0; i < primitives.size(); ++i) {
        if (primitives[i].first_terminal) terminal_nodes.push_back(primitive_from[i]);
        if (primitives[i].last_terminal) terminal_nodes.push_back(primitive_to[i]);
    }

    // 压缩面内交点合并的端点，保证骨干线与后续连接线共用同一节点。
    const size_t invalid_index = std::numeric_limits<size_t>::max();
    std::vector<size_t> compact_index(node_points.size(), invalid_index);
    std::vector<Vec3d> compact_points;
    std::vector<NonMotorBackboneNodeType> compact_types;
    std::vector<std::string> compact_sources;
    const auto compactNode = [&](size_t node) {
        const size_t root = rootOf(node);
        if (compact_index[root] == invalid_index) {
            compact_index[root] = compact_points.size();
            compact_points.push_back(node_points[root]);
            compact_types.push_back(node_types[root]);
            compact_sources.push_back(node_sources[root]);
        }
        return compact_index[root];
    };
    for (size_t i = 0; i < primitives.size(); ++i) {
        primitive_from[i] = "n" + std::to_string(compactNode(
            static_cast<size_t>(std::stoul(primitive_from[i].substr(1)))));
        primitive_to[i] = "n" + std::to_string(compactNode(
            static_cast<size_t>(std::stoul(primitive_to[i].substr(1)))));
    }
    for (GapConnection& gap : isolated_gaps) {
        gap.from_node = compactNode(gap.from_node);
        gap.to_node = compactNode(gap.to_node);
    }
    std::vector<std::string> compact_terminals;
    for (const std::string& id : terminal_nodes) {
        compact_terminals.push_back("n" + std::to_string(compactNode(
            static_cast<size_t>(std::stoul(id.substr(1))))));
    }
    std::set<size_t> compact_unresolved_isolated_nodes;
    for (const size_t unresolved_node : unresolved_isolated_nodes)
        compact_unresolved_isolated_nodes.insert(compactNode(unresolved_node));
    node_points.swap(compact_points);
    node_types.swap(compact_types);
    node_sources.swap(compact_sources);
    terminal_nodes.swap(compact_terminals);
    sortTerminalNodes();
    std::vector<NonMotorBackboneEdge> backbone_edges;
    for (size_t i = 0; i < primitives.size(); ++i) {
        NonMotorBackboneEdge edge;
        edge.id = primitives[i].id;
        edge.from = primitive_from[i];
        edge.to = primitive_to[i];
        edge.type = primitives[i].type;
        // Composite沿用参与合并的来源面，后续挂接偏移按原面约束，不能靠长ID反推归属。
        std::istringstream source_stream(primitives[i].source_id);
        std::string source_id;
        while (std::getline(source_stream, source_id, '+')) {
            if (std::find(edge.source_surface_ids.begin(), edge.source_surface_ids.end(), source_id) ==
                edge.source_surface_ids.end()) edge.source_surface_ids.push_back(source_id);
        }
        edge.ring = primitives[i].first_terminal && primitives[i].last_terminal;
        edge.branch = !edge.ring;
        edge.geometry.points = {node_points[static_cast<size_t>(std::stoul(edge.from.substr(1)))],
                                node_points[static_cast<size_t>(std::stoul(edge.to.substr(1)))]};
        backbone_edges.push_back(edge);
    }

    // 近邻端点未合并时，以来源线方向确定首尾，添加大环之外的局部 Gap 连接。
    std::set<std::pair<size_t, size_t>> explicit_gap_pairs;
    std::vector<NonMotorBackboneEdge> local_gap_edges;
    for (const GapConnection& gap : isolated_gaps) {
        if (gap.from_node == gap.to_node) continue;
        const std::pair<size_t, size_t> key(std::min(gap.from_node, gap.to_node),
                                             std::max(gap.from_node, gap.to_node));
        if (!explicit_gap_pairs.insert(key).second) continue;
        NonMotorBackboneEdge edge;
        edge.id = gap.id;
        edge.from = "n" + std::to_string(gap.from_node);
        edge.to = "n" + std::to_string(gap.to_node);
        edge.type = NonMotorBackboneEdgeType::Gap;
        edge.branch = true;
        edge.geometry.points = gap.points.empty()
            ? std::vector<Vec3d>{node_points[gap.from_node], node_points[gap.to_node]} : gap.points;
        edge.geometry.points.front() = node_points[gap.from_node];
        edge.geometry.points.back() = node_points[gap.to_node];
        local_gap_edges.push_back(edge);
    }

    // 按骨干对象中心的逆时针顺序排主环，再在每条骨干段之间连接最近的一对端点。
    std::vector<size_t> ring_spines;
    Vec2d route_center = Vec2d::Zero();
    for (size_t i = 0; i < primitives.size(); ++i) {
        if (!primitives[i].first_terminal || !primitives[i].last_terminal) continue;
        ring_spines.push_back(i);
        route_center += 0.5 * (primitives[i].first.xy() + primitives[i].last.xy());
    }
    if (!ring_spines.empty()) route_center /= static_cast<double>(ring_spines.size());
    std::sort(ring_spines.begin(), ring_spines.end(), [&](size_t a, size_t b) {
        const Vec2d a_center = 0.5 * (primitives[a].first.xy() + primitives[a].last.xy());
        const Vec2d b_center = 0.5 * (primitives[b].first.xy() + primitives[b].last.xy());
        const double a_angle = std::atan2(a_center.y() - route_center.y(), a_center.x() - route_center.x());
        const double b_angle = std::atan2(b_center.y() - route_center.y(), b_center.x() - route_center.x());
        if (a_angle != b_angle) return a_angle < b_angle;
        return primitives[a].id < primitives[b].id;
    });

    // 保留孤立骨干在原有环序中的位置，不压缩相邻对象；否则会跨过孤立骨干重新补出错误 Gap。
    // 若所有骨干都孤立，则回退原有全量大环逻辑，兼容没有可判定近邻的简单输入。
    size_t eligible_ring_spines = 0;
    for (const size_t spine_index : ring_spines) {
        const size_t first = static_cast<size_t>(std::stoul(primitive_from[spine_index].substr(1)));
        const size_t last = static_cast<size_t>(std::stoul(primitive_to[spine_index].substr(1)));
        if (compact_unresolved_isolated_nodes.count(first) == 0 &&
            compact_unresolved_isolated_nodes.count(last) == 0) ++eligible_ring_spines;
    }
    const bool suppress_isolated_ring_gaps = eligible_ring_spines >= 2;

    std::vector<NonMotorBackboneEdge> route_edges;
    std::vector<std::string> ordered_ring_nodes;
    std::set<std::pair<size_t, size_t>> perimeter_gap_pairs;
    std::vector<bool> ring_spine_written(primitives.size(), false);
    // 大环补 Gap 前建立骨干线邻接表。检查某个相邻端点对时屏蔽“下一条骨干线”本身，
    // 只有存在其他骨干线组成的替代路径才跳过 Gap；这样不会误删两条骨干线首尾相接时
    // 用于闭环的必要 Gap，同时可过滤通过安全岛间接连通的重复 Gap。
    std::vector<std::vector<std::pair<size_t, size_t>>> ring_neighbors(node_points.size());
    for (size_t edge_index = 0; edge_index < backbone_edges.size(); ++edge_index) {
        const NonMotorBackboneEdge& edge = backbone_edges[edge_index];
        const size_t from = static_cast<size_t>(std::stoul(edge.from.substr(1)));
        const size_t to = static_cast<size_t>(std::stoul(edge.to.substr(1)));
        ring_neighbors[from].push_back({to, edge_index});
        ring_neighbors[to].push_back({from, edge_index});
    }
    const auto hasAlternateRingPath = [&](size_t start, size_t target,
                                           size_t blocked_edge_index) {
        std::vector<bool> visited(node_points.size(), false);
        std::vector<size_t> pending(1, start);
        visited[start] = true;
        for (size_t cursor = 0; cursor < pending.size(); ++cursor) {
            const size_t node = pending[cursor];
            if (node == target) return true;
            for (const auto& neighbor : ring_neighbors[node]) {
                if (neighbor.second == blocked_edge_index || visited[neighbor.first]) continue;
                visited[neighbor.first] = true;
                pending.push_back(neighbor.first);
            }
        }
        return false;
    };
    bool ring_complete = true;
    for (size_t i = 0; ring_spines.size() >= 2 && i < ring_spines.size(); ++i) {
        const size_t spine_index = ring_spines[i];
        NonMotorBackboneEdge spine_edge = backbone_edges[spine_index];
        const std::string entry_id = primitive_from[spine_index];
        const std::string exit_id = primitive_to[spine_index];
        if (ordered_ring_nodes.empty() || ordered_ring_nodes.back() != entry_id)
            ordered_ring_nodes.push_back(entry_id);
        if (ordered_ring_nodes.back() != exit_id) ordered_ring_nodes.push_back(exit_id);
        route_edges.push_back(spine_edge);
        ring_spine_written[spine_index] = true;

        const size_t next = (i + 1) % ring_spines.size();
        const size_t exit_node = static_cast<size_t>(std::stoul(exit_id.substr(1)));
        const size_t next_spine = ring_spines[next];
        const std::string next_entry_id = primitive_from[next_spine];
        const size_t next_entry = static_cast<size_t>(std::stoul(next_entry_id.substr(1)));
        const std::pair<size_t, size_t> key(std::min(exit_node, next_entry),
                                             std::max(exit_node, next_entry));
        if (exit_node != next_entry) {
            const size_t spine_first = static_cast<size_t>(std::stoul(entry_id.substr(1)));
            const size_t spine_last = static_cast<size_t>(std::stoul(exit_id.substr(1)));
            const bool current_spine_is_unresolved =
                compact_unresolved_isolated_nodes.count(spine_first) != 0 ||
                compact_unresolved_isolated_nodes.count(spine_last) != 0;
            const std::string next_exit_id = primitive_to[next_spine];
            const size_t next_exit = static_cast<size_t>(std::stoul(next_exit_id.substr(1)));
            const bool next_spine_is_unresolved =
                compact_unresolved_isolated_nodes.count(next_entry) != 0 ||
                compact_unresolved_isolated_nodes.count(next_exit) != 0;
            if (suppress_isolated_ring_gaps &&
                (current_spine_is_unresolved || next_spine_is_unresolved)) {
                // 保持涉及孤立骨干的两侧端点断开，避免绕过孤立对象重新生成环 Gap。
                ring_complete = false;
                continue;
            }
            if (hasAlternateRingPath(exit_node, next_entry, next_spine)) {
                // 端点已有不依赖下一条骨干线的替代路径；不再用大环 Gap 重复连接。
                ring_complete = false;
                continue;
            }
            perimeter_gap_pairs.insert(key);
            NonMotorBackboneEdge gap;
            const auto local_gap = std::find_if(local_gap_edges.begin(), local_gap_edges.end(),
                [&](const NonMotorBackboneEdge& candidate) {
                    const size_t from = static_cast<size_t>(std::stoul(candidate.from.substr(1)));
                    const size_t to = static_cast<size_t>(std::stoul(candidate.to.substr(1)));
                    return std::min(from, to) == key.first && std::max(from, to) == key.second;
                });
            gap.id = local_gap == local_gap_edges.end()
                         ? "gap_ring_" + std::to_string(i) : local_gap->id;
            gap.from = local_gap == local_gap_edges.end() ? exit_id : local_gap->from;
            gap.to = local_gap == local_gap_edges.end() ? next_entry_id : local_gap->to;
            gap.type = NonMotorBackboneEdgeType::Gap;
            gap.ring = true;
            gap.geometry.points = {
                node_points[static_cast<size_t>(std::stoul(gap.from.substr(1)))],
                node_points[static_cast<size_t>(std::stoul(gap.to.substr(1)))]};
            if (local_gap != local_gap_edges.end()) gap.geometry = local_gap->geometry;
            route_edges.push_back(gap);
            if (i + 1 < ring_spines.size() && ordered_ring_nodes.back() != next_entry_id)
                ordered_ring_nodes.push_back(next_entry_id);
        }
    }
    if (ordered_ring_nodes.size() > 1 && ordered_ring_nodes.front() == ordered_ring_nodes.back())
        ordered_ring_nodes.pop_back();
    if (!ring_complete) {
        // 跳过重复 Gap 后不再存在一条连续的大环，避免输出不连续的 ring_nodes/ring_edges。
        ordered_ring_nodes.clear();
        for (NonMotorBackboneEdge& edge : route_edges) {
            if (!edge.ring) continue;
            edge.ring = false;
            edge.branch = true;
        }
    }
    for (size_t i = 0; i < backbone_edges.size(); ++i)
        if (!ring_spine_written[i]) {
            NonMotorBackboneEdge edge = backbone_edges[i];
            if (edge.ring) {
                edge.ring = false;
                edge.branch = false;
            }
            route_edges.push_back(edge);
        }

    // 同端点对已由大环连接时复用大环 Gap，避免平行重复边。
    for (const NonMotorBackboneEdge& gap : local_gap_edges) {
        const size_t from = static_cast<size_t>(std::stoul(gap.from.substr(1)));
        const size_t to = static_cast<size_t>(std::stoul(gap.to.substr(1)));
        if (perimeter_gap_pairs.count({std::min(from, to), std::max(from, to)}) == 0)
            route_edges.push_back(gap);
    }

    // 最终按连通分量优先补道路Gap，失败后尝试安全直线Gap，弥补旧规则留下的孤立骨干。
    // 每次仅连接不同分量，不重复加闭环边；受搜索范围和空间约束限制。
    std::vector<size_t> network_parent(node_points.size());
    for (size_t i = 0; i < network_parent.size(); ++i) network_parent[i] = i;
    const auto network_root = [&](size_t node) {
        while (network_parent[node] != node) node = network_parent[node];
        return node;
    };
    for (const auto& edge : route_edges) {
        const size_t a = static_cast<size_t>(std::stoul(edge.from.substr(1)));
        const size_t b = static_cast<size_t>(std::stoul(edge.to.substr(1)));
        network_parent[network_root(a)] = network_root(b);
    }
    // 原始来源面按合并边source_id回溯，选择靠近候选端点的来源面供道路裁剪。
    const auto repair_surface = [&](const Primitive& primitive, const Vec3d& endpoint) -> const Polygon2d* {
        const Polygon2d* best = nullptr;
        double distance = std::numeric_limits<double>::infinity();
        const std::string sources = "+" + primitive.source_id + "+";
        const auto consider = [&](const std::string& id, const Polygon2d& polygon) {
            if (sources.find("+" + id + "+") == std::string::npos) return;
            double d = pointInPolygonSurface(endpoint.xy(), polygon) ? 0 : islandBoundaryDistance(endpoint.xy(), polygon);
            if (d < distance) { distance = d; best = &polygon; }
        };
        for (const auto& crosswalk : input.crosswalks) consider(crosswalk.id, crosswalk.geometry);
        for (const auto& island : input.safety_islands) consider(island.id, island.geometry);
        return best;
    };
    // 每个分量最多增加一条当前最近的连接后重算度，至多执行C-1轮，无全局路径穷举。
    size_t repair_number = 0;
    for (;;) {
        std::vector<size_t> degree(node_points.size(), 0);
        for (const auto& edge : route_edges) {
            ++degree[static_cast<size_t>(std::stoul(edge.from.substr(1)))];
            ++degree[static_cast<size_t>(std::stoul(edge.to.substr(1)))];
        }
        struct RepairCandidate { double distance; size_t from, to; };
        std::vector<RepairCandidate> repairs;
        for (size_t a = 0; a < node_points.size(); ++a) for (size_t b = a + 1; b < node_points.size(); ++b) {
            const double distance = dist(node_points[a], node_points[b]);
            const bool road_supported = road_extension_enabled &&
                distance > gap_options_.search_radius + kPointTolerance &&
                road_edge_near(node_points[a]) <= road_edge_radius + kPointTolerance &&
                road_edge_near(node_points[b]) <= road_edge_radius + kPointTolerance;
            if (gap_options_.search_radius > 0 && degree[a] && degree[b] &&
                network_root(a) != network_root(b) &&
                (distance <= road_edge_radius + kPointTolerance || road_supported))
                repairs.push_back({distance, a, b});
        }
        std::sort(repairs.begin(), repairs.end(), [](const RepairCandidate& a, const RepairCandidate& b) {
            if (a.distance != b.distance) return a.distance < b.distance;
            if (a.from != b.from) return a.from < b.from;
            return a.to < b.to;
        });
        bool joined = false;
        for (const auto& candidate : repairs) {
            // 道路优先：从断头端接入另一分量，可在目标骨干内部拆出公共节点，主环原端点不移动。
            size_t from = candidate.from, to = candidate.to;
            if (degree[from] != 1 && degree[to] == 1) std::swap(from, to);
            std::vector<Vec3d> selected;
            size_t source_index = route_edges.size(), target_index = route_edges.size();
            for (size_t i = 0; degree[from] == 1 && i < route_edges.size() && selected.empty(); ++i) {
                if (route_edges[i].from != "n" + std::to_string(from) && route_edges[i].to != "n" + std::to_string(from)) continue;
                const auto source = std::find_if(primitives.begin(), primitives.end(), [&](const Primitive& p) { return p.id == route_edges[i].id; });
                if (source == primitives.end()) continue;
                for (size_t j = 0; j < route_edges.size() && selected.empty(); ++j) {
                    if (route_edges[j].from != "n" + std::to_string(to) && route_edges[j].to != "n" + std::to_string(to)) continue;
                    const auto target = std::find_if(primitives.begin(), primitives.end(), [&](const Primitive& p) { return p.id == route_edges[j].id; });
                    if (target == primitives.end()) continue;
                    const Polygon2d* first_surface = repair_surface(*source, node_points[from]);
                    const Polygon2d* second_surface = repair_surface(*target, node_points[to]);
                    if (!first_surface || !second_surface) continue;
                    Primitive first = *source, second = *target;
                    first.first = route_edges[i].geometry.points.front(); first.last = route_edges[i].geometry.points.back();
                    second.first = route_edges[j].geometry.points.front(); second.last = route_edges[j].geometry.points.back();
                    const bool first_start = route_edges[i].from == "n" + std::to_string(from);
                    const bool second_start = route_edges[j].from == "n" + std::to_string(to);
                    auto points = road_network.path(node_points[from], node_points[to],
                        roadGapExpandedSurface(*first_surface, gap_options_.surface_expansion),
                        roadGapExpandedSurface(*second_surface, gap_options_.surface_expansion), road_edge_radius);
                    NonMotorGapOptions road_options = gap_options_;
                    if (candidate.distance > gap_options_.search_radius + kPointTolerance)
                        road_options.endpoint_adjustment += 2.0 * gap_options_.inward_offset;
                    if (!shiftedRoadGap(points, gap_center, first, first_start, second, second_start, road_options)) continue;
                    if (degree[to] > 1 && pointToSegment(points.back(), second.first, second.last).first > kPointTolerance) continue;
                    std::vector<Vec3d> whole{first_start ? first.last : first.first};
                    whole.insert(whole.end(), points.begin(), points.end());
                    whole.push_back(second_start ? second.last : second.first);
                    const std::vector<Vec3d> allowed_contacts{
                        first.first, first.last, second.first, second.last};
                    if (!safeRoadGapPath(whole, input, allowed_contacts)) continue;
                    bool crosses = false;
                    for (size_t k = 0; k < route_edges.size(); ++k) {
                        if (k == i || k == j) continue;
                        const auto& edge = route_edges[k];
                        for (size_t x = 1; x < points.size(); ++x) for (size_t y = 1; y < edge.geometry.points.size(); ++y) {
                            std::vector<double> hits;
                            roadGapCuts(points[x - 1], points[x], edge.geometry.points[y - 1], edge.geometry.points[y], hits);
                            for (double t : hits) {
                                const Vec3d hit = points[x - 1] + t * (points[x] - points[x - 1]);
                                if (dist(hit, node_points[to]) > kPointTolerance ||
                                    (edge.from != "n" + std::to_string(to) && edge.to != "n" + std::to_string(to))) crosses = true;
                            }
                        }
                    }
                    if (crosses) continue;
                    selected = points; source_index = i; target_index = j;
                }
            }
            if (!selected.empty()) {
                auto& source = route_edges[source_index];
                node_points[from] = selected.front();
                (source.from == "n" + std::to_string(from) ? source.geometry.points.front() : source.geometry.points.back()) = selected.front();
                if (degree[to] == 1) {
                    node_points[to] = selected.back();
                    auto& target = route_edges[target_index];
                    (target.from == "n" + std::to_string(to) ? target.geometry.points.front() : target.geometry.points.back()) = selected.back();
                } else if (dist(selected.back(), node_points[to]) > kPointTolerance) {
                    // 接入主环/共享骨干内部：保留原边两端，拆边后在公共点连接，避免删掉有效范围。
                    NonMotorBackboneEdge tail = route_edges[target_index];
                    const size_t new_node = node_points.size();
                    node_points.push_back(selected.back()); node_types.push_back(NonMotorBackboneNodeType::Gap); node_sources.push_back(tail.id);
                    network_parent.push_back(network_root(to));
                    route_edges[target_index].to = "n" + std::to_string(new_node);
                    route_edges[target_index].geometry.points.back() = selected.back();
                    tail.from = "n" + std::to_string(new_node); tail.geometry.points.front() = selected.back();
                    tail.id += "_road_split_" + std::to_string(repair_number);
                    route_edges.push_back(tail); to = new_node;
                }
            } else {
                // 超出普通Gap范围的候选只有在道路边缘路径成功时才允许使用，
                // 避免独立放大的道路搜索范围退化成超长直线补连。
                if (candidate.distance > gap_options_.search_radius + kPointTolerance) continue;
                selected = {node_points[from], node_points[to]};
                if (!safeRoadGapPath(selected, input)) continue;
                bool crosses = false;
                for (const auto& edge : route_edges) for (size_t i = 1; i < edge.geometry.points.size(); ++i) {
                    std::vector<double> hits;
                    roadGapCuts(selected.front(), selected.back(), edge.geometry.points[i - 1], edge.geometry.points[i], hits);
                    for (double t : hits) {
                        const std::string node = "n" + std::to_string(t < 0.5 ? from : to);
                        if ((t > 1e-7 && t < 1 - 1e-7) || (node != edge.from && node != edge.to)) crosses = true;
                    }
                }
                if (crosses) continue;
            }
            NonMotorBackboneEdge edge;
            edge.id = (source_index < route_edges.size() ? "gap_road_component_" : "gap_component_") + std::to_string(repair_number++);
            edge.from = "n" + std::to_string(from); edge.to = "n" + std::to_string(to);
            edge.type = NonMotorBackboneEdgeType::Gap; edge.branch = true; edge.geometry.points = selected;
            route_edges.push_back(edge); network_parent[network_root(from)] = network_root(to);
            joined = true; break;
        }
        if (!joined) break;
    }
    std::set<size_t> active_components;
    for (const auto& edge : route_edges)
        active_components.insert(network_root(static_cast<size_t>(std::stoul(edge.from.substr(1)))));
    bb.connected_components = active_components.size();
    if (bb.connected_components > 1) {
        bb.diagnostics.push_back("unconnected backbone components=" + std::to_string(bb.connected_components) +
            "; no safe gap within search radius=" + std::to_string(gap_options_.search_radius) +
            ", road edge radius=" + std::to_string(road_edge_radius));
        std::cerr << "[nonmotor] " << bb.diagnostics.back() << std::endl;
    }

    // 从最终骨干图提取真实闭环：反复剥离度为1的叶节点，剩余2-core中的简单环
    // 才能进入 ring_edges/ring_geometry；挂在环上的安全岛/横道分支不会被误纳入主环。
    const size_t invalid_edge = std::numeric_limits<size_t>::max();
    const auto isPerimeterGap = [&](size_t edge_index) {
        return route_edges[edge_index].id.find("gap_ring_") == 0;
    };
    const auto computeCycleCore = [&](bool include_perimeter_gaps) {
        std::vector<std::vector<std::pair<size_t, size_t>>> neighbors(node_points.size());
        std::vector<size_t> degree(node_points.size(), 0);
        for (size_t edge_index = 0; edge_index < route_edges.size(); ++edge_index) {
            if (!include_perimeter_gaps && isPerimeterGap(edge_index)) continue;
            const size_t from = static_cast<size_t>(std::stoul(route_edges[edge_index].from.substr(1)));
            const size_t to = static_cast<size_t>(std::stoul(route_edges[edge_index].to.substr(1)));
            if (from == to) continue;
            neighbors[from].push_back({to, edge_index});
            neighbors[to].push_back({from, edge_index});
            ++degree[from];
            ++degree[to];
        }
        std::vector<bool> active(route_edges.size(), false);
        for (size_t edge_index = 0; edge_index < route_edges.size(); ++edge_index)
            active[edge_index] = include_perimeter_gaps || !isPerimeterGap(edge_index);
        std::queue<size_t> leaves;
        for (size_t node = 0; node < degree.size(); ++node)
            if (degree[node] <= 1) leaves.push(node);
        while (!leaves.empty()) {
            const size_t node = leaves.front();
            leaves.pop();
            if (degree[node] > 1) continue;
            for (const auto& neighbor : neighbors[node]) {
                if (!active[neighbor.second]) continue;
                active[neighbor.second] = false;
                if (degree[node] > 0) --degree[node];
                if (degree[neighbor.first] > 0 && --degree[neighbor.first] <= 1)
                    leaves.push(neighbor.first);
            }
        }
        return std::make_pair(active, neighbors);
    };
    const auto isSimpleCycleCore = [&](const std::pair<std::vector<bool>,
                                                       std::vector<std::vector<std::pair<size_t, size_t>>>>& core) {
        size_t active_edges = 0;
        for (const bool active : core.first) if (active) ++active_edges;
        if (active_edges < 3) return false;
        for (size_t node = 0; node < core.second.size(); ++node) {
            size_t active_degree = 0;
            for (const auto& neighbor : core.second[node])
                if (core.first[neighbor.second]) ++active_degree;
            if (active_degree != 0 && active_degree != 2) return false;
        }
        return true;
    };
    std::pair<std::vector<bool>, std::vector<std::vector<std::pair<size_t, size_t>>>> cycle_core =
        computeCycleCore(false);
    const bool has_non_gap_cycle = isSimpleCycleCore(cycle_core);
    bool has_semantic_branch = false;
    for (size_t edge_index = 0; edge_index < route_edges.size(); ++edge_index)
        if (route_edges[edge_index].branch && !isPerimeterGap(edge_index) &&
            route_edges[edge_index].id.find("gap_road_") != 0 &&
            route_edges[edge_index].id.find("gap_component_") != 0) {
            has_semantic_branch = true;
            break;
        }
    // 只有没有真实骨干闭环、且不存在明确分支时才允许大环 Gap 参与闭环，兼容简单输入。
    if (!has_non_gap_cycle && !has_semantic_branch)
        cycle_core = computeCycleCore(true);
    bool has_simple_cycle = isSimpleCycleCore(cycle_core);
    std::vector<std::string> cycle_ring_nodes;
    std::vector<size_t> ordered_cycle_edges;
    if (has_simple_cycle) {
        std::vector<size_t> cycle_edges;
        for (size_t edge_index = 0; edge_index < route_edges.size(); ++edge_index)
            if (cycle_core.first[edge_index]) cycle_edges.push_back(edge_index);
        const size_t start_edge = cycle_edges.front();
        const size_t start_node = static_cast<size_t>(std::stoul(
            route_edges[start_edge].from.substr(1)));
        size_t current_node = start_node;
        size_t previous_edge = invalid_edge;
        while (ordered_cycle_edges.size() < cycle_edges.size()) {
            cycle_ring_nodes.push_back("n" + std::to_string(current_node));
            size_t selected_edge = invalid_edge;
            size_t next_node = current_node;
            for (const auto& neighbor : cycle_core.second[current_node]) {
                if (!cycle_core.first[neighbor.second] || neighbor.second == previous_edge) continue;
                selected_edge = neighbor.second;
                next_node = neighbor.first;
                break;
            }
            if (selected_edge == invalid_edge) break;
            ordered_cycle_edges.push_back(selected_edge);
            previous_edge = selected_edge;
            current_node = next_node;
            if (current_node == start_node) break;
        }
        if (current_node != start_node || ordered_cycle_edges.size() != cycle_edges.size()) {
            cycle_ring_nodes.clear();
            ordered_cycle_edges.clear();
            has_simple_cycle = false;
        }
    }
    // 遍历起边的输入方向不稳定；按完整折线面积统一为逆时针，保证左转截取方向。
    double cycle_area_twice = 0.0;
    for (size_t i = 0; i < ordered_cycle_edges.size(); ++i) {
        const auto& edge = route_edges[ordered_cycle_edges[i]];
        const double direction = edge.from == cycle_ring_nodes[i] ? 1.0 : -1.0;
        for (size_t j = 1; j < edge.geometry.points.size(); ++j)
            cycle_area_twice += direction * cross2d(edge.geometry.points[j - 1].xy(), edge.geometry.points[j].xy());
    }
    if (cycle_area_twice < -1e-9) {
        std::reverse(cycle_ring_nodes.begin() + 1, cycle_ring_nodes.end());
        std::reverse(ordered_cycle_edges.begin(), ordered_cycle_edges.end());
    }
    for (NonMotorBackboneEdge& edge : route_edges) {
        edge.ring = false;
        edge.branch = true;
    }
    if (has_simple_cycle) {
        std::vector<NonMotorBackboneEdge> reordered_edges;
        reordered_edges.reserve(route_edges.size());
        for (const size_t edge_index : ordered_cycle_edges) {
            route_edges[edge_index].ring = true;
            route_edges[edge_index].branch = false;
            reordered_edges.push_back(route_edges[edge_index]);
        }
        for (size_t edge_index = 0; edge_index < route_edges.size(); ++edge_index)
            if (!cycle_core.first[edge_index]) reordered_edges.push_back(route_edges[edge_index]);
        route_edges.swap(reordered_edges);
        ordered_ring_nodes = cycle_ring_nodes;
    } else {
        ordered_ring_nodes.clear();
    }

    std::set<std::string> ring_node_set(ordered_ring_nodes.begin(), ordered_ring_nodes.end());
    terminal_nodes.assign(ring_node_set.begin(), ring_node_set.end());

    // 去掉折点产生的孤立节点，重新连续编号骨干图节点。
    std::vector<bool> used_nodes(node_points.size(), false);
    for (const NonMotorBackboneEdge& edge : route_edges) {
        used_nodes[static_cast<size_t>(std::stoul(edge.from.substr(1)))] = true;
        used_nodes[static_cast<size_t>(std::stoul(edge.to.substr(1)))] = true;
    }
    for (const std::string& id : terminal_nodes)
        used_nodes[static_cast<size_t>(std::stoul(id.substr(1)))] = true;
    std::vector<size_t> remap(node_points.size(), invalid_index);
    std::vector<Vec3d> final_points;
    std::vector<NonMotorBackboneNodeType> final_types;
    std::vector<std::string> final_sources;
    for (size_t i = 0; i < node_points.size(); ++i) {
        if (!used_nodes[i]) continue;
        remap[i] = final_points.size();
        final_points.push_back(node_points[i]);
        final_types.push_back(node_types[i]);
        final_sources.push_back(node_sources[i]);
    }
    const auto remapId = [&](const std::string& id) {
        const size_t old_index = static_cast<size_t>(std::stoul(id.substr(1)));
        return std::string("n") + std::to_string(remap[old_index]);
    };
    for (NonMotorBackboneEdge& edge : route_edges) {
        edge.from = remapId(edge.from);
        edge.to = remapId(edge.to);
        edge.geometry.points.front() = final_points[static_cast<size_t>(std::stoul(edge.from.substr(1)))];
        edge.geometry.points.back() = final_points[static_cast<size_t>(std::stoul(edge.to.substr(1)))];
    }
    for (std::string& id : terminal_nodes) id = remapId(id);
    for (std::string& id : ordered_ring_nodes) id = remapId(id);
    node_points.swap(final_points);
    node_types.swap(final_types);
    node_sources.swap(final_sources);

    // 环提取完成后以最终 ring 标记统一派生 branch 标记，保证两类边严格互斥；
    // 没有进入闭合主环的骨干段（包括局部 Gap）全部作为分支保留。
    for (NonMotorBackboneEdge& edge : route_edges)
        edge.branch = !edge.ring;

    // 主环边保持“骨干段、Gap、下一骨干段”的路径顺序，分枝按角度排在主环之后。
    const auto first_branch = std::find_if(route_edges.begin(), route_edges.end(),
        [](const NonMotorBackboneEdge& edge) { return !edge.ring; });
    route_center = Vec2d::Zero();
    for (const Vec3d& point : node_points) route_center += point.xy();
    if (!node_points.empty()) route_center /= static_cast<double>(node_points.size());
    std::sort(first_branch, route_edges.end(), [&](const NonMotorBackboneEdge& a,
                                                   const NonMotorBackboneEdge& b) {
        const Vec2d a_midpoint = 0.5 * (a.geometry.points.front().xy() + a.geometry.points.back().xy());
        const Vec2d b_midpoint = 0.5 * (b.geometry.points.front().xy() + b.geometry.points.back().xy());
        const double a_angle = std::atan2(a_midpoint.y() - route_center.y(), a_midpoint.x() - route_center.x());
        const double b_angle = std::atan2(b_midpoint.y() - route_center.y(), b_midpoint.x() - route_center.x());
        if (a_angle != b_angle) return a_angle < b_angle;
        return a.id < b.id;
    });
    for (const NonMotorBackboneEdge& edge : route_edges) {
        bb.edges.push_back(edge);
        if (edge.ring) bb.ring_edges.push_back(edge.id);
        if (!edge.ring) bb.branch_edges.push_back(edge.id);
    }
    bb.ring_nodes = ordered_ring_nodes;
    // 主环沿拓扑遍历方向拼接完整边几何，不丢失RoadEdge Gap的折点。
    for (size_t i = 0; i < bb.ring_edges.size(); ++i) {
        const auto& edge = bb.edges[i];
        auto points = edge.geometry.points;
        if (edge.from != bb.ring_nodes[i]) std::reverse(points.begin(), points.end());
        for (const auto& p : points)
            if (bb.ring_geometry.points.empty() || dist(bb.ring_geometry.points.back(), p) > kPointTolerance)
                bb.ring_geometry.points.push_back(p);
    }
    // 几何填充完成后再确认闭环，避免有效的 2-core 环被误标为未闭合。
    bb.closed = bb.ring_nodes.size() >= 3 &&
                bb.ring_edges.size() == bb.ring_nodes.size() &&
                bb.ring_geometry.points.size() >= bb.ring_nodes.size() + 1 &&
                dist(bb.ring_geometry.points.front(), bb.ring_geometry.points.back()) <= kPointTolerance;

    for (size_t i = 0; i < node_points.size(); ++i) {
        NonMotorBackboneNode node;
        node.id = "n" + std::to_string(i);
        node.type = node_types[i];
        node.point = node_points[i];
        node.source_id = node_sources[i];
        bb.nodes.push_back(node);
    }
    return bb;
}

bool NonMotorBackboneExtractor::output(const NonMotorBackbone& backbone, std::string out_dir) const {
    return save(backbone, out_dir, "");
}

}  // namespace isg
