#pragma once

#include "types.h"
#include "utils.h"
#include <algorithm>
#include <cmath>
#include <limits>

namespace isg {
namespace nonmotor_detail {

/// 面边界净距包含外环和洞；外环边界属于覆盖范围，洞内部不属于覆盖范围。
inline double surfaceBoundaryDistance(const Vec2d& point, const Polygon2d& polygon) {
    double distance = std::numeric_limits<double>::infinity();
    const auto visit = [&](const std::vector<Vec3d>& ring) {
        for (size_t i = 0; i < ring.size(); ++i)
            distance = std::min(distance, pointToSegment(point, ring[i], ring[(i+1)%ring.size()]).first);
    };
    visit(polygon.outer);
    for (const auto& hole : polygon.holes) visit(hole);
    return distance;
}

inline bool insideSourceSurface(const Vec2d& point, const Polygon2d& polygon) {
    if (!pointInPolygon(point, polygon.outer) && surfaceBoundaryDistance(point, polygon) > 1e-7) return false;
    for (const auto& hole : polygon.holes) if (pointInPolygon(point, hole)) return false;
    return true;
}

/// 无限骨干支撑线与来源面边界精确求交，返回连续覆盖的站距区间；洞和凹口会切开区间。
/// 平行边的顶点也加入切分，避免线与边界共线时丢失交点。
inline std::vector<std::pair<double,double>> sourceSurfaceIntervals(
    const Vec2d& origin, const Vec2d& tangent, const Polygon2d& polygon) {
    std::vector<double> cuts;
    const auto visit = [&](const std::vector<Vec3d>& ring) {
        for (size_t i = 0; i < ring.size(); ++i) {
            const Vec2d a = ring[i].xy(), edge = ring[(i+1)%ring.size()].xy()-a;
            const double det = cross2d(tangent, edge);
            if (std::abs(det) < 1e-10) {
                if (std::abs(cross2d(tangent,a-origin)) < 1e-8) {
                    cuts.push_back((a-origin).dot(tangent));
                    cuts.push_back((a+edge-origin).dot(tangent));
                }
            } else {
                const double u = cross2d(a-origin,tangent)/det;
                if (u >= -1e-9 && u <= 1+1e-9) cuts.push_back(cross2d(a-origin,edge)/det);
            }
        }
    };
    visit(polygon.outer);
    for (const auto& hole : polygon.holes) visit(hole);
    std::sort(cuts.begin(),cuts.end());
    cuts.erase(std::unique(cuts.begin(),cuts.end(),[](double a,double b){return std::abs(a-b)<1e-8;}),cuts.end());
    std::vector<std::pair<double,double>> result;
    for (size_t i = 1; i < cuts.size(); ++i) {
        if (!insideSourceSurface(origin+tangent*((cuts[i-1]+cuts[i])*0.5),polygon)) continue;
        if (!result.empty() && std::abs(result.back().second-cuts[i-1]) < 1e-8) result.back().second = cuts[i];
        else result.push_back({cuts[i-1],cuts[i]});
    }
    return result;
}

/// 原骨干在来源面中的有效连续范围，端部允许沿支撑线重新求交；中部保留覆盖和净距。
struct SourceSurfaceSpan {
    const Polygon2d* polygon = nullptr;
    std::string source_id;
    Vec2d origin, tangent;
    double first = 0, last = 0, clearance = 0;
};

/// 精确检查中部线段覆盖和净距，同时禁止平移后在凹面/洞之间离面再重入。
inline bool sourceOffsetFits(const SourceSurfaceSpan& span, double offset) {
    const Vec2d origin = span.origin+rotLeft(span.tangent)*offset;
    const auto intervals = sourceSurfaceIntervals(origin,span.tangent,*span.polygon);
    if (intervals.size() != 1 || intervals[0].first > span.first+1e-8 ||
        intervals[0].second < span.last-1e-8) return false;
    const Vec3d a(origin+span.tangent*span.first), b(origin+span.tangent*span.last);
    double nearest = std::numeric_limits<double>::infinity();
    const auto visit = [&](const std::vector<Vec3d>& ring) {
        for (size_t i = 0; i < ring.size(); ++i) {
            const auto& c = ring[i]; const auto& d = ring[(i+1)%ring.size()];
            nearest = std::min(nearest,std::min({pointToSegment(a.xy(),c,d).first,
                pointToSegment(b.xy(),c,d).first,pointToSegment(c.xy(),a,b).first,pointToSegment(d.xy(),a,b).first}));
        }
    };
    visit(span.polygon->outer);
    for (const auto& hole : span.polygon->holes) visit(hole);
    return nearest+1e-7 >= span.clearance;
}

} // namespace nonmotor_detail
} // namespace isg
