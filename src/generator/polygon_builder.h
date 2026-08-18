#pragma once

#include "types.h"
#include "utils.h"

#include <algorithm>
#include <string>
#include <vector>

namespace isg {

struct EdgeEndpoint {
    Vec3d pt;
    int roadEdgeIdx = -1;
    bool isStart = true;
    EdgeEndpoint(Vec3d point, int edge_index, bool start)
        : pt(point), roadEdgeIdx(edge_index), isStart(start) {}
};

enum class AreaSourceKind { Lane, LaneEdge, Boundary };

struct BoundaryLine {
    std::vector<Vec3d> pts;
    std::vector<std::string> source_ids;
    std::vector<Vec2d> segment_dirs;
    std::vector<std::pair<Vec3d, Vec3d>> source_endpoints;
    AreaSourceKind source_kind = AreaSourceKind::Boundary;
};

/// 根据车道组切线和道路边界构建精细路口面。
class IntersectionAreaBuilder {
public:
    explicit IntersectionAreaBuilder(
        double group_cut_offset,
        double cut_line_extension = 7.0,
        double snap_tolerance = 0.5,
        std::string winding = "clockwise");

    IntersectionArea build(
        const IntersectionInput& input,
        const std::vector<ConnectivityCurve>& centerlines,
        const std::vector<ConnectivityLaneEdge>& edge_lines);

    /// 几何测试和旧调用方使用的兼容工具。
    template <typename P>
    static std::vector<std::vector<P>> MergeBoundaries(
        std::vector<std::vector<P>> boundaries) {
        struct Chain { std::vector<P> points; };
        const auto same = [](const P& a, const P& b) { return dist(a, b) <= 1e-9; };
        const auto append = [](std::vector<P>& dst, const std::vector<P>& src,
                               int begin, int end, int step) {
            for (int i = begin; i != end; i += step)
                if (dst.empty() || dist(dst.back(), src[i]) > 1e-9) dst.push_back(src[i]);
        };
        const auto merge = [&](const std::vector<P>& a, const std::vector<P>& b,
                               std::vector<P>& out) {
            if (a.empty() || b.empty()) return false;
            out.clear();
            if (same(a.back(), b.front())) { out = a; append(out, b, 1, (int)b.size(), 1); return true; }
            if (same(a.back(), b.back())) { out = a; append(out, b, (int)b.size() - 2, -1, -1); return true; }
            if (same(a.front(), b.back())) { out = b; append(out, a, 1, (int)a.size(), 1); return true; }
            if (same(a.front(), b.front())) {
                append(out, b, (int)b.size() - 1, -1, -1);
                append(out, a, 1, (int)a.size(), 1);
                return true;
            }
            return false;
        };
        std::vector<Chain> chains;
        for (auto& boundary : boundaries)
            if (boundary.size() >= 2) { Chain chain; chain.points = std::move(boundary); chains.push_back(std::move(chain)); }
        bool changed = true;
        while (changed) {
            changed = false;
            for (int i = 0; i < (int)chains.size() && !changed; ++i) {
                for (int j = i + 1; j < (int)chains.size(); ++j) {
                    std::vector<P> joined;
                    if (!merge(chains[i].points, chains[j].points, joined)) continue;
                    chains[i].points = std::move(joined);
                    chains.erase(chains.begin() + j);
                    changed = true;
                    break;
                }
            }
        }
        std::vector<std::vector<P>> result;
        for (auto& chain : chains) result.push_back(std::move(chain.points));
        return result;
    }

    static std::vector<BoundaryLine> MergeBoundaryLines(
        std::vector<BoundaryLine> boundaries);

private:
    double group_cut_offset_;
    double cut_line_extension_;
    double snap_tolerance_;
    std::string winding_;
};

}  // 命名空间 isg
