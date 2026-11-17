// 仅由nonmotor_backbone_extractor.cpp在isg命名空间内包含。
// RoadEdge按真实交点分段建图；查询只做局部多源Dijkstra，不跨道路断缝。
namespace {

/// 将线段与另一线段的真实相交位置加入参数数组；共线重叠只插入端点。
void roadGapCuts(const Vec3d& a, const Vec3d& b, const Vec3d& c, const Vec3d& d,
                 std::vector<double>& cuts) {
    const Vec2d r = b.xy() - a.xy(), s = d.xy() - c.xy();
    if (r.squaredNorm() < 1e-16) return;
    const double denominator = cross2d(r, s);
    if (std::abs(denominator) > 1e-10) {
        const double t = cross2d(c.xy() - a.xy(), s) / denominator;
        const double u = cross2d(c.xy() - a.xy(), r) / denominator;
        if (t >= -1e-8 && t <= 1 + 1e-8 && u >= -1e-8 && u <= 1 + 1e-8)
            cuts.push_back(std::max(0.0, std::min(1.0, t)));
    } else if (std::abs(cross2d(c.xy() - a.xy(), r)) < 1e-8 * r.norm()) {
        const double t0 = (c.xy() - a.xy()).dot(r) / r.squaredNorm();
        const double t1 = (d.xy() - a.xy()).dot(r) / r.squaredNorm();
        const double lo = std::max(0.0, std::min(t0, t1)), hi = std::min(1.0, std::max(t0, t1));
        if (lo <= hi + 1e-8) { cuts.push_back(lo); cuts.push_back(hi); }
    }
}

/// 排序并合并相同线段参数，避免零长图边。
void uniqueRoadGapCuts(std::vector<double>& cuts) {
    std::sort(cuts.begin(), cuts.end());
    cuts.erase(std::unique(cuts.begin(), cuts.end(), [](double a, double b) {
        return std::abs(a - b) < 1e-9;
    }), cuts.end());
}

/// 使用原生Clipper按米制外扩面；外环/洞先统一方向，避免包装器的delta单位差异。
std::vector<Polygon2d> roadGapExpandedSurface(const Polygon2d& surface, double expansion) {
    constexpr double scale = 10000.0;
    ClipperLib::ClipperOffset offset(2.0, 0.001 * scale);
    const auto add_ring = [&](const std::vector<Vec3d>& ring, bool outer) {
        if (ring.size() < 3) return;
        ClipperLib::Path path;
        for (const Vec3d& p : ring) path.emplace_back(std::llround(p.x() * scale), std::llround(p.y() * scale));
        if (ClipperLib::Orientation(path) != outer) std::reverse(path.begin(), path.end());
        offset.AddPath(path, ClipperLib::jtRound, ClipperLib::etClosedPolygon);
    };
    add_ring(surface.outer, true);
    for (const auto& hole : surface.holes) add_ring(hole, false);
    ClipperLib::Paths paths;
    offset.Execute(paths, std::max(0.0, expansion) * scale);
    std::vector<Polygon2d> result;
    for (const auto& path : paths) {
        Polygon2d polygon;
        for (const auto& p : path) polygon.outer.emplace_back(p.X / scale, p.Y / scale, 0);
        result.push_back(polygon);
    }
    return result;
}

/// 外扩结果使用奇偶包含规则，保留洞；边界本身算相交。
bool roadGapInExpanded(const Vec3d& p, const std::vector<Polygon2d>& polygons) {
    bool inside = false;
    for (const auto& polygon : polygons) {
        for (size_t i = 0; i < polygon.outer.size(); ++i)
            if (pointToSegment(p, polygon.outer[i], polygon.outer[(i + 1) % polygon.outer.size()]).first < 1e-6) return true;
        if (pointInPolygon(p.xy(), polygon.outer)) inside = !inside;
    }
    return inside;
}

/// 新增连接的空间门禁：完整线段不越fence/洞、不碰Boundary/障碍，且折线不自交。
/// fence通过边界切分后逐区间检查，避免稀疏采样漏掉窄洞。
/// @param allowed_contacts 允许与已有骨干公共端点相切的坐标；这些接触是道路Gap的合法
/// 公共节点，不应被当作新的Boundary穿越。其它边界交点仍全部拒绝。
bool safeRoadGapPath(const std::vector<Vec3d>& points, const IntersectionInput& input,
                     const std::vector<Vec3d>& allowed_contacts = {}) {
    if (points.size() < 2) return false;
    for (size_t i = 1; i < points.size(); ++i) {
        const Vec3d& a = points[i - 1];
        const Vec3d& b = points[i];
        if (!std::isfinite(a.x()) || !std::isfinite(a.y()) || !std::isfinite(a.z()) ||
            !std::isfinite(b.x()) || !std::isfinite(b.y()) || !std::isfinite(b.z()) ||
            dist(a, b) < 1e-7) return false;
        for (const auto& boundary : input.boundaries)
            for (size_t j = 1; j < boundary.geometry.points.size(); ++j) {
                std::vector<double> hits;
                roadGapCuts(a, b, boundary.geometry.points[j - 1], boundary.geometry.points[j], hits);
                for (const double t : hits) {
                    const Vec3d hit = a + t * (b - a);
                    bool allowed = false;
                    for (const auto& contact : allowed_contacts)
                        if (dist(hit, contact) <= 1e-5) { allowed = true; break; }
                    if (!allowed) return false;
                }
            }
        for (const auto& obstacle : input.obstacles) {
            const Polygon2d& polygon = obstacle.buffered_geometry.outer.empty() ? obstacle.geometry : obstacle.buffered_geometry;
            if (pointInPolygonSurface(a.xy(), polygon) || pointInPolygonSurface(b.xy(), polygon) ||
                !segmentPolygonBoundaryEdges(a, b, polygon).empty()) return false;
        }
        if (!input.area.geometry.outer.empty()) {
            std::vector<double> cuts{0, 1};
            const auto add_ring = [&](const std::vector<Vec3d>& ring) {
                for (size_t j = 0; j < ring.size(); ++j) roadGapCuts(a, b, ring[j], ring[(j + 1) % ring.size()], cuts);
            };
            add_ring(input.area.geometry.outer);
            for (const auto& hole : input.area.geometry.holes) add_ring(hole);
            uniqueRoadGapCuts(cuts);
            for (size_t j = 0; j < cuts.size(); ++j) {
                const double t = j == 0 ? cuts[j] : 0.5 * (cuts[j - 1] + cuts[j]);
                if (!pointInPolygonSurface((a + t * (b - a)).xy(), input.area.geometry)) return false;
            }
            if (!pointInPolygonSurface(b.xy(), input.area.geometry)) return false;
        }
        for (size_t j = 1; j + 1 < i; ++j) {
            std::vector<double> hits;
            roadGapCuts(a, b, points[j - 1], points[j], hits);
            if (!hits.empty()) return false;
        }
        if (i >= 2) {
            const Vec2d previous = a.xy() - points[i - 2].xy(), next = b.xy() - a.xy();
            if (std::abs(cross2d(previous, next)) < 1e-8 && previous.dot(next) < 0) return false;
        }
    }
    return true;
}

/// 与其他骨干只能共享整条候选路径的端点；正长度重叠和内部交叉均须拒绝。
bool roadGapCrossesOther(const std::vector<Vec3d>& path, const std::vector<Vec3d>& other) {
    for (size_t i = 1; i < path.size(); ++i) for (size_t j = 1; j < other.size(); ++j) {
        std::vector<double> hits;
        roadGapCuts(path[i - 1], path[i], other[j - 1], other[j], hits);
        for (double t : hits) {
            const Vec3d p = path[i - 1] + t * (path[i] - path[i - 1]);
            if ((dist(p, path.front()) > 1e-6 && dist(p, path.back()) > 1e-6) ||
                (dist(p, other.front()) > 1e-6 && dist(p, other.back()) > 1e-6)) return true;
        }
    }
    return false;
}

/// 连续道路边缘网络。真实交点拆边，跨层交点不连通；所有段保留线性插值高程。
struct RoadGapNetwork {
    struct Segment { Vec3d a, b; };
    std::vector<Segment> segments;

    /// 每次路口只建立一次RoadEdge拓扑，包围盒预筛后才求精确交点。
    explicit RoadGapNetwork(const IntersectionInput& input) {
        std::vector<Segment> raw;
        for (const auto& boundary : input.boundaries) {
            if (boundary.type != Boundary::Type::RoadEdge) continue;
            for (size_t j = 1; j < boundary.geometry.points.size(); ++j)
                if (dist(boundary.geometry.points[j - 1], boundary.geometry.points[j]) > 1e-6)
                    raw.push_back({boundary.geometry.points[j - 1], boundary.geometry.points[j]});
        }
        std::vector<std::vector<double>> cuts(raw.size(), std::vector<double>{0, 1});
        for (size_t i = 0; i < raw.size(); ++i) for (size_t j = i + 1; j < raw.size(); ++j) {
            if ((raw[i].a.xy().cwiseMin(raw[i].b.xy()).array() >
                 raw[j].a.xy().cwiseMax(raw[j].b.xy()).array() + 1e-6).any() ||
                (raw[j].a.xy().cwiseMin(raw[j].b.xy()).array() >
                 raw[i].a.xy().cwiseMax(raw[i].b.xy()).array() + 1e-6).any()) continue;
            std::vector<double> first, second;
            roadGapCuts(raw[i].a, raw[i].b, raw[j].a, raw[j].b, first);
            roadGapCuts(raw[j].a, raw[j].b, raw[i].a, raw[i].b, second);
            for (double t : first) for (double u : second) {
                const Vec3d p = raw[i].a + t * (raw[i].b - raw[i].a);
                const Vec3d q = raw[j].a + u * (raw[j].b - raw[j].a);
                if ((p.xy() - q.xy()).norm() < 1e-6 && std::abs(p.z() - q.z()) < 0.1) {
                    cuts[i].push_back(t); cuts[j].push_back(u);
                }
            }
        }
        for (size_t i = 0; i < raw.size(); ++i) {
            uniqueRoadGapCuts(cuts[i]);
            for (size_t j = 1; j < cuts[i].size(); ++j)
                segments.push_back({raw[i].a + cuts[i][j - 1] * (raw[i].b - raw[i].a),
                                    raw[i].a + cuts[i][j] * (raw[i].b - raw[i].a)});
        }
    }

    /// 对两个外扩面裁剪RoadEdge，局部多源最短路按道路长度和到骨干端点距离选取有效区间。
    std::vector<Vec3d> path(const Vec3d& from, const Vec3d& to,
                           const std::vector<Polygon2d>& first, const std::vector<Polygon2d>& second,
                           double radius) const {
        std::vector<Vec3d> nodes;
        std::vector<std::vector<std::pair<size_t, double>>> neighbors;
        // 半径约束端点与道路网络的局部走廊，而不是要求每一条中间边同时贴近
        // 两个端点；否则长但连续的RoadEdge会在中段被剪断，无法形成有效Gap。
        const Vec2d low = from.xy().cwiseMin(to.xy()) - Vec2d(radius, radius);
        const Vec2d high = from.xy().cwiseMax(to.xy()) + Vec2d(radius, radius);
        const auto node_for = [&](const Vec3d& p) {
            for (size_t i = 0; i < nodes.size(); ++i)
                if ((nodes[i].xy() - p.xy()).norm() < 1e-6 && std::abs(nodes[i].z() - p.z()) < 0.1) return i;
            nodes.push_back(p); neighbors.emplace_back(); return nodes.size() - 1;
        };
        for (const Segment& segment : segments) {
            if ((segment.a.xy().cwiseMin(segment.b.xy()).array() > high.array()).any() ||
                (segment.a.xy().cwiseMax(segment.b.xy()).array() < low.array()).any()) continue;
            std::vector<double> cuts{0, 1};
            for (const auto* polygons : {&first, &second})
                for (const auto& polygon : *polygons)
                    for (size_t i = 0; i < polygon.outer.size(); ++i)
                        roadGapCuts(segment.a, segment.b, polygon.outer[i], polygon.outer[(i + 1) % polygon.outer.size()], cuts);
            uniqueRoadGapCuts(cuts);
            for (size_t i = 1; i < cuts.size(); ++i) {
                const size_t a = node_for(segment.a + cuts[i - 1] * (segment.b - segment.a));
                const size_t b = node_for(segment.a + cuts[i] * (segment.b - segment.a));
                if (a == b) continue;
                const double length = (nodes[a].xy() - nodes[b].xy()).norm();
                neighbors[a].push_back({b, length}); neighbors[b].push_back({a, length});
            }
        }
        const size_t invalid = nodes.size();
        std::vector<double> distances(nodes.size(), std::numeric_limits<double>::infinity());
        std::vector<size_t> previous(nodes.size(), invalid);
        std::priority_queue<std::pair<double, size_t>, std::vector<std::pair<double, size_t>>,
                            std::greater<std::pair<double, size_t>>> queue;
        for (size_t i = 0; i < nodes.size(); ++i) if (roadGapInExpanded(nodes[i], first)) {
            distances[i] = (nodes[i].xy() - from.xy()).norm(); queue.push({distances[i], i});
        }
        size_t best = invalid; double best_cost = std::numeric_limits<double>::infinity();
        while (!queue.empty()) {
            const auto current = queue.top(); queue.pop();
            if (current.first > distances[current.second] + 1e-8 || current.first > 4 * radius) continue;
            if (roadGapInExpanded(nodes[current.second], second)) {
                const double cost = current.first + (nodes[current.second].xy() - to.xy()).norm();
                if (cost < best_cost && previous[current.second] != invalid) { best = current.second; best_cost = cost; }
            }
            for (const auto& neighbor : neighbors[current.second]) {
                const double cost = current.first + neighbor.second;
                if (cost + 1e-8 >= distances[neighbor.first]) continue;
                distances[neighbor.first] = cost; previous[neighbor.first] = current.second;
                queue.push({cost, neighbor.first});
            }
        }
        std::vector<Vec3d> result;
        if (best == invalid) return result;
        for (size_t node = best; node != invalid; node = previous[node]) result.push_back(nodes[node]);
        std::reverse(result.begin(), result.end());
        // 去除RoadEdge拆分/裁剪引入的共线内部点，真实转折及Z变化均保留。
        for (size_t i = 1; i + 1 < result.size();) {
            const Vec3d direction = result[i + 1] - result[i - 1];
            const double t = (result[i].xy() - result[i - 1].xy()).dot(direction.xy()) / direction.xy().squaredNorm();
            if (pointToSegment(result[i], result[i - 1], result[i + 1]).first < 1e-7 &&
                std::abs(result[i].z() - (result[i - 1].z() + t * direction.z())) < 1e-6) result.erase(result.begin() + i);
            else ++i;
        }
        return result;
    }
};

/// 沿一致的内侧法向偏移道路折线；用支撑线交点保持平行形态，过长尖角不接受。
/// 首尾与原骨干无限线求交但禁止吞掉/反转骨干，交点Z沿原骨干插值。
bool shiftedRoadGap(std::vector<Vec3d>& points, const Vec2d& center,
                    const Primitive& first, bool first_start, const Primitive& second, bool second_start,
                    const NonMotorGapOptions& options) {
    if (points.size() < 2 || options.inward_offset <= 0) return false;
    std::vector<Vec2d> directions, normals;
    double side_score = 0;
    for (size_t i = 1; i < points.size(); ++i) {
        const Vec2d direction = (points[i].xy() - points[i - 1].xy()).normalized();
        directions.push_back(direction); normals.emplace_back(-direction.y(), direction.x());
        side_score += normals.back().dot(center - 0.5 * (points[i].xy() + points[i - 1].xy()));
    }
    if (std::abs(side_score) < 1e-8) return false;
    for (size_t i = 0; i < normals.size(); ++i) {
        normals[i] *= side_score > 0 ? options.inward_offset : -options.inward_offset;
        // 整条连续道路统一选侧；局部若转到路口中心的反侧则拒绝，不输出局部向外偏移。
        if (normals[i].dot(center - 0.5 * (points[i].xy() + points[i + 1].xy())) < -1e-8) return false;
    }
    std::vector<Vec3d> shifted = points;
    shifted.front() = Vec3d(points.front().xy() + normals.front(), points.front().z());
    shifted.back() = Vec3d(points.back().xy() + normals.back(), points.back().z());
    for (size_t i = 1; i + 1 < points.size(); ++i) {
        const double denominator = cross2d(directions[i - 1], directions[i]);
        if (std::abs(denominator) < 1e-8) {
            if (directions[i - 1].dot(directions[i]) < 0) return false;
            shifted[i] = Vec3d(points[i].xy() + normals[i], points[i].z());
        } else {
            const double t = cross2d(normals[i] - normals[i - 1], directions[i]) / denominator;
            const Vec2d offset = normals[i - 1] + t * directions[i - 1];
            if (offset.norm() > 4 * options.inward_offset) return false;
            shifted[i] = Vec3d(points[i].xy() + offset, points[i].z());
        }
    }
    const auto common_point = [&](const Primitive& spine, bool start, size_t index, const Vec2d& direction) {
        const Vec3d endpoint = start ? spine.first : spine.last;
        const Vec3d other = start ? spine.last : spine.first;
        const Vec2d along = endpoint.xy() - other.xy();
        const double denominator = cross2d(along, direction);
        if (std::abs(denominator) < 1e-9) return false;
        const double t = cross2d(shifted[index].xy() - other.xy(), direction) / denominator;
        const Vec3d p(other.xy() + t * along, other.z() + t * (endpoint.z() - other.z()));
        if (t <= 1e-6 || dist(p, endpoint) > options.endpoint_adjustment) return false;
        shifted[index] = p; return true;
    };
    if (!common_point(first, first_start, 0, directions.front()) ||
        !common_point(second, second_start, shifted.size() - 1, directions.back())) return false;
    for (size_t i = 1; i < shifted.size(); ++i)
        if ((shifted[i].xy() - shifted[i - 1].xy()).dot(directions[i - 1]) <= 1e-6) return false;
    points.swap(shifted); return true;
}

}  // namespace
