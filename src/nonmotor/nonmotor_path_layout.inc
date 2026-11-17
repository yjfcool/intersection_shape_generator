// 私有实现：实际有向几何段与占用区间一一对应，不通过点下标猜测来源边。
struct LayoutUse { size_t path, segment, corridor; double first, last; };
struct LayoutPath {
    const Connectivity* connection;
    Vec3d entry, exit;
    Vec2d entry_dir, exit_dir;
    std::vector<Vec3d> raw;
    std::vector<size_t> uses;
    LineString2d body, geometry;
    double radius = 0.0;
    bool fixed = false, connected = true;
    size_t rank = 0, count = 1;
};
struct LayoutCorridor {
    Vec2d origin, tangent;
    std::vector<size_t> uses;
    std::vector<double> stations;
    std::vector<nonmotor_detail::SourceSurfaceSpan> surfaces;
    double min_offset = -std::numeric_limits<double>::infinity();
    double max_offset = std::numeric_limits<double>::infinity();
    bool allow_surface_shift = false;
};

/// 按明确来源面及共线几何匹配通道。复合骨干每个来源面分别裁剪；Gap不借用相邻面。
/// 横向搜索只在零偏移所在的连续可行区间内进行，结合面顶点事件/32等分及16次二分；
/// 极窄洞的横向事件也必须检查，不能靠均匀扫描越过一段不可行区间。
void constrainSourceSurfaceCorridors(std::vector<LayoutCorridor>& corridors,
    const NonMotorBackbone& backbone, const IntersectionInput& input, double clearance) {
    for (auto& corridor : corridors) for (const auto& edge : backbone.edges) {
        for (size_t j = 1; j < edge.geometry.points.size(); ++j) {
            const Vec2d a = edge.geometry.points[j-1].xy(), b = edge.geometry.points[j].xy();
            if ((b-a).norm() < 1e-7 || std::abs(cross2d(corridor.tangent,(b-a).normalized())) > 1e-7 ||
                std::abs(cross2d(corridor.tangent,a-corridor.origin)) > 1e-6) continue;
            for (const auto& id : edge.source_surface_ids) {
                const Polygon2d* polygon = nullptr;
                for (const auto& cw : input.crosswalks) if (cw.id == id) polygon = &cw.geometry;
                for (const auto& si : input.safety_islands) if (si.id == id) polygon = &si.geometry;
                if (!polygon || polygon->outer.size() < 3) continue;
                const auto intervals = nonmotor_detail::sourceSurfaceIntervals(corridor.origin,corridor.tangent,*polygon);
                for (const auto& interval : intervals) {
                    const double lo = std::max(interval.first,std::min((a-corridor.origin).dot(corridor.tangent),(b-corridor.origin).dot(corridor.tangent)));
                    const double hi = std::min(interval.second,std::max((a-corridor.origin).dot(corridor.tangent),(b-corridor.origin).dot(corridor.tangent)));
                    if (hi-lo < 0.1) continue;
                    nonmotor_detail::SourceSurfaceSpan span;
                    span.polygon = polygon; span.source_id = id;
                    span.origin = corridor.origin; span.tangent = corridor.tangent;
                    // 端部重新求交，仅从中部约束净距；短面按20%留过渡区，长面每端最多2米。
                    const double end_trim = std::min(2.0,(hi-lo)*0.2);
                    span.first = lo+end_trim; span.last = hi-end_trim;
                    span.clearance = std::max(0.0,clearance);
                    // 原骨干已靠边时不强求不可满足的净距，但仍禁止偏移进一步越界。
                    if (!nonmotor_detail::sourceOffsetFits(span,0)) {
                        double low = 0, high = span.clearance;
                        for (int k = 0; k < 16; ++k) {
                            span.clearance = (low+high)*0.5;
                            if (nonmotor_detail::sourceOffsetFits(span,0)) low = span.clearance;
                            else high = span.clearance;
                        }
                        span.clearance = low;
                    }
                    corridor.surfaces.push_back(span);
                    double extent = 0;
                    for (const auto& point : polygon->outer)
                        extent = std::max(extent,std::abs(cross2d(corridor.tangent,point.xy()-corridor.origin)));
                    for (int side : {-1,1}) {
                        double safe = 0, blocked = extent+1;
                        std::vector<double> events{0,extent+1};
                        const auto add_events = [&](const std::vector<Vec3d>& ring) {
                            for (const auto& point : ring) {
                                const double position = side*cross2d(corridor.tangent,point.xy()-corridor.origin);
                                for (double event : {position-span.clearance,position,position+span.clearance})
                                    if (event > 0 && event < extent+1) events.push_back(event);
                            }
                        };
                        add_events(polygon->outer);
                        for (const auto& hole : polygon->holes) add_events(hole);
                        for (int step = 1; step < 32; ++step) events.push_back((extent+1)*step/32);
                        std::sort(events.begin(),events.end());
                        std::vector<double> probes;
                        for (size_t k = 1; k < events.size(); ++k) {
                            probes.push_back((events[k-1]+events[k])*0.5); probes.push_back(events[k]);
                        }
                        for (double trial : probes) {
                            if (!nonmotor_detail::sourceOffsetFits(span,side*trial)) { blocked = trial; break; }
                            safe = trial;
                        }
                        for (int k = 0; k < 16; ++k) {
                            const double trial = (safe+blocked)*0.5;
                            if (nonmotor_detail::sourceOffsetFits(span,side*trial)) safe = trial;
                            else blocked = trial;
                        }
                        if (side < 0) corridor.min_offset = std::max(corridor.min_offset,-safe);
                        else corridor.max_offset = std::min(corridor.max_offset,safe);
                    }
                }
            }
        }
    }
}

// 对称共线投影检查；只将正长度交集视为重叠，点接触另行检查。
double overlapLength(const Vec2d& a, const Vec2d& b, const Vec2d& c, const Vec2d& d) {
    const Vec2d ab = b - a, cd = d - c;
    const double la = ab.norm(), lb = cd.norm();
    if (la < 1e-8 || lb < 1e-8) return 0;
    const Vec2d t = ab / la;
    if (std::abs(cross2d(t, cd / lb)) > 1e-8 ||
        std::abs(cross2d(t, c - a)) > 1e-8 || std::abs(cross2d(t, d - a)) > 1e-8) return 0;
    const double x = (c - a).dot(t), y = (d - a).dot(t);
    return std::max(0.0, std::min(la, std::max(x, y)) - std::max(0.0, std::min(x, y)));
}

size_t layoutRoot(std::vector<size_t>& roots, size_t i) {
    while (roots[i] != i) { roots[i] = roots[roots[i]]; i = roots[i]; }
    return i;
}

// 不同来源的共线几何归入同一物理通道；站距切分记录完整/部分重叠。
// 通道左法向统一朝向几何中心，仅作为轨位排序参考，不替代实际空间约束验收。
void buildOccupancy(std::vector<LayoutPath>& paths, std::vector<LayoutUse>& uses,
                    std::vector<LayoutCorridor>& corridors, std::vector<size_t>& roots, const Vec2d& center) {
    for (size_t i = 0; i < paths.size(); ++i) {
        for (size_t j = 1; j < paths[i].raw.size(); ++j) {
            const Vec2d a = paths[i].raw[j - 1].xy(), b = paths[i].raw[j].xy();
            if ((b - a).norm() < 1e-7) continue;
            size_t k = 0;
            for (; k < corridors.size(); ++k) {
                const auto& c = corridors[k];
                if (std::abs(cross2d(c.tangent, (b - a).normalized())) < 1e-8 &&
                    std::abs(cross2d(c.tangent, a - c.origin)) < 1e-7) break;
            }
            if (k == corridors.size()) {
                LayoutCorridor c; c.origin = a; c.tangent = (b - a).normalized();
                const double side = cross2d(c.tangent, center-a);
                if (side < -1e-9 || (std::abs(side) <= 1e-9 &&
                    (c.tangent.x() < -1e-9 || (std::abs(c.tangent.x()) <= 1e-9 && c.tangent.y() < 0))))
                    c.tangent = -c.tangent;
                corridors.push_back(c);
            }
            paths[i].uses.push_back(uses.size());
            corridors[k].uses.push_back(uses.size());
            uses.push_back({i, j - 1, k, 0, 0});
        }
    }
    for (auto& u : uses) {
        auto& c = corridors[u.corridor]; const auto& p = paths[u.path];
        u.first = (p.raw[u.segment].xy() - c.origin).dot(c.tangent);
        u.last = (p.raw[u.segment + 1].xy() - c.origin).dot(c.tangent);
        c.stations.push_back(u.first); c.stations.push_back(u.last);
    }
    for (auto& c : corridors) {
        std::sort(c.stations.begin(), c.stations.end());
        c.stations.erase(std::unique(c.stations.begin(), c.stations.end(),
            [](double a, double b) { return std::abs(a - b) < 1e-7; }), c.stations.end());
        for (size_t i = 0; i < c.uses.size(); ++i) {
            for (size_t j = i + 1; j < c.uses.size(); ++j) {
                const auto& a = uses[c.uses[i]]; const auto& b = uses[c.uses[j]];
                if (std::min(std::max(a.first, a.last), std::max(b.first, b.last)) -
                    std::max(std::min(a.first, a.last), std::min(b.first, b.last)) > 1e-7)
                    roots[layoutRoot(roots, b.path)] = layoutRoot(roots, a.path);
            }
        }
    }
    std::map<size_t, std::vector<size_t>> components;
    for (size_t i = 0; i < paths.size(); ++i) components[layoutRoot(roots, i)].push_back(i);
    for (const auto& component : components) {
        for (size_t j = 0; j < component.second.size(); ++j) {
            paths[component.second[j]].rank = j;
            paths[component.second[j]].count = component.second.size();
        }
    }
}

// 掉头按进出车道端点中位所在侧选轨，独立于骨干/路径行驶方向；端点恰在线上时不猜测侧别。
int uturnTrackSide(const LayoutPath& path, const LayoutCorridor& corridor) {
    if (path.fixed || !isUTurn(path.connection->turn_type)) return 0;
    const double side = ((path.entry.xy()+path.exit.xy())*0.5-corridor.origin).dot(rotLeft(corridor.tangent));
    return side < -1e-7 ? -1 : side > 1e-7 ? 1 : 0;
}

// 只对当前站位实际占用者排序和去重；掉头分别放在其车道组侧的外缘，其他路径保稳定几何序。
std::vector<size_t> stationPaths(const LayoutCorridor& c, double station,
    const std::vector<LayoutPath>& paths, const std::vector<LayoutUse>& uses) {
    std::vector<size_t> active;
    for (size_t index : c.uses) {
        const auto& u = uses[index];
        if (station < std::min(u.first, u.last) - 1e-7 || station > std::max(u.first, u.last) + 1e-7) continue;
        active.push_back(u.path);
    }
    std::sort(active.begin(), active.end(), [&](size_t a, size_t b) {
        return std::make_tuple(uturnTrackSide(paths[a], c), paths[a].rank, a) <
               std::make_tuple(uturnTrackSide(paths[b], c), paths[b].rank, b);
    });
    active.erase(std::unique(active.begin(), active.end()), active.end());
    return active;
}

// 局部占用者居中分配等间距轨位；单独掉头也允许向车道组侧偏移，固定形态继续占零轨位。
double unconstrainedStationTrackOffset(size_t index, const std::vector<size_t>& active, const LayoutCorridor& c,
    const std::vector<LayoutPath>& paths, double spacing, double shift) {
    if (spacing <= 0 || paths[index].fixed || active.empty()) return 0;
    const auto found = std::find(active.begin(), active.end(), index);
    if (found == active.end()) return 0;
    const size_t position = static_cast<size_t>(found-active.begin());
    // 单独普通路径默认零偏移；来源面内首轮受阻后允许在已验证的面宽内避让。
    if (active.size() == 1 && uturnTrackSide(paths[index], c) == 0)
        return c.allow_surface_shift && !c.surfaces.empty() ? shift : 0;
    size_t lo = kInvalidIndex, hi = 0, negative = 0, positive = 0;
    for (size_t i = 0; i < active.size(); ++i) {
        const auto& path = paths[active[i]];
        if (path.fixed) {
            lo = std::min(lo, i); hi = std::max(hi, i);
        }
        negative += uturnTrackSide(path, c) < 0 ? 1 : 0;
        positive += uturnTrackSide(path, c) > 0 ? 1 : 0;
    }
    if (lo != kInvalidIndex) {
        if (position < lo) return (static_cast<double>(position)-lo)*spacing;
        if (position > hi) return (static_cast<double>(position)-hi)*spacing;
        return 0;
    }
    const double count = static_cast<double>(active.size());
    // 多数占用者都是同侧掉头时仅作满足其侧别的最小平移；普通共享区间保持左右居中。
    const double lower = static_cast<double>(positive)-count*0.5;
    const double upper = count*0.5-static_cast<double>(negative);
    const double center_shift = std::max(lower, std::min(upper, 0.0));
    return (static_cast<double>(position)-0.5*(count-1)+center_shift)*spacing+shift;
}

/// 同一共享区间所有轨位执行同一个仿射变换。面宽不足时等比例压缩间距，
/// 禁止逐条clamp到边界形成重叠；固定曲线继续占零轨位，不参与平移。
double stationTrackOffset(size_t index, const std::vector<size_t>& active, const LayoutCorridor& c,
    const std::vector<LayoutPath>& paths, double spacing, double shift) {
    // 整组为无归属Gap避障所作的平移不能把面内骨干一起推向面边界；面内轨位独立靠近原骨干。
    const double local_shift = c.surfaces.empty() || c.allow_surface_shift ? shift : 0;
    const double offset = unconstrainedStationTrackOffset(index,active,c,paths,spacing,local_shift);
    if (spacing <= 0 || paths[index].fixed || active.empty() || c.surfaces.empty()) return offset;
    double lower = std::numeric_limits<double>::infinity(), upper = -lower;
    bool fixed = false;
    for (size_t other : active) {
        const double value = unconstrainedStationTrackOffset(other,active,c,paths,spacing,local_shift);
        lower = std::min(lower,value); upper = std::max(upper,value);
        fixed = fixed || paths[other].fixed;
    }
    const double min_offset = c.min_offset, max_offset = c.max_offset;
    if (fixed) {
        double scale = 1;
        if (lower < -1e-9) scale = std::min(scale,min_offset/lower);
        if (upper > 1e-9) scale = std::min(scale,max_offset/upper);
        return offset*std::max(0.0,scale);
    }
    const double width = std::max(0.0,max_offset-min_offset);
    const double scale = upper-lower > 1e-9 ? std::min(1.0,width/(upper-lower)) : 1;
    const double half = (upper-lower)*scale*0.5;
    const double center = std::max(min_offset+half,std::min(max_offset-half,(lower+upper)*0.5));
    return center+(offset-(lower+upper)*0.5)*scale;
}

// 事件站位和中间采样使用同一局部轨位规则，进出骨干交点不再被单独锁定。
double stationOffset(const LayoutPath& path, const LayoutCorridor& c, double station,
    const std::vector<LayoutPath>& paths, const std::vector<LayoutUse>& uses, double spacing, double shift) {
    const auto active = stationPaths(c, station, paths, uses);
    const auto found = std::find_if(active.begin(), active.end(), [&](size_t i) { return &paths[i] == &path; });
    return found == active.end() ? 0 : stationTrackOffset(*found, active, c, paths, spacing, shift);
}

// 同等约束下按占用长度比较掉头侧别、两侧均衡和偏移幅度，防止首个可用的单侧平移候选直接胜出。
std::array<double, 3> distributionPreference(const std::vector<size_t>& component,
    const std::vector<LayoutPath>& paths, const std::vector<LayoutUse>& uses,
    const std::vector<LayoutCorridor>& corridors, double spacing, double shift) {
    std::array<double, 3> score{{0, 0, 0}};
    if (spacing <= 0) return score;
    for (const auto& c : corridors) for (size_t j = 1; j < c.stations.size(); ++j) {
        const double length = c.stations[j]-c.stations[j-1];
        const auto active = stationPaths(c, (c.stations[j]+c.stations[j-1])*0.5, paths, uses);
        if (active.empty() || std::find(component.begin(), component.end(), active.front()) == component.end()) continue;
        double total = 0;
        int negative = 0, positive = 0;
        for (size_t index : active) {
            const double offset = stationTrackOffset(index, active, c, paths, spacing, shift)/spacing;
            const int side = uturnTrackSide(paths[index], c);
            if (side) score[0] += length*std::max(0.0, 0.5-side*offset);
            total += offset;
            negative += offset < -1e-7 ? 1 : 0;
            positive += offset > 1e-7 ? 1 : 0;
            score[2] += length*std::abs(offset);
        }
        score[1] += length*(std::abs(total/active.size()) +
            std::max(0, std::abs(positive-negative)-static_cast<int>(active.size()%2)));
    }
    return score;
}

// 对占用事件区间加密；所有骨干点（含进出交点）均使用对应轨位，不在首尾恢复原交点。
std::vector<Vec3d> displacedKnots(const LayoutPath& path, const std::vector<LayoutPath>& paths,
    const std::vector<LayoutUse>& uses, const std::vector<LayoutCorridor>& corridors,
    double spacing, double shift, double corner_radius) {
    std::vector<Vec3d> knots;
    for (size_t index : path.uses) {
        const auto& u = uses[index]; const auto& c = corridors[u.corridor];
        std::vector<double> stations{u.first, u.last};
        for (double s : c.stations) {
            if (s <= std::min(u.first, u.last) + 1e-7 || s >= std::max(u.first, u.last) - 1e-7) continue;
            // 圆角保护带内的占用事件不能再制造短折线；该事件的轨位变化
            // 延后到折点外侧连续过渡，避免在真实骨干转角前后形成二次折角。
            const double guard = std::max(0.2, corner_radius);
            if (std::min(std::abs(s-u.first), std::abs(s-u.last)) < guard) continue;
            stations.push_back(s);
        }
        std::sort(stations.begin(), stations.end());
        if (u.last < u.first) std::reverse(stations.begin(), stations.end());
        std::vector<Vec3d> segment;
        for (size_t j = 1; j < stations.size(); ++j) {
            const double a = stations[j - 1], b = stations[j];
            const int pieces = std::max(4, static_cast<int>(std::ceil(std::abs(b - a) / 2.0)));
            const double da = stationOffset(path, c, a, paths, uses, spacing, shift);
            const double db = stationOffset(path, c, b, paths, uses, spacing, shift);
            const double mid = stationOffset(path, c, 0.5 * (a + b), paths, uses, spacing, shift);
            for (int k = j == 1 ? 0 : 1; k <= pieces; ++k) {
                const double t = static_cast<double>(k) / pieces;
                const double v = t <= 0.5 ? t * 2 : (t - 0.5) * 2, w = v * v * (3 - 2 * v);
                const double offset = t <= 0.5 ? da + (mid - da) * w : mid + (db - mid) * w;
                const double f = (a + (b - a) * t - u.first) / (u.last - u.first);
                appendUnique(segment, offsetXY(lerp3(path.raw[u.segment], path.raw[u.segment + 1], f), rotLeft(c.tangent) * offset));
            }
        }
        if (!knots.empty() && !segment.empty()) {
            // 偏移支撑线求公共点，不能用两端均值让相邻直段的轨位在转角前后渐变回原骨干。
            knots.back() = knots.size() > 1 && segment.size() > 1
                ? nonmotor_detail::joinOffsetSegments(knots[knots.size()-2], knots.back(),
                    segment.front(), segment[1], path.raw[u.segment], true)
                : (knots.back()+segment.front())*0.5;
            size_t first_append = 1;
            const double guard = std::max(0.2, corner_radius);
            while (first_append + 1 < segment.size() &&
                (segment[first_append].xy()-segment.front().xy()).norm() < guard) ++first_append;
            // 保留该段末点，避免短骨干被保护带全部吞掉；长段从保护带外继续采样。
            for (size_t j = first_append; j < segment.size(); ++j) appendUnique(knots, segment[j]);
        } else knots = segment;
        appendUnique(knots, segment.back());
    }
    return knots;
}

/// 清理偏移后产生的局部回折点。真实非掉头路径不允许相邻有效段反向；
/// 删除短保护带内的反向伪折点，或删除接近180度的明显回折，保留道路骨干的正常大转角交给贝塞尔圆角处理。
std::vector<Vec3d> sanitizeDisplacedKnots(const std::vector<Vec3d>& source, double corner_radius,
                                          bool uturn) {
    if (uturn || source.size() < 3) return source;
    std::vector<Vec3d> points = source;
    const double short_limit = std::max(0.35, std::min(0.8, 0.4 * std::max(corner_radius, 1.0)));
    bool changed = true;
    while (changed && points.size() >= 3) {
        changed = false;
        for (size_t i = 1; i + 1 < points.size(); ++i) {
            const Vec2d previous = points[i].xy() - points[i - 1].xy();
            const Vec2d next = points[i + 1].xy() - points[i].xy();
            const double first_length = previous.norm(), last_length = next.norm();
            if (first_length < 1e-9 || last_length < 1e-9) {
                points.erase(points.begin() + i); changed = true; break;
            }
            const double dot = previous.dot(next) / (first_length * last_length);
            // 短段反向是轨位事件/支撑线回退产生的伪折点；删除后由相邻真实段一次圆角连接。
            if ((std::min(first_length, last_length) <= short_limit && dot < -0.15) || dot < -0.82) {
                points.erase(points.begin() + i); changed = true; break;
            }
        }
    }
    return points;
}

// 最终输出先统一处理端部与骨干折点，再以受限小半径修正残留尖折；有效骨干审计
// 排除车道射线上的直线进出段，保留全部有效骨干转角。
struct LayoutGeometry { LineString2d body, full; };
LayoutGeometry makeLayoutGeometry(const LayoutPath& path, const std::vector<Vec3d>& knots, double radius) {
    LayoutGeometry result;
    if (path.fixed) { result.body.points = result.full.points = path.raw; return result; }
    result.full = nonmotor_detail::roundAttachedPolyline(knots, path.entry, path.entry_dir,
        path.exit, path.exit_dir, radius);
    if (!isUTurn(path.connection->turn_type) && radius > 0.0 && result.full.points.size() > 3) {
        // 一次端部过渡后，偏移轨位交点仍可能留下短的内部miter；再以较小截距
        // 只圆化内部真实折点，避免把车道端点或掉头最短路径重新拉偏。
        // 端部窗口合并后，内部偏移事件仍可能留下采样尺度的短折角；
        // 采用不超过目标截距的有界稳定半径，避免在长骨干中保留针尖式二次转弯。
        const double stabilizer_radius = std::min(1.5, 0.75 * radius);
        result.full = nonmotor_detail::roundPolylineCorners(result.full.points, stabilizer_radius);
        result.full = nonmotor_detail::smoothSharpCorners(result.full, stabilizer_radius);
    }
    if (result.full.points.size() < 2) { result.body = result.full; return result; }
    size_t first = 0, last = result.full.points.size()-1;
    const auto along_lane = [](const Vec3d& point, const Vec3d& origin, const Vec2d& direction) {
        return direction.norm() > 1e-9 && std::abs(cross2d(point.xy()-origin.xy(), direction.normalized())) < 1e-7;
    };
    // 只剔除到原始截取位置附近的共线进出段，不把同方向的真实骨干误删。
    const double entry_limit = path.raw.empty() ? 0 : (path.raw.front().xy()-path.entry.xy()).norm()+4*radius;
    const double exit_limit = path.raw.empty() ? 0 : (path.raw.back().xy()-path.exit.xy()).norm()+4*radius;
    while (first+1 < last && along_lane(result.full.points[first+1], path.entry, path.entry_dir) &&
        (result.full.points[first+1].xy()-path.entry.xy()).norm() <= entry_limit+1e-7) ++first;
    while (last > first+1 && along_lane(result.full.points[last-1], path.exit, path.exit_dir) &&
        (result.full.points[last-1].xy()-path.exit.xy()).norm() <= exit_limit+1e-7) --last;
    result.body.points.assign(result.full.points.begin()+first, result.full.points.begin()+last+1);
    return result;
}

// X扫描和包围盒剪枝；仅对空间接近的输出线段执行精确验收。
struct AuditSegment { Vec3d a, b; size_t path, index; double minx, maxx, miny, maxy; };
template<class Visitor>
void visitLayoutPairs(const std::vector<LineString2d>& lines, Visitor visitor) {
    std::vector<AuditSegment> segments;
    for (size_t p = 0; p < lines.size(); ++p) {
        for (size_t i = 1; i < lines[p].points.size(); ++i) {
            const auto& a = lines[p].points[i - 1]; const auto& b = lines[p].points[i];
            if ((a.xy() - b.xy()).norm() < 1e-8) continue;
            segments.push_back({a, b, p, i - 1, std::min(a.x(), b.x()), std::max(a.x(), b.x()), std::min(a.y(), b.y()), std::max(a.y(), b.y())});
        }
    }
    std::sort(segments.begin(), segments.end(), [](const AuditSegment& a, const AuditSegment& b) { return a.minx < b.minx; });
    for (size_t i = 0; i < segments.size(); ++i) {
        for (size_t j = i + 1; j < segments.size() && segments[j].minx <= segments[i].maxx + 1e-8; ++j) {
            const auto& a = segments[i]; const auto& b = segments[j];
            if (b.miny > a.maxy + 1e-8 || b.maxy < a.miny - 1e-8) continue;
            visitor(a, b);
        }
    }
}

// 进出段现可弯曲，单条曲线的自交必须检查完整输出；不同曲线间仍不检查同簇非端点相交。
bool layoutSelfIntersects(const LineString2d& line) {
    bool hit = false;
    visitLayoutPairs({line}, [&](const AuditSegment& a, const AuditSegment& b) {
        if (std::max(a.index, b.index)-std::min(a.index, b.index) > 1 &&
            segmentsIntersect(a.a.xy(), a.b.xy(), b.a.xy(), b.b.xy(), nullptr)) hit = true;
    });
    return hit;
}

/// 最终圆角输出仍须覆盖来源面：按最近原始路径段归属，裁去面端部重新求交区，
/// 对中部实际线段做精确面裁剪，洞和凹口不能被采样间距漏检；不同车道连接段不套用此约束。
bool layoutLeavesSourceSurface(const LineString2d& line, const LayoutPath& path,
    const std::vector<LayoutUse>& uses, const std::vector<LayoutCorridor>& corridors,
    std::string* source_id = nullptr) {
    if (path.fixed) return false;
    for (size_t index : path.uses) {
        const auto& use = uses[index]; const auto& corridor = corridors[use.corridor];
        for (const auto& span : corridor.surfaces) for (size_t i = 1; i < line.points.size(); ++i) {
            // 挂接端和相邻骨干折点的圆角仍可从面外进入，验收仅覆盖实际截取段中部。
            const double trim = std::min(2.0,std::abs(use.last-use.first)*0.2);
            const double middle_first = std::max(span.first,std::min(use.first,use.last)+trim);
            const double middle_last = std::min(span.last,std::max(use.first,use.last)-trim);
            if (middle_last <= middle_first+1e-8) continue;
            const Vec2d a = line.points[i-1].xy(), b = line.points[i].xy(), d = b-a;
            const double first = (a-span.origin).dot(span.tangent), delta = d.dot(span.tangent);
            if (std::abs(delta) < 1e-9) continue;
            double lo = (middle_first-first)/delta, hi = (middle_last-first)/delta;
            if (lo > hi) std::swap(lo,hi);
            lo = std::max(0.0,lo); hi = std::min(1.0,hi);
            if (hi-lo < 1e-8) continue;
            const Vec2d midpoint = a+d*((lo+hi)*0.5);
            size_t nearest = kInvalidIndex; double distance = std::numeric_limits<double>::infinity();
            for (size_t j = 1; j < path.raw.size(); ++j) {
                const double value = pointToSegment(midpoint,path.raw[j-1],path.raw[j]).first;
                if (value < distance) { nearest = j-1; distance = value; }
            }
            if (nearest != use.segment) continue;
            // 车道端点到骨干的连接不属于面内骨干。必须与进出段共同比较最近来源，
            // 否则远处入口直线的法向投影也会被错误归入第一条横道骨干。
            if (pointToSegment(midpoint,path.entry,path.raw.front()).first+1e-7 < distance ||
                pointToSegment(midpoint,path.raw.back(),path.exit).first+1e-7 < distance) continue;
            const Vec2d start = a+d*lo, end = a+d*hi;
            const double length = (end-start).norm();
            if (length < 1e-8) continue;
            const auto intervals = nonmotor_detail::sourceSurfaceIntervals(start,(end-start)/length,*span.polygon);
            bool covered = false;
            for (const auto& interval : intervals)
                if (interval.first <= 1e-7 && interval.second >= length-1e-7) covered = true;
            if (!covered) {
                if (source_id) *source_id = span.source_id + " near (" +
                    std::to_string(start.x()) + "," + std::to_string(start.y()) + ")";
                return true;
            }
        }
    }
    return false;
}

// 分开计数穿越类约束和净距约束，避免用减少一项净距告警换取新增穿越障碍或边界。
std::pair<int, int> layoutSpatialScore(const LineString2d& full, const LayoutPath& path,
    const IntersectionInput& input, const std::vector<ConnectivityCurve>& motors, const NonmotorConifg& config) {
    int crossing = (pathHitsObstacle(full, input.obstacles) ? 1 : 0) +
        (pathOutsideFence(full, input.area.geometry) ? 1 : 0) +
        (pathHitsBoundary(full, input.boundaries) ? 1 : 0);
    int clearance = minimumRoadEdgeDistance(full, input.boundaries) + 1e-9 < config.road_edge_clearance ? 1 : 0;
    const Lane* en = input.findLane(path.connection->entry_lane_id);
    const Lane* ex = input.findLane(path.connection->exit_lane_id);
    const double motor_clearance = en && ex && en->motoway && !ex->motoway ? config.motor_to_nm_clearance : 0;
    for (const auto& motor : motors) {
        if (minimumPolylineDistance(full, motor.geometry) + 1e-9 < motor_clearance) ++clearance;
        bool overlap = false;
        for (size_t a = 1; a < full.points.size() && !overlap; ++a)
            for (size_t b = 1; b < motor.geometry.points.size() && !overlap; ++b)
                overlap = overlapLength(full.points[a-1].xy(), full.points[a].xy(),
                                        motor.geometry.points[b-1].xy(), motor.geometry.points[b].xy()) > 1e-6;
        if (overlap) ++crossing;
    }
    return {crossing, clearance};
}

#include "nonmotor/nonmotor_boundary_avoidance.inc"
#include "nonmotor/nonmotor_attachment_windows.inc"

std::vector<ConnectivityCurve> layoutNonMotorPaths(const IntersectionInput& input,
    const std::vector<ConnectivityCurve>& motors, NonMotorBackbone& backbone, const NonmotorConifg& config) {
    recoverClosedRingFromEdgeGraph(backbone);
    std::vector<LayoutPath> paths;
    for (const auto& connection : input.connectivities) {
        LayoutPath p; p.connection = &connection;
        const auto entry = input.entryPtDir(connection.entry_lane_id), exit = input.exitPtDir(connection.exit_lane_id);
        p.entry = Vec3d(entry.first, laneEndpointZ(input, connection.entry_lane_id, true));
        p.exit = Vec3d(exit.first, laneEndpointZ(input, connection.exit_lane_id, false));
        p.entry_dir = entry.second; p.exit_dir = exit.second;
        p.fixed = connection.fixed_shape && connection.geometry.points.size() >= 2;
        if (p.fixed) p.raw = connection.geometry.points;
        else {
            auto from = attachToBackbone(backbone, entry.first, entry.second, true);
            auto to = attachToBackbone(backbone, exit.first, -exit.second, true);
            if (std::isfinite(from.distance) && std::isfinite(to.distance)) {
                if (!isUTurn(connection.turn_type)) {
                    collapseShortSameSegmentAttachments(backbone, from, to,
                        config.short_same_segment_max_length);
                }
                p.raw = backbonePath(backbone, from, to, isUTurn(connection.turn_type)).points;
                if (p.raw.empty()) { p.raw = {from.point, to.point}; p.connected = false; }
            } else p.connected = false;
        }
        paths.push_back(p);
    }
    // 内部稳定几何排序，最后恢复输入顺序；输入数组下标不再决定偏移方向。
    std::sort(paths.begin(), paths.end(), [](const LayoutPath& a, const LayoutPath& b) {
        return std::make_tuple(a.entry.x(), a.entry.y(), a.exit.x(), a.exit.y(), a.connection->id) <
               std::make_tuple(b.entry.x(), b.entry.y(), b.exit.x(), b.exit.y(), b.connection->id);
    });
    std::vector<size_t> roots(paths.size());
    for (size_t i = 0; i < roots.size(); ++i) roots[i] = i;
    std::vector<LayoutUse> uses; std::vector<LayoutCorridor> corridors;
    // 主环（无主环时用路口面）顶点均值只用于统一左右参考，实际可行性仍按原几何精确验收。
    const auto& reference = backbone.ring_geometry.points.empty() ? input.area.geometry.outer : backbone.ring_geometry.points;
    size_t reference_count = reference.size();
    if (reference_count > 1 && (reference.front().xy()-reference.back().xy()).norm() < 1e-7) --reference_count;
    Vec2d center(0, 0);
    for (size_t i = 0; i < reference_count; ++i) center += reference[i].xy();
    if (reference_count) center /= static_cast<double>(reference_count);
    buildOccupancy(paths, uses, corridors, roots, center);
    constrainSourceSurfaceCorridors(corridors,backbone,input,config.source_surface_clearance);
    std::map<size_t, std::vector<size_t>> components;
    for (size_t i = 0; i < paths.size(); ++i) components[layoutRoot(roots, i)].push_back(i);
    for (const auto& component : components) {
        // 非机动不验证同簇非端点相交；约束优先，随后比较掉头侧别、左右均衡与偏移幅度。
        using Score = std::tuple<int, int, int, int, int, double, double, double>;
        Score best_score(1000000, 1000000, 1000000, 1000000, 1000000, 1e100, 1e100, 1e100);
        std::vector<LineString2d> best, best_full;
        std::vector<double> best_radii;
        // 首轮面内轨位靠近原骨干；若仍有空间/重叠/自交冲突，第二轮才允许在安全面宽内避让。
        // 两轮各14个保序候选，候选数有界，不为单纯的分布偏好把骨干推到面边界。
        for (int source_trial = 0; source_trial < 2; ++source_trial) {
        if (source_trial) {
            if (std::get<0>(best_score) == 0 && std::get<1>(best_score) == 0 &&
                std::get<2>(best_score) == 0 && std::get<3>(best_score) == 0 && std::get<4>(best_score) == 0) break;
            for (auto& corridor : corridors) corridor.allow_surface_shift = true;
        }
        for (int mirror = 0; mirror < 2; ++mirror) {
        if (mirror) for (size_t index : component.second) paths[index].rank = paths[index].count - 1 - paths[index].rank;
        // 两种保序方向各最多七次整体平移，候选数固定，不做排列穷举。
        for (double shift : {0.0, -1.0, 1.0, -2.0, 2.0, -4.0, 4.0}) {
            std::vector<LineString2d> candidate, candidate_full;
            std::vector<double> radii;
            int surface_crossing = 0, crossing = 0, clearance = 0, overlaps = 0, self = 0;
            for (size_t index : component.second) {
                auto& p = paths[index];
                p.radius = config.corner_radius;
                const auto raw_knots = p.fixed || p.uses.empty() ? p.raw :
                    displacedKnots(p, paths, uses, corridors, config.overlap_offset,
                        shift * config.overlap_offset, p.radius);
                const auto has_backtracking = [&]() {
                    for (size_t k = 1; k + 1 < raw_knots.size(); ++k) {
                        const Vec2d a = raw_knots[k].xy() - raw_knots[k-1].xy();
                        const Vec2d b = raw_knots[k+1].xy() - raw_knots[k].xy();
                        if (a.norm() < 1e-9 || b.norm() < 1e-9) return true;
                        const double cosine = a.dot(b) / (a.norm() * b.norm());
                        if (cosine < -0.15 && std::min(a.norm(), b.norm()) <= 0.8) return true;
                        if (cosine < -0.82) return true;
                    }
                    return false;
                }();
                const auto knots = has_backtracking && !isUTurn(p.connection->turn_type)
                    ? sanitizeDisplacedKnots(raw_knots, p.radius, false) : raw_knots;
                auto layout = makeLayoutGeometry(p, knots, p.radius);
                auto full = layout.full;
                auto spatial = layoutSpatialScore(full, p, input, motors, config);
                // 同一轨位只在圆角新增空间违规时缩小截距，最多6次，不回退轨位或移动车道端点。
                if (!p.fixed && !p.uses.empty() && config.corner_radius > 0 &&
                    (spatial != std::make_pair(0, 0) || layoutSelfIntersects(full))) {
                    const auto unsmoothed = layoutSpatialScore(makeLayoutGeometry(p, knots, 0).full, p, input, motors, config);
                    const bool current_self = layoutSelfIntersects(full);
                    if (unsmoothed < spatial || current_self) {
                        for (int step = 1; step <= 6; ++step) {
                            const double radius = std::ldexp(config.corner_radius, -step);
                            const auto trial = makeLayoutGeometry(p, knots, radius);
                            const auto score = layoutSpatialScore(trial.full, p, input, motors, config);
                            if ((score < spatial || current_self) && !layoutSelfIntersects(trial.full)) {
                                spatial = score; full = trial.full; layout = trial; p.radius = radius;
                            }
                            if (spatial <= unsmoothed) break;
                        }
                    }
                }
                self += layoutSelfIntersects(full) ? 1 : 0;
                surface_crossing += layoutLeavesSourceSurface(full,p,uses,corridors) ? 1 : 0;
                crossing += spatial.first; clearance += spatial.second;
                candidate.push_back(layout.body);
                candidate_full.push_back(full);
                radii.push_back(p.radius);
                
            }
            visitLayoutPairs(candidate, [&](const AuditSegment& a, const AuditSegment& b) {
                if (overlapLength(a.a.xy(), a.b.xy(), b.a.xy(), b.b.xy()) > 1e-6) ++overlaps;
            });
            const auto preference = distributionPreference(component.second, paths, uses, corridors,
                config.overlap_offset, shift*config.overlap_offset);
            const Score score(crossing, surface_crossing, clearance, overlaps, self, preference[0], preference[1], preference[2]);
            if (score < best_score) { best_score = score; best = candidate; best_full = candidate_full; best_radii = radii; }
            if (config.overlap_offset <= 0) break;
        }
        if (mirror) for (size_t index : component.second) paths[index].rank = paths[index].count - 1 - paths[index].rank;
        if (config.overlap_offset <= 0) break;
        }
        }
        for (auto& corridor : corridors) corridor.allow_surface_shift = false;
        for (size_t i = 0; i < component.second.size(); ++i) {
            paths[component.second[i]].body = best[i];
            paths[component.second[i]].geometry = best_full[i];
            paths[component.second[i]].radius = best_radii[i];
        }
    }
    // 轨位/圆角只能改变骨干附近的形态，不能绕开进出段上的有限边缘；在最终几何上局部绕行。
    const auto barriers = nonmotorBoundarySegments(input,config);
    // 先在原始边界绕行前合并端部短窗口，避免绕行阶段把已平滑的挂接再次折断；
    // 绕行随后仍按完整边界/障碍/来源面约束复核并保留固定端点。
    for (auto& path : paths) repairAttachmentWindows(path,input,motors,config,barriers,uses,corridors);
    for (auto& path : paths) avoidNonmotorBoundaries(path,input,motors,config,barriers,uses,corridors,paths);
    std::vector<ConnectivityCurve> result(paths.size()); std::vector<LineString2d> bodies;
    const auto issue = [&](size_t i, const std::string& reason, bool hard) {
        auto& out = result[i];
        if (hard) { out.status = CurveStatus::Infeasible; out.violation.type = ViolationInfo::InfeasibilityType::TopologicalBlock; }
        else if (out.status == CurveStatus::OK) { out.status = CurveStatus::Degraded; out.violation.type = ViolationInfo::InfeasibilityType::NarrowPassage; }
        if (out.violation.reason.find(reason) == std::string::npos) {
            if (!out.violation.reason.empty()) out.violation.reason += "; ";
            out.violation.reason += reason;
        }
    };
    for (size_t i = 0; i < paths.size(); ++i) {
        const auto& p = paths[i]; const auto& c = *p.connection; auto& out = result[i];
        out.id = c.id; out.entry_lane_id = c.entry_lane_id; out.exit_lane_id = c.exit_lane_id; out.turn_type = c.turn_type;
        out.lane_type = MotorLaneType::NonMotorLane; out.logic_type = MotorLogicType::NonMotorLogic;
        out.fixed_shape = c.fixed_shape; out.update_shape = false; out.geometry = p.geometry;
        bodies.push_back(p.body); 
        if (pathHitsObstacle(out.geometry, input.obstacles)) {
            out.status = CurveStatus::WarnA2; out.violation.type = ViolationInfo::InfeasibilityType::ForcedCross;
            out.violation.reason = "non-motor backbone intersects an obstacle";
        }
        if (!p.connected) issue(i, "no connected backbone route between lane attachments", true);
        if (layoutSelfIntersects(out.geometry)) issue(i, "non-motor path self-intersection", true);
        std::string source_id;
        if (layoutLeavesSourceSurface(out.geometry,p,uses,corridors,&source_id))
            issue(i, "non-motor path leaves source surface: " + source_id, true);
        if (pathOutsideFence(out.geometry, input.area.geometry)) issue(i, "non-motor path leaves intersection fence", false);
        if (pathHitsBoundary(out.geometry, input.boundaries)) issue(i, "non-motor path intersects Boundary/RoadEdge", true);
        // 平滑或窗口拼接也不能使绕行净距失效；固定端点/固定形态不满足时保留诊断并硬失败。
        if (std::isfinite(nonmotorFirstBoundaryConflict(out.geometry,barriers)))
            issue(i, "non-motor path violates Boundary avoidance clearance", true);
        if (minimumRoadEdgeDistance(out.geometry, input.boundaries) + 1e-9 < config.road_edge_clearance) issue(i, "non-motor path is too close to RoadEdge", false);
        for (const auto& motor : motors) {
            const Lane* en = input.findLane(c.entry_lane_id); const Lane* ex = input.findLane(c.exit_lane_id);
            const double clearance = en && ex && en->motoway && !ex->motoway ? config.motor_to_nm_clearance : 0;
            if (clearance > 0 && minimumPolylineDistance(out.geometry, motor.geometry) + 1e-9 < clearance)
                issue(i, "motor clearance conflict: " + motor.id, false);
            bool overlap = false;
            for (size_t a = 1; a < out.geometry.points.size() && !overlap; ++a)
                for (size_t b = 1; b < motor.geometry.points.size() && !overlap; ++b)
                    overlap = overlapLength(out.geometry.points[a - 1].xy(), out.geometry.points[a].xy(),
                                            motor.geometry.points[b - 1].xy(), motor.geometry.points[b].xy()) > 1e-6;
            if (overlap) issue(i, "motor overlap conflict: " + motor.id, true);
        }
    }
    visitLayoutPairs(bodies, [&](const AuditSegment& a, const AuditSegment& b) {
        const double overlap = overlapLength(a.a.xy(), a.b.xy(), b.b.xy(), b.a.xy());
        if (overlap > 1e-6) {
            const std::string reason = "effective backbone overlap: " + paths[a.path].connection->id + "/" + paths[b.path].connection->id +
                " segments " + std::to_string(a.index) + "/" + std::to_string(b.index) + " length=" + std::to_string(overlap);
            issue(a.path, reason, true); issue(b.path, reason, true);
        }
    });
    std::map<const Connectivity*, size_t> positions;
    for (size_t i = 0; i < paths.size(); ++i) positions[paths[i].connection] = i;
    std::vector<ConnectivityCurve> ordered;
    for (const auto& c : input.connectivities) ordered.push_back(result[positions[&c]]);
    return ordered;
}
