// 私有实现：在最终挂接路径的冲突邻域绕行有限Boundary，所有类型都禁止穿越。
// 每次仅替换一个弧长窗口，窗口外折线形态和线性高程保持不变，仅合并共线冗余采样；
// 图规模、窗口扩张和修复次数均有上限。
struct NonmotorBoundarySegment {
    Vec3d a, b;
    double clearance;
};

/// 每路口复用边界线段；所有类型统一应用绕行净距，RoadEdge同时满足原有净距。
std::vector<NonmotorBoundarySegment> nonmotorBoundarySegments(const IntersectionInput& input,
    const NonmotorConifg& config) {
    std::vector<NonmotorBoundarySegment> result;
    for (const auto& boundary : input.boundaries)
        for (size_t i = 1; i < boundary.geometry.points.size(); ++i)
            if ((boundary.geometry.points[i].xy()-boundary.geometry.points[i-1].xy()).norm() > 1e-8)
                result.push_back({boundary.geometry.points[i-1],boundary.geometry.points[i],
                    std::max(0.02,std::max(config.boundary_avoidance_clearance,
                        boundary.type == Boundary::Type::RoadEdge ? config.road_edge_clearance : 0.0))});
    return result;
}

/// 净距和求交均使用完整线段；包围盒预筛不能替代精确检查。
bool nonmotorBoundaryVisible(const Vec3d& a, const Vec3d& b,
    const std::vector<NonmotorBoundarySegment>& barriers, const IntersectionInput& input) {
    for (const auto& edge : barriers) {
        const double clearance = edge.clearance;
        if ((a.xy().cwiseMin(b.xy()).array() > edge.a.xy().cwiseMax(edge.b.xy()).array()+clearance).any() ||
            (edge.a.xy().cwiseMin(edge.b.xy()).array() > a.xy().cwiseMax(b.xy()).array()+clearance).any()) continue;
        if (segmentDistance(a,b,edge.a,edge.b) < clearance-1e-7) return false;
    }
    LineString2d line; line.points = {a,b};
    if (pathHitsObstacle(line,input.obstacles)) return false;
    // 精确裁剪fence，避免长可见边跨过凹口或细洞而稀疏采样漏检。
    if (!input.area.geometry.outer.empty()) {
        const Vec2d delta = b.xy()-a.xy(); const double length = delta.norm();
        if (length < 1e-8) return false;
        bool covered = false;
        for (const auto& interval : nonmotor_detail::sourceSurfaceIntervals(a.xy(),delta/length,input.area.geometry))
            if (interval.first <= 1e-7 && interval.second >= length-1e-7) covered = true;
        if (!covered) return false;
    }
    return true;
}

/// 有限边界按净距与圆角余量作方端外扩，取联合轮廓节点，绕真实端点而非穿过边界或无限延伸。
/// 至多512个节点的可见图用Dijkstra求最短路径，最多检查局部12米外扩范围。
std::vector<Vec3d> nonmotorBoundaryDetour(const Vec3d& start, const Vec3d& end,
    const std::vector<NonmotorBoundarySegment>& barriers, const IntersectionInput& input,
    double padding, double search_radius, const LayoutPath& path,
    const std::vector<LayoutUse>& uses, const std::vector<LayoutCorridor>& corridors) {
    constexpr double scale = 10000;
    const Vec2d low = start.xy().cwiseMin(end.xy())-Vec2d(search_radius,search_radius);
    const Vec2d high = start.xy().cwiseMax(end.xy())+Vec2d(search_radius,search_radius);
    ClipperLib::Paths expanded;
    for (const auto& edge : barriers) {
        if ((edge.a.xy().cwiseMin(edge.b.xy()).array() > high.array()).any() ||
            (edge.a.xy().cwiseMax(edge.b.xy()).array() < low.array()).any()) continue;
        ClipperLib::Path segment;
        for (const auto& point : {edge.a,edge.b})
            segment.emplace_back(std::llround(point.x()*scale),std::llround(point.y()*scale));
        ClipperLib::ClipperOffset offset;
        offset.AddPath(segment,ClipperLib::jtMiter,ClipperLib::etOpenSquare);
        ClipperLib::Paths polygons; offset.Execute(polygons,(edge.clearance+padding)*scale);
        expanded.insert(expanded.end(),polygons.begin(),polygons.end());
    }
    ClipperLib::Clipper clipper; clipper.AddPaths(expanded,ClipperLib::ptSubject,true);
    ClipperLib::Paths contours; clipper.Execute(ClipperLib::ctUnion,contours,ClipperLib::pftNonZero,ClipperLib::pftNonZero);
    std::vector<Vec3d> nodes{start,end};
    for (const auto& contour : contours) for (const auto& point : contour) {
        Vec3d candidate(point.X/scale,point.Y/scale,0);
        if ((candidate.xy().array() < low.array()).any() || (candidate.xy().array() > high.array()).any()) continue;
        bool duplicate = false;
        for (const auto& other : nodes) if ((candidate.xy()-other.xy()).norm() < 1e-5) { duplicate = true; break; }
        if (!duplicate) nodes.push_back(candidate);
    }
    // 来源面也是路由约束，须在搜索时纳入，而不能选出面外最短路后只在末尾拒绝。
    // 加入面内轮廓及已有偏移路径的稀疏节点，给狭窄面内通道提供可见转向节点。
    const auto add_node = [&](const Vec3d& point) {
        if ((point.xy().array()<low.array()).any() || (point.xy().array()>high.array()).any()) return;
        for (const auto& other:nodes) if ((point.xy()-other.xy()).norm()<1e-5) return;
        nodes.push_back(point);
    };
    Vec3d last_sample = path.geometry.points.front();
    for (const auto& point:path.geometry.points) if ((point.xy()-last_sample.xy()).norm()>=2) {
        add_node(point); last_sample=point;
    }
    for (size_t use:path.uses) for (const auto& surface:corridors[uses[use].corridor].surfaces) {
        ClipperLib::Path polygon;
        for (const auto& point:surface.polygon->outer) polygon.emplace_back(std::llround(point.x()*scale),std::llround(point.y()*scale));
        if (!ClipperLib::Orientation(polygon)) std::reverse(polygon.begin(),polygon.end());
        ClipperLib::ClipperOffset inner; inner.AddPath(polygon,ClipperLib::jtMiter,ClipperLib::etClosedPolygon);
        ClipperLib::Paths rings; inner.Execute(rings,-std::max(0.02,surface.clearance+padding*0.1)*scale);
        for (const auto& ring:rings) for (const auto& point:ring) add_node(Vec3d(point.X/scale,point.Y/scale,0));
    }
    if (nodes.size() > 512) return {};
    const size_t count = nodes.size();
    std::vector<double> cost(count,std::numeric_limits<double>::infinity());
    std::vector<size_t> previous(count,kInvalidIndex); std::vector<bool> done(count,false);
    cost[0] = 0;
    for (size_t step = 0; step < count; ++step) {
        size_t current = kInvalidIndex;
        for (size_t i = 0; i < count; ++i) if (!done[i] &&
            (current == kInvalidIndex || cost[i] < cost[current])) current = i;
        if (current == kInvalidIndex || !std::isfinite(cost[current])) return {};
        if (current == 1) break;
        done[current] = true;
        for (size_t next = 0; next < count; ++next) {
            if (done[next]) continue;
            const double trial = cost[current]+(nodes[next].xy()-nodes[current].xy()).norm();
            if (trial >= cost[next] || !nonmotorBoundaryVisible(nodes[current],nodes[next],barriers,input)) continue;
            LineString2d candidate; candidate.points={nodes[current],nodes[next]};
            if (layoutLeavesSourceSurface(candidate,path,uses,corridors)) continue;
            cost[next] = trial; previous[next] = current;
        }
    }
    if (!std::isfinite(cost[1])) return {};
    std::vector<Vec3d> result;
    for (size_t index = 1; index != kInvalidIndex; index = previous[index]) result.push_back(nodes[index]);
    std::reverse(result.begin(),result.end());
    double length = 0;
    for (size_t i = 1; i < result.size(); ++i) length += (result[i].xy()-result[i-1].xy()).norm();
    double station = 0;
    for (size_t i = 0; i < result.size(); ++i) {
        if (i) station += (result[i].xy()-result[i-1].xy()).norm();
        result[i].z() = start.z()+(end.z()-start.z())*station/std::max(1e-9,length);
    }
    return result;
}

/// 找首个穿越或净距不足段的最近冲突站位。检查完整线段，不能漏掉顶点之间的贴边。
/// 固定车道端点也不豁免净距；端点本身空间不足时最终明确返回不可行。
double nonmotorFirstBoundaryConflict(const LineString2d& line,
    const std::vector<NonmotorBoundarySegment>& barriers) {
    double station = 0;
    for (size_t i = 1; i < line.points.size(); ++i) {
        const auto& a = line.points[i-1]; const auto& b = line.points[i];
        double first = std::numeric_limits<double>::infinity();
        for (const auto& edge : barriers) {
            if (segmentDistance(a,b,edge.a,edge.b) >= edge.clearance-1e-7) continue;
            Vec2d point;
            if (segmentsIntersect(a.xy(),b.xy(),edge.a.xy(),edge.b.xy(),&point)) {
                first = std::min(first,(point-a.xy()).norm());
            } else {
                // 两端点到对方线段的四个投影包含有限线段间的最近点。
                double nearest = pointToSegment(a.xy(),edge.a,edge.b).first, fraction = 0;
                const double last = pointToSegment(b.xy(),edge.a,edge.b).first;
                if (last < nearest) { nearest=last; fraction=1; }
                for (const auto& endpoint : {edge.a,edge.b}) {
                    const auto hit=pointToSegment(endpoint.xy(),a,b);
                    if (hit.first < nearest) { nearest=hit.first; fraction=hit.second; }
                }
                first=std::min(first,fraction*(b.xy()-a.xy()).norm());
            }
        }
        if (std::isfinite(first)) return station+first;
        station += (b.xy()-a.xy()).norm();
    }
    return std::numeric_limits<double>::infinity();
}

/// 只修复冲突窗口。局部最短折线先圆角并按0.2米采样，验收失败则缩小圆角，
/// 禁止将检测失败当成可用曲线；不可行时保留诊断几何，由最终硬约束报告失败。
void avoidNonmotorBoundaries(LayoutPath& path, const IntersectionInput& input,
    const std::vector<ConnectivityCurve>& motors, const NonmotorConifg& config,
    const std::vector<NonmotorBoundarySegment>& barriers, const std::vector<LayoutUse>& uses,
    const std::vector<LayoutCorridor>& corridors, const std::vector<LayoutPath>& all_paths) {
    if (path.fixed || path.geometry.points.size() < 2) return;
    if (!std::isfinite(nonmotorFirstBoundaryConflict(path.geometry,barriers))) return;
    for (int repair = 0; repair < 12; ++repair) {
        const double crossing = nonmotorFirstBoundaryConflict(path.geometry,barriers);
        if (!std::isfinite(crossing)) break;
        const auto original = path.geometry;
        const auto original_spatial = layoutSpatialScore(original,path,input,motors,config);
        const bool original_obstacle = pathHitsObstacle(original,input.obstacles);
        const bool original_outside = pathOutsideFence(original,input.area.geometry);
        std::vector<double> stations(original.points.size(),0);
        for (size_t i = 1; i < stations.size(); ++i)
            stations[i] = stations[i-1]+(original.points[i].xy()-original.points[i-1].xy()).norm();
        const auto at = [&](double station) {
            const auto upper = std::upper_bound(stations.begin(),stations.end(),station);
            size_t i = std::min(stations.size()-1,std::max(size_t(1),size_t(upper-stations.begin())));
            return lerp3(original.points[i-1],original.points[i],
                (station-stations[i-1])/std::max(1e-9,stations[i]-stations[i-1]));
        };
        // 进出段上的绕行必须接回各自偏移骨干，不能先回到公共车道延长线再分离，
        // 否则多条挂接会在绕行后新增一段完全重叠的中间直线。
        const auto nearest_station = [&](const Vec3d& point) {
            double nearest = 1e100, result = 0;
            for (size_t j = 1; j < original.points.size(); ++j) {
                const auto hit = pointToSegment(point.xy(),original.points[j-1],original.points[j]);
                if (hit.first < nearest) { nearest = hit.first; result = stations[j-1]+hit.second*(stations[j]-stations[j-1]); }
            }
            return result;
        };
        const double entry_station = path.raw.empty() ? 0 : nearest_station(path.raw.front());
        const double exit_station = path.raw.empty() ? stations.back() : nearest_station(path.raw.back());
        bool accepted = false;
        // 每条曲线独立的微小安全余量防止同一边缘绕行轨迹完全重叠，不限制不同曲线交叉。
        const double margin = 0.15+0.004*(&path-&all_paths.front());
        for (double window : {3.0,6.0,10.0,16.0}) {
            const double lo = std::max(0.0,std::min(crossing-window,crossing >= exit_station ? exit_station-1.0 : crossing-window));
            const double hi = std::min(stations.back(),std::max(crossing+window,crossing <= entry_station ? entry_station+1.0 : crossing+window));
            const Vec3d a = at(lo), b = at(hi);
            const double guard = std::min(1.2,(hi-lo)*0.15);
            const Vec2d da = lo < 1e-8 ? path.entry_dir.normalized() : (at(lo+0.05).xy()-a.xy()).normalized();
            const Vec2d db = hi > stations.back()-1e-8 ? path.exit_dir.normalized() : (b.xy()-at(hi-0.05).xy()).normalized();
            Vec3d ga = offsetXY(a,da*guard), gb = offsetXY(b,-db*guard);
            ga.z() = at(lo+guard).z(); gb.z() = at(hi-guard).z();
            if (!nonmotorBoundaryVisible(a,ga,barriers,input) || !nonmotorBoundaryVisible(gb,b,barriers,input)) {
                continue;
            }
            auto route = nonmotorBoundaryDetour(ga,gb,barriers,input,margin,12,path,uses,corridors);
            if (route.empty()) {
                continue;
            }
            route.insert(route.begin(),a); route.push_back(b);
            // corner_radius<=0遵循关闭圆角的配置；默认最多7个递减截距候选。
            for (int step = 0; step < (config.corner_radius > 0 ? 7 : 1); ++step) {
                const double radius = config.corner_radius > 0 ? std::ldexp(config.corner_radius,-step) : 0;
                auto local = nonmotor_detail::roundPolylineCorners(route,radius,&da,&db);
                bool safe = true;
                for (size_t j = 1; j < local.points.size(); ++j)
                    if (!nonmotorBoundaryVisible(local.points[j-1],local.points[j],barriers,input)) { safe = false; break; }
                if (!safe) {
                    continue;
                }
                for (size_t j = 1; radius > 0 && j+1 < local.points.size(); ++j) {
                    const Vec2d before = local.points[j].xy()-local.points[j-1].xy();
                    const Vec2d after = local.points[j+1].xy()-local.points[j].xy();
                    if (before.norm() > 0.05 && after.norm() > 0.05 &&
                        before.normalized().dot(after.normalized()) < std::cos(45*DEG2RAD)) { safe = false; break; }
                }
                if (!safe) continue;
                LineString2d trial;
                for (size_t j = 0; j < original.points.size() && stations[j] < lo-1e-8; ++j) appendUnique(trial.points,original.points[j]);
                for (const auto& point : local.points) appendUnique(trial.points,point);
                for (size_t j = 0; j < original.points.size(); ++j) if (stations[j] > hi+1e-8) appendUnique(trial.points,original.points[j]);
                trial.points.front() = original.points.front(); trial.points.back() = original.points.back();
                // 窗口端点不能把原有的一条共线车道直段切成可被误认作骨干的多个段。
                std::vector<Vec3d> compact;
                for (const auto& point : trial.points) {
                    appendUnique(compact,point);
                    while (compact.size() >= 3) {
                        const size_t n = compact.size();
                        const Vec2d first = compact[n-2].xy()-compact[n-3].xy(), last = compact[n-1].xy()-compact[n-2].xy();
                        if (first.dot(last) <= 0 || std::abs(cross2d(first,last)) > 1e-8*first.norm()*last.norm()) break;
                        const double f = first.norm()/(first.norm()+last.norm());
                        if (std::abs(compact[n-2].z()-(compact[n-3].z()+(compact[n-1].z()-compact[n-3].z())*f)) > 1e-5) break;
                        compact[n-2] = compact.back(); compact.pop_back();
                    }
                }
                trial.points = std::move(compact);
                if (layoutSelfIntersects(trial) || layoutLeavesSourceSurface(trial,path,uses,corridors) ||
                    (!original_obstacle && pathHitsObstacle(trial,input.obstacles)) ||
                    (!original_outside && pathOutsideFence(trial,input.area.geometry)) ||
                    layoutSpatialScore(trial,path,input,motors,config) > original_spatial) {
                    continue;
                }
                // 分项保护机动车约束，不能用消除边界穿越交换新增机动车重叠或净距亏欠。
                const Lane* entry_lane = input.findLane(path.connection->entry_lane_id);
                const Lane* exit_lane = input.findLane(path.connection->exit_lane_id);
                const double motor_clearance = entry_lane && exit_lane && entry_lane->motoway && !exit_lane->motoway
                    ? config.motor_to_nm_clearance : 0;
                for (const auto& motor : motors) {
                    if (minimumPolylineDistance(trial,motor.geometry)+1e-8 <
                        std::min(motor_clearance,minimumPolylineDistance(original,motor.geometry))) { safe=false; break; }
                    const auto overlaps = [&](const LineString2d& candidate) {
                        for (size_t j=1;j<candidate.points.size();++j)
                            for (size_t k=1;k<motor.geometry.points.size();++k)
                                if (overlapLength(candidate.points[j-1].xy(),candidate.points[j].xy(),
                                    motor.geometry.points[k-1].xy(),motor.geometry.points[k].xy())>1e-6) return true;
                        return false;
                    };
                    if (overlaps(trial) && !overlaps(original)) { safe=false; break; }
                }
                if (!safe) continue;
                path.geometry = std::move(trial); accepted = true; break;
            }
            if (accepted) break;
        }
        if (!accepted) {
            break;
        }
    }
    // 重新提取有效骨干范围：进入/退出直段可以保留共线部分，曲线绕行段参与重叠验收。
    path.body = path.geometry;
    while (path.body.points.size() > 2 && std::abs(cross2d(path.body.points[1].xy()-path.entry.xy(),path.entry_dir)) < 1e-7)
        path.body.points.erase(path.body.points.begin());
    while (path.body.points.size() > 2 && std::abs(cross2d(path.body.points[path.body.points.size()-2].xy()-path.exit.xy(),path.exit_dir)) < 1e-7)
        path.body.points.pop_back();
}
