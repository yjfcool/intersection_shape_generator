// 私有实现：共享轨位确定后合并端部近邻转角。只替换有曲率缺陷的有界窗口，
// 窗口外坐标/顺序不变；近反向入口与短Gap整体构造，不再逐采样点缩短圆角。
struct AttachmentWindowCandidate {
    nonmotor_detail::CornerBezier curve;
    double start, end, quality;
};

/// 独立于采样点距的三点圆曲率；用于发现窗口，候选自身另验解析一/二阶导数。
double attachmentPolylineCurvature(const Vec3d& a,const Vec3d& b,const Vec3d& c) {
    const Vec2d u=b.xy()-a.xy(),v=c.xy()-b.xy();
    const double denominator=u.norm()*v.norm()*(u+v).norm();
    return denominator>1e-12 ? 2*std::abs(cross2d(u,v))/denominator : 0;
}

/// 每端最多4种接回站距×3种车道站位×25组独立把手；仅最优32个进入空间审计。
/// 候选数不乘进共享轨位的28次布局。固定端点已有净距不足只允许保留，不能扩大违规。
void repairAttachmentWindows(LayoutPath& path,const IntersectionInput& input,
    const std::vector<ConnectivityCurve>& motors,const NonmotorConifg& config,
    const std::vector<NonmotorBoundarySegment>& barriers,
    const std::vector<LayoutUse>& uses,const std::vector<LayoutCorridor>& corridors) {
    if (path.fixed || config.corner_radius<=0 || path.geometry.points.size()<3) return;
    for (int side=0;side<2;++side) {
        LineString2d source=path.geometry;
        if (side) std::reverse(source.points.begin(),source.points.end());
        // 仅对真实挂接点紧邻短骨干段的窗口启用合并；长直骨干的普通圆角
        // 不应被这个局部修复器重新分配轨位（保持历史共享直线的确定性）。
        if (path.raw.size() < 2) continue;
        const size_t raw_index = side ? path.raw.size()-1 : 0;
        const size_t raw_neighbor = side ? raw_index-1 : raw_index+1;
        if ((path.raw[raw_neighbor].xy()-path.raw[raw_index].xy()).norm() > 6.0)
            continue;
        const Vec2d lane=side ? Vec2d(-path.exit_dir) : path.entry_dir;
        if (lane.norm()<1e-8) continue;
        const Vec2d direction=lane.normalized();
        std::vector<double> stations(source.points.size(),0);
        for (size_t i=1;i<stations.size();++i)
            stations[i]=stations[i-1]+(source.points[i].xy()-source.points[i-1].xy()).norm();
        const Vec3d raw_end=side ? path.raw.back() : path.raw.front();
        const double attachment_distance=(raw_end.xy()-source.points.front().xy()).norm();
        const double limit=std::min(stations.back()*0.5,std::max(12.0,attachment_distance+8.0));
        double last_bad=-1;
        for (size_t i=1;i+1<stations.size() && stations[i]<limit;++i)
            if (attachmentPolylineCurvature(source.points[i-1],source.points[i],source.points[i+1])>1.0)
                last_bad=stations[i];
        if (last_bad<0) continue;
        const auto sample=[&](double s,Vec2d& tangent) {
            size_t next=std::upper_bound(stations.begin(),stations.end(),s)-stations.begin();
            next=std::max(size_t(1),std::min(next,stations.size()-1));
            const double length=stations[next]-stations[next-1];
            tangent=(source.points[next].xy()-source.points[next-1].xy())/std::max(1e-9,length);
            return Vec3d(source.points[next-1]+(source.points[next]-source.points[next-1])*((s-stations[next-1])/std::max(1e-9,length)));
        };
        std::vector<AttachmentWindowCandidate> candidates;
        for (double extra : {1.5,3.0,5.0,8.0}) {
            const double end=std::min(stations.back()-0.5,last_bad+extra);
            Vec2d tangent; const Vec3d b=sample(end,tangent);
            if (tangent.norm()<1e-8) continue;
            for (double start : {0.0,1.0,2.0}) {
                if (start>std::max(0.0,attachment_distance-2.0)) continue;
                const Vec3d a(source.points.front().xy()+direction*start,source.points.front().z());
                const double chord=(b.xy()-a.xy()).norm();
                if (chord<0.5) continue;
                for (double first : {0.25,0.5,0.75,1.0,1.5})
                    for (double last : {0.25,0.5,0.75,1.0,1.5}) {
                        const nonmotor_detail::CornerBezier curve{{{a,
                            Vec3d(a.xy()+direction*(chord*first),a.z()+(b.z()-a.z())/3),
                            Vec3d(b.xy()-tangent*(chord*last),b.z()-(b.z()-a.z())/3),b}}};
                        const double curvature=nonmotor_detail::attachmentBezierCurvature(curve);
                        if (curvature>1.0) continue;
                        // 优先短窗口和小偏离；达标后无需为追求极小曲率吞掉更多有效骨干。
                        candidates.push_back({curve,start,end,end+2*curvature+0.1*chord*(first+last)});
                    }
            }
        }
        std::stable_sort(candidates.begin(),candidates.end(),[](const AttachmentWindowCandidate& a,const AttachmentWindowCandidate& b) {return a.quality<b.quality;});
        const auto original_spatial=layoutSpatialScore(path.geometry,path,input,motors,config);
        const auto boundary_min=[&](const LineString2d& line,const NonmotorBoundarySegment& edge) {
            double minimum=std::numeric_limits<double>::infinity();
            for (size_t i=1;i<line.points.size();++i)
                minimum=std::min(minimum,segmentDistance(line.points[i-1],line.points[i],edge.a,edge.b));
            return minimum;
        };
        std::vector<double> original_clearance;
        for (const auto& edge:barriers) original_clearance.push_back(boundary_min(path.geometry,edge));
        for (size_t index=0;index<std::min(size_t(32),candidates.size());++index) {
            const auto& candidate=candidates[index];
            LineString2d trial;
            nonmotor_detail::appendCornerPoint(trial.points,source.points.front());
            for (const auto& p:nonmotor_detail::sampleCornerBezier(candidate.curve)) nonmotor_detail::appendCornerPoint(trial.points,p);
            for (size_t j=0;j<source.points.size();++j)
                if (stations[j]>candidate.end+1e-8) nonmotor_detail::appendCornerPoint(trial.points,source.points[j]);
            if (side) std::reverse(trial.points.begin(),trial.points.end());
            if (layoutSelfIntersects(trial) || layoutLeavesSourceSurface(trial,path,uses,corridors)) continue;
            const auto spatial=layoutSpatialScore(trial,path,input,motors,config);
            if (spatial.first>original_spatial.first || spatial.second>original_spatial.second) continue;
            bool safe=true;
            for (size_t j=0;j<barriers.size();++j)
                if (boundary_min(trial,barriers[j])+1e-5<std::min(barriers[j].clearance,original_clearance[j])) {safe=false;break;}
            if (!safe) continue;
            path.geometry=std::move(trial); break;
        }
    }
    // 有效骨干重新从最终输出提取，曲线过渡参与重叠验收，只有车道共线进出段可剔除。
    path.body=path.geometry;
    while (path.body.points.size()>2 && std::abs(cross2d(path.body.points[1].xy()-path.entry.xy(),path.entry_dir))<1e-7)
        path.body.points.erase(path.body.points.begin());
    while (path.body.points.size()>2 && std::abs(cross2d(path.body.points[path.body.points.size()-2].xy()-path.exit.xy(),path.exit_dir))<1e-7)
        path.body.points.pop_back();
}
