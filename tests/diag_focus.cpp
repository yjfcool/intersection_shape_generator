// Focused diagnostic for conn 74|76 (100000643.json) and conn 26 (intersection_jd.json)
// - For 74|76: dump full geometry, find ALL interior intersections, check crosswalk position
// - For 26: dump geometry, check G1, self-intersection, fence overflow
#include "constraints/cluster_order.h"
#include "constraints/fence_check.h"
#include "curve/curve_utils.h"
#include "intersection_shape_generator.h"
#include "io/iodata_json.h"
#include "optimizer/sdf_field.h"
#include "utils.h"
#include <cstdio>
#include <unordered_map>
#include <unordered_set>

using namespace isg;

static std::unordered_map<ConnId, const ConnectivityCurve*> curveMap(const IntersectionOutput& output) {
    std::unordered_map<ConnId, const ConnectivityCurve*> curves;
    for (const auto& cc : output.connectivity_curves) curves[cc.id] = &cc;
    return curves;
}

static void dumpConn(const char* id, const IntersectionInput& input, const IntersectionOutput& output) {
    auto cmap = curveMap(output);
    auto it = cmap.find(id);
    if (it == cmap.end() || !it->second->curve) {
        printf("  %s: missing\n", id);
        return;
    }
    const auto& c = *it->second->curve;
    auto entry = input.entryPtDir(it->second->entry_lane_id);
    auto exit_ = input.exitPtDir(it->second->exit_lane_id);
    Vec2d p0 = entry.first; Vec2d T0 = entry.second.norm() > 1e-8 ? entry.second.normalized() : Vec2d(1,0);
    Vec2d p1 = exit_.first;  Vec2d T1 = exit_.second.norm() > 1e-8 ? exit_.second.normalized() : Vec2d(1,0);
    printf("\n==== %s ====\n", id);
    printf("  entry_lane=%s exit_lane=%s\n", it->second->entry_lane_id.c_str(), it->second->exit_lane_id.c_str());
    printf("  p0=(%.2f,%.2f) T0=(%.3f,%.3f)\n", p0.x(), p0.y(), T0.x(), T0.y());
    printf("  p1=(%.2f,%.2f) T1=(%.3f,%.3f)\n", p1.x(), p1.y(), T1.x(), T1.y());
    printf("  chord=%.3f arc=%.3f arc/chord=%.3f maxk=%.3f\n",
           (p1-p0).norm(), c.arcLength(),
           (p1-p0).norm()>1e-6 ? c.arcLength()/(p1-p0).norm() : 0.0,
           c.maxCurvature(40));
    Vec2d st = c.startTan().norm() > 1e-8 ? c.startTan().normalized() : Vec2d(1,0);
    Vec2d et = c.endTan().norm()   > 1e-8 ? c.endTan().normalized()   : Vec2d(1,0);
    printf("  G1 cos0=%.4f cos1=%.4f\n", st.dot(T0), et.dot(T1));
    printf("  status=%d reason=%s\n", (int)it->second->status, it->second->violation.reason.c_str());
    printf("  max_obstacle_pen=%.3f max_fence_overflow=%.3f exempt_crosses=%zu\n",
           it->second->violation.max_obstacle_penetration,
           it->second->violation.max_fence_overflow,
           it->second->violation.exempt_crosses.size());
    // Control points
    for (size_t s = 0; s < c.segs.size(); ++s) {
        printf("  seg[%zu] ctrl: ", s);
        for (int k = 0; k < 4; ++k)
            printf("(%.2f,%.2f)%s", c.segs[s].ctrl[k].x(), c.segs[s].ctrl[k].y(), k<3?" ":"");
        printf("\n");
    }
    // Self-intersection check
    bool self_x = curveSelfIntersectsBusiness(c, 1.0);
    printf("  self_intersects(ep=1.0): %s\n", self_x ? "YES" : "no");
    // Fence overflow
    if (!input.area.geometry.outer.empty()) {
        auto pts = c.sampleByArcLength(80);
        double max_out = 0; int n_out = 0;
        for (auto& pt : pts) {
            if (!polygonContains(input.area.geometry, pt)) {
                double d = pointToPolygonDist(pt, input.area.geometry);
                if (d > 0.10) { max_out = std::max(max_out, d); ++n_out; }
            }
        }
        printf("  fence_overflow: max=%.3f n_out_pts=%d/%zu\n", max_out, n_out, pts.size());
    }
}

static void findAllIntersections(const char* a_id, const char* b_id,
                                  const IntersectionInput& input, const IntersectionOutput& output) {
    auto cmap = curveMap(output);
    auto ia = cmap.find(a_id);
    auto ib = cmap.find(b_id);
    if (ia == cmap.end() || ib == cmap.end() || !ia->second->curve || !ib->second->curve) {
        printf("  %s-%s: missing\n", a_id, b_id);
        return;
    }
    const auto& a = *ia->second->curve;
    const auto& b = *ib->second->curve;
    printf("\n  ---- Intersections %s <-> %s ----\n", a_id, b_id);
    printf("  a: start=(%.2f,%.2f) end=(%.2f,%.2f)\n", a.startPt().x(), a.startPt().y(), a.endPt().x(), a.endPt().y());
    printf("  b: start=(%.2f,%.2f) end=(%.2f,%.2f)\n", b.startPt().x(), b.startPt().y(), b.endPt().x(), b.endPt().y());
    auto pa = a.sampleByArcLength(120);
    auto pb = b.sampleByArcLength(120);
    int n_interior = 0, n_endpoint = 0;
    double min_interior_d = 1e18, max_interior_d = 0;
    for (size_t i = 0; i+1 < pa.size(); ++i) {
        for (size_t j = 0; j+1 < pb.size(); ++j) {
            Vec2d ipt;
            if (!segmentsIntersect(pa[i], pa[i+1], pb[j], pb[j+1], &ipt)) continue;
            double da = std::min((ipt - a.startPt()).norm(), (ipt - a.endPt()).norm());
            double db = std::min((ipt - b.startPt()).norm(), (ipt - b.endPt()).norm());
            double dmin = std::min(da, db);
            if (dmin > 1.5) {
                ++n_interior;
                min_interior_d = std::min(min_interior_d, dmin);
                max_interior_d = std::max(max_interior_d, dmin);
                if (n_interior <= 6)
                    printf("    interior isect at (%.2f,%.2f)  d_to_a=%.2f d_to_b=%.2f\n",
                           ipt.x(), ipt.y(), da, db);
            } else {
                ++n_endpoint;
            }
        }
    }
    printf("    total: %d interior (d_min=%.2f, d_max=%.2f), %d endpoint-zone\n",
           n_interior, min_interior_d==1e18?0.0:min_interior_d, max_interior_d, n_endpoint);
    // Cluster pair info
    ClusterOrderSolver solver;
    solver.build(input.connectivities, input.lanes, input.lane_groups, input.crosswalks);
    const CurvePair* pair = nullptr;
    for (auto& p : solver.pairs()) {
        if ((p.id_a==a_id && p.id_b==b_id) || (p.id_a==b_id && p.id_b==a_id)) { pair = &p; break; }
    }
    if (pair) {
        printf("    pair: exempt=%d shared_ep=%d expected_side=%d ref_perp=(%.3f,%.3f)\n",
               (int)pair->exempt, (int)pair->shared_endpoint, pair->expected_side,
               pair->ref_perp.x(), pair->ref_perp.y());
        printf("    pair_ids: id_a=%s id_b=%s\n",
               pair->id_a.c_str(), pair->id_b.c_str());
    }
}

static void dumpCrosswalks(const IntersectionInput& input) {
    printf("\n  ---- Crosswalks (count=%zu) ----\n", input.crosswalks.size());
    for (auto& cw : input.crosswalks) {
        printf("  cw %s: outer_pts=%zu\n", cw.id.c_str(), cw.geometry.outer.size());
        for (auto& p : cw.geometry.outer)
            printf("    (%.2f,%.2f)\n", p.x(), p.y());
    }
}

int main(int argc, char** argv) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <data.json>\n", argv[0]);
        return 1;
    }
    std::string path = argv[1];
    IntersectionInput input = IntersectionIO::loadFromFile(path);
    fprintf(stderr, "Loaded %s: %zu conns, %zu lanes, %zu crosswalks\n",
            path.c_str(), input.connectivities.size(), input.lanes.size(), input.crosswalks.size());

    IntersectionShapeGenerator gen;
    IntersectionOutput output;
    if (!gen.generate(input, output)) {
        fprintf(stderr, "Generation failed\n");
        return 1;
    }
    fprintf(stderr, "Generated %zu curves\n", output.connectivity_curves.size());

    dumpCrosswalks(input);

    // Pick which conns to inspect from argv[2..]
    if (argc >= 3) {
        for (int i = 2; i < argc; ++i) {
            const char* target = argv[i];
            // If it contains a |, treat as a pair
            std::string s(target);
            auto bar = s.find('|');
            if (bar != std::string::npos) {
                std::string a = s.substr(0, bar);
                std::string b = s.substr(bar+1);
                findAllIntersections(a.c_str(), b.c_str(), input, output);
            } else {
                dumpConn(target, input, output);
            }
        }
    }
    return 0;
}
