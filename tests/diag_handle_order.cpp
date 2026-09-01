// Diagnostic: same-entry single-segment cubic handle ordering / middle-control
// polyline crossing for a shared-entry cluster family.
//
// Usage: diag_handle_order [json_path] [conn_id ...]
#include "constraints/cluster_order.h"
#include "curve/curve_utils.h"
#include "intersection_shape_generator.h"
#include "io/iodata_json.h"
#include <cstdio>
#include <string>
#include <unordered_map>
#include <vector>

using namespace isg;

int main(int argc, char** argv) {
    std::string path = (argc > 1) ? argv[1]
        : (std::string(PROJECT_ROOT_DIR) + "/datas/110003449.json");
    std::vector<ConnId> ids;
    for (int i = 2; i < argc; ++i) ids.push_back(argv[i]);
    if (ids.empty()) ids = {"4", "5", "6"};

    IntersectionInput input = IntersectionIO::loadFromFile(path);
    IntersectionShapeGenerator gen;
    IntersectionOutput output;
    if (!gen.generate(input, output)) { fprintf(stderr, "GEN FAILED\n"); return 1; }

    std::unordered_map<ConnId, const ConnectivityCurve*> cmap;
    for (const auto& cc : output.connectivity_curves) cmap[cc.id] = &cc;

    printf("==== %s ====\n", path.c_str());
    for (const ConnId& id : ids) {
        auto it = cmap.find(id);
        if (it == cmap.end() || !it->second->curve) { printf("conn %s: missing\n", id.c_str()); continue; }
        const BezierCurve& c = *it->second->curve;
        printf("conn %-3s nseg=%d  arc=%.3f  maxkappa=%.4f\n",
               id.c_str(), c.numSegments(), c.arcLength(), c.maxCurvature(60));
        for (int s = 0; s < c.numSegments(); ++s) {
            const auto& g = c.segs[s].ctrl;
            printf("   seg%d P0=(%.3f,%.3f) P1=(%.3f,%.3f) P2=(%.3f,%.3f) P3=(%.3f,%.3f)"
                   "  h0=%.3f h1=%.3f chord=%.3f\n",
                   s, g[0].x(), g[0].y(), g[1].x(), g[1].y(), g[2].x(), g[2].y(),
                   g[3].x(), g[3].y(), (g[1]-g[0]).norm(), (g[3]-g[2]).norm(),
                   (g[3]-g[0]).norm());
        }
    }

    printf("\n---- pairwise ----\n");
    for (size_t i = 0; i < ids.size(); ++i) {
        for (size_t j = i + 1; j < ids.size(); ++j) {
            auto ia = cmap.find(ids[i]), ib = cmap.find(ids[j]);
            if (ia == cmap.end() || ib == cmap.end()) continue;
            if (!ia->second->curve || !ib->second->curve) continue;
            const BezierCurve& a = *ia->second->curve;
            const BezierCurve& b = *ib->second->curve;
            bool hit = curvesIntersectBusiness(a, b, 1.5);
            printf("%s|%s curve-cross(ep=1.5)=%s", ids[i].c_str(), ids[j].c_str(),
                   hit ? "TRUE" : "false");
            if (a.numSegments() == 1 && b.numSegments() == 1) {
                const auto& ga = a.segs[0].ctrl;
                const auto& gb = b.segs[0].ctrl;
                Vec2d ipt;
                bool mid_hit = segmentsIntersect(ga[1], ga[2], gb[1], gb[2], &ipt);
                printf("  mid-poly-cross=%s", mid_hit ? "TRUE" : "false");
                if (mid_hit) printf(" at (%.3f,%.3f)", ipt.x(), ipt.y());
                printf("  h0: %s=%.3f %s=%.3f  h1: %s=%.3f %s=%.3f",
                       ids[i].c_str(), (ga[1]-ga[0]).norm(),
                       ids[j].c_str(), (gb[1]-gb[0]).norm(),
                       ids[i].c_str(), (ga[3]-ga[2]).norm(),
                       ids[j].c_str(), (gb[3]-gb[2]).norm());
                bool same_entry = (ga[0] - gb[0]).norm() < 1e-6;
                bool same_exit = (ga[3] - gb[3]).norm() < 1e-6;
                printf("  same_entry=%d same_exit=%d", (int)same_entry, (int)same_exit);
            }
            printf("\n");
            if (hit) {
                auto pts = curveCrossings(a, b, 0.02);
                for (const auto& p : pts) {
                    double d = distToAllEndpoints(p, a, b);
                    printf("    crossing (%.3f,%.3f) d_endpoint=%.3f\n", p.x(), p.y(), d);
                }
            }
        }
    }

    printf("\n---- cluster pair table ----\n");
    ClusterOrderSolver solver;
    solver.build(input.connectivities, input.lanes, input.lane_groups, input.crosswalks);
    for (const auto& p : solver.pairs()) {
        bool in_a = false, in_b = false;
        for (const ConnId& id : ids) { if (p.id_a == id) in_a = true; if (p.id_b == id) in_b = true; }
        if (!in_a || !in_b) continue;
        printf("pair %s|%s exempt=%d shared_ep=%d\n", p.id_a.c_str(), p.id_b.c_str(),
               (int)p.exempt, (int)p.shared_endpoint);
    }
    return 0;
}
