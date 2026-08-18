// Focused test: does curvesIntersectBusiness report 67-68 as crossing?
#include "constraints/cluster_order.h"
#include "curve/curve_utils.h"
#include "intersection_shape_generator.h"
#include "io/iodata_json.h"
#include <cstdio>

using namespace isg;

int main(int argc, char** argv) {
    std::string path = (argc > 1) ? argv[1]
        : (std::string(PROJECT_ROOT_DIR) + "/datas/100000643-1.json");
    IntersectionInput input = IntersectionIO::loadFromFile(path);

    IntersectionShapeGenerator gen;
    IntersectionOutput output;
    if (!gen.generate(input, output)) return 1;

    std::unordered_map<ConnId, const ConnectivityCurve*> cmap;
    for (const auto& cc : output.connectivity_curves) cmap[cc.id] = &cc;

    auto test_pair = [&](const char* a_id, const char* b_id) {
        auto ia = cmap.find(a_id);
        auto ib = cmap.find(b_id);
        if (ia == cmap.end() || ib == cmap.end() || !ia->second->curve || !ib->second->curve) {
            printf("  %s-%s: missing\n", a_id, b_id);
            return;
        }
        const auto& a = *ia->second->curve;
        const auto& b = *ib->second->curve;
        printf("  %s-%s: a.start=(%.2f,%.2f) a.end=(%.2f,%.2f)  b.start=(%.2f,%.2f) b.end=(%.2f,%.2f)\n",
               a_id, b_id,
               a.startPt().x(), a.startPt().y(), a.endPt().x(), a.endPt().y(),
               b.startPt().x(), b.startPt().y(), b.endPt().x(), b.endPt().y());
        bool hit = curvesIntersectBusiness(a, b, 1.5);
        printf("    curvesIntersectBusiness ep=1.5: %s\n", hit ? "TRUE" : "false");

        // Manual check: sample both and find all segment intersections
        auto pa = a.sampleByArcLength(80);
        auto pb = b.sampleByArcLength(80);
        int n_interior = 0, n_endpoint = 0;
        for (size_t i = 0; i+1 < pa.size(); ++i) {
            for (size_t j = 0; j+1 < pb.size(); ++j) {
                Vec2d ipt;
                if (!segmentsIntersect(pa[i], pa[i+1], pb[j], pb[j+1], &ipt)) continue;
                double da = std::min((ipt - a.startPt()).norm(), (ipt - a.endPt()).norm());
                double db = std::min((ipt - b.startPt()).norm(), (ipt - b.endPt()).norm());
                double dmin = std::min(da, db);
                if (dmin > 1.5) {
                    ++n_interior;
                    if (n_interior <= 3)
                        printf("    interior isect at (%.2f,%.2f)  d_to_a=%.2f d_to_b=%.2f\n",
                               ipt.x(), ipt.y(), da, db);
                } else {
                    ++n_endpoint;
                }
            }
        }
        printf("    total: %d interior, %d endpoint-zone\n", n_interior, n_endpoint);
    };

    printf("==== Pair-by-pair intersection analysis ====\n");
    test_pair("67", "68");
    test_pair("67", "107");
    test_pair("94", "99");
    test_pair("94", "101");
    test_pair("88", "99");
    test_pair("90", "99");
    test_pair("116", "91");
    return 0;
}
