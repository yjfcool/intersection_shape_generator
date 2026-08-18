// Check all data files for same-cluster non-endpoint intersection violations.
#include "constraints/cluster_order.h"
#include "curve/curve_utils.h"
#include "intersection_shape_generator.h"
#include "io/iodata_json.h"
#include "optimizer/sdf_field.h"
#include <cstdio>
#include <dirent.h>
#include <sys/stat.h>
#include <unordered_map>
#include <vector>

using namespace isg;

static int checkFile(const std::string& path) {
    IntersectionInput input = IntersectionIO::loadFromFile(path);
    if (input.connectivities.empty()) return 0;

    IntersectionShapeGenerator gen;
    IntersectionOutput output;
    if (!gen.generate(input, output)) {
        fprintf(stderr, "  GENERATION FAILED\n");
        return -1;
    }

    std::unordered_map<ConnId, const ConnectivityCurve*> cmap;
    for (const auto& cc : output.connectivity_curves) cmap[cc.id] = &cc;

    ClusterOrderSolver solver;
    solver.build(input.connectivities, input.lanes, input.lane_groups, input.crosswalks);

    int violations = 0;
    int total_constrained_pairs = 0;
    int total_exempt_pairs = 0;

    for (const auto& p : solver.pairs()) {
        if (p.exempt == CrossExemption::StructuralCross) {
            ++total_exempt_pairs;
            continue;
        }
        ++total_constrained_pairs;
        auto ia = cmap.find(p.id_a);
        auto ib = cmap.find(p.id_b);
        if (ia == cmap.end() || ib == cmap.end()) continue;
        if (!ia->second->curve || !ib->second->curve) continue;
        if (curvesIntersectBusiness(*ia->second->curve, *ib->second->curve, 1.5)) {
            ++violations;
            printf("  VIOLATION: %s <-> %s (exempt=%d, shared_ep=%d)\n",
                   p.id_a.c_str(), p.id_b.c_str(),
                   (int)p.exempt, (int)p.shared_endpoint);
        }
    }
    printf("  constrained_pairs=%d, exempt_pairs=%d, violations=%d\n",
           total_constrained_pairs, total_exempt_pairs, violations);
    return violations;
}

int main(int argc, char** argv) {
    std::string dir = (argc > 1) ? argv[1]
        : (std::string(PROJECT_ROOT_DIR) + "/datas");
    DIR* d = opendir(dir.c_str());
    if (!d) {
        fprintf(stderr, "Cannot open %s\n", dir.c_str());
        return 1;
    }
    int total_violations = 0;
    struct dirent* ent;
    while ((ent = readdir(d)) != nullptr) {
        std::string name = ent->d_name;
        if (name.size() < 6) continue;
        if (name.substr(name.size() - 5) != ".json") continue;
        std::string path = dir + "/" + name;
        printf("\n==== %s ====\n", name.c_str());
        int v = checkFile(path);
        if (v > 0) total_violations += v;
    }
    closedir(d);
    printf("\n==== TOTAL VIOLATIONS: %d ====\n", total_violations);
    return total_violations > 0 ? 1 : 0;
}
