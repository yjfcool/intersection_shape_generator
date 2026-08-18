// Dump geometry of U-turns and their crossing siblings as JSON for inspection.
#include "constraints/cluster_order.h"
#include "curve/curve_utils.h"
#include "intersection_shape_generator.h"
#include "io/iodata_json.h"
#include "optimizer/sdf_field.h"
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <unordered_map>

using namespace isg;

int main(int argc, char** argv) {
    std::string path = (argc > 1) ? argv[1]
        : (std::string(PROJECT_ROOT_DIR) + "/datas/100000643-1.json");
    IntersectionInput input = IntersectionIO::loadFromFile(path);
    if (std::getenv("ISG_FOCUS_UTURN_FAMILY")) {
        input.connectivities.erase(
            std::remove_if(
                input.connectivities.begin(), input.connectivities.end(),
                [](const Connectivity& conn) {
                    return conn.id != "48" && conn.id != "49" &&
                           conn.id != "50" && conn.id != "51" &&
                           conn.id != "52";
                }),
            input.connectivities.end());
    }

    IntersectionShapeGenerator gen;
    IntersectionOutput output;
    if (!gen.generate(input, output)) return 1;

    std::unordered_map<ConnId, const ConnectivityCurve*> cmap;
    for (const auto& cc : output.connectivity_curves) cmap[cc.id] = &cc;

    // Dump a few problem cases.
    auto dump = [&](const char* id) {
        auto it = cmap.find(id);
        if (it == cmap.end() || !it->second->curve) return;
        const auto& c = *it->second->curve;
        auto entry = input.entryPtDir(it->second->entry_lane_id);
        auto exit_ = input.exitPtDir(it->second->exit_lane_id);
        printf("\n==== %s ====\n", id);
        printf("  status=%d reason=%s segments=%d\n",
               static_cast<int>(it->second->status),
               it->second->violation.reason.c_str(), c.numSegments());
        printf("  p0=(%.2f,%.2f) T0=(%.3f,%.3f)\n", entry.first.x(), entry.first.y(),
               entry.second.x(), entry.second.y());
        printf("  p1=(%.2f,%.2f) T1=(%.3f,%.3f)\n", exit_.first.x(), exit_.first.y(),
               exit_.second.x(), exit_.second.y());
        printf("  arc=%.2f  maxk=%.3f  pts={", c.arcLength(), c.maxCurvature(40));
        auto pts = c.sampleByArcLength(20);
        for (size_t i = 0; i < pts.size(); ++i) {
            if (i) printf(",");
            printf("(%.1f,%.1f)", pts[i].x(), pts[i].y());
        }
        printf("}\n");
        // Control points
        for (size_t s = 0; s < c.segs.size(); ++s) {
            printf("  seg[%zu] ctrl: ", s);
            for (int k = 0; k < 4; ++k)
                printf("(%.2f,%.2f)%s", c.segs[s].ctrl[k].x(), c.segs[s].ctrl[k].y(),
                       k<3?" ":"");
            printf("\n");
        }
    };
    const std::vector<const char*> ids =
        input.id == "100000012"
            ? std::vector<const char*>{"19", "20", "21", "22"}
        : input.id == "110003285"
            ? std::vector<const char*>{"48", "49", "50", "51", "52"}
            : std::vector<const char*>{
                "67", "68", "69", "107", "108", "109", "88", "90",
                "94", "99", "100", "101", "91", "116", "64"};
    for (const char* id : ids)
        dump(id);
    return 0;
}
