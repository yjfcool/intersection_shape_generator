// Diagnostic: print cluster order details for a specific exit group
#include "constraints/cluster_order.h"
#include "intersection_shape_generator.h"
#include "io/iodata_json.h"
#include <cstdio>
#include <unordered_map>

using namespace isg;

int main(int argc, char** argv) {
    std::string path = (argc > 1) ? argv[1] : (std::string(PROJECT_ROOT_DIR) + "/datas/intersection_jd.json");
    IntersectionInput input = IntersectionIO::loadFromFile(path);

    ClusterOrderSolver solver;
    solver.build(input.connectivities, input.lanes, input.lane_groups, input.crosswalks);

    printf("==== Exit group orders ====\n");
    for (auto& kv : solver.exitGroupOrder()) {
        printf("Exit group %s: [", kv.first.c_str());
        for (size_t i = 0; i < kv.second.size(); ++i) {
            if (i) printf(", ");
            printf("%s(r%zu)", kv.second[i].c_str(), i);
        }
        printf("]\n");
    }

    printf("\n==== Entry group orders ====\n");
    for (auto& kv : solver.entryGroupOrder()) {
        printf("Entry group %s: [", kv.first.c_str());
        for (size_t i = 0; i < kv.second.size(); ++i) {
            if (i) printf(", ");
            printf("%s(r%zu)", kv.second[i].c_str(), i);
        }
        printf("]\n");
    }

    printf("\n==== Pairs involving conn 74 ====\n");
    for (auto& p : solver.pairs()) {
        if (p.id_a == "74" || p.id_b == "74") {
            printf("  pair (%s, %s): exempt=%d shared_ep=%d expected_side=%d ref_perp=(%.3f,%.3f)\n",
                   p.id_a.c_str(), p.id_b.c_str(), (int)p.exempt, (int)p.shared_endpoint,
                   p.expected_side, p.ref_perp.x(), p.ref_perp.y());
        }
    }

    printf("\n==== All pairs in exit group 43100280 (same_exit) ====\n");
    // Find conns in exit group 43100280
    std::unordered_map<ConnId, std::string> conn_exit_group;
    for (auto& c : input.connectivities) {
        conn_exit_group[c.id] = c.exitGroupId;
    }
    for (auto& p : solver.pairs()) {
        if (conn_exit_group[p.id_a] == "43100280" && conn_exit_group[p.id_b] == "43100280" &&
            p.id_a != p.id_b) {
            // Check if same exit group but different exit lanes
            const Connectivity* ca = nullptr, *cb = nullptr;
            for (auto& c : input.connectivities) {
                if (c.id == p.id_a) ca = &c;
                if (c.id == p.id_b) cb = &c;
            }
            if (ca && cb && ca->exit_lane_id != cb->exit_lane_id) {
                printf("  pair (%s, %s): exempt=%d shared_ep=%d expected_side=%d  entry_lanes=(%s,%s) exit_lanes=(%s,%s)\n",
                       p.id_a.c_str(), p.id_b.c_str(), (int)p.exempt, (int)p.shared_endpoint,
                       p.expected_side,
                       ca->entry_lane_id.c_str(), cb->entry_lane_id.c_str(),
                       ca->exit_lane_id.c_str(), cb->exit_lane_id.c_str());
            }
        }
    }
    return 0;
}
