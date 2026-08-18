// Print all cluster pairs involving U-turns
#include "constraints/cluster_order.h"
#include "intersection_shape_generator.h"
#include "io/iodata_json.h"
#include <cstdio>

using namespace isg;

int main(int argc, char** argv) {
    std::string path = (argc > 1) ? argv[1]
        : (std::string(PROJECT_ROOT_DIR) + "/datas/100000643-1.json");
    IntersectionInput input = IntersectionIO::loadFromFile(path);

    ClusterOrderSolver solver;
    solver.build(input.connectivities, input.lanes, input.lane_groups, input.crosswalks);

    printf("==== Cluster pairs involving U-turns ====\n");
    printf("%-6s %-6s %-12s %-12s %-6s %-12s %-12s %-8s %-6s\n",
           "A", "B", "egrp_A", "egrp_B", "same_e", "xgrp_A", "xgrp_B", "exempt", "side");
    auto findC = [&](const ConnId& id) -> const Connectivity* {
        for (auto& c : input.connectivities) if (c.id == id) return &c;
        return nullptr;
    };
    auto isUTurn = [](const Connectivity& c) {
        return c.turn_type == ConnTurnType::UTurnLeft || c.turn_type == ConnTurnType::UTurnRight;
    };
    auto exempt_str = [](CrossExemption e) -> const char* {
        switch (e) {
            case CrossExemption::None: return "None";
            case CrossExemption::StructuralCross: return "Struct";
            case CrossExemption::ObstacleCross: return "Obst";
        }
        return "?";
    };
    for (const auto& p : solver.pairs()) {
        const Connectivity* ca = findC(p.id_a);
        const Connectivity* cb = findC(p.id_b);
        if (!ca || !cb) continue;
        if (!isUTurn(*ca) && !isUTurn(*cb)) continue;
        bool same_e = (ca->enterGroupId == cb->enterGroupId);
        bool same_x = (ca->exitGroupId == cb->exitGroupId);
        if (!same_e && !same_x) continue;  // only cluster pairs
        printf("%-6s %-6s %-12s %-12s %-6s %-12s %-12s %-8s %-6d\n",
               p.id_a.c_str(), p.id_b.c_str(),
               ca->enterGroupId.c_str(), cb->enterGroupId.c_str(),
               same_e ? "Y" : "N",
               ca->exitGroupId.c_str(), cb->exitGroupId.c_str(),
               exempt_str(p.exempt), p.expected_side);
    }
    return 0;
}
