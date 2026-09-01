// 家族横向阶梯诊断：直接调用 UTurnFamilyBuilder::familyLateralLadder，
// 打印指定数据集里每个掉头家族的成员次序、走廊宽度与最终分档站位。
//
// 用法: diag_family_ladder <input.json> [conn_id ...]
// 不给 conn_id 时遍历全部几何掉头。
//
// 输入必须与生成流程同源（InputNormalizer + ConnectivityDirectionNormalizer），
// 否则端点/切向不在同一坐标系里，走廊宽度与站位不可与会话日志直接比对。
// 名义步长口径与 connectivity_generation_session.cpp 的
// uturnSharedEndpointStagger 一致：物理场景且家族规模 >= 3 时为 0.25，否则 0.01。

#include "constraints/cluster_order.h"
#include "io/iodata_json.h"
#include "preprocessing/uturn_family_builder.h"
#include "toolkits/toolkits.h"
#include "types.h"

#include <cstdio>
#include <set>
#include <string>
#include <vector>

using namespace isg;

namespace {

bool pureGeometricInput(const IntersectionInput& input) {
    return input.obstacles.empty() && input.boundaries.empty() &&
           input.crosswalks.empty() && input.stop_lines.empty();
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        std::fprintf(stderr, "用法: %s <input.json> [conn_id ...]\n", argv[0]);
        return 2;
    }
    IntersectionInput input = IntersectionIO::loadFromFile(argv[1]);
    input = InputNormalizer(input);
    ConnectivityDirectionNormalizer(input, ConnectivityDirectionConfig());

    std::set<ConnId> wanted;
    for (int i = 2; i < argc; ++i)
        wanted.insert(ConnId(argv[i]));

    const UTurnFamilyBuilder builder;
    ClusterOrderSolver solver;
    solver.build(input.connectivities, input.lanes, input.lane_groups,
                 input.crosswalks);
    const UTurnAlignmentScope scope = UTurnAlignmentScope::LaneEndpoint;
    const bool require_shared_pair = true;
    const bool physical = !pureGeometricInput(input);

    std::set<ConnId> reported;
    for (const auto& conn : input.connectivities) {
        if (!builder.isGeometricUTurn(conn, input))
            continue;
        if (!wanted.empty() && !wanted.count(conn.id))
            continue;
        if (reported.count(conn.id))
            continue;
        const std::vector<const Connectivity*> component =
            builder.alignmentComponent(conn, input, scope, &solver,
                                       require_shared_pair);
        const double nominal_step =
            physical && component.size() >= 3 ? 0.25 : 0.01;
        std::printf("==== family of conn %s (size=%zu, step=%.3f) ====\n",
                    conn.id.c_str(), component.size(), nominal_step);
        for (const Connectivity* member : component) {
            reported.insert(member->id);
            const UTurnFamilyLadder ladder = builder.familyLateralLadder(
                *member, input, scope, nominal_step, &solver,
                require_shared_pair);
            const UTurnFamilyRank rank = builder.radiusRank(
                *member, input, scope, &solver, require_shared_pair);
            const double radius = builder.radiusKey(*member, input);
            std::printf(
                "  conn %4s  radius=%8.4f  rank=%zu/%zu  nominal=%.4f  "
                "corridor_cap=%.4f  entry=%.4f  exit=%.4f\n",
                member->id.c_str(), radius, rank.reverse_radius_rank,
                rank.family_size, ladder.nominal_stagger, 0.25 * radius,
                ladder.entry_stagger, ladder.exit_stagger);
        }
    }
    return 0;
}
