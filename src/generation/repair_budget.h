#pragma once

namespace isg {

/// 全局修复循环的兼容预算；默认值保持现有候选和迭代上限。
struct RepairBudget {
    int mode2_uturn_pair_passes;
    int mode2_uturn_pairs_per_pass;
    int shared_endpoint_pair_repairs;
    int early_straight_uturn_repairs;
    int final_straight_uturn_repairs;
    int pure_topology_passes;
    int pure_topology_repairs_per_pass;
    int shape_restores;
    int base_expression_passes;
    int pair_safe_passes;
    int pair_safe_phases;
    int regeneration_passes;
    int regenerations_per_pass;

    RepairBudget()
        : mode2_uturn_pair_passes(2), mode2_uturn_pairs_per_pass(4),
          shared_endpoint_pair_repairs(20),
          early_straight_uturn_repairs(4),
          final_straight_uturn_repairs(6), pure_topology_passes(3),
          pure_topology_repairs_per_pass(6), shape_restores(8),
          base_expression_passes(1), pair_safe_passes(3),
          // 0~3 为共享端点转向阶段，4 为非共享普通转弯，5 为共享近直行
          // 曲线的有界修复。
          pair_safe_phases(6), regeneration_passes(1),
          regenerations_per_pass(8) {}
};

}  // 命名空间 isg
