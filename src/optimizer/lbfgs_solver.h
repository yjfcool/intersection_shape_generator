#pragma once
#include "types.h"
#include <functional>

namespace isg {

/// 求解结果
struct SolveResult {
    VecXd x;                    ///< 最优解
    double final_cost = 0;      ///< 最终代价
    int iterations = 0;         ///< 迭代次数
    bool converged = false;     ///< 是否收敛
    int fn_evals = 0;           ///< 代价函数总调用次数（诊断用）
};

/// 代价函数类型: 返回f，并填充grad
/// 函数必须正确计算梯度——lineSearch不再重新评估
using CostFn = std::function<double(const VecXd&, VecXd&)>;

/// LBFGS拟牛顿求解器：适用于大规模优化问题
class LBFGSSolver {
public:
    explicit LBFGSSolver(const LBFGSConfig& cfg = {}) : cfg_(cfg) {}

    /// 冷启动求解
    SolveResult solve(CostFn fn, const VecXd& x0);

    /// 热启动求解：从已有解继续优化
    SolveResult solveWarm(CostFn fn, const VecXd& x_warm);

private:
    LBFGSConfig cfg_;

    /// Strong Wolfe线搜索
    /// 返回接受的alpha，存储接受的(xn, fn_val, gn)
    /// 不在接受点重新调用fn——复用最后一次评估
    double lineSearch(CostFn& fn, const VecXd& x, const VecXd& dir,
                      double f0, const VecXd& g0, VecXd& xn, double& fn_val, VecXd& gn, int& evals);
};

}
