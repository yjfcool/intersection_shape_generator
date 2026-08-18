#pragma once

#include "types.h"
#include "sdf_field.h"
#include "lbfgs_solver.h"
#include <vector>

namespace isg {

/// 同簇相邻曲线（用于簇约束）
struct SiblingCurve {
    ConnId id;
    BezierCurve curve;
    bool exempt_a1 = false;           ///< true=结构性交叉或不同arm交叉交通→跳过惩罚
    double exempt_a2_radius = 0.0;    ///< >0=障碍物强制交叉（弱化权重0.4×，非零）
    int expected_side = 0;            ///< +1=本曲线在兄弟左侧; -1=右侧; 0=未知
    Vec2d ref_perp{0, 0};             ///< 配对特定的固定横向参考轴
    bool shared_endpoint = false;     ///< 共享首/尾连接点，仅用于连接点识别和侧向排序
    bool fixed_shape = false;         ///< 兄弟曲线来自保留的输入固有形态
};

/// 惩罚项权重
struct PenaltyWeights {
    double smooth = 1.0;
    double obstacle = 30.0;   ///< 最高优先级——避障
    double boundary = 10.0;
    double fence = 6.0;
    double cluster = 20.0;    ///< 拓扑序——低于避障
    double crosswalk = 0.5;

    /// 根据违规量自适应放大权重
    void update(double op, double bp, double fp, double cp) {
        if (op > 1e-3) obstacle = std::min(obstacle * 1.8, 240.0);
        if (bp > 1e-3) boundary = std::min(boundary * 2.0, 160.0);
        if (fp > 1e-3) fence = std::min(fence * 2.0, 96.0);
        if (cp > 1e-3) cluster = std::min(cluster * 1.5, 96.0);
    }
};

struct BoundarySegment {
    Vec2d a, b;
    double len = 0.0;
    double center_signed = 0.0;
    bool road_edge = false;
    bool center_side_reliable = true;
    bool a_is_boundary_endpoint = false;
    bool b_is_boundary_endpoint = false;
};

/// 惩罚项缓存：避免每次评估重复采样
struct PenaltyCostCache {
    struct SiblingPoly {
        std::vector<Vec2d> pts;
        bool exempt_a1;
        double a2_radius;
        int expected_side;
        Vec2d ref_perp;
        bool shared_endpoint = false;  ///< 共享首/尾连接点，仅用于连接点识别和侧向排序
        bool fixed_shape = false;      ///< 兄弟曲线来自保留的输入固有形态
    };

    std::vector<BoundarySegment> bnd_segs;
    std::vector<SiblingPoly> sib_polys;
    bool fence_empty = true;

    void rebuild(const std::vector<Boundary>& bounds, const std::vector<SiblingCurve>& sibls,
                 const Polygon2d& fence, const Vec2d& center);
};

/// 惩罚代价函数：包含平滑、避障、边界、围栏、同簇序、人行横道等项
class PenaltyCost {
public:
    BezierCurve proto;
    PenaltyWeights weights;
    const SDFField* sdf = nullptr;
    std::vector<Boundary> boundaries;
    Vec2d road_center{0, 0};
    Polygon2d fence;
    std::vector<SiblingCurve> siblings;
    std::vector<Crosswalk> crosswalks;
    double obstacle_clearance = 0.0;
    double road_edge_clearance = 0.0;
    bool full_param_mode = false;
    bool constrain_single_cubic_axes = true;
    Vec2d start_tan_dir{0, 0};  ///< 进入车道末点切向
    Vec2d end_tan_dir{0, 0};    ///< 退出车道首点切向

    /// 重建缓存（每次外层迭代前调用）
    void buildCache();

    /// 评估总代价并填充梯度
    double operator()(const VecXd& params, VecXd& grad);

    double evalSmooth(const BezierCurve& c) const;     ///< 平滑项（曲率平方积分）
    double evalObstacle(const BezierCurve& c) const;   ///< 障碍物穿透惩罚
    double evalBoundary(const BezierCurve& c) const;   ///< 边界线穿越惩罚
    double evalFence(const BezierCurve& c) const;      ///< 围栏越界惩罚
    double evalCluster(const BezierCurve& c) const;    ///< 同簇非交约束惩罚
    double evalCrosswalk(const BezierCurve& c) const;  ///< 人行横道斜穿惩罚

private:
    PenaltyCostCache cache_;

    /// 障碍物项解析梯度: d/dp_j [ Σ_i max(0, cl - d(c(t_i)))² ]
    void addObstacleGrad(const BezierCurve& c, const VecXd& params, double w, VecXd& grad) const;

    /// 边界+围栏+簇+人行横道项的中心差分数值梯度
    void addNumericGrad(const VecXd& params, double w_bnd, double w_fence, double w_cluster, double w_xwalk, VecXd& grad);

    /// 平滑项解析梯度（基于曲率导数）
    void addSmoothGrad(const BezierCurve& c, const VecXd& params, double w, VecXd& grad) const;
};

/// 用LBFGS求解器优化曲线
BezierCurve optimiseCurve(PenaltyCost& cost, LBFGSSolver& solver, const BezierCurve& initial, int outer_iters = 5);

}
