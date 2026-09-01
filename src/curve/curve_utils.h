#pragma once
#include "bezier.h"

namespace isg {

struct SDFField;

/// 普通单段三次 Bezier 的端点方向轴及把手有效范围。
/// direction_intersection_* 仅在两条前向方向射线存在稳定交点时有效。
struct OrdinarySingleCubicHandleBounds {
    Vec2d start_dir{1, 0};
    Vec2d end_dir{1, 0};
    double start_min = 0.05;
    double end_min = 0.05;
    double start_max = 0.05;
    double end_max = 0.05;
    bool has_direction_intersection = false;
    double direction_intersection_start = 0.0;
    double direction_intersection_end = 0.0;
};

/// 三点局部曲率（基于外接圆半径倒数）
double localCurvature(const Vec2d& a, const Vec2d& b, const Vec2d& c);

/// 三点外接圆圆心
Vec2d circumcenter(const Vec2d& a, const Vec2d& b, const Vec2d& c);

/// 弹性带平滑：将点列视为弹性带，在曲率超限处沿外接圆径向移动
/// 顶点以降低曲率，受SDF和围栏约束
std::vector<Vec2d> elasticBandSmooth(const std::vector<Vec2d>&, const SDFField&,
                                     const Polygon2d&, double kappa_max = 0.25, double move_step = 0.03,
                                     int max_iter = 50, double min_sdf = 0.1);

/// 自适应版本弹性带平滑：根据曲线长度与曲率变化自适应调整采样数
std::vector<Vec2d> elasticBandSmoothAdaptive(const std::vector<Vec2d>&, const SDFField&,
                                     const Polygon2d&, double kappa_max = 0.25, double move_step = 0.03,
                                     int max_iter = 50, double min_sdf = 0.1);

/// 抽稀点列至目标数量
std::vector<Vec2d> downsamplePoints(const std::vector<Vec2d>&, int target_count);

/// 由平滑后点列重建Bezier曲线
BezierCurve rebuildFromSmoothedPts(const std::vector<Vec2d>&, const Vec2d&, const Vec2d&);

/// 沿两条曲线弧长方向均匀采样中点
std::vector<Vec2d> midlineSampleByArcLength(const BezierCurve&, const BezierCurve&, int n = 20);

/// 点到直线的有符号距离
double signedDistToLine(const Vec2d& pt, const Vec2d& p0, const Vec2d& p1);

/// 两线段是否相交，out返回交点
bool segmentsIntersect(const Vec2d& a, const Vec2d& b,
                       const Vec2d& c, const Vec2d& d, Vec2d* out = nullptr);

/// 曲线沿采样点的最小SDF值
double minSDFAlongCurve(const BezierCurve&, const SDFField&, int sps = 20);

/// 交点到两条曲线所有端点的最短距离
double distToAllEndpoints(const Vec2d&, const BezierCurve&, const BezierCurve&);

/// 返回两条曲线的所有交点（容差tol）
std::vector<Vec2d> curveCrossings(const BezierCurve&, const BezierCurve&, double tol = 0.01);

/// 两条曲线包围盒是否重叠
bool bboxOverlap(const BezierCurve&, const BezierCurve&);

/// 业务相交判定：两条曲线是否存在非端点附近相交（ep=端点排除半径，米）
bool curvesIntersectBusiness(const BezierCurve&, const BezierCurve&, double ep = 0.01);

/// 曲线自身相交判定：排除首尾相接的情况
bool curveSelfIntersectsBusiness(const BezierCurve&, double ep = 0.01);

/// 计算普通单段曲线的端点方向和把手上下界。
OrdinarySingleCubicHandleBounds ordinarySingleCubicHandleBounds(
    const Vec2d& p0, const Vec2d& start_tan,
    const Vec2d& p1, const Vec2d& end_tan,
    bool cap_at_direction_intersection = true);

/// 将普通单段三次Bezier的两个内部控制点限制在端点切向轴的有效范围内。
/// U型调头和物理避让候选由调用方显式跳过本约束。
void constrainOrdinarySingleCubicControls(
    BezierCurve& curve, const Vec2d& p0, const Vec2d& start_tan,
    const Vec2d& p1, const Vec2d& end_tan,
    bool cap_at_direction_intersection = true);

/// 检查普通单段三次Bezier是否满足端点方向轴和有效范围约束。
bool ordinarySingleCubicControlsValid(
    const BezierCurve& curve, const Vec2d& p0, const Vec2d& start_tan,
    const Vec2d& p1, const Vec2d& end_tan, double tol = 1e-6,
    bool cap_at_direction_intersection = true);

/// 共享进入或退出端点的两条单段三次Bezier，其中间两个控制点的连线
/// (P1→P2) 是否相交。
///
/// 同入/同出扇出族的控制多边形必须保持嵌套。两条曲线共享一个端点、
/// 在该端点切向相同时，中段形态完全由两个中间控制点决定：一旦两条
/// P1→P2 连线相交，控制多边形互穿，曲线必然在非端点处交叉。因此该
/// 判定是"同簇非端点不相交"在单段贝塞尔表达下的几何前兆，可直接用于
/// 候选剔除，比采样相交判定更早、更稳定。
///
/// 任一曲线不是单段时返回 false（多段表达不适用该充要前兆）。
bool sharedEndpointMidControlSegmentsCross(
    const BezierCurve& a, const BezierCurve& b);

/// 共享进入或退出端点的两条单段三次Bezier，其控制多边形折线
/// (P0→P1→P2→P3) 是否存在"非重合、非控制点"的交叉。
///
/// `sharedEndpointMidControlSegmentsCross` 只比较两条 P1→P2 连线，
/// 无法覆盖一条曲线的中间连线穿过另一条曲线首/尾把手连线的情形
/// （例如 A 的 P1→P2 与 B 的 P2→P3 相交）。控制多边形是曲线的凸包骨架，
/// 任意一处互穿都意味着两族控制多边形不再嵌套，曲线在非端点处交叉。
/// 因此本判定是同簇非端点不相交的完整几何前兆。
///
/// 排除两类合法接触：
/// - 重合：两段共线（同入侧首把手沿同一进入切向天然共线重叠）；
/// - 控制点处相接：交点落在任一控制点的 `ctrl_tol` 邻域内（共享端点
///   本身允许在端点容差内重合）。
///
/// 任一曲线不是单段时返回 false。
bool sharedEndpointControlPolylinesCross(
    const BezierCurve& a, const BezierCurve& b, double ctrl_tol = 0.10);

/// 单段三次Bezier在首端点(`at_start=true`)或尾端点处的有符号曲率，即
/// 共享端点扇出族的"转入率"。
///
/// 共享同一端点与同一端点切向的一族单段三次Bezier，在该端点附近的横向
/// 偏离满足 `lat(s) ≈ kappa0 * s^2 / 2`（s 为自端点起算的弧长），因此这一族
/// 在端点邻域的左右次序**完全**由该有符号曲率决定：值越大越偏向切向左侧。
///
/// 展开后 `kappa0 = (2/3) * (cross(T0, chord) - h1 * cross(T0, T1)) / h0^2`，
/// 其中 h0/h1 为首/尾把手长度。可见首把手 h0 单调并不能保证转入率单调：
/// 尾把手 h1 独立选取时同样改变次序。这正是"逐对比较共享侧把手投影"无法
/// 表达扇出总序的原因（110003285 的左转 80/81：h0 16.67 < 17.57 已单调，
/// h1 23.34 > 17.57 却把转入率压成 0.0325 < 0.0428，次序反转并在距端点
/// 18.4m 处真实相交）。
///
/// 曲线为空或对应把手退化（长度 < 1e-9）时返回 0。多段曲线只取对应端点
/// 所在的那一段，因为端点邻域的次序只由该段决定。
double singleCubicSignedEndCurvature(const BezierCurve& c, bool at_start);

/// 共享端点"汇入/分流漏斗"半径：两条曲线在共享端点附近彼此始终贴近的那一段
/// 的作用半径，`0` 表示两条曲线没有共享端点或一离开端点就分开。
///
/// 两条连接汇入同一条车道（共享出口端点）或自同一条车道分流（共享入口端点）时，
/// 它们在端点附近必然收敛到同一点同一切向。收敛段内的贴近与至多一次穿越是
/// 车道拓扑本身决定的，不是任何一条曲线的形态缺陷，任何重新造型都消不掉它——
/// 例如 110000703-u 的掉头 41 必须在跨过人行横道（9.32m）之后才起拱，而汇入
/// 同一条出口车道的近直行 13 在这 9.3m 里与 41 的入口直段横向只差 0.005~0.31m，
/// 并且 13 的终点就落在 41 走廊内部，因此"贴行 + 一次穿越"在几何上无法避免。
/// 固定的 0.15/1.5m 端点容差覆盖不到这么长的收敛段，会把结构性汇入误判成违约，
/// 逼迫候选搜索退回单段拱形而破坏掉头的人行横道分段约束。
///
/// 半径定义为：自共享端点起算、两条曲线到对方的距离**始终**小于 `band` 的那一段
/// 弧长（两条曲线各自计算后取小），再截断到 `cap`。一旦间距超过 `band`，收敛段
/// 即告结束，之后的冲突照常按原容差判定，因此本半径不会豁免真正远离端点的冲突。
///
/// 端点是否共享按 `endpoint_tol` 的坐标重合判定；两条曲线同时共享两个端点时
/// 取较大的半径。
double sharedEndpointMergeFunnelRadius(const BezierCurve& a,
                                       const BezierCurve& b,
                                       double band = 0.50,
                                       double cap = 12.0,
                                       double endpoint_tol = 0.30);

/// 两条曲线坐标重合的端点集合（`endpoint_tol` 内视为同一点）。
std::vector<Vec2d> sharedEndpointsOf(const BezierCurve& a, const BezierCurve& b,
                                     double endpoint_tol = 0.30);

/// 与 `curvesIntersectBusiness` 完全相同的判定，但额外豁免落在 `centers` 任一点
/// `radius` 邻域内的交点。用于结构性汇入/分流：只放行收敛段内的穿越，远离共享
/// 端点的真实穿越仍然禁止，因此比直接放大 `ep` 精确得多——放大 `ep` 会把四个
/// 端点的邻域一起放大。
bool curvesIntersectBusinessOutsideBalls(const BezierCurve& a,
                                         const BezierCurve& b, double ep,
                                         const std::vector<Vec2d>& centers,
                                         double radius);

/// 曲线切向的"单向绕转跨度"（弧度）：沿曲线累加**有符号**切向转角得到running
/// 累加曲线 `S(s)`，返回 `max S - min S`，即任意子区间内同一方向连续绕过的
/// 最大角度。
///
/// 为什么用它判"折返/绕圈"，而不用「总转角 - 首尾切向夹角」的超量：
///   * 换道弧是合法的 S 形（先左后右），总转角是两段弯的和、首尾切向却几乎
///     不变，超量因此等于 `2·min(两段弯)`——110003285 的换道弧 14 实测超量
///     120.4°，与真正的病态折返同量级，用超量设阈值必然误伤。S 形的绕转跨度
///     只等于其中较大的**单侧**弯角，仍然很小。
///   * 病态曲线是同一方向绕满一整圈（100000385-u 的 14 修复候选实测总转角
///     448°、首尾切向差只有 88°，局部曲率半径 0.5m），绕转跨度就是 448°，
///     与正常转向（≤ ~150°）拉开整整一个量级。
///
/// 它同时补上了曲率符号翻转判据的盲区：后者每段只取有限个内部采样点、并用
/// 死区过滤小抖动，因此绕圈过程中曲率始终同号时完全看不到。
///
/// 曲线为空、`samples_per_seg < 1`、或有效切向不足两个时返回 0。
double curveTurningSpan(const BezierCurve& curve, int samples_per_seg = 48);

}
