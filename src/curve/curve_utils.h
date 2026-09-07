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
    /// 弦方向联合预算：`start_chord_proj * λ + end_chord_proj * μ <= chord_len`
    /// （λ/μ 为首/尾把手长度）。`has_chord_budget=false` 表示该几何不适用
    /// （弦退化，或某侧方向轴不朝弦前方，例如掉头）。
    bool has_chord_budget = false;
    double chord_len = 0.0;
    double start_chord_proj = 0.0;
    double end_chord_proj = 0.0;
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

/// 单段三次 Bezier 的控制多边形是否沿自身弦单调，即 `proj(P1) <= proj(P2)`
/// （投影取到弦方向上）。
///
/// 这是需求 1.3.1「不得越过端点弦长的有效范围」的**联合**读法。原有实现只给
/// 两个把手各自的上界（各自可吃满整条弦），两侧同时吃满时控制多边形沿弦倒序，
/// 曲线被折成 Z 形或 L 形：110000385-u 的直行 55 弦长 53.96m、两个把手各
/// 43.17m，`proj(P1)=41.97 > proj(P2)=12.00`，弦向净退 29.96m；两段避让曲线
/// 21 的两段同样各超支 7.91m 和 3.11m。它们都逐个满足旧的单侧上界。
///
/// 展开成把手长度：`λ·(T0·u) + μ·(T1·u) <= |chord|`（u 为弦方向单位向量），
/// 即两个把手在弦方向上的投影之和不得超过弦长——一个天然的联合预算。
///
/// 它与方向交点上限相容而非冲突：设方向交点距两端点为 s、e，则恒有
/// `s·(T0·u) + e·(T1·u) = |chord|`（把交点分别沿两条射线投影到弦上即得）。
/// 于是 `λ<=s 且 μ<=e` 必然满足本预算，交点上限生效时本预算自动成立；本预算
/// 只在直行、无前向交点、或调用方显式关闭交点上限（密集 Boundary 链需要沿轴
/// 延长把手）时才起作用——正是那些原本完全没有联合约束的场景。同时它比单侧
/// 交点上限更宽松也更准确：一侧把手越过交点、另一侧相应缩短仍可保持单调，
/// 避让需要的轴向延长因此不被误杀。
///
/// 某侧方向轴不朝弦前方时（掉头及极端反向几何）返回 true——该场景下弦向单调
/// 没有意义，由掉头专用约束负责。
///
/// 判定带弦长比例松量（实现里的 kChordBudgetSlackFraction = 5%）：等号形态就是
/// 把手正好落在方向交点，零容差会把只越界百分之几的候选一并淘汰、连带改变候选
/// 集（实测使 100000443 的连接 5 退到粗围栏外）。真正折形的越界量是弦长的
/// 16%~56%，与 5% 松量相差一个量级。
bool cubicControlPolygonMonotone(const BezierSegment& segment, double tol = 1e-6);

/// 曲线折形的**幅度**：各段控制多边形沿自身弦的倒退量占该段弦长的比例，取最大值；
/// `0` 表示全部段都沿弦单调（含预算松量内的临界形态）。
///
/// `cubicControlPolygonMonotone` 只回答"是否折形"，物理修复搜索需要的是"折得多狠"：
/// Boundary / Obstacle / Fence 修复的候选里，折形候选有时是唯一物理安全的选择
/// （100000385-u 的左转 21 只有折形两段候选能绕开 Boundary），一律淘汰会把形态问题
/// 换成更严重的物理违规。因此修复搜索按本函数的返回值给折形候选加打分惩罚，
/// 让"能不折就不折"，而不是"一折就否"。
///
/// 比例口径与 kChordProjectionBudgetSlack 一致：返回值已扣除松量，因此
/// `overshoot > 0` 与形态闸门判定的"折形"完全同义。段的某侧把手方向不朝弦前方时
/// （掉头轴向段等）该段不计入。
double curveChordBudgetOvershoot(const BezierCurve& curve);

/// 将普通单段三次Bezier的两个内部控制点限制在端点切向轴的有效范围内。
/// 除单侧上下界外还施加弦方向联合预算（见 cubicControlPolygonMonotone）：
/// 超支时在最小把手之上等比缩放两侧，保持两侧把手比例即拱顶位置不变。
/// U型调头和物理避让候选由调用方显式跳过本约束。
void constrainOrdinarySingleCubicControls(
    BezierCurve& curve, const Vec2d& p0, const Vec2d& start_tan,
    const Vec2d& p1, const Vec2d& end_tan,
    bool cap_at_direction_intersection = true);

/// 检查普通单段三次Bezier是否满足端点方向轴、有效范围和弦方向联合预算约束。
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

/// 普通转弯曲线 arc/chord 比例的下限（同侧单拱形态的"够弧"口径）。
///
/// `turn_angle_rad` 是进入切向与退出切向之间的夹角，`chord_len` 是端点弦长。
///
/// 需求 1.4② 把下限写成常数："长距离普通转弯 arc/chord 原则上不低于 1.08，
/// 短距离不低于 1.03"。但 arc/chord 的下界由转角唯一决定、与"长短"无关：给定
/// 弦长 c 和首尾切向夹角 θ，弧长最短的 G1 单拱就是圆弧，此时
/// arc/chord = θ / (2 sin(θ/2))。θ=90° 时该值是 1.1107，1.08 恰好是它的
/// 97%，1.03 恰好是它的 93%——常数阈值真正表达的是"不低于圆弧的 97%/93%"，
/// 只是被按 90° 折算成了数字。
///
/// 一旦转角变浅，常数就变成不可达的：θ=52° 的圆弧只有 1.036，θ=66° 只有
/// 1.058，任何满足 1.08 的曲线都必须比圆弧更长，也就是必须鼓包或走 S 形，
/// 与同一条需求"不能被压成 S 弯"直接矛盾。100000385-u 的 69/70/71/83/84/85/86
/// 就是这样被误判的：它们的转角只有 52°~67°，实测比例 1.037~1.049 已经达到
/// 甚至超过同转角圆弧（83/84/85/86 分别为圆弧的 100.1%~100.7%），却因为
/// 弦长大于 12m 而被拿去和 1.06 比较。
///
/// 因此这里以"按转角折算的圆弧基准比例"为准：长弦取圆弧的 97%、短弦（<12m，
/// 起停段占比大、几何余量小）取 93%，θ=90° 时与需求写明的 1.08 / 1.03 完全
/// 一致。同时把需求原有的业务常数（长弦 1.06、短弦 1.02）保留为**上界**：转角
/// 足够大时仍按常数裁定，只有圆弧本身都达不到常数的浅转才让位给圆弧基准。
/// 于是这条口径只做一件事——把不可达的阈值降到可达，不放松任何原本可达的场景。
/// 另有 1.005 的绝对地板，使近乎直行的浅转不会退化成"任意比例合法"。
double ordinaryTurnArcChordFloor(double turn_angle_rad, double chord_len);

/// 形态恢复目标的 arc/chord 下限，口径同 ordinaryTurnArcChordFloor，但业务常数
/// 上界取需求里更严的一档（长弦 1.08、短弦 1.03）。它用于"这条曲线是否值得再
/// 试一次纯几何恢复"的判定：目标不可达时恢复永远不会被接受，既白花修复预算，
/// 又会把已经合规的浅转当成待修对象反复改形。
double ordinaryTurnArcChordRestoreFloor(double turn_angle_rad,
                                        double chord_len);

/// 曲线首尾切向之间的夹角(弧度, [0, π])。用于把 arc/chord 下限按转角折算。
/// 端点 G1 由独立约束保证，因此对合规曲线它与车道端点切向夹角一致。
double curveEndpointTurnAngle(const BezierCurve& curve);

}
