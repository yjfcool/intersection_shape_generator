#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "optimizer/lbfgs_solver.h"
#include "optimizer/penalty_cost.h"
#include "curve/curve_utils.h"
#include "curve/hermite_init.h"
#include "optimizer/sdf_field.h"

#include "io/iodata_shapefile.h"

using Catch::Matchers::WithinAbs;

using namespace isg;

// ──────────────────────────────────────────────────────────────
//  L-BFGS: minimise a simple quadratic  f = ||x - x*||^2
// ──────────────────────────────────────────────────────────────
TEST_CASE("L-BFGS converges on quadratic", "[lbfgs]") {
    VecXd target(4);
    target << 1.0, 2.0, -1.0, 3.0;

    LBFGSSolver solver;
    VecXd x0 = VecXd::Zero(4);

    auto res = solver.solve([&](const VecXd& x, VecXd& g) {
        g = 2.0*(x - target);
        return (x-target).squaredNorm();
    }, x0);

    REQUIRE(res.converged);
    REQUIRE((res.x - target).norm() < 1e-4);
}

TEST_CASE("L-BFGS warm start produces same minimum", "[lbfgs]") {
    VecXd target(2);
    target << 3.0, -2.0;

    LBFGSSolver solver;
    VecXd x0(2); x0 << 0,0;

    auto r1 = solver.solve([&](const VecXd& x, VecXd& g){
        g = 2*(x-target); return (x-target).squaredNorm();
    }, x0);

    auto r2 = solver.solveWarm([&](const VecXd& x, VecXd& g){
        g = 2*(x-target); return (x-target).squaredNorm();
    }, r1.x);

    REQUIRE((r2.x - r1.x).norm() < 1e-6);
}

// ──────────────────────────────────────────────────────────────
//  PenaltyCost: smooth term non-negative
// ──────────────────────────────────────────────────────────────
TEST_CASE("PenaltyCost smooth term >= 0", "[penalty]") {
    Vec2d p0(0,0), t0(1,0), p1(5,0), t1(1,0);

    SDFField sdf;
    BoundingBox2d roi; roi.min_pt={-5,-5}; roi.max_pt={10,10};
    sdf.build(roi, {}, 0.2, 0.0);

    BezierCurve init = buildInitialCurve(p0,t0,p1,t1,sdf,Polygon2d{});
    PenaltyCost cost;
    cost.proto = init;
    cost.sdf   = &sdf;
    cost.fence.outer = {{-10,-10},{10,-10},{10,10},{-10,10}};

    VecXd params = curveToParams(init);
    VecXd grad(params.size());
    double f = cost(params, grad);

    REQUIRE(f >= 0.0);
}

// ──────────────────────────────────────────────────────────────
//  Full optimisation: straight-line without obstacles
//  Curve should stay near the straight line and remain G1
// ──────────────────────────────────────────────────────────────
TEST_CASE("optimiseCurve straight no obstacles: endpoints preserved", "[optimizer]") {
    Vec2d p0(0,0), t0(1,0), p1(8,0), t1(1,0);

    SDFField sdf;
    BoundingBox2d roi; roi.min_pt={-2,-5}; roi.max_pt={12,5};
    sdf.build(roi, {}, 0.2, 0.0);

    Polygon2d fence;
    fence.outer = {{-2,-5},{12,-5},{12,5},{-2,5}};

    BezierCurve init = buildInitialCurve(p0,t0,p1,t1,sdf,fence);

    PenaltyCost cost;
    cost.proto = init;
    cost.sdf   = &sdf;
    cost.fence = fence;

    LBFGSSolver solver;
    BezierCurve result = optimiseCurve(cost, solver, init, 3);

    // Endpoints must match
    REQUIRE_THAT((result.startPt()-p0).norm(), WithinAbs(0.0, 0.1));
    REQUIRE_THAT((result.endPt()  -p1).norm(), WithinAbs(0.0, 0.1));

    // Start tangent must be close to t0
    double cos_a0 = result.startTan().dot(t0.normalized());
    REQUIRE(cos_a0 > 0.95);

    double cos_a1 = result.endTan().dot(t1.normalized());
    REQUIRE(cos_a1 > 0.95);
}

TEST_CASE("optimiseCurve enforces ordinary single cubic control-axis bounds",
          "[optimizer][constraints]") {
    Vec2d p0(0, 0), p1(10, 10), t0(1, 0), t1(0, 1);
    BezierCurve initial;
    initial.segs.push_back(makeCubicG1(p0, t0, p1, t1, 0.4));
    initial.segs.front().ctrl[1] = Vec2d(12, 2);
    initial.segs.front().ctrl[2] = Vec2d(8, -2);

    SDFField sdf;
    BoundingBox2d roi;
    roi.min_pt = {-5, -5};
    roi.max_pt = {15, 15};
    sdf.build(roi, {}, 0.2, 0.0);

    PenaltyCost cost;
    cost.sdf = &sdf;
    cost.start_tan_dir = t0;
    cost.end_tan_dir = t1;
    cost.constrain_single_cubic_axes = true;
    LBFGSSolver solver;
    BezierCurve result = optimiseCurve(cost, solver, initial, 1);

    CHECK(ordinarySingleCubicControlsValid(result, p0, t0, p1, t1));
    CHECK_THAT(result.segs.front().ctrl[1].x(), WithinAbs(10.0, 1e-6));
    CHECK_THAT(result.segs.front().ctrl[1].y(), WithinAbs(0.0, 1e-6));
    CHECK_THAT(result.segs.front().ctrl[2].x(), WithinAbs(10.0, 1e-6));
    CHECK_THAT(result.segs.front().ctrl[2].y(), WithinAbs(0.0, 1e-6));
}

// ──────────────────────────────────────────────────────────────
//  Optimisation with obstacle: final curve must not penetrate
// ──────────────────────────────────────────────────────────────
TEST_CASE("optimiseCurve avoids central obstacle", "[optimizer]") {
    // Entry: left, exit: right, obstacle in the middle (below straight line)
    Vec2d p0(0,0), t0(1,0), p1(8,0), t1(1,0);

    // Obstacle blocking the straight path
    Polygon2d obs_poly;
    obs_poly.outer = {{3,-1},{5,-1},{5,1},{3,1}};
    Obstacle obs; obs.geometry = obs_poly;

    BoundingBox2d roi; roi.min_pt={-2,-6}; roi.max_pt={12,6};
    SDFField sdf;
    sdf.build(roi, {obs}, 0.2, 0.3);

    Polygon2d fence;
    fence.outer = {{-2,-6},{12,-6},{12,6},{-2,6}};

    BezierCurve init = buildInitialCurve(p0,t0,p1,t1,sdf,fence);

    PenaltyCost cost;
    cost.proto = init;
    cost.sdf   = &sdf;
    cost.fence = fence;
    cost.obstacle_clearance = 0.1;
    cost.weights.obstacle = 20.0;

    LBFGSSolver solver;
    BezierCurve result = optimiseCurve(cost, solver, init, 5);

    // Check no sample point penetrates the raw obstacle
    for (auto& pt : result.sampleByArcLength(40)) {
        std::pair<double,Vec2d> _q = sdf.queryWithGrad(pt);
        // Allow small numerical tolerance
        REQUIRE(_q.first > -0.15);
    }
}

// ─────────────────────────────────────────────────────────────────────────────
//  Two-level bypass strategy tests
// ─────────────────────────────────────────────────────────────────────────────
#include "curve/hermite_init.h"
#include "constraints/fence_check.h"

static Polygon2d squareFence(Vec2d c,double h){
    Polygon2d p;
    p.outer={{c[0]-h,c[1]-h},{c[0]+h,c[1]-h},
             {c[0]+h,c[1]+h},{c[0]-h,c[1]+h}};
    return p;
}
static Obstacle makeObs(Vec2d c,double hw,double hh){
    Obstacle o;
    o.geometry.outer={{c[0]-hw,c[1]-hh},{c[0]+hw,c[1]-hh},
                      {c[0]+hw,c[1]+hh},{c[0]-hw,c[1]+hh}};
    return o;
}
static double angleDeg(const Vec2d& a,const Vec2d& b){
    double c=a.normalized().dot(b.normalized());
    return std::acos(std::max(-1.0,std::min(1.0,c)))*180.0/M_PI;
}

// ── Level-1: geometric bypass produces smooth arch, not staircase ─────────
TEST_CASE("Level-1: buildInitialCurve returns arch for single obstacle", "[bypass][L1]")
{
    Vec2d p0(0,0),t0(1,0),p1(8,0),t1(1,0);
    Obstacle obs=makeObs({4,0},1.0,1.0);
    BoundingBox2d roi;roi.min_pt={-2,-6};roi.max_pt={12,6};
    SDFField sdf;sdf.build(roi,{obs},0.2,0.3);
    Polygon2d fence=squareFence({4,0},8);

    BezierCurve init=buildInitialCurve(p0,t0,p1,t1,sdf,fence);
    REQUIRE(!init.empty());

    // Must have ≥2 segments (Level-1 returns 2-seg arch)
    REQUIRE(init.numSegments() >= 2);

    // Arch: all sample points on the same side of the direct line (no zigzag)
    auto pts=init.sampleByArcLength(30);
    Vec2d along=(p1-p0).normalized();
    Vec2d perp{-along[1],along[0]};
    // Find which side the arch goes
    double max_lat=0;
    for(auto& pt:pts) max_lat=std::max(max_lat,std::abs((pt-p0).dot(perp)));
    // All interior points should be on one consistent side (no staircase zigzag)
    double first_side=((pts[pts.size()/2]-p0).dot(perp)>0)?1:-1;
    int wrong=0;
    for(int i=2;i<(int)pts.size()-2;++i){
        double lat=(pts[i]-p0).dot(perp);
        if(lat*first_side<-0.05) wrong++;   // 5cm tolerance
    }
    REQUIRE(wrong == 0);  // all on same side → arch, not staircase
}

// ── Level-1: arch clears the obstacle ────────────────────────────────────
TEST_CASE("Level-1: initial arch has no obstacle penetration", "[bypass][L1][clear]")
{
    Vec2d p0(0,0),t0(1,0),p1(8,0),t1(1,0);
    Obstacle obs=makeObs({4,0},1.0,1.0);
    BoundingBox2d roi;roi.min_pt={-2,-6};roi.max_pt={12,6};
    SDFField sdf;sdf.build(roi,{obs},0.2,0.3);
    Polygon2d fence=squareFence({4,0},8);

    BezierCurve init=buildInitialCurve(p0,t0,p1,t1,sdf,fence);
    REQUIRE(!init.empty());

    SDFField sdf_raw;sdf_raw.build(roi,{obs},0.2,0.0);
    double min_d=1e18;
    for(auto& pt:init.sampleByArcLength(50)){
        std::pair<double,Vec2d> _q=sdf_raw.queryWithGrad(pt);
        min_d=std::min(min_d,_q.first);
    }
    REQUIRE(min_d > -0.15);  // ≤15cm penetration in initial (before optim)
}

// ── Level-1: endpoint tangents G1 in initial curve ────────────────────────
TEST_CASE("Level-1: initial arch preserves G1 endpoint tangents", "[bypass][L1][g1]")
{
    Vec2d p0(0,0),t0(1,0),p1(8,0),t1(1,0);
    Obstacle obs=makeObs({4,0},1.0,1.0);
    BoundingBox2d roi;roi.min_pt={-2,-6};roi.max_pt={12,6};
    SDFField sdf;sdf.build(roi,{obs},0.2,0.3);
    Polygon2d fence=squareFence({4,0},8);

    BezierCurve init=buildInitialCurve(p0,t0,p1,t1,sdf,fence);
    REQUIRE(!init.empty());

    REQUIRE_THAT((init.startPt()-p0).norm(), WithinAbs(0.0,0.01));
    REQUIRE_THAT((init.endPt()  -p1).norm(), WithinAbs(0.0,0.01));
    REQUIRE(angleDeg(init.startTan(),t0) < 5.0);
    REQUIRE(angleDeg(init.endTan(),  t1) < 5.0);
}

// ── Level-2: full-param optimisation improves clearance ──────────────────
TEST_CASE("Level-2: full-param optimiseCurve avoids obstacle", "[bypass][L2][clear]")
{
    Vec2d p0(0,0),t0(1,0),p1(8,0),t1(1,0);
    Obstacle obs=makeObs({4,0},1.0,1.0);
    BoundingBox2d roi;roi.min_pt={-2,-6};roi.max_pt={12,6};
    SDFField sdf;sdf.build(roi,{obs},0.2,0.3);
    Polygon2d fence=squareFence({4,0},8);

    BezierCurve init=buildInitialCurve(p0,t0,p1,t1,sdf,fence);

    PenaltyCost cost;
    cost.proto=init;cost.sdf=&sdf;cost.fence=fence;
    cost.obstacle_clearance=0.1;
    cost.start_tan_dir=t0.normalized();
    cost.end_tan_dir  =t1.normalized();
    cost.full_param_mode=(init.numSegments()>1);
    cost.buildCache();

    LBFGSSolver solver;
    BezierCurve result=optimiseCurve(cost,solver,init,4);
    REQUIRE(!result.empty());

    // Endpoints preserved
    REQUIRE_THAT((result.startPt()-p0).norm(),WithinAbs(0.0,0.15));
    REQUIRE_THAT((result.endPt()  -p1).norm(),WithinAbs(0.0,0.15));

    // No deep obstacle penetration
    SDFField sdf_raw;sdf_raw.build(roi,{obs},0.2,0.0);
    double min_d=1e18;
    for(auto& pt:result.sampleByArcLength(60)){
        std::pair<double,Vec2d> _q=sdf_raw.queryWithGrad(pt);
        min_d=std::min(min_d,_q.first);
    }
    REQUIRE(min_d > -0.25);
}

// ── Level-2: G1 at endpoints after full-param optimisation ────────────────
TEST_CASE("Level-2: G1 endpoint tangents preserved after optimisation", "[bypass][L2][g1]")
{
    Vec2d p0(0,0),t0(1,0),p1(8,0),t1(1,0);
    Obstacle obs=makeObs({4,0},1.0,1.0);
    BoundingBox2d roi;roi.min_pt={-2,-6};roi.max_pt={12,6};
    SDFField sdf;sdf.build(roi,{obs},0.2,0.3);
    Polygon2d fence=squareFence({4,0},8);

    BezierCurve init=buildInitialCurve(p0,t0,p1,t1,sdf,fence);

    PenaltyCost cost;
    cost.proto=init;cost.sdf=&sdf;cost.fence=fence;
    cost.obstacle_clearance=0.1;
    cost.start_tan_dir=t0.normalized();
    cost.end_tan_dir  =t1.normalized();
    cost.full_param_mode=(init.numSegments()>1);
    cost.buildCache();

    LBFGSSolver solver;
    BezierCurve result=optimiseCurve(cost,solver,init,4);

    // Endpoint tangent angle < 5 degrees from required direction
    REQUIRE(angleDeg(result.startTan(),t0) < 5.0);
    REQUIRE(angleDeg(result.endTan(),  t1) < 5.0);
}

// ── Arch smoothness: no direction reversal (no U-turns in the arch) ───────
TEST_CASE("Level-1: arch has monotone longitudinal progress (no backtracking)", "[bypass][L1][mono]")
{
    Vec2d p0(0,0),t0(1,0),p1(8,0),t1(1,0);
    Obstacle obs=makeObs({4,0},1.0,1.0);
    BoundingBox2d roi;roi.min_pt={-2,-6};roi.max_pt={12,6};
    SDFField sdf;sdf.build(roi,{obs},0.2,0.3);
    Polygon2d fence=squareFence({4,0},8);

    BezierCurve init=buildInitialCurve(p0,t0,p1,t1,sdf,fence);
    REQUIRE(!init.empty());

    Vec2d along=(p1-p0).normalized();
    auto pts=init.sampleByArcLength(40);
    int backward=0;
    for(int i=1;i<(int)pts.size();++i){
        Vec2d step=pts[i]-pts[i-1];
        if(step.dot(along)<-0.05) backward++;  // allow tiny numerical noise
    }
    REQUIRE(backward==0);  // no backtracking along travel direction
}

// ── Two obstacles: Level-2 fallback handles complex case ──────────────────
TEST_CASE("Level-2: two obstacles — arch initialisation + full-param optimisation", "[bypass][L2][two_obs]")
{
    Vec2d p0(0,0),t0(1,0),p1(10,0),t1(1,0);
    Obstacle obs1=makeObs({3,0},0.8,0.8);
    Obstacle obs2=makeObs({7,0},0.8,0.8);
    BoundingBox2d roi;roi.min_pt={-2,-6};roi.max_pt={14,6};
    SDFField sdf;sdf.build(roi,{obs1,obs2},0.2,0.3);
    Polygon2d fence=squareFence({5,0},9);

    BezierCurve init=buildInitialCurve(p0,t0,p1,t1,sdf,fence);
    REQUIRE(!init.empty());

    PenaltyCost cost;
    cost.proto=init;cost.sdf=&sdf;cost.fence=fence;
    cost.obstacle_clearance=0.1;
    cost.start_tan_dir=t0.normalized();
    cost.end_tan_dir  =t1.normalized();
    cost.full_param_mode=(init.numSegments()>1);
    cost.buildCache();

    LBFGSSolver solver;
    BezierCurve result=optimiseCurve(cost,solver,init,4);
    REQUIRE(!result.empty());

    // Endpoints must be close to required points
    REQUIRE_THAT((result.startPt()-p0).norm(),WithinAbs(0.0,0.2));
    REQUIRE_THAT((result.endPt()  -p1).norm(),WithinAbs(0.0,0.2));
}

// ── No obstacle: Level-0 straight returns single segment ──────────────────
TEST_CASE("Level-0: straight path returns single-segment Bezier", "[bypass][L0]")
{
    Vec2d p0(0,0),t0(1,0),p1(8,0),t1(1,0);
    BoundingBox2d roi;roi.min_pt={-2,-6};roi.max_pt={12,6};
    SDFField sdf;sdf.build(roi,{},0.2,0.0);   // no obstacles
    Polygon2d fence=squareFence({4,0},8);

    BezierCurve init=buildInitialCurve(p0,t0,p1,t1,sdf,fence);
    REQUIRE(!init.empty());
    REQUIRE(init.numSegments()==1);   // straight line → single segment

    REQUIRE_THAT((init.startPt()-p0).norm(),WithinAbs(0.0,0.01));
    REQUIRE_THAT((init.endPt()  -p1).norm(),WithinAbs(0.0,0.01));
    REQUIRE(angleDeg(init.startTan(),t0)<1.0);
    REQUIRE(angleDeg(init.endTan(),  t1)<1.0);
}

// ── curveToParamsFull / curveFromParamsFull round-trip ────────────────────
TEST_CASE("Level-2: full-param encode/decode preserves control points", "[bypass][L2][params]")
{
    std::vector<Vec2d> pts={Vec2d(0,0),Vec2d(3,2),Vec2d(6,2),Vec2d(9,0)};
    std::vector<Vec2d> tans={Vec2d(1,0),Vec2d(1,0),Vec2d(1,0),Vec2d(1,0)};
    BezierCurve c=makeCurveFromKnots(pts,tans);
    REQUIRE(c.numSegments()==3);

    VecXd params=curveToParamsFull(c);
    // Expected size: 4 + 6*(3-1) = 16
    REQUIRE(params.size()==16);

    BezierCurve c2=curveFromParamsFull(params,c);
    REQUIRE(c2.numSegments()==3);

    // Endpoints fixed
    REQUIRE_THAT((c2.startPt()-c.startPt()).norm(),WithinAbs(0.0,1e-6));
    REQUIRE_THAT((c2.endPt()  -c.endPt()  ).norm(),WithinAbs(0.0,1e-6));

    // G1 at join points: tangent directions must be continuous
    for(int k=0;k<2;++k){
        Vec2d t_end   =(c2.segs[k].ctrl[3]-c2.segs[k].ctrl[2]).normalized();
        Vec2d t_start =(c2.segs[k+1].ctrl[1]-c2.segs[k+1].ctrl[0]).normalized();
        REQUIRE(angleDeg(t_end,t_start)<1.0);
    }
}
