// 普通曲线单段性专项诊断：找出“非掉头曲线最终不是单段 cubic”的连接，并给出
// “把它换成首选单段 cubic 会违反哪些约束”的逐条判据。
//
// 用途：`shape.ordinary.single_segment` 报违约时，需要区分两种成因——
//   1) 单段 cubic 本身不可行（穿障 / 越 Boundary / 出围栏 / 同簇相交），多段是
//      修复的代价；
//   2) 单段 cubic 完全可行，多段只是搜索顺序或分支条件的副产物。
// 只有第 2 类才应该在生成侧收敛回单段。
//
// 用法：
//   diag_single_cubic <数据文件> [--id 连接id]... [--all]
//     --id   只看指定连接（可重复）；缺省看全部非单段的普通曲线
//     --all  连同单段的普通曲线一起列出（用于确认基线未被改坏）
#include "constraints/boundary_safety.h"
#include "constraints/cluster_order.h"
#include "constraints/constraint_evaluator.h"
#include "curve/curve_utils.h"
#include "domain/scene_context.h"
#include "generation/connectivity_generation_context.h"
#include "initialization/ordinary_curve_initializer.h"
#include "intersection_shape_generator.h"
#include "io/iodata_json.h"
#include "optimizer/sdf_field.h"
#include "preprocessing/uturn_family_builder.h"
#include "toolkits/toolkits.h"
#include "utils.h"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using namespace isg;

namespace {

bool isGeometricUTurnLocal(const std::pair<Vec2d, Vec2d>& entry,
                           const std::pair<Vec2d, Vec2d>& exit) {
    return entry.second.norm() > 1e-8 && exit.second.norm() > 1e-8 &&
           entry.second.normalized().dot(exit.second.normalized()) < -0.5;
}

/// 与 diag_all_violations 的 auditProfile 同口径，另外打开同簇检查，因为这里要
/// 回答“单段 cubic 是否会引入新的同簇相交”。
ConstraintProfile diagProfile(int mode) {
    ConstraintProfile p;
    p.enforce_fence = true;
    p.check_self_intersection = true;
    p.check_obstacle = true;
    p.check_boundary = true;
    p.check_g1 = true;
    p.check_curvature = false;
    p.check_ordinary_shape = true;
    p.check_uturn_shape = true;
    p.check_crosswalk = true;
    p.check_cluster = true;
    p.obstacle_clearance = 0.0;
    p.road_edge_clearance = mode == 2 ? 1.0 : 0.0;
    p.cluster_endpoint_tolerance = 0.30;
    p.samples = 32;
    return p;
}

/// 与 IntersectionShapeGenerator::generate 同源的 SDF 取景框。
BoundingBox2d outerBoundingBox(const IntersectionInput& input) {
    BoundingBox2d roi;
    if (!input.area.geometry.outer.empty()) {
        roi = input.area.geometry.bbox();
    } else {
        for (size_t i = 0; i < input.lanes.size(); ++i)
            for (size_t j = 0; j < input.lanes[i].geometry.points.size(); ++j)
                roi.expand(input.lanes[i].geometry.points[j]);
        roi.min_pt -= Vec2d(20, 20);
        roi.max_pt += Vec2d(20, 20);
    }
    return roi;
}

std::string violatedIds(const ConstraintReport& report) {
    std::string out;
    for (const auto& r : report.results) {
        if (r.state != ConstraintState::Violated) continue;
        if (!out.empty()) out += " ";
        out += r.id;
    }
    return out.empty() ? std::string("-") : out;
}

const char* boundaryTypeName(Boundary::Type t) {
    switch (t) {
        case Boundary::Type::RoadEdge: return "RoadEdge";
        case Boundary::Type::MedianStrip: return "MedianStrip";
        case Boundary::Type::GreenBelt: return "GreenBelt";
        default: return "Other";
    }
}

/// 点到折线段的最短距离。
double distPointToSegment(const Vec2d& p, const Vec2d& a, const Vec2d& b) {
    const Vec2d ab = b - a;
    const double L2 = ab.squaredNorm();
    if (L2 < 1e-18) return (p - a).norm();
    double t = (p - a).dot(ab) / L2;
    t = t < 0.0 ? 0.0 : (t > 1.0 ? 1.0 : t);
    return (p - (a + ab * t)).norm();
}

/// 从 boundary 的某个端点（tip_at_front=true 表示首点）出发，沿折线走出“与车道切向
/// 基本共线”的连续段区间 [lo, hi]（段序号闭区间）。这是真正贴着车道走的那一截，
/// 而不是整条折线——发夹形路缘的另一条腿必须留在检查范围内。
void alignedRun(const std::vector<Vec2d>& bpts, bool tip_at_front, const Vec2d& tangent,
                double min_dot, size_t* lo, size_t* hi) {
    const size_t nseg = bpts.size() - 1;
    auto aligned = [&](size_t s) {
        const Vec2d d = bpts[s + 1] - bpts[s];
        return d.norm() > 1e-9 && std::abs(d.normalized().dot(tangent)) >= min_dot;
    };
    if (tip_at_front) {
        *lo = 0;
        size_t s = 0;
        while (s < nseg && aligned(s)) ++s;
        *hi = s == 0 ? 0 : s - 1;
    } else {
        *hi = nseg - 1;
        size_t s = nseg;
        while (s > 0 && aligned(s - 1)) --s;
        *lo = s >= nseg ? nseg - 1 : s;
    }
}

/// 逐条 Boundary 复算，定位到底是哪一条、以哪种方式（穿越 / 越到路缘外侧）判违。
void reportBoundaryDetail(const char* tag, const BezierCurve& curve,
                          const IntersectionInput& input,
                          const std::pair<Vec2d, Vec2d>& entry,
                          const std::pair<Vec2d, Vec2d>& exit) {
    if (curve.empty() || input.boundaries.empty()) return;
    const Vec2d center = boundarySafetyCenter(input);
    const Vec2d T0 = entry.second.norm() > 1e-8 ? entry.second.normalized() : Vec2d(1, 0);
    const Vec2d T1 = exit.second.norm() > 1e-8 ? exit.second.normalized() : Vec2d(1, 0);
    const std::vector<Vec2d> samples = curve.sampleByArcLength(
        std::max(64, std::min(240, (int)std::ceil(curve.arcLength() / 0.18) + 1)));
    for (const auto& bnd : input.boundaries) {
        const std::vector<Boundary> one{bnd};
        const BoundarySafetyResult bs =
            curveBoundarySafety(curve, one, center, 64, 0.75, 0.10, 0.05);
        if (!bs.intersects && !bs.outside_road_edge) continue;
        std::printf("      %s boundary id=%s type=%s intersects=%d outside_road_edge=%d "
                    "penalty=%.5f pts=%zu\n",
                    tag, bnd.id.c_str(), boundaryTypeName(bnd.type),
                    (int)bs.intersects, (int)bs.outside_road_edge, bs.outside_penalty,
                    bnd.geometry.points.size());
        const std::vector<Vec2d> bpts = toVec2dArray(bnd.geometry.points);
        if (bpts.size() < 2) continue;
        const Vec2d chord_dir = (bpts.back() - bpts.front()).norm() > 1e-8
                                    ? (bpts.back() - bpts.front()).normalized()
                                    : Vec2d(1, 0);
        // 整条 boundary 的首尾弦方向（现有 filterEndpointAdherentBoundaries 用的量）
        std::printf("            chord_dir=(%.3f,%.3f) |dot T0|=%.3f |dot T1|=%.3f "
                    "d(first,p0)=%.2f d(last,p0)=%.2f d(first,p1)=%.2f d(last,p1)=%.2f\n",
                    chord_dir.x(), chord_dir.y(), std::abs(chord_dir.dot(T0)),
                    std::abs(chord_dir.dot(T1)),
                    (bpts.front() - entry.first).norm(), (bpts.back() - entry.first).norm(),
                    (bpts.front() - exit.first).norm(), (bpts.back() - exit.first).norm());
        // 终端腿方向：真正与车道方向对齐的那一段
        const Vec2d first_leg = (bpts[1] - bpts[0]).norm() > 1e-8
                                    ? (bpts[1] - bpts[0]).normalized() : Vec2d(1, 0);
        const Vec2d last_leg =
            (bpts[bpts.size() - 1] - bpts[bpts.size() - 2]).norm() > 1e-8
                ? (bpts[bpts.size() - 1] - bpts[bpts.size() - 2]).normalized() : Vec2d(1, 0);
        std::printf("            first_leg=(%.3f,%.3f) |dot T0|=%.3f |dot T1|=%.3f  "
                    "last_leg=(%.3f,%.3f) |dot T0|=%.3f |dot T1|=%.3f\n",
                    first_leg.x(), first_leg.y(), std::abs(first_leg.dot(T0)),
                    std::abs(first_leg.dot(T1)), last_leg.x(), last_leg.y(),
                    std::abs(last_leg.dot(T0)), std::abs(last_leg.dot(T1)));
        // 逐个真实交点：与曲线端点的弧长/直线距离、命中的 boundary 段序号
        for (size_t i = 0; i + 1 < samples.size(); ++i)
            for (size_t j = 0; j + 1 < bpts.size(); ++j) {
                Vec2d isect;
                if (!segmentsIntersect(samples[i], samples[i + 1], bpts[j], bpts[j + 1],
                                       &isect))
                    continue;
                std::printf("            cross @ (%.2f,%.2f) bseg=%zu/%zu "
                            "d(start)=%.2f d(end)=%.2f\n",
                            isect.x(), isect.y(), j, bpts.size() - 2,
                            (isect - samples.front()).norm(),
                            (isect - samples.back()).norm());
            }
        // 贴合端点判据：boundary 的某个端点与曲线端点重合时，从该端点走出的“共线走向段”
        // 上的相交只是共享端点造成的数值擦碰，还是真的切进了路缘？用曲线相对该走向段
        // 的最大横向偏移来量化。
        for (int side = 0; side < 2; ++side) {
            const Vec2d cp = side == 0 ? entry.first : exit.first;
            const Vec2d ct = side == 0 ? T0 : T1;
            const bool tip_front = (bpts.front() - cp).norm() <= 0.60;
            const bool tip_back = (bpts.back() - cp).norm() <= 0.60;
            if (!tip_front && !tip_back) continue;
            size_t lo = 0, hi = 0;
            alignedRun(bpts, tip_front, ct, 0.70, &lo, &hi);
            // 曲线上从贴合端点出发、直到最远那个落在 run 内的交点为止的那一段，
            // 相对 run 折线的最大横向偏移 = 擦碰深度。
            double graze = 0.0;
            int far_idx = -1;
            for (size_t i = 0; i + 1 < samples.size(); ++i)
                for (size_t j = lo; j <= hi; ++j) {
                    if (segmentsIntersect(samples[i], samples[i + 1], bpts[j], bpts[j + 1],
                                          nullptr))
                        far_idx = (int)(side == 0 ? i + 1 : i);
                }
            if (far_idx >= 0) {
                const int a = side == 0 ? 0 : far_idx;
                const int b = side == 0 ? far_idx : (int)samples.size() - 1;
                for (int i = a; i <= b; ++i) {
                    double d = 1e18;
                    for (size_t j = lo; j <= hi; ++j)
                        d = std::min(d, distPointToSegment(samples[i], bpts[j], bpts[j + 1]));
                    graze = std::max(graze, d);
                }
            }
            std::printf("            adhere@%s run=[%zu,%zu] cross_in_run=%d "
                        "graze_depth=%.3fm\n",
                        side == 0 ? "p0" : "p1", lo, hi, far_idx >= 0, graze);
        }
    }
}

void printCtrl(const char* tag, const BezierCurve& c) {
    std::printf("      %s segs=%d", tag, c.numSegments());
    for (int i = 0; i < c.numSegments(); ++i) {
        std::printf(" | ");
        for (int k = 0; k < 4; ++k)
            std::printf("(%.2f,%.2f)", c.segs[i].ctrl[k].x(), c.segs[i].ctrl[k].y());
    }
    std::printf("\n");
}

}  // namespace

int main(int argc, char** argv) {
    std::string path;
    std::unordered_set<std::string> only_ids;
    bool show_all = false;
    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "--id") == 0 && i + 1 < argc)
            only_ids.insert(argv[++i]);
        else if (std::strcmp(argv[i], "--all") == 0)
            show_all = true;
        else if (path.empty())
            path = argv[i];
    }
    if (path.empty()) {
        std::fprintf(stderr, "用法: diag_single_cubic <数据文件> [--id 连接id]... [--all]\n");
        return 2;
    }

    IntersectionInput raw = IntersectionIO::loadFromFile(path);
    if (raw.connectivities.empty()) {
        std::fprintf(stderr, "数据集无连通关系\n");
        return 2;
    }

    IntersectionShapeGenerator gen;
    IntersectionOutput output;
    if (!gen.generate(raw, output)) {
        std::fprintf(stderr, "generate() 返回 false\n");
        return 1;
    }

    IntersectionInput input = InputNormalizer(raw);
    const ConnectivityDirectionConfig dir_cfg;
    ConnectivityDirectionNormalizer(input, dir_cfg);

    SDFField sdf;
    if (!input.obstacles.empty())
        sdf.build(outerBoundingBox(input), input.obstacles, 0.2,
                  obstacleAvoidanceClearanceForMode(input.mode));

    ClusterOrderSolver solver;
    solver.build(input.connectivities, input.lanes, input.lane_groups, input.crosswalks);
    const ClusterTopology topology = solver.topology();

    SceneContext scene(input);
    scene.bindSdf(&sdf);
    scene.bindClusterTopology(&topology);

    std::unordered_map<ConnId, const ConnectivityCurve*> cmap;
    GenerationState state;
    for (const auto& cc : output.connectivity_curves) {
        cmap[cc.id] = &cc;
        if (cc.curve) state.accepted_curves[cc.id] = *cc.curve;
    }

    const UTurnFamilyBuilder family_builder;
    const ConstraintEvaluator evaluator;
    const OrdinaryCurveInitializer ordinary;
    const ConstraintProfile profile = diagProfile(input.mode);

    int n_multi = 0;
    for (const auto& conn : input.connectivities) {
        auto it = cmap.find(conn.id);
        if (it == cmap.end() || !it->second->curve || it->second->curve->empty())
            continue;
        if (!only_ids.empty() && only_ids.count(conn.id) == 0)
            continue;
        const ConnectivityCurve& cc = *it->second;
        const BezierCurve& final_curve = *cc.curve;

        const UTurnFamilyInfo family = family_builder.build(
            conn, input, dir_cfg.uturn_alignment_scope, &solver, false);
        CurveGenerationContext context = CurveGenerationContextBuilder().build(
            scene, conn, dir_cfg.uturn_alignment_scope, &solver, &family);
        context.profile = profile;
        if (isGeometricUTurnLocal(context.entry, context.exit))
            continue;
        if (!show_all && final_curve.numSegments() == 1)
            continue;
        if (final_curve.numSegments() != 1) ++n_multi;

        const Vec2d chord = context.exit.first - context.entry.first;
        const double chord_len = chord.norm();
        const Vec2d t0 = context.entry.second.norm() > 1e-8
                             ? context.entry.second.normalized()
                             : chord.normalized();
        const double turn_strength =
            chord_len > 1e-8 ? std::abs(cross2d(t0, chord.normalized())) : 0.0;

        // 排除自身，避免“与自己相交”的伪报。
        GenerationState others = state;
        others.accepted_curves.erase(conn.id);

        const ConstraintReport final_report =
            evaluator.evaluate(final_curve, context, others);

        const BezierCurve single = ordinary.buildPreferredSingleCubic(
            context.entry.first, context.entry.second, context.exit.first,
            context.exit.second);
        const ConstraintReport single_report =
            evaluator.evaluate(single, context, others);

        const BezierCurve alpha_single = ordinary.buildSingleCubic(
            context.entry.first, context.entry.second, context.exit.first,
            context.exit.second, 0.4);
        const ConstraintReport alpha_report =
            evaluator.evaluate(alpha_single, context, others);

        std::printf("== %s  turn=%d fixed=%d lane_type=%d in_geom_pts=%zu\n",
                    conn.id.c_str(), (int)context.turn, (int)conn.fixed_shape,
                    (int)conn.lane_type, conn.geometry.points.size());
        std::printf("   chord=%.2fm turn_strength=%.3f status=%d reason=\"%s\"\n",
                    chord_len, turn_strength, (int)cc.status,
                    cc.violation.reason.c_str());
        std::printf("   FINAL  segs=%d arc/chord=%.3f maxk=%.3f  violated: %s\n",
                    final_curve.numSegments(),
                    chord_len > 1e-8 ? final_curve.arcLength() / chord_len : 0.0,
                    final_curve.maxCurvature(40), violatedIds(final_report).c_str());
        std::printf("   PREF1  segs=%d arc/chord=%.3f maxk=%.3f  violated: %s\n",
                    single.numSegments(),
                    chord_len > 1e-8 ? single.arcLength() / chord_len : 0.0,
                    single.maxCurvature(40), violatedIds(single_report).c_str());
        std::printf("   A0.4   segs=%d arc/chord=%.3f maxk=%.3f  violated: %s\n",
                    alpha_single.numSegments(),
                    chord_len > 1e-8 ? alpha_single.arcLength() / chord_len : 0.0,
                    alpha_single.maxCurvature(40), violatedIds(alpha_report).c_str());
        printCtrl("FINAL", final_curve);
        printCtrl("PREF1", single);
        reportBoundaryDetail("FINAL", final_curve, input, context.entry, context.exit);
        reportBoundaryDetail("PREF1", single, input, context.entry, context.exit);
    }
    std::printf("\n非单段普通曲线数量: %d\n", n_multi);
    return 0;
}
