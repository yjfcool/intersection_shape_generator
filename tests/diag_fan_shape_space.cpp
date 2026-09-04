// Diagnostic: for one member of a shared-endpoint fan, enumerate its whole
// single-cubic handle space and report which handle pairs (h0, h1) clear every
// *already generated* sibling curve.
//
// Motivation: 110003285's entry fan 43103892 (straight 14, lefts 15/16/18/19/24)
// keeps two same-cluster crossings (14|18, 14|24) after generation. The natural
// shapes cross too, so the question is not "who broke the order" but "does a
// non-crossing shape exist at all, and where in the handle space is it".
//
// For the target connectivity it prints, per feasible cell, the handle pair,
// kappa0, arc/chord, max curvature, and the set of siblings it still crosses.
//
// Usage: diag_fan_shape_space <data.json> <target_id> <sibling_id> [...]
#include "constraints/shape_constraint.h"
#include "constraints/constraint_evaluator.h"
#include "curve/curve_utils.h"
#include "domain/scene_context.h"
#include "generation/connectivity_generation_context.h"
#include "intersection_shape_generator.h"
#include "io/iodata_json.h"
#include "toolkits/toolkits.h"
#include "utils.h"

#include <cstdio>
#include <string>
#include <unordered_map>
#include <vector>

using namespace isg;

namespace {

BezierCurve makeCubic(const Vec2d& p0, const Vec2d& t0, const Vec2d& p1,
                      const Vec2d& t1, double h0, double h1) {
    BezierSegment seg;
    seg.ctrl[0] = p0;
    seg.ctrl[1] = p0 + t0.normalized() * h0;
    seg.ctrl[2] = p1 - t1.normalized() * h1;
    seg.ctrl[3] = p1;
    BezierCurve curve;
    curve.segs.push_back(seg);
    return curve;
}

bool hasPhysicalViolation(const BezierCurve& curve,
                          const SceneContext& scene,
                          const Connectivity& conn) {
    CurveGenerationContext context = CurveGenerationContextBuilder().build(
        scene, conn);
    context.profile.enforce_fence = true;
    // Match the generator's strict candidate gate rather than the evaluator's
    // lower default sample count; narrow fence notches must not be skipped.
    context.profile.samples = 64;
    context.profile.check_self_intersection = true;
    context.profile.check_obstacle = true;
    context.profile.check_boundary = true;
    context.profile.check_g1 = false;
    context.profile.check_curvature = false;
    context.profile.check_ordinary_shape = false;
    context.profile.check_uturn_shape = false;
    context.profile.check_crosswalk = false;
    context.profile.check_cluster = false;
    context.profile.road_edge_clearance = 0.0;
    const ConstraintReport report = ConstraintEvaluator().evaluate(
        curve, context, GenerationState());
    for (const auto& result : report.results)
        if (result.state == ConstraintState::Violated)
            return true;
    return false;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 4) {
        fprintf(stderr,
                "usage: diag_fan_shape_space <data.json> <target> <sibling>...\n");
        return 2;
    }
    const std::string path = argv[1];
    const ConnId target_id = argv[2];
    std::vector<ConnId> sibling_ids;
    for (int i = 3; i < argc; ++i)
        sibling_ids.push_back(argv[i]);

    // The generator owns the direction-normalization pass. Keep its input raw
    // for the formal output, and build a separate once-normalized copy for the
    // candidate-space inspection below; passing an already normalized input to
    // generate() would apply the pass twice and change exit tangents.
    const IntersectionInput raw_input = IntersectionIO::loadFromFile(path);
    IntersectionInput input = InputNormalizer(raw_input);
    ConnectivityDirectionNormalizer(input, ConnectivityDirectionConfig{});
    IntersectionShapeGenerator gen;
    IntersectionOutput output;
    if (!gen.generate(raw_input, output)) {
        fprintf(stderr, "GEN FAILED\n");
        return 1;
    }
    std::unordered_map<ConnId, const BezierCurve*> finals;
    for (const auto& cc : output.connectivity_curves) {
        if (cc.curve)
            finals[cc.id] = cc.curve.get();
    }
    const Connectivity* target = nullptr;
    for (const auto& c : input.connectivities) {
        if (c.id == target_id)
            target = &c;
    }
    if (!target) {
        fprintf(stderr, "target %s not in input\n", target_id.c_str());
        return 1;
    }

    SceneContext scene(input);
    const CurveGenerationContextBuilder builder;
    const CurveGenerationContext context = builder.build(scene, *target);

    const std::pair<Vec2d, Vec2d> entry = input.entryPtDir(target->entry_lane_id);
    const std::pair<Vec2d, Vec2d> exit = input.exitPtDir(target->exit_lane_id);
    const OrdinarySingleCubicHandleBounds bounds = ordinarySingleCubicHandleBounds(
        entry.first, entry.second, exit.first, exit.second);
    const double chord = (exit.first - entry.first).norm();
    printf("==== %s target=%s chord=%.3f h0=[%.3f,%.3f] h1=[%.3f,%.3f] ====\n",
           path.c_str(), target_id.c_str(), chord, bounds.start_min,
           bounds.start_max, bounds.end_min, bounds.end_max);
    if (finals.count(target_id)) {
        const BezierCurve& c = *finals[target_id];
        if (c.numSegments() == 1) {
            const auto& g = c.segs.front().ctrl;
            printf("  final: h0=%.3f h1=%.3f kappa0=%+.6f arc/chord=%.4f maxk=%.4f\n",
                   (g[1] - g[0]).norm(), (g[3] - g[2]).norm(),
                   singleCubicSignedEndCurvature(c, true),
                   c.arcLength() / std::max(1e-9, chord), c.maxCurvature(60));
        } else {
            printf("  final: %d segments\n", c.numSegments());
        }
    }
    for (const ConnId& sibling_id : sibling_ids) {
        auto sibling_it = finals.find(sibling_id);
        if (sibling_it == finals.end() || sibling_it->second->numSegments() != 1)
            continue;
        const auto& g = sibling_it->second->segs.front().ctrl;
        printf("  sibling %s: h0=%.3f h1=%.3f\n", sibling_id.c_str(),
               (g[1] - g[0]).norm(), (g[3] - g[2]).norm());
    }

    const int kSteps = 16;
    int feasible = 0;
    int clean = 0;
    for (int i = 0; i <= kSteps; ++i) {
        const double h0 = bounds.start_min +
            (bounds.start_max - bounds.start_min) * i / kSteps;
        for (int j = 0; j <= kSteps; ++j) {
            const double h1 = bounds.end_min +
                (bounds.end_max - bounds.end_min) * j / kSteps;
            const BezierCurve candidate = makeCubic(
                entry.first, entry.second, exit.first, exit.second, h0, h1);
            if (!ordinarySingleCubicControlsValid(candidate, entry.first,
                                                  entry.second, exit.first,
                                                  exit.second))
                continue;
            const ConstraintResult shape = evaluateOrdinaryShape(candidate, context);
            if (shape.state == ConstraintState::Violated) {
                if (std::getenv("ISG_DIAG_SHOW_REJECT"))
                    printf("  h0=%8.3f h1=%8.3f REJECT %s\n", h0, h1,
                           shape.reason.c_str());
                continue;
            }
            ++feasible;
            // 同时按两套口径统计：1.5m 是最终审计口径（diag_all_violations），
            // 0.15m 是配对修复内部的候选硬门口径。两者不一致时，审计已经
            // 干净的候选仍可能在修复阶段被 0.15m 判为穿越而剔除，因此把
            // tight 列表单独打出来，便于定位这类口径落差。
            std::string crossed;
            std::string crossed_tight;
            for (const ConnId& sid : sibling_ids) {
                if (!finals.count(sid))
                    continue;
                if (curvesIntersectBusiness(candidate, *finals[sid], 1.5)) {
                    crossed += sid;
                    crossed += " ";
                }
                if (curvesIntersectBusiness(candidate, *finals[sid], 0.15)) {
                    crossed_tight += sid;
                    crossed_tight += " ";
                }
            }
            if (crossed.empty())
                ++clean;
            printf("  h0=%8.3f h1=%8.3f kappa0=%+.6f arc/chord=%.4f maxk=%.4f  %s | tight: %s\n",
                   h0, h1, singleCubicSignedEndCurvature(candidate, true),
                   candidate.arcLength() / std::max(1e-9, chord),
                   candidate.maxCurvature(60),
                   crossed.empty() ? "CLEAN" : ("cross: " + crossed).c_str(),
                   crossed_tight.empty() ? "clean" : crossed_tight.c_str());
        }
    }
    printf("  shape-feasible cells=%d, of which clean=%d\n", feasible, clean);

    if (std::getenv("ISG_DIAG_PHYSICAL") && sibling_ids.size() >= 1) {
        std::vector<const Connectivity*> members;
        auto find_conn = [&](const ConnId& id) -> const Connectivity* {
            for (const auto& conn : input.connectivities)
                if (conn.id == id)
                    return &conn;
            return nullptr;
        };
        members.push_back(target);
        for (const ConnId& id : sibling_ids) {
            const Connectivity* conn = find_conn(id);
            if (conn)
                members.push_back(conn);
        }
        SceneContext physical_scene(input);
        std::vector<std::vector<BezierCurve>> safe(members.size());
        const int physical_steps = 32;
        for (std::size_t mi = 0; mi < members.size(); ++mi) {
            const auto entry_m = input.entryPtDir(members[mi]->entry_lane_id);
            const auto exit_m = input.exitPtDir(members[mi]->exit_lane_id);
            const auto bounds_m = ordinarySingleCubicHandleBounds(
                entry_m.first, entry_m.second, exit_m.first, exit_m.second,
                false);
            const double chord_m = (exit_m.first - entry_m.first).norm();
            const double turn_m = chord_m > 1e-8
                ? std::abs(cross2d(entry_m.second.normalized(),
                                   (exit_m.first - entry_m.first).normalized()))
                : 0.0;
            for (int i0 = 0; i0 <= physical_steps; ++i0) {
                const double f0 = 0.02 + 0.98 * i0 / physical_steps;
                const double h0 = std::max(bounds_m.start_min,
                                           f0 * bounds_m.start_max);
                for (int i1 = 0; i1 <= physical_steps; ++i1) {
                    const double f1 = 0.02 + 0.98 * i1 / physical_steps;
                    const double h1 = std::max(bounds_m.end_min,
                                               f1 * bounds_m.end_max);
                    const BezierCurve candidate = makeCubic(
                        entry_m.first, entry_m.second, exit_m.first,
                        exit_m.second, h0, h1);
                    if (!ordinarySingleCubicControlsValid(
                            candidate, entry_m.first, entry_m.second,
                            exit_m.first, exit_m.second, 1e-5, false) ||
                        curveSelfIntersectsBusiness(candidate, 1.0))
                        continue;
                    const double ratio = candidate.arcLength() / chord_m;
                    const double min_ratio = chord_m < 12.0 ? 1.02 : 1.06;
                    const double max_ratio = chord_m <= 15.0 ? 1.40 : 1.35;
                    if ((turn_m > 0.25 &&
                         (ratio < min_ratio || ratio > max_ratio ||
                          candidate.maxCurvature(20) > 2.5)) ||
                        (turn_m <= 0.25 &&
                         (ratio > 1.08 || candidate.maxCurvature(40) > 3.0)) ||
                        hasPhysicalViolation(candidate, physical_scene, *members[mi]))
                        continue;
                    safe[mi].push_back(candidate);
                }
            }
            printf("  physical-safe %s: %zu candidates\n",
                   members[mi]->id.c_str(), safe[mi].size());
            if (std::getenv("ISG_DIAG_PRINT_SAFE")) {
                printf("    safe handles %s:", members[mi]->id.c_str());
                for (const auto& candidate : safe[mi]) {
                    printf(" %.3f/%.3f",
                           (candidate.segs.front().ctrl[1] -
                            candidate.segs.front().ctrl[0]).norm(),
                           (candidate.segs.front().ctrl[3] -
                            candidate.segs.front().ctrl[2]).norm());
                }
                printf("\n");
            }
        }
        if (members.size() >= 2) {
            int pair_ok = 0;
            for (const auto& a : safe[0])
                for (const auto& b : safe[1])
                    if (!curvesIntersectBusiness(a, b, 1.5) &&
                        !sharedEndpointControlPolylinesCross(a, b, 0.30)) {
                        ++pair_ok;
                        if (pair_ok <= 8)
                            printf("    pair h=%.3f/%.3f %.3f/%.3f\n",
                                   (a.segs[0].ctrl[1] - a.segs[0].ctrl[0]).norm(),
                                   (a.segs[0].ctrl[3] - a.segs[0].ctrl[2]).norm(),
                                   (b.segs[0].ctrl[1] - b.segs[0].ctrl[0]).norm(),
                                   (b.segs[0].ctrl[3] - b.segs[0].ctrl[2]).norm());
                    }
            printf("  physical-safe pair %s|%s combinations=%d\n",
                   members[0]->id.c_str(), members[1]->id.c_str(), pair_ok);
        }
        if (members.size() == 3) {
            int triples = 0;
            for (const auto& a : safe[0]) {
                for (const auto& b : safe[1]) {
                    if (curvesIntersectBusiness(a, b, 1.5) ||
                        sharedEndpointControlPolylinesCross(a, b, 0.30))
                        continue;
                    for (const auto& c : safe[2]) {
                        if (!curvesIntersectBusiness(a, c, 1.5) &&
                            !sharedEndpointControlPolylinesCross(a, c, 0.30) &&
                            !curvesIntersectBusiness(b, c, 1.5) &&
                            !sharedEndpointControlPolylinesCross(b, c, 0.30)) {
                            ++triples;
                            if (triples <= 5)
                                printf("    triple h=%.3f/%.3f %.3f/%.3f %.3f/%.3f\n",
                                       (a.segs[0].ctrl[1] - a.segs[0].ctrl[0]).norm(),
                                       (a.segs[0].ctrl[3] - a.segs[0].ctrl[2]).norm(),
                                       (b.segs[0].ctrl[1] - b.segs[0].ctrl[0]).norm(),
                                       (b.segs[0].ctrl[3] - b.segs[0].ctrl[2]).norm(),
                                       (c.segs[0].ctrl[1] - c.segs[0].ctrl[0]).norm(),
                                       (c.segs[0].ctrl[3] - c.segs[0].ctrl[2]).norm());
                        }
                    }
                }
            }
            printf("  physical-safe triples=%d\n", triples);
        }
    }
    return 0;
}
