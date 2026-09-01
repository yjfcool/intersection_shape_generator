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
#include "curve/curve_utils.h"
#include "domain/scene_context.h"
#include "generation/connectivity_generation_context.h"
#include "intersection_shape_generator.h"
#include "io/iodata_json.h"
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

    IntersectionInput input = IntersectionIO::loadFromFile(path);
    IntersectionShapeGenerator gen;
    IntersectionOutput output;
    if (!gen.generate(input, output)) {
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
    return 0;
}
