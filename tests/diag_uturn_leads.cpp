// Standalone diagnostic: explain why a geometric U-turn fails to keep the
// straight-arc-straight (three-segment) form.
//
// For each requested connectivity id it prints
//   * the crosswalk clearance floors measured ahead of the entry and behind
//     the exit (the physical "must fully cross the crosswalk" requirement),
//   * the family-aggregated lead floors / aligned station used by generation,
//   * the generated curve as produced by the full pipeline, and
//   * a gate-by-gate verdict for every (lead_extra0, lead_extra1, arc_alpha)
//     candidate SegmentedUTurnCandidateSearch would enumerate, so the gate
//     that empties the candidate set is visible directly.
//
// Usage: diag_uturn_leads <data.json> <conn_id> [conn_id...]
#include "constraints/cluster_order.h"
#include "curve/curve_utils.h"
#include "generation/uturn_shape.h"
#include "initialization/uturn_curve_initializer.h"
#include "intersection_shape_generator.h"
#include "io/iodata_json.h"
#include "preprocessing/crosswalk_clearance_calculator.h"
#include "preprocessing/uturn_family_builder.h"
#include "toolkits/toolkits.h"
#include "utils.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

using namespace isg;

namespace {

const Connectivity* findConn(const IntersectionInput& input, const ConnId& id) {
    for (const auto& conn : input.connectivities)
        if (conn.id == id)
            return &conn;
    return nullptr;
}

// Optional overrides so the candidate sweep can be replayed with exactly the
// (min_lead0, min_lead1, aligned_station, stagger) quadruple that generation
// logged via ISG_DEBUG_UTURN, instead of the builder's own recomputation.
double envOverride(const char* name, double fallback) {
    const char* raw = std::getenv(name);
    if (!raw || !*raw)
        return fallback;
    return std::atof(raw);
}

void printClearance(const char* label, const CrosswalkClearanceResult& r) {
    if (!r.found) {
        printf("  %s: no crosswalk on ray\n", label);
        return;
    }
    printf("  %s: cw=%s near=%.3f far=%.3f clearance(far+0.30)=%.3f\n",
           label, r.crosswalk_id.c_str(), r.near, r.far, r.clearance);
}

void printCurve(const char* label, const BezierCurve& curve) {
    printf("  %s: segments=%d arc=%.3f maxk=%.3f\n", label,
           curve.numSegments(), curve.arcLength(), curve.maxCurvature(40));
    for (std::size_t s = 0; s < curve.segs.size(); ++s) {
        const auto& g = curve.segs[s].ctrl;
        printf("    seg[%zu] len=%.3f ctrl:", s, (g[3] - g[0]).norm());
        for (int k = 0; k < 4; ++k)
            printf(" (%.2f,%.2f)", g[k].x(), g[k].y());
        printf("\n");
    }
}

// Mirrors the gate order inside SegmentedUTurnCandidateSearch::search so the
// first failing gate for each candidate is directly comparable with it.
std::string firstFailingGate(
    const BezierCurve& candidate, double min_lead0, double min_lead1,
    const Vec2d& axis, double chord_len, double alignment_station,
    double max_curvature, const std::vector<Crosswalk>& clearance_crosswalks) {
    char buf[160];
    if (candidate.empty())
        return "empty";
    if (candidate.numSegments() != 3) {
        snprintf(buf, sizeof(buf), "segments=%d", candidate.numSegments());
        return buf;
    }
    if (!segmentedUTurnHasMinimumStraightLeads(candidate, min_lead0, min_lead1)) {
        snprintf(buf, sizeof(buf), "min_lead lead0=%.3f/%.3f lead1=%.3f/%.3f",
                 (candidate.segs.front().ctrl[3] -
                  candidate.segs.front().ctrl[0]).norm(), min_lead0,
                 (candidate.segs.back().ctrl[3] -
                  candidate.segs.back().ctrl[0]).norm(), min_lead1);
        return buf;
    }
    if (curveSelfIntersectsBusiness(candidate, 1.0))
        return "self_intersect";
    const double curvature = candidate.maxCurvature(40);
    if (curvature >= max_curvature) {
        snprintf(buf, sizeof(buf), "curvature %.3f>=%.3f", curvature, max_curvature);
        return buf;
    }
    const Vec2d q0 = candidate.segs.front().ctrl[3];
    const Vec2d q1 = candidate.segs.back().ctrl[0];
    if (std::abs((q0 - q1).dot(axis)) > 0.05) {
        snprintf(buf, sizeof(buf), "q0/q1 axis skew %.4f", (q0 - q1).dot(axis));
        return buf;
    }
    if (std::isfinite(alignment_station) &&
        (std::abs(q0.dot(axis) - alignment_station) > 0.05 ||
         std::abs(q1.dot(axis) - alignment_station) > 0.05)) {
        snprintf(buf, sizeof(buf), "station q0=%.3f q1=%.3f target=%.3f",
                 q0.dot(axis), q1.dot(axis), alignment_station);
        return buf;
    }
    const double arc_chord = candidate.arcLength() / chord_len;
    if (arc_chord < 1.35) {
        snprintf(buf, sizeof(buf), "arc/chord %.3f<1.35", arc_chord);
        return buf;
    }
    if (chord_len >= 1.0 && !segmentedUTurnMiddleArcLooksRound(candidate, axis)) {
        const BezierSegment& arc = candidate.segs[1];
        const double gap = (arc.ctrl[3] - arc.ctrl[0]).norm();
        snprintf(buf, sizeof(buf), "arc_not_round gap=%.3f arc/gap=%.3f",
                 gap, gap > 1e-9 ? arc.arcLength(32) / gap : 0.0);
        return buf;
    }
    if (!segmentedUTurnMiddleArcClearsCrosswalks(candidate, clearance_crosswalks))
        return "arc_in_crosswalk";
    return "PASS";
}

void diagnose(const IntersectionInput& input, const IntersectionOutput& output,
              const ClusterOrderSolver& solver, const ConnId& id) {
    const Connectivity* conn = findConn(input, id);
    if (!conn) {
        printf("\n==== %s: not found ====\n", id.c_str());
        return;
    }
    printf("\n==== %s ====\n", id.c_str());
    printf("  entry_lane=%s exit_lane=%s turn_type=%d\n",
           conn->entry_lane_id.c_str(), conn->exit_lane_id.c_str(),
           conn->turn_type);

    const std::pair<Vec2d, Vec2d> entry = input.entryPtDir(conn->entry_lane_id);
    const std::pair<Vec2d, Vec2d> exit = input.exitPtDir(conn->exit_lane_id);
    printf("  p0=(%.3f,%.3f) T0=(%.4f,%.4f)\n", entry.first.x(), entry.first.y(),
           entry.second.x(), entry.second.y());
    printf("  p1=(%.3f,%.3f) T1=(%.4f,%.4f)\n", exit.first.x(), exit.first.y(),
           exit.second.x(), exit.second.y());

    const UTurnFamilyBuilder builder;
    if (!builder.isGeometricUTurn(*conn, input)) {
        printf("  not a geometric U-turn\n");
        return;
    }

    const CrosswalkClearanceCalculator calculator;
    printClearance("entry ray ahead",
                   calculator.ahead(entry.first, entry.second, input));
    printClearance("exit ray behind",
                   calculator.behind(exit.first, exit.second, input));

    const UTurnLeadFloors own = builder.leadFloors(
        entry.first, entry.second, exit.first, exit.second, input);
    printf("  own lead floors: lead0=%.3f lead1=%.3f cw0=%d cw1=%d\n",
           own.lead0, own.lead1, (int)own.crosswalk0, (int)own.crosswalk1);

    const UTurnFamilyInfo family = builder.build(
        *conn, input, UTurnAlignmentScope::LaneEndpoint, &solver, true);
    printf("  family: size=%zu rank=%zu radius=%.3f axis=(%.4f,%.4f)\n",
           family.rank.family_size, family.rank.reverse_radius_rank,
           family.radius, family.axis.x(), family.axis.y());
    printf("  family lead floors: lead0=%.3f lead1=%.3f aligned_station=%.3f\n",
           family.lead0, family.lead1, family.aligned_station);
    printf("  clearance crosswalks: %zu of %zu\n",
           family.clearance_crosswalks.size(), input.crosswalks.size());

    const std::vector<const Connectivity*> component = builder.alignmentComponent(
        *conn, input, UTurnAlignmentScope::LaneEndpoint, &solver, false);
    printf("  alignment component:");
    for (const Connectivity* m : component)
        printf(" %s", m->id.c_str());
    printf("\n");

    const double max_curvature = segmentedUTurnMaxCurvatureLimit(
        entry.first, entry.second, exit.first, exit.second);
    const double chord_len = (exit.first - entry.first).norm();
    printf("  chord=%.3f max_curvature_limit=%.3f\n", chord_len, max_curvature);

    // Alignment station is only enforced for multi-member families, matching
    // ConnectivityGenerationSession.
    const double alignment_station = envOverride(
        "ISG_DIAG_STATION",
        family.rank.family_size > 1
            ? family.aligned_station
            : std::numeric_limits<double>::quiet_NaN());
    const double min_lead0 = envOverride("ISG_DIAG_MIN_LEAD0", family.lead0);
    const double min_lead1 = envOverride("ISG_DIAG_MIN_LEAD1", family.lead1);
    const double stagger = envOverride("ISG_DIAG_STAGGER", 0.0);
    printf("  sweep params: min_lead0=%.4f min_lead1=%.4f station=%.4f "
           "stagger=%.4f\n", min_lead0, min_lead1, alignment_station, stagger);
    printf("  p0.axis=%.4f p1.axis=%.4f => station needs lead0=%.4f lead1=%.4f\n",
           entry.first.dot(family.axis), exit.first.dot(family.axis),
           alignment_station - entry.first.dot(family.axis),
           alignment_station - exit.first.dot(family.axis));

    for (const auto& cc : output.connectivity_curves) {
        if (cc.id != id || !cc.curve)
            continue;
        printCurve("generated", *cc.curve);
        printf("    status=%d reason=%s\n", (int)cc.status,
               cc.violation.reason.c_str());
    }

    const std::vector<double> lead_extras =
        {0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0};
    const std::vector<double> arc_alphas =
        {2.0 / 3.0, 0.75, 0.85, 1.0, 1.25, 1.50, 1.75, 2.0,
         0.58, 0.50, 0.38, 0.28, 0.16};
    std::map<std::string, int> verdicts;
    std::string first_pass;
    int total = 0;
    for (double extra0 : lead_extras) {
        for (double extra1 : lead_extras) {
            for (double alpha : arc_alphas) {
                const BezierCurve candidate =
                    UTurnCurveInitializer().buildSegmented(
                        entry.first, entry.second, exit.first, exit.second,
                        min_lead0, min_lead1, alpha, stagger,
                        extra0, extra1, stagger, stagger);
                std::string gate = firstFailingGate(
                    candidate, min_lead0, min_lead1, family.axis,
                    chord_len, alignment_station, max_curvature,
                    family.clearance_crosswalks);
                ++total;
                if (extra0 == 0.0 && extra1 == 0.0 &&
                    std::abs(alpha - 2.0 / 3.0) < 1e-9 &&
                    candidate.numSegments() == 3) {
                    const Vec2d q0 = candidate.segs.front().ctrl[3];
                    const Vec2d q1 = candidate.segs.back().ctrl[0];
                    printf("  base candidate (extra 0/0 alpha=2/3): "
                           "lead0=%.4f lead1=%.4f q0.axis=%.4f q1.axis=%.4f "
                           "d_station=%.4f/%.4f\n",
                           (q0 - candidate.segs.front().ctrl[0]).norm(),
                           (candidate.segs.back().ctrl[3] - q1).norm(),
                           q0.dot(family.axis), q1.dot(family.axis),
                           q0.dot(family.axis) - alignment_station,
                           q1.dot(family.axis) - alignment_station);
                }
                // Collapse the numeric detail so the histogram stays readable.
                const std::size_t space = gate.find(' ');
                verdicts[space == std::string::npos ? gate : gate.substr(0, space)] += 1;
                if (gate == "PASS" && first_pass.empty()) {
                    char buf[64];
                    snprintf(buf, sizeof(buf), "extra0=%.2f extra1=%.2f alpha=%.3f",
                             extra0, extra1, alpha);
                    first_pass = buf;
                    printCurve("first passing candidate", candidate);
                }
                if (gate != "PASS" && extra0 == 0.0 && extra1 == 0.0)
                    printf("    [extra 0/0 alpha=%.3f] %s\n", alpha, gate.c_str());
            }
        }
    }
    printf("  candidate gate histogram (%d candidates):\n", total);
    for (const auto& kv : verdicts)
        printf("    %-24s %d\n", kv.first.c_str(), kv.second);
    if (!first_pass.empty())
        printf("  first pass at %s\n", first_pass.c_str());
    else
        printf("  NO CANDIDATE PASSES THE SHAPE GATES\n");
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 3) {
        fprintf(stderr, "Usage: %s <data.json> <conn_id> [conn_id...]\n", argv[0]);
        return 1;
    }
    const std::string path = argv[1];
    // 生成流程内部先做 InputNormalizer 规范化，再由
    // ConnectivityDirectionNormalizer 统一首尾方向；车道端点/切向与人行横道
    // 都在这两步之后的坐标系里。诊断必须用同一份输入，否则站位、lead 下限
    // 都不可与 ISG_DEBUG_UTURN 的日志直接比对。
    IntersectionInput input = IntersectionIO::loadFromFile(path);
    input = InputNormalizer(input);
    ConnectivityDirectionNormalizer(input, ConnectivityDirectionConfig());
    fprintf(stderr, "Loaded %s: %zu conns, %zu crosswalks\n", path.c_str(),
            input.connectivities.size(), input.crosswalks.size());

    IntersectionShapeGenerator gen;
    IntersectionOutput output;
    if (!gen.generate(input, output)) {
        fprintf(stderr, "Generation failed\n");
        return 1;
    }
    ClusterOrderSolver solver;
    solver.build(input.connectivities, input.lanes, input.lane_groups,
                 input.crosswalks);

    for (int i = 2; i < argc; ++i)
        diagnose(input, output, solver, ConnId(argv[i]));
    return 0;
}
