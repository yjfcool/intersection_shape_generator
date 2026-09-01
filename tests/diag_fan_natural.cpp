// Diagnostic: does a shared-endpoint fan-out family come out of *initialization*
// already ordered, and if not, which stage broke it?
//
// For each member it prints two shapes side by side:
//   - natural: OrdinaryCurveInitializer::buildPreferredSingleCubic from the raw
//     entry/exit point+tangent, i.e. before any candidate search or pair repair
//   - final:   the curve the generator actually emitted
// plus the shared-endpoint turn-in rate kappa0 for both, and the pairwise
// business intersection for both. If natural is monotone/non-crossing and final
// is not, the ordering was broken downstream, not at initialization.
//
// Usage: diag_fan_natural <data.json> <conn_id> <conn_id> [...]
#include "curve/curve_utils.h"
#include "initialization/ordinary_curve_initializer.h"
#include "intersection_shape_generator.h"
#include "io/iodata_json.h"
#include "utils.h"

#include <cstdio>
#include <string>
#include <unordered_map>
#include <vector>

using namespace isg;

namespace {

struct Row {
    ConnId id;
    BezierCurve natural;
    const BezierCurve* final_curve = nullptr;
    double far_lat = 0.0;
    double turn_deg = 0.0;
};

void printShape(const char* tag, const ConnId& id, const BezierCurve& c,
                bool at_start) {
    if (c.empty()) {
        printf("  %-4s %-8s (empty)\n", id.c_str(), tag);
        return;
    }
    const auto& g = c.segs.front().ctrl;
    const auto& t = c.segs.back().ctrl;
    printf("  %-4s %-8s segs=%d h0=%8.3f h1=%8.3f arc=%8.3f maxk=%.4f "
           "kappa0=%+.6f\n",
           id.c_str(), tag, c.numSegments(), (g[1] - g[0]).norm(),
           (t[3] - t[2]).norm(), c.arcLength(), c.maxCurvature(60),
           singleCubicSignedEndCurvature(c, at_start));
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 4) {
        fprintf(stderr, "usage: diag_fan_natural <data.json> <id> <id> [...]\n");
        return 2;
    }
    const std::string path = argv[1];
    std::vector<ConnId> ids;
    for (int i = 2; i < argc; ++i)
        ids.push_back(argv[i]);

    IntersectionInput input = IntersectionIO::loadFromFile(path);
    IntersectionShapeGenerator gen;
    IntersectionOutput output;
    if (!gen.generate(input, output)) {
        fprintf(stderr, "GEN FAILED\n");
        return 1;
    }
    std::unordered_map<ConnId, const ConnectivityCurve*> cmap;
    for (const auto& cc : output.connectivity_curves)
        cmap[cc.id] = &cc;
    std::unordered_map<ConnId, const Connectivity*> conns;
    for (const auto& c : input.connectivities)
        conns[c.id] = &c;

    const OrdinaryCurveInitializer initializer;
    std::vector<Row> rows;
    for (const ConnId& id : ids) {
        if (!conns.count(id)) {
            printf("conn %s: missing in input\n", id.c_str());
            continue;
        }
        const Connectivity& conn = *conns[id];
        const std::pair<Vec2d, Vec2d> entry = input.entryPtDir(conn.entry_lane_id);
        const std::pair<Vec2d, Vec2d> exit = input.exitPtDir(conn.exit_lane_id);
        Row row;
        row.id = id;
        row.natural = initializer.buildPreferredSingleCubic(
            entry.first, entry.second, exit.first, exit.second);
        row.final_curve = cmap.count(id) ? cmap[id]->curve.get() : nullptr;
        row.far_lat = cross2d(exit.second.normalized(), exit.first - entry.first);
        const double dot = entry.second.normalized().dot(exit.second.normalized());
        row.turn_deg = std::acos(std::max(-1.0, std::min(1.0, dot))) * 180.0 / M_PI;
        rows.push_back(row);
    }

    printf("==== %s : natural vs final (shared entry endpoint) ====\n",
           path.c_str());
    for (const Row& r : rows) {
        printf("  %-4s turn=%7.2fdeg far_lat=%+9.3f\n", r.id.c_str(), r.turn_deg,
               r.far_lat);
        printShape("natural", r.id, r.natural, true);
        if (r.final_curve)
            printShape("final", r.id, *r.final_curve, true);
    }

    printf("\n---- pairwise (natural | final) ----\n");
    for (std::size_t i = 0; i < rows.size(); ++i) {
        for (std::size_t j = i + 1; j < rows.size(); ++j) {
            const Row& a = rows[i];
            const Row& b = rows[j];
            const bool nat_cross = !a.natural.empty() && !b.natural.empty() &&
                curvesIntersectBusiness(a.natural, b.natural, 1.5);
            const bool nat_poly = !a.natural.empty() && !b.natural.empty() &&
                sharedEndpointControlPolylinesCross(a.natural, b.natural, 0.30);
            const bool fin_cross = a.final_curve && b.final_curve &&
                curvesIntersectBusiness(*a.final_curve, *b.final_curve, 1.5);
            const bool fin_poly = a.final_curve && b.final_curve &&
                sharedEndpointControlPolylinesCross(*a.final_curve,
                                                    *b.final_curve, 0.30);
            printf("  %-4s|%-4s natural cross=%d poly=%d   final cross=%d poly=%d%s\n",
                   a.id.c_str(), b.id.c_str(), nat_cross ? 1 : 0,
                   nat_poly ? 1 : 0, fin_cross ? 1 : 0, fin_poly ? 1 : 0,
                   (!nat_cross && fin_cross) ? "   <== broken downstream" : "");
        }
    }
    return 0;
}
