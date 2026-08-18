// Check if U-turn curves clear the crosswalk after optimization
#include "constraints/cluster_order.h"
#include "constraints/fence_check.h"
#include "curve/curve_utils.h"
#include "intersection_shape_generator.h"
#include "io/iodata_json.h"
#include "utils.h"
#include <cstdio>
#include <unordered_map>

using namespace isg;

int main(int argc, char** argv) {
    std::string path = (argc > 1) ? argv[1]
        : (std::string(PROJECT_ROOT_DIR) + "/datas/100000643.json");
    IntersectionInput input = IntersectionIO::loadFromFile(path);
    IntersectionShapeGenerator gen;
    IntersectionOutput output;
    if (!gen.generate(input, output)) return 1;

    std::unordered_map<ConnId, const ConnectivityCurve*> cmap;
    for (const auto& cc : output.connectivity_curves) cmap[cc.id] = &cc;

    // For each U-turn, check if the arc region is inside the crosswalk
    for (const auto& cc : output.connectivity_curves) {
        if (!cc.curve) continue;
        auto entry = input.entryPtDir(cc.entry_lane_id);
        auto exit_ = input.exitPtDir(cc.exit_lane_id);
        Vec2d p0 = entry.first;
        Vec2d t0 = entry.second.norm() > 1e-8 ? entry.second.normalized() : Vec2d(1,0);
        Vec2d p1 = exit_.first;
        Vec2d t1 = exit_.second.norm() > 1e-8 ? exit_.second.normalized() : Vec2d(1,0);
        // Only check U-turns
        if (t0.dot(t1) >= -0.5) continue;
        const auto& c = *cc.curve;
        printf("\n==== U-turn %s ====\n", cc.id.c_str());
        printf("  p0=(%.2f,%.2f) p1=(%.2f,%.2f) arc=%.2f segs=%zu\n",
               p0.x(), p0.y(), p1.x(), p1.y(), c.arcLength(), c.segs.size());

        // Check crosswalk clearance
        // Requirement (4.2): entry/exit straights MUST cross the crosswalk;
        // only the ARC (between q0 and q1) must be OUTSIDE the crosswalk.
        // The arc's westernmost point (apex) is the deepest part of the arc.
        // If the apex is outside the crosswalk, the entire arc is outside.
        auto pts = c.sampleByArcLength(120);
        // Find the apex (furthest point along U-turn axis = westernmost for west-facing U-turn)
        Vec2d T0_axis = (t0 - t1).normalized();
        if (T0_axis.dot(t0) < 0) T0_axis = -T0_axis;  // ensure axis points "into" the U-turn
        Vec2d apex_pt = p0;
        double max_proj = -1e18;
        for (auto& pt : pts) {
            double proj = (pt - p0).dot(T0_axis);
            if (proj > max_proj) { max_proj = proj; apex_pt = pt; }
        }
        // Check if apex is inside any crosswalk
        bool apex_inside = false;
        double apex_pen = 0;
        for (const auto& cw : input.crosswalks) {
            if (polygonContains(cw.geometry, apex_pt)) {
                apex_inside = true;
                apex_pen = std::max(apex_pen, pointToPolygonDist(apex_pt, cw.geometry));
            }
        }
        // Also check points within 3m of the apex (the true arc region)
        int n = pts.size();
        int arc_inside = 0;
        int arc_total = 0;
        double arc_max_pen = 0;
        for (int i = 0; i < n; ++i) {
            if ((pts[i] - apex_pt).norm() > 3.0) continue;  // only check near apex
            ++arc_total;
            for (const auto& cw : input.crosswalks) {
                if (polygonContains(cw.geometry, pts[i])) {
                    ++arc_inside;
                    arc_max_pen = std::max(arc_max_pen, pointToPolygonDist(pts[i], cw.geometry));
                    break;
                }
            }
        }
        printf("  apex=(%.2f,%.2f) apex_inside=%s apex_pen=%.3f\n",
               apex_pt.x(), apex_pt.y(), apex_inside ? "YES" : "no", apex_pen);
        printf("  arc_region (within 3m of apex): %d/%d pts inside crosswalk, max_pen=%.3f\n",
               arc_inside, arc_total, arc_max_pen);
        if (apex_inside || arc_inside > 2) {
            printf("  *** (4.2) VIOLATION: U-turn arc enters crosswalk! ***\n");
        }
    }
    return 0;
}
