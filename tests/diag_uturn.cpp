// Standalone diagnostic: run generator on 100000643-1.json and emit per-curve
// metrics for U-turns (G1 cos at p0/p1, maxk, arc/chord, same-cluster crossing
// pairs, obstacle penetration, fence overflow, sibling that crosses).
#include "constraints/cluster_order.h"
#include "curve/curve_utils.h"
#include "intersection_shape_generator.h"
#include "io/iodata_json.h"
#include "optimizer/sdf_field.h"
#include <cstdio>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using namespace isg;

static std::unordered_map<ConnId, const ConnectivityCurve*> curveMap(const IntersectionOutput& output) {
    std::unordered_map<ConnId, const ConnectivityCurve*> curves;
    for (const auto& cc : output.connectivity_curves)
        curves[cc.id] = &cc;
    return curves;
}

static bool isUTurn(const ConnectivityCurve& cc) {
    return cc.turn_type == ConnTurnType::UTurnLeft || cc.turn_type == ConnTurnType::UTurnRight;
}

int main(int argc, char** argv) {
    std::string path = (argc > 1) ? argv[1]
        : (std::string(PROJECT_ROOT_DIR) + "/datas/100000643-1.json");
    IntersectionInput input = IntersectionIO::loadFromFile(path);
    fprintf(stderr, "Loaded %s: %zu conns, %zu lanes\n",
            path.c_str(), input.connectivities.size(), input.lanes.size());

    IntersectionShapeGenerator gen;
    IntersectionOutput output;
    if (!gen.generate(input, output)) {
        fprintf(stderr, "Generation failed\n");
        return 1;
    }
    fprintf(stderr, "Generated %zu connectivity curves\n", output.connectivity_curves.size());

    ClusterOrderSolver solver;
    solver.build(input.connectivities, input.lanes, input.lane_groups, input.crosswalks);

    auto cmap = curveMap(output);

    // Build SDF (mirror main.cpp's coarse SDF) for obstacle penetration.
    SDFField sdf;
    auto roi = input.area.geometry.empty() ? input.area.geometry.bbox()
                                            : input.area.geometry.bbox();
    if (roi.width() < 1) {
        for (auto& l : input.lanes)
            for (auto& p : l.geometry.points) roi.expand(p);
        roi.min_pt -= Vec2d(20, 20);
        roi.max_pt += Vec2d(20, 20);
    }
    if (!input.obstacles.empty())
        sdf.build(roi, input.obstacles, 0.2, 0.0);

    // For each U-turn curve, gather metrics.
    printf("==== U-turn detailed report ====\n");
    printf("%-4s %8s %8s %8s %8s %8s %6s %6s %6s %6s %s\n",
           "id", "G1cos0", "G1cos1", "maxk", "arcLen", "chord", "a/c", "tg(m)", "obst", "sibX", "crossed_with [isect_pt]");

    std::vector<std::string> problem_ids;
    for (const auto& cc : output.connectivity_curves) {
        if (!cc.curve) continue;
        if (!isUTurn(cc)) continue;
        const auto& c = *cc.curve;

        auto entry = input.entryPtDir(cc.entry_lane_id);
        auto exit_ = input.exitPtDir(cc.exit_lane_id);
        Vec2d p0 = entry.first;
        Vec2d T0 = entry.second.norm() > 1e-8 ? entry.second.normalized() : Vec2d(1, 0);
        Vec2d p1 = exit_.first;
        Vec2d T1 = exit_.second.norm() > 1e-8 ? exit_.second.normalized() : Vec2d(1, 0);

        Vec2d st = c.startTan().norm() > 1e-8 ? c.startTan().normalized() : Vec2d(1, 0);
        Vec2d et = c.endTan().norm()   > 1e-8 ? c.endTan().normalized()   : Vec2d(1, 0);
        double g1_0 = st.dot(T0);
        double g1_1 = et.dot(T1);

        double maxk = c.maxCurvature(40);
        double arc = c.arcLength();
        double chord = (p1 - p0).norm();
        double ac = (chord > 1e-6) ? arc / chord : 0.0;

        // turn_gap (lateral separation in U-turn axis)
        Vec2d axis = T0 - T1;
        if (axis.norm() < 1e-8) axis = T0;
        axis.normalize();
        Vec2d lat_dir{-axis[1], axis[0]};
        double tg = std::abs((p1 - p0).dot(lat_dir));

        // obstacle penetration (min SDF along curve)
        double obst_pen = 0.0;
        if (sdf.valid()) {
            double ms = minSDFAlongCurve(c, sdf, 60);
            obst_pen = std::max(0.0, -ms);
        }

        // same-cluster sibling crossings — also report intersection point
        std::vector<std::string> crossed_with;
        for (const auto& p : solver.pairs()) {
            if (p.exempt == CrossExemption::StructuralCross) continue;
            ConnId other;
            if (p.id_a == cc.id) other = p.id_b;
            else if (p.id_b == cc.id) other = p.id_a;
            else continue;
            auto it = cmap.find(other);
            if (it == cmap.end() || !it->second->curve) continue;
            Vec2d isect_pt(0,0);
            bool hit = curvesIntersectBusiness(c, *it->second->curve, 1.5);
            if (hit) {
                char buf[64];
                // find approx isect point by sampling
                auto my_pts = c.sampleByArcLength(80);
                auto& other_c = *it->second->curve;
                auto oth_pts = other_c.sampleByArcLength(80);
                Vec2d best_pt(0,0); double best_d = 1e18;
                bool found_pt = false;
                for (size_t i = 0; i+1 < my_pts.size(); ++i) {
                    for (size_t j = 0; j+1 < oth_pts.size(); ++j) {
                        Vec2d ipt;
                        if (segmentsIntersect(my_pts[i], my_pts[i+1],
                                              oth_pts[j], oth_pts[j+1], &ipt)) {
                            double d = std::min((ipt - c.startPt()).norm(),
                                                (ipt - c.endPt()).norm());
                            if (d < best_d) { best_d = d; best_pt = ipt; found_pt = true; }
                        }
                    }
                }
                if (found_pt)
                    snprintf(buf, sizeof(buf), "%s@(%.1f,%.1f,d=%.1f)",
                             other.c_str(), best_pt.x(), best_pt.y(), best_d);
                else
                    snprintf(buf, sizeof(buf), "%s@(?)", other.c_str());
                crossed_with.push_back(buf);
            }
        }
        int sib_x = (int)crossed_with.size();

        printf("%-4s %8.3f %8.3f %8.3f %8.2f %8.2f %6.2f %6.2f %6.2f %6d  ",
               cc.id.c_str(), g1_0, g1_1, maxk, arc, chord, ac, tg, obst_pen, sib_x);
        for (size_t i = 0; i < crossed_with.size(); ++i) {
            if (i) printf(",");
            printf("%s", crossed_with[i].c_str());
        }
        printf("\n");

        // Flag problems: G1 < 0.99, or sibling crossing > 0, or obst_pen > 0.05
        if (g1_0 < 0.99 || g1_1 < 0.99 || sib_x > 0 || obst_pen > 0.05)
            problem_ids.push_back(cc.id);
    }

    printf("\n==== Problem U-turns (G1<0.99 OR sib_x>0 OR obst>0.05) ====\n");
    for (const auto& id : problem_ids) printf("  %s\n", id.c_str());

    // Overall summary: total same-cluster crossings among U-turns
    int total_uturn_x = 0;
    for (const auto& p : solver.pairs()) {
        if (p.exempt == CrossExemption::StructuralCross) continue;
        auto ia = cmap.find(p.id_a);
        auto ib = cmap.find(p.id_b);
        if (ia == cmap.end() || ib == cmap.end()) continue;
        if (!ia->second->curve || !ib->second->curve) continue;
        if (!curvesIntersectBusiness(*ia->second->curve, *ib->second->curve, 1.5)) continue;
        bool a_u = isUTurn(*ia->second);
        bool b_u = isUTurn(*ib->second);
        if (a_u || b_u) {
            ++total_uturn_x;
            printf("  U-involved crossing: %s(%s) - %s(%s)\n",
                   p.id_a.c_str(), a_u ? "U" : "-",
                   p.id_b.c_str(), b_u ? "U" : "-");
        }
    }
    printf("Total U-involved same-cluster crossings: %d\n", total_uturn_x);

    return 0;
}
