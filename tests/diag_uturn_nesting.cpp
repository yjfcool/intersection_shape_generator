// 同入/同出掉头家族的"嵌套"诊断：小径掉头必须整体落在大径掉头的控制面内。
//
// 打印内容
//   1. 家族每个成员的 familyLateralLadder 分档量与**实测**横向站位
//      （首/尾直段末端相对共享端点、沿 U 轴法线的有符号偏移，正方向指向
//      出口侧即"更靠里"）。阶梯反号在这里一眼可见。
//   2. 指定成员两两之间：小径中段（首中间控制点、尾中间控制点）是否落在
//      大径全部控制点围成的多边形内，以及小径中段把手折线与大径全部把手
//      折线是否相交。
//   3. 曲线级 curvesIntersectBusiness 结论，用于把"控制面前兆"与真实违约对齐。
//
// 用法: diag_uturn_nesting <data.json> [conn_id ...]
//       不给 conn_id 时遍历全部几何掉头家族。
#include "constraints/cluster_order.h"
#include "curve/curve_utils.h"
#include "intersection_shape_generator.h"
#include "io/iodata_json.h"
#include "preprocessing/uturn_family_builder.h"
#include "toolkits/toolkits.h"
#include "types.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

using namespace isg;

namespace {

bool pureGeometricInput(const IntersectionInput& input) {
    return input.obstacles.empty() && input.boundaries.empty() &&
           input.crosswalks.empty() && input.stop_lines.empty();
}

std::vector<Vec2d> allControlPoints(const BezierCurve& curve) {
    std::vector<Vec2d> pts;
    for (const auto& seg : curve.segs) {
        for (int k = 0; k < 4; ++k) {
            if (!pts.empty() && (pts.back() - seg.ctrl[k]).norm() < 1e-9)
                continue;
            pts.push_back(seg.ctrl[k]);
        }
    }
    return pts;
}

// 控制点按序首尾相连围成的多边形（奇偶判定）。
bool pointInPolygon(const Vec2d& p, const std::vector<Vec2d>& poly) {
    bool inside = false;
    const std::size_t n = poly.size();
    for (std::size_t i = 0, j = n - 1; i < n; j = i++) {
        const double yi = poly[i].y(), yj = poly[j].y();
        if ((yi > p.y()) == (yj > p.y()))
            continue;
        const double x = poly[j].x() +
            (p.y() - yj) / (yi - yj) * (poly[i].x() - poly[j].x());
        if (p.x() < x)
            inside = !inside;
    }
    return inside;
}

double distToPolygonEdges(const Vec2d& p, const std::vector<Vec2d>& poly) {
    double best = 1e18;
    const std::size_t n = poly.size();
    for (std::size_t i = 0, j = n - 1; i < n; j = i++) {
        const Vec2d a = poly[j], b = poly[i];
        const Vec2d ab = b - a;
        const double len2 = ab.squaredNorm();
        double t = len2 > 1e-18 ? (p - a).dot(ab) / len2 : 0.0;
        t = std::max(0.0, std::min(1.0, t));
        best = std::min(best, (p - (a + t * ab)).norm());
    }
    return best;
}

struct Handle {
    Vec2d a, b;
    std::string label;
};

std::vector<Handle> handlesOf(const BezierCurve& curve) {
    std::vector<Handle> out;
    for (std::size_t s = 0; s < curve.segs.size(); ++s) {
        const auto& g = curve.segs[s].ctrl;
        for (int k = 0; k < 3; ++k) {
            Handle h;
            h.a = g[k];
            h.b = g[k + 1];
            char buf[32];
            snprintf(buf, sizeof(buf), "seg%zu.P%d-P%d", s, k, k + 1);
            h.label = buf;
            out.push_back(h);
        }
    }
    return out;
}

// 中段：三段式取中弧，退化成单段时取该单段本身。用户口径里的"首/尾中间控制点"
// 就是这一段的 P1/P2，"中段把手"就是这一段的三条把手。
std::size_t middleSegmentIndex(const BezierCurve& curve) {
    return static_cast<std::size_t>(curve.numSegments()) / 2;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        fprintf(stderr, "用法: %s <data.json> [conn_id ...]\n", argv[0]);
        return 2;
    }
    IntersectionInput input = IntersectionIO::loadFromFile(argv[1]);
    input = InputNormalizer(input);
    ConnectivityDirectionNormalizer(input, ConnectivityDirectionConfig());

    std::set<ConnId> wanted;
    for (int i = 2; i < argc; ++i)
        wanted.insert(ConnId(argv[i]));

    IntersectionShapeGenerator gen;
    IntersectionOutput output;
    if (!gen.generate(input, output)) {
        fprintf(stderr, "Generation failed\n");
        return 1;
    }
    std::unordered_map<ConnId, const BezierCurve*> curves;
    for (const auto& cc : output.connectivity_curves)
        if (cc.curve)
            curves[cc.id] = cc.curve.get();

    ClusterOrderSolver solver;
    solver.build(input.connectivities, input.lanes, input.lane_groups,
                 input.crosswalks);
    const UTurnFamilyBuilder builder;
    const UTurnAlignmentScope scope = UTurnAlignmentScope::LaneEndpoint;
    const bool physical = !pureGeometricInput(input);

    std::set<ConnId> reported;
    for (const auto& conn : input.connectivities) {
        if (!builder.isGeometricUTurn(conn, input))
            continue;
        if (!wanted.empty() && !wanted.count(conn.id))
            continue;
        if (reported.count(conn.id))
            continue;
        const std::vector<const Connectivity*> component =
            builder.alignmentComponent(conn, input, scope, &solver, true);
        const double step = physical && component.size() >= 3 ? 0.25 : 0.01;
        printf("==== family of %s (size=%zu step=%.3f) ====\n",
               conn.id.c_str(), component.size(), step);

        // 家族参考系：取最大半径成员的入口切向为轴，法线正方向指向它的出口。
        const Connectivity* ref = nullptr;
        double ref_radius = -1.0;
        for (const Connectivity* m : component) {
            const double r = builder.radiusKey(*m, input);
            if (r > ref_radius) { ref_radius = r; ref = m; }
        }
        const std::pair<Vec2d, Vec2d> ref_entry =
            input.entryPtDir(ref->entry_lane_id);
        const std::pair<Vec2d, Vec2d> ref_exit =
            input.exitPtDir(ref->exit_lane_id);
        Vec2d T0 = ref_entry.second.normalized();
        Vec2d axis = (T0 - ref_exit.second.normalized()).normalized();
        if (axis.dot(T0) < 0.0) axis = -axis;
        Vec2d lateral{-axis.y(), axis.x()};
        if ((ref_exit.first - ref_entry.first).dot(lateral) < 0.0)
            lateral = -lateral;
        printf("  frame: axis=(%.5f,%.5f) lateral=(%.5f,%.5f) (正=指向出口侧/更靠里)\n",
               axis.x(), axis.y(), lateral.x(), lateral.y());

        struct Row { ConnId id; double radius; };
        std::vector<Row> rows;
        for (const Connectivity* m : component) {
            reported.insert(m->id);
            rows.push_back({m->id, builder.radiusKey(*m, input)});
        }
        std::sort(rows.begin(), rows.end(),
                  [](const Row& a, const Row& b) { return a.radius > b.radius; });

        for (const Row& row : rows) {
            const Connectivity* m = nullptr;
            for (const Connectivity* c : component)
                if (c->id == row.id) m = c;
            const UTurnFamilyLadder ladder = builder.familyLateralLadder(
                *m, input, scope, step, &solver, true);
            const std::pair<Vec2d, Vec2d> e = input.entryPtDir(m->entry_lane_id);
            const std::pair<Vec2d, Vec2d> x = input.exitPtDir(m->exit_lane_id);
            printf("  %-4s radius=%8.4f corridor_cap=%.4f ladder(entry=%.4f exit=%.4f)",
                   row.id.c_str(), row.radius, 0.25 * row.radius,
                   ladder.entry_stagger, ladder.exit_stagger);
            auto it = curves.find(row.id);
            if (it == curves.end()) { printf("  <no curve>\n"); continue; }
            const BezierCurve& c = *it->second;
            printf("  nseg=%d maxk=%.4f\n", c.numSegments(), c.maxCurvature(60));
            if (c.numSegments() == 3) {
                // 实测横向站位：直段末端相对本侧车道端点、沿家族法线的偏移。
                const Vec2d q0 = c.segs.front().ctrl[3];
                const Vec2d q1 = c.segs.back().ctrl[0];
                const double lead0 = (q0 - c.segs.front().ctrl[0]).norm();
                const double lead1 = (c.segs.back().ctrl[3] - q1).norm();
                const double lat0 =
                    (q0 - (e.first + lead0 * e.second.normalized())).dot(lateral);
                const double lat1 =
                    (q1 - (x.first - lead1 * x.second.normalized())).dot(lateral);
                printf("       lead0=%.4f lead1=%.4f  实测横向: entry=%+.4f exit=%+.4f"
                       "  (期望 entry=%+.4f exit=%+.4f)\n",
                       lead0, lead1, lat0, lat1,
                       ladder.entry_stagger, -ladder.exit_stagger);
                printf("       q0.lat=%+.4f q1.lat=%+.4f corridor=%.4f\n",
                       (q0 - e.first).dot(lateral), (q1 - e.first).dot(lateral),
                       (q1 - q0).dot(lateral));
            }
            for (std::size_t s = 0; s < c.segs.size(); ++s) {
                const auto& g = c.segs[s].ctrl;
                printf("       seg%zu", s);
                for (int k = 0; k < 4; ++k)
                    printf(" (%.4f,%.4f)", g[k].x(), g[k].y());
                printf("\n");
            }
        }

        // 两两嵌套判定：半径小的一方必须落在半径大的一方控制面内。
        printf("  ---- 嵌套判定（内=小径, 外=大径）----\n");
        for (std::size_t i = 0; i < rows.size(); ++i) {
            for (std::size_t j = i + 1; j < rows.size(); ++j) {
                const Row& outer = rows[i];
                const Row& inner = rows[j];
                auto io = curves.find(outer.id), ii = curves.find(inner.id);
                if (io == curves.end() || ii == curves.end()) continue;
                const BezierCurve& co = *io->second;
                const BezierCurve& ci = *ii->second;
                const std::vector<Vec2d> poly = allControlPoints(co);
                printf("  内 %-4s / 外 %-4s : curve_cross=%s",
                       inner.id.c_str(), outer.id.c_str(),
                       curvesIntersectBusiness(ci, co, 1.5) ? "TRUE" : "false");
                if (ci.numSegments() >= 1) {
                    const auto& mg = ci.segs[middleSegmentIndex(ci)].ctrl;
                    const Vec2d head = mg[1];
                    const Vec2d tail = mg[2];
                    const bool in_head = pointInPolygon(head, poly);
                    const bool in_tail = pointInPolygon(tail, poly);
                    printf("  首中间控制点%s(d=%.4f) 尾中间控制点%s(d=%.4f)",
                           in_head ? "在面内" : "**在面外**",
                           distToPolygonEdges(head, poly),
                           in_tail ? "在面内" : "**在面外**",
                           distToPolygonEdges(tail, poly));
                }
                printf("\n");
                // 小径中段把手 × 大径全部把手
                if (ci.numSegments() >= 1) {
                    const auto mid = handlesOf(ci);
                    const auto all_outer = handlesOf(co);
                    const std::size_t base = 3 * middleSegmentIndex(ci);
                    for (std::size_t a = base; a < base + 3 && a < mid.size();
                         ++a) {
                        for (const Handle& hb : all_outer) {
                            Vec2d ipt;
                            if (!segmentsIntersect(mid[a].a, mid[a].b,
                                                   hb.a, hb.b, &ipt))
                                continue;
                            // 共享端点处两条曲线的首/末把手必然相交于端点本身，
                            // 这不是控制面互穿，单独标注以免掩盖真实相交。
                            const bool at_shared_ep =
                                (ipt - ci.segs.front().ctrl[0]).norm() < 1e-6 ||
                                (ipt - ci.segs.back().ctrl[3]).norm() < 1e-6;
                            printf("      把手相交%s: 内 %s × 外 %s at (%.4f,%.4f)\n",
                                   at_shared_ep ? "(共享端点,可忽略)" : "",
                                   mid[a].label.c_str(), hb.label.c_str(),
                                   ipt.x(), ipt.y());
                        }
                    }
                }
            }
        }
    }
    return 0;
}
