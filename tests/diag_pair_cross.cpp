// Standalone diagnostic: locate exactly where two generated connectivity
// curves cross, and how far each crossing sits from the shared endpoint.
//
// Motivation: "same-cluster non-endpoint non-intersection" is reported as a
// single boolean by diag_all_violations, which cannot distinguish a real
// mid-curve crossing from a near-tangential overlap that hugs a shared
// entry/exit endpoint. When a straight-through and a three-segment U-turn
// share an exit lane, the U-turn's straight lead-out lies almost on top of the
// straight-through's approach, so the crossing station relative to the shared
// endpoint is the deciding measurement.
//
// Usage: diag_pair_cross <data.json> <conn_a> <conn_b> [<conn_a> <conn_b> ...]
#include "curve/curve_utils.h"
#include "intersection_shape_generator.h"
#include "io/iodata_json.h"
#include "toolkits/toolkits.h"
#include "utils.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <string>
#include <vector>

using namespace isg;

namespace {

const BezierCurve* findCurve(const IntersectionOutput& output, const ConnId& id) {
    for (const auto& cc : output.connectivity_curves)
        if (cc.id == id && cc.curve)
            return &*cc.curve;
    return nullptr;
}

const Connectivity* findConn(const IntersectionInput& input, const ConnId& id) {
    for (const auto& conn : input.connectivities)
        if (conn.id == id)
            return &conn;
    return nullptr;
}

// Cumulative arc length of a polyline up to vertex i, used to express a
// crossing as a station along each curve instead of a raw sample index.
std::vector<double> cumulativeLength(const std::vector<Vec2d>& pts) {
    std::vector<double> s(pts.size(), 0.0);
    for (std::size_t i = 1; i < pts.size(); ++i)
        s[i] = s[i - 1] + (pts[i] - pts[i - 1]).norm();
    return s;
}

struct SegHit {
    Vec2d point;
    double s_a = 0.0;
    double s_b = 0.0;
    // 交点处两条曲线的切向夹角余弦绝对值：接近 1 表示近平行的汇流贴行，
    // 明显小于 1 才是真正的横穿。
    double abs_cos_tan = 0.0;
};

// Proper segment-segment intersection; parallel/collinear pairs are skipped
// because an exactly overlapping straight lead is reported by the adherence
// metric below rather than as a crossing.
bool segmentsCross(const Vec2d& a0, const Vec2d& a1, const Vec2d& b0,
                   const Vec2d& b1, double& t, double& u) {
    const Vec2d r = a1 - a0;
    const Vec2d s = b1 - b0;
    const double denom = cross2d(r, s);
    if (std::abs(denom) < 1e-12)
        return false;
    t = cross2d(b0 - a0, s) / denom;
    u = cross2d(b0 - a0, r) / denom;
    return t >= 0.0 && t <= 1.0 && u >= 0.0 && u <= 1.0;
}

void report(const IntersectionInput& input, const IntersectionOutput& output,
            const ConnId& id_a, const ConnId& id_b) {
    printf("\n==== %s <-> %s ====\n", id_a.c_str(), id_b.c_str());
    const Connectivity* ca = findConn(input, id_a);
    const Connectivity* cb = findConn(input, id_b);
    const BezierCurve* a = findCurve(output, id_a);
    const BezierCurve* b = findCurve(output, id_b);
    if (!ca || !cb || !a || !b) {
        printf("  missing connectivity or curve\n");
        return;
    }
    printf("  %s: entry=%s exit=%s segs=%d arc=%.3f\n", id_a.c_str(),
           ca->entry_lane_id.c_str(), ca->exit_lane_id.c_str(),
           a->numSegments(), a->arcLength());
    printf("  %s: entry=%s exit=%s segs=%d arc=%.3f\n", id_b.c_str(),
           cb->entry_lane_id.c_str(), cb->exit_lane_id.c_str(),
           b->numSegments(), b->arcLength());
    const bool shared_entry = ca->entry_lane_id == cb->entry_lane_id;
    const bool shared_exit = ca->exit_lane_id == cb->exit_lane_id;
    printf("  shared_entry=%d shared_exit=%d\n", (int)shared_entry,
           (int)shared_exit);

    // 控制点与端点切向：判断"弯曲直行 vs 掉头"这类交叉时，必须能看到直行
    // 曲线到底把控制点甩到了掉头弧的哪一侧。
    auto dumpCtrl = [](const ConnId& id, const BezierCurve& c) {
        printf("  %s ctrl:", id.c_str());
        for (int s = 0; s < c.numSegments(); ++s) {
            const BezierSegment& seg = c.segs[(std::size_t)s];
            for (int k = (s == 0 ? 0 : 1); k < 4; ++k)
                printf(" (%.3f,%.3f)", seg.ctrl[k].x(), seg.ctrl[k].y());
        }
        printf("\n");
        const Vec2d t_start = c.startTan().normalized();
        const Vec2d t_end = c.endTan().normalized();
        printf("  %s startTan=(%.4f,%.4f) endTan=(%.4f,%.4f) maxK=%.4f\n",
               id.c_str(), t_start.x(), t_start.y(), t_end.x(), t_end.y(),
               c.maxCurvature(60));
    };
    dumpCtrl(id_a, *a);
    dumpCtrl(id_b, *b);

    for (double ep : {0.01, 0.30, 1.00, 1.50, 3.00}) {
        printf("  curvesIntersectBusiness(ep=%.2f) = %d\n", ep,
               (int)curvesIntersectBusiness(*a, *b, ep));
    }

    const std::vector<Vec2d> pa = a->sampleByArcLength(400);
    const std::vector<Vec2d> pb = b->sampleByArcLength(400);
    const std::vector<double> sa = cumulativeLength(pa);
    const std::vector<double> sb = cumulativeLength(pb);
    const double len_a = sa.empty() ? 0.0 : sa.back();
    const double len_b = sb.empty() ? 0.0 : sb.back();

    std::vector<SegHit> hits;
    for (std::size_t i = 0; i + 1 < pa.size(); ++i) {
        for (std::size_t j = 0; j + 1 < pb.size(); ++j) {
            double t = 0.0, u = 0.0;
            if (!segmentsCross(pa[i], pa[i + 1], pb[j], pb[j + 1], t, u))
                continue;
            SegHit hit;
            hit.point = pa[i] + t * (pa[i + 1] - pa[i]);
            hit.s_a = sa[i] + t * (sa[i + 1] - sa[i]);
            hit.s_b = sb[j] + u * (sb[j + 1] - sb[j]);
            const Vec2d da = pa[i + 1] - pa[i];
            const Vec2d db = pb[j + 1] - pb[j];
            if (da.norm() > 1e-12 && db.norm() > 1e-12) {
                hit.abs_cos_tan =
                    std::abs(da.normalized().dot(db.normalized()));
            }
            // Collapse the chain of sample-level hits produced by a single
            // geometric crossing (or by a long near-tangential overlap).
            if (!hits.empty() && (hit.point - hits.back().point).norm() < 0.05)
                continue;
            hits.push_back(hit);
        }
    }
    printf("  raw crossings: %zu\n", hits.size());
    for (const SegHit& hit : hits) {
        const double d_ends = std::min({
            (hit.point - a->startPt()).norm(),
            (hit.point - a->endPt()).norm(),
            (hit.point - b->startPt()).norm(),
            (hit.point - b->endPt()).norm()});
        printf("    at (%.3f,%.3f) s_%s=%.3f/%.3f s_%s=%.3f/%.3f "
               "dist_to_nearest_endpoint=%.3f |cos_tan|=%.6f\n",
               hit.point.x(), hit.point.y(), id_a.c_str(), hit.s_a, len_a,
               id_b.c_str(), hit.s_b, len_b, d_ends, hit.abs_cos_tan);
    }
    if (!hits.empty()) {
        double min_d = 1e30, max_d = 0.0;
        for (const SegHit& hit : hits) {
            const double d = std::min({
                (hit.point - a->startPt()).norm(),
                (hit.point - a->endPt()).norm(),
                (hit.point - b->startPt()).norm(),
                (hit.point - b->endPt()).norm()});
            min_d = std::min(min_d, d);
            max_d = std::max(max_d, d);
        }
        printf("  endpoint distance range: %.3f .. %.3f\n", min_d, max_d);
    }

    // 共享端点附近的横向分离剖面：以共享端点为原点、该端点切向为轴，打印
    // 两条曲线各自的横向偏离。两条汇入同一出口车道的曲线在端点必然一阶贴
    // 合，只有横向分离随站位增长才谈得上"仅端点相交"；若分离长期停留在厘米
    // 级，则是事实上的重合贴行，交点只是这条窄缝内符号翻转的位置。
    if (shared_exit) {
        const Vec2d origin = a->endPt();
        Vec2d tan = a->endTan();
        if (tan.norm() > 1e-8) {
            tan.normalize();
            const Vec2d lat(-tan.y(), tan.x());
            printf("  lateral separation profile about shared exit endpoint "
                   "(axis = -exit tangent):\n");
            printf("    %8s %10s %10s %10s\n", "station", id_a.c_str(),
                   id_b.c_str(), "gap");
            for (double station : {0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0,
                                   7.93, 9.0, 10.0, 12.0}) {
                // 按"自终点回溯的弧长"定位，而不是按轴向投影：U 型调头会多次
                // 经过同一轴向站位，投影匹配会在中弧与出口直段之间乱跳。
                auto lateralAt = [&](const std::vector<Vec2d>& pts,
                                     const std::vector<double>& cum) {
                    const double total = cum.empty() ? 0.0 : cum.back();
                    if (total < station)
                        return std::numeric_limits<double>::quiet_NaN();
                    double best_err = 1e30;
                    double best_lat = 0.0;
                    for (std::size_t i = 0; i < pts.size(); ++i) {
                        const double from_end = total - cum[i];
                        const double err = std::abs(from_end - station);
                        if (err < best_err) {
                            best_err = err;
                            best_lat = (pts[i] - origin).dot(lat);
                        }
                    }
                    return best_lat;
                };
                const double la = lateralAt(pa, sa);
                const double lb = lateralAt(pb, sb);
                printf("    %8.2f %10.4f %10.4f %10.4f\n", station, la, lb,
                       la - lb);
            }
        }
    }
    // 共享进入端点附近的横向分离剖面：与出口侧对称，用于"弯曲直行 + 同入
    // 掉头"这一类。轴取进入切向，站位自起点向前量取；掉头的横向偏离会随中弧
    // 迅速拉开，而弯曲直行的横向偏离方向决定它会不会追进掉头弧的内侧。
    if (shared_entry) {
        const Vec2d origin = a->startPt();
        Vec2d tan = a->startTan();
        if (tan.norm() > 1e-8) {
            tan.normalize();
            const Vec2d lat(-tan.y(), tan.x());
            printf("  lateral separation profile about shared entry endpoint "
                   "(axis = +entry tangent):\n");
            printf("    %8s %10s %10s %10s\n", "station", id_a.c_str(),
                   id_b.c_str(), "gap");
            for (double station : {0.5, 1.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0,
                                   14.0, 16.0, 18.0, 20.0}) {
                auto lateralAt = [&](const std::vector<Vec2d>& pts,
                                     const std::vector<double>& cum) {
                    const double total = cum.empty() ? 0.0 : cum.back();
                    if (total < station)
                        return std::numeric_limits<double>::quiet_NaN();
                    double best_err = 1e30;
                    double best_lat = 0.0;
                    for (std::size_t i = 0; i < pts.size(); ++i) {
                        const double err = std::abs(cum[i] - station);
                        if (err < best_err) {
                            best_err = err;
                            best_lat = (pts[i] - origin).dot(lat);
                        }
                    }
                    return best_lat;
                };
                const double la = lateralAt(pa, sa);
                const double lb = lateralAt(pb, sb);
                printf("    %8.2f %10.4f %10.4f %10.4f\n", station, la, lb,
                       la - lb);
            }
        }
    }
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 4 || (argc % 2) != 0) {
        fprintf(stderr, "Usage: %s <data.json> <conn_a> <conn_b> [...]\n",
                argv[0]);
        return 1;
    }
    const std::string path = argv[1];
    IntersectionInput input = IntersectionIO::loadFromFile(path);
    IntersectionShapeGenerator gen;
    IntersectionOutput output;
    if (!gen.generate(input, output)) {
        fprintf(stderr, "Generation failed\n");
        return 1;
    }
    // Generation normalizes its own copy; re-apply the same two preprocessing
    // steps here so lane ids/endpoints line up with the generated curves.
    input = InputNormalizer(input);
    ConnectivityDirectionNormalizer(input, ConnectivityDirectionConfig());
    fprintf(stderr, "Loaded %s: %zu conns, %zu curves\n", path.c_str(),
            input.connectivities.size(), output.connectivity_curves.size());
    for (int i = 2; i + 1 < argc; i += 2)
        report(input, output, ConnId(argv[i]), ConnId(argv[i + 1]));
    return 0;
}
