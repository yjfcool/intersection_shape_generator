// Standalone diagnostic: shared-endpoint fan-out lateral ordering.
//
// Motivation: a single entry lane can feed 10+ exit lanes. All those curves
// share P0 and the entry tangent, so near the shared endpoint they separate
// only by their turn-in rate (initial curvature). "Same-cluster non-endpoint
// non-intersection" for such a family is a *total order* property: at every
// station the members must appear in the same left-to-right order as their
// final turn angle. Pairwise handle-length comparisons cannot express that,
// which is why pairwise repair oscillates.
//
// This tool reports, for one shared-endpoint family:
//   - the canonical fan order key (signed total turn angle, then turn radius)
//   - each member's handles, initial curvature and lateral profile
//   - the actual lateral order at a ladder of stations vs. the canonical one
//   - a per-member (h0,h1) grid sweep answering "does a monotone, non-crossing
//     assignment exist inside the shape gates at all?"
//
// Usage: diag_fanout_order <data.json> <conn_id> <conn_id> [...]
#include "curve/bezier.h"
#include "curve/curve_utils.h"
#include "intersection_shape_generator.h"
#include "io/iodata_json.h"
#include "utils.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

using namespace isg;

namespace {

struct Member {
    ConnId id;
    Vec2d p0{0, 0};
    Vec2d t0{1, 0};
    Vec2d p1{0, 0};
    Vec2d t1{1, 0};
    double h0 = 0.0;
    double h1 = 0.0;
    double turn_deg = 0.0;      ///< 自进入切向到退出切向的有符号转角
    double far_lat = 0.0;       ///< 远端点在远端切向坐标系里的有符号横向位置
    double kappa0 = 0.0;        ///< 共享端点处初始曲率
    BezierCurve curve;
};

const ConnectivityCurve* findCurve(const IntersectionOutput& out, const ConnId& id) {
    for (const auto& cc : out.connectivity_curves)
        if (cc.id == id && cc.curve)
            return &cc;
    return nullptr;
}

double signedTurnDeg(const Vec2d& t0, const Vec2d& t1) {
    return std::atan2(cross2d(t0, t1), t0.dot(t1)) * 180.0 / M_PI;
}

// 单段三次 Bezier 在 t=0 处的有符号曲率：lat(s) ≈ kappa0*s^2/2，
// 因此共享端点附近的左右顺序完全由该量决定。
double signedKappa0(const BezierCurve& c) {
    if (c.empty())
        return 0.0;
    const auto& g = c.segs.front().ctrl;
    const Vec2d d1 = g[1] - g[0];
    const Vec2d d2 = g[2] - g[1];
    const double h = d1.norm();
    if (h < 1e-9)
        return 0.0;
    return (2.0 / 3.0) * cross2d(d1, d2) / (h * h * h);
}

BezierCurve makeSingleCubic(const Member& m, double h0, double h1) {
    BezierSegment seg;
    seg.ctrl[0] = m.p0;
    seg.ctrl[1] = m.p0 + m.t0 * h0;
    seg.ctrl[2] = m.p1 - m.t1 * h1;
    seg.ctrl[3] = m.p1;
    BezierCurve c;
    c.segs.push_back(seg);
    return c;
}

// 复刻 isNonUTurnTurnShapeAcceptable 的关键门禁（弧弦比与曲率上限），
// 供可行性搜索使用；诊断工具不链接会话内部的静态函数。
bool shapeGateOk(const BezierCurve& c, double chord_len, double turn_strength) {
    if (chord_len < 1e-6 || turn_strength <= 0.35)
        return true;
    const double arc_chord = c.arcLength() / chord_len;
    const double min_arc = chord_len < 12.0 ? 1.02 : 1.06;
    if (arc_chord < min_arc || arc_chord > 1.35)
        return false;
    return c.maxCurvature(20) <= 2.5;
}

double lateralAt(const BezierCurve& c, const Vec2d& origin, const Vec2d& lat,
                 double station) {
    const std::vector<Vec2d> pts = c.sampleByArcLength(400);
    double cum = 0.0;
    for (std::size_t i = 1; i < pts.size(); ++i) {
        const double step = (pts[i] - pts[i - 1]).norm();
        if (cum + step >= station) {
            const double f = step > 1e-12 ? (station - cum) / step : 0.0;
            const Vec2d p = pts[i - 1] + f * (pts[i] - pts[i - 1]);
            return (p - origin).dot(lat);
        }
        cum += step;
    }
    return std::numeric_limits<double>::quiet_NaN();
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 4) {
        fprintf(stderr,
                "Usage: diag_fanout_order <data.json> <conn_id> <conn_id> [...]\n");
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
        fprintf(stderr, "GENERATION FAILED\n");
        return 1;
    }

    std::vector<Member> ms;
    for (const ConnId& id : ids) {
        const ConnectivityCurve* cc = findCurve(output, id);
        if (!cc) {
            printf("conn %s: missing curve\n", id.c_str());
            continue;
        }
        Member m;
        m.id = id;
        m.curve = *cc->curve;
        m.p0 = m.curve.startPt();
        m.p1 = m.curve.endPt();
        m.t0 = m.curve.startTan().normalized();
        m.t1 = m.curve.endTan().normalized();
        const auto& gf = m.curve.segs.front().ctrl;
        const auto& gb = m.curve.segs.back().ctrl;
        m.h0 = (gf[1] - gf[0]).norm();
        m.h1 = (gb[3] - gb[2]).norm();
        m.turn_deg = signedTurnDeg(m.t0, m.t1);
        m.far_lat = cross2d(m.t1, m.p1 - m.p0);
        m.kappa0 = signedKappa0(m.curve);
        ms.push_back(m);
    }
    if (ms.size() < 2) {
        fprintf(stderr, "need at least 2 resolvable members\n");
        return 2;
    }

    // 共享端点与其切向：以第一个成员的进入端为参考系。
    const Vec2d origin = ms.front().p0;
    const Vec2d axis = ms.front().t0;
    const Vec2d lat(-axis.y(), axis.x());
    printf("==== %s : shared-entry fan-out family ====\n", path.c_str());
    printf("shared endpoint (%.3f,%.3f) tangent (%.4f,%.4f)\n", origin.x(),
           origin.y(), axis.x(), axis.y());
    for (const Member& m : ms) {
        printf("  %-3s segs=%d h0=%8.3f h1=%8.3f arc=%7.3f chord=%7.3f "
               "turn=%+8.2fdeg far_lat=%+8.3f kappa0=%+.6f startPt_off=%.3f\n",
               m.id.c_str(), m.curve.numSegments(), m.h0, m.h1,
               m.curve.arcLength(), (m.p1 - m.p0).norm(), m.turn_deg,
               m.far_lat, m.kappa0, (m.p0 - origin).norm());
    }

    // 规范扇出顺序：先按有符号转角（左转角越大越靠左）；转角差在 10 度内
    // 视为同 arm，再按"远端点在远端切向坐标系里的有符号横向位置"降序——远端
    // 更靠左的成员在共享端点侧也必须更靠左，否则左右次序在两端之间必然翻转
    // 一次，强制产生一个非端点交点。
    std::vector<const Member*> canonical;
    for (const Member& m : ms)
        canonical.push_back(&m);
    std::sort(canonical.begin(), canonical.end(),
              [](const Member* a, const Member* b) {
                  if (std::abs(a->turn_deg - b->turn_deg) > 10.0)
                      return a->turn_deg > b->turn_deg;
                  return a->far_lat > b->far_lat;
              });
    printf("\ncanonical left-to-right order (turn angle desc, then far lateral "
           "desc):\n  ");
    for (const Member* m : canonical)
        printf("%s ", m->id.c_str());
    printf("\n");

    printf("\nactual lateral offsets about the shared endpoint "
           "(+ = left of entry tangent):\n    %8s", "station");
    for (const Member& m : ms)
        printf(" %9s", m.id.c_str());
    printf("   order(left->right)\n");
    for (double station : {1.0, 2.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 40.0}) {
        printf("    %8.2f", station);
        std::vector<std::pair<double, ConnId>> row;
        for (const Member& m : ms) {
            const double v = lateralAt(m.curve, origin, lat, station);
            printf(" %9.4f", v);
            if (std::isfinite(v))
                row.emplace_back(-v, m.id);
        }
        std::sort(row.begin(), row.end());
        printf("   ");
        for (const auto& r : row)
            printf("%s ", r.second.c_str());
        printf("\n");
    }

    printf("\npairwise: canonical relation vs actual crossing\n");
    for (std::size_t i = 0; i < ms.size(); ++i) {
        for (std::size_t j = i + 1; j < ms.size(); ++j) {
            const Member& a = ms[i];
            const Member& b = ms[j];
            const bool cross = curvesIntersectBusiness(a.curve, b.curve, 1.5);
            const bool poly = sharedEndpointControlPolylinesCross(
                a.curve, b.curve, 0.30);
            std::size_t ra = 0, rb = 0;
            for (std::size_t k = 0; k < canonical.size(); ++k) {
                if (canonical[k]->id == a.id) ra = k;
                if (canonical[k]->id == b.id) rb = k;
            }
            const ConnId& should_be_left = ra < rb ? a.id : b.id;
            printf("  %-3s|%-3s cross=%d polyCross=%d canonical_left=%s "
                   "kappa0=%+.6f/%+.6f\n",
                   a.id.c_str(), b.id.c_str(), (int)cross, (int)poly,
                   should_be_left.c_str(), a.kappa0, b.kappa0);
        }
    }

    // 可行性搜索：对每个成员在 (h0,h1) 网格内枚举形态门禁内的候选，再按
    // 规范顺序逐个贪心挑选与已选成员全部不相交的候选。若能全部落位，说明
    // 这一族的非交叉解在单段表达内确实存在，问题出在排序而非几何不可能。
    printf("\nfeasibility sweep (single-cubic, shape gates enforced)\n");
    std::vector<BezierCurve> chosen;
    std::vector<ConnId> chosen_ids;
    bool all_ok = true;
    for (const Member* mp : canonical) {
        const Member& m = *mp;
        const Vec2d chord = m.p1 - m.p0;
        const double chord_len = chord.norm();
        const double turn_strength =
            std::abs(cross2d(m.t0, chord.normalized()));
        const OrdinarySingleCubicHandleBounds bounds =
            ordinarySingleCubicHandleBounds(m.p0, m.t0, m.p1, m.t1, true);
        bool placed = false;
        double best_h0 = 0.0, best_h1 = 0.0, best_score = 1e30;
        for (int i0 = 0; i0 <= 48 && !placed; ++i0) {
            const double f0 = 0.02 + 0.98 * i0 / 48.0;
            const double h0 = std::max(bounds.start_min, f0 * bounds.start_max);
            for (int i1 = 0; i1 <= 48; ++i1) {
                const double f1 = 0.02 + 0.98 * i1 / 48.0;
                const double h1 = std::max(bounds.end_min, f1 * bounds.end_max);
                BezierCurve cand = makeSingleCubic(m, h0, h1);
                if (!ordinarySingleCubicControlsValid(
                        cand, m.p0, m.t0, m.p1, m.t1, 1e-5, true))
                    continue;
                if (curveSelfIntersectsBusiness(cand, 1.0))
                    continue;
                if (!shapeGateOk(cand, chord_len, turn_strength))
                    continue;
                bool ok = true;
                for (const BezierCurve& c : chosen) {
                    if (curvesIntersectBusiness(cand, c, 1.5) ||
                        sharedEndpointControlPolylinesCross(cand, c, 0.30)) {
                        ok = false;
                        break;
                    }
                }
                if (!ok)
                    continue;
                const double score =
                    std::abs(cand.arcLength() / chord_len - 1.10) +
                    cand.maxCurvature(40);
                if (score < best_score) {
                    best_score = score;
                    best_h0 = h0;
                    best_h1 = h1;
                    placed = true;
                }
            }
        }
        if (!placed) {
            printf("  %-3s : NO feasible (h0,h1) against already-placed {",
                   m.id.c_str());
            for (const ConnId& cid : chosen_ids)
                printf("%s ", cid.c_str());
            printf("}\n");
            all_ok = false;
            continue;
        }
        printf("  %-3s : h0=%8.3f h1=%8.3f (bounds start<=%.3f end<=%.3f) "
               "score=%.4f\n",
               m.id.c_str(), best_h0, best_h1, bounds.start_max,
               bounds.end_max, best_score);
        chosen.push_back(makeSingleCubic(m, best_h0, best_h1));
        chosen_ids.push_back(m.id);
    }
    printf("  => monotone non-crossing assignment %s\n",
           all_ok ? "EXISTS" : "NOT FOUND by this greedy sweep");
    return 0;
}
