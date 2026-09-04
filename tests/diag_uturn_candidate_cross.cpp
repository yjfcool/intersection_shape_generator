// 掉头三段式候选的同簇交叉归因诊断。
//
// SegmentedUTurnCandidateSearch 只统计 rejected_sibling_cross 的个数，不记录被
// 谁挡住。本工具用真实家族参数（家族 lead 下限、平齐站位、familyLateralLadder
// 分档量）重建 extra=0/0 的各 arc_alpha 候选，逐个点名它与哪些同簇兄弟发生
// curvesHaveForbiddenSameClusterIntersection，并打印候选自身的形态读数，用来
// 区分"形态不可行"与"形态可行但被某条兄弟曲线挡住"。
//
// 用法: diag_uturn_candidate_cross <data.json> <conn_id> [entry_bias exit_bias]
#include "constraints/cluster_order.h"
#include "constraints/boundary_safety.h"
#include "curve/curve_utils.h"
#include "generation/uturn_shape.h"
#include "initialization/uturn_curve_initializer.h"
#include "intersection_shape_generator.h"
#include "io/iodata_json.h"
#include "preprocessing/uturn_family_builder.h"
#include "toolkits/toolkits.h"
#include "types.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <unordered_map>
#include <vector>

using namespace isg;

namespace {

constexpr double kClusterEndpointTol = 1.5;

double envOverride(const char* name, double fallback) {
    const char* value = std::getenv(name);
    if (!value || *value == '\0')
        return fallback;
    char* end = nullptr;
    const double parsed = std::strtod(value, &end);
    return end != value && *end == '\0' && std::isfinite(parsed)
        ? parsed : fallback;
}

std::vector<std::string> rawBoundaryHits(
    const BezierCurve& curve, const IntersectionInput& input) {
    std::vector<std::string> hits;
    const bool print_contacts = std::getenv("ISG_DIAG_CONTACT") != nullptr;
    if (curve.empty())
        return hits;
    const std::vector<Vec2d> sampled = curve.sampleByArcLength(std::max(
        64, std::min(240, (int)std::ceil(curve.arcLength() / 0.18) + 1)));
    for (const auto& boundary : input.boundaries) {
        const auto& points = boundary.geometry.points;
        bool hit = false;
        for (int i = 0; i + 1 < (int)sampled.size() && !hit; ++i) {
            for (int j = 0; j + 1 < (int)points.size(); ++j) {
                Vec2d intersection;
                if (segmentsIntersect(sampled[i], sampled[i + 1],
                                      points[j], points[j + 1], &intersection) &&
                    segmentHasForbiddenBoundaryContact(
                        sampled[i], sampled[i + 1], points[j], points[j + 1],
                        sampled.front(), sampled.back())) {
                    hit = true;
                    if (print_contacts) {
                        printf("\n    contact %s: (%.6f,%.6f) curve_end_dist=%.6f/%.6f\n",
                               boundary.id.c_str(), intersection.x(), intersection.y(),
                               (intersection - sampled.front()).norm(),
                               (intersection - sampled.back()).norm());
                    }
                    break;
                }
            }
        }
        if (hit)
            hits.push_back(boundary.id);
    }
    return hits;
}

bool pureGeometricInput(const IntersectionInput& input) {
    return input.obstacles.empty() && input.boundaries.empty() &&
           input.crosswalks.empty() && input.stop_lines.empty();
}

// 打印两条曲线采样折线的所有交点，坐标折算到 (沿轴回退距离, 有符号横向) 口径，
// 便于直接对照 buildSegmented 里的 q0/q1 站位读数。
void printCrossings(const BezierCurve& a, const BezierCurve& b,
                    const Vec2d& origin, const Vec2d& axis,
                    const Vec2d& lateral, double side) {
    const std::vector<Vec2d> pa = a.sample(240);
    const std::vector<Vec2d> pb = b.sample(240);
    for (std::size_t i = 0; i + 1 < pa.size(); ++i) {
        for (std::size_t j = 0; j + 1 < pb.size(); ++j) {
            Vec2d hit;
            if (!segmentsIntersect(pa[i], pa[i + 1], pb[j], pb[j + 1], &hit))
                continue;
            const double back = (hit - origin).dot(axis);
            const double lat = (hit - origin).dot(lateral) * side;
            printf("      hit at back_axis=%+.3f lat=%+.4f  (a_t=%.3f b_t=%.3f)"
                   " dist_to_a_ends=%.3f/%.3f dist_to_b_ends=%.3f/%.3f\n",
                   back, lat,
                   static_cast<double>(i) / (pa.size() - 1),
                   static_cast<double>(j) / (pb.size() - 1),
                   (hit - pa.front()).norm(), (hit - pa.back()).norm(),
                   (hit - pb.front()).norm(), (hit - pb.back()).norm());
            return;
        }
    }
}

}  // namespace

namespace {

// 复刻 connectivity_generation_session.cpp 里的静态判定，用来区分"相交被禁止"
// 究竟出自结构性豁免不成立，还是出自非端点贴行否决了豁免。
double cross2dLocal(const Vec2d& u, const Vec2d& v) {
    return u.x() * v.y() - u.y() * v.x();
}

bool looksLeftRightTurn(const BezierCurve& curve) {
    if (curve.empty() || curveLooksUTurnForClusterExemption(curve))
        return false;
    const Vec2d chord = curve.endPt() - curve.startPt();
    const Vec2d st = curve.startTan();
    if (chord.norm() < 1e-8 || st.norm() < 1e-8)
        return false;
    return std::abs(cross2dLocal(st.normalized(), chord.normalized())) > 0.35;
}

// 返回最远离全部端点的贴行见证点距离；<0 表示不存在贴行。
double adherenceWitnessDistance(const BezierCurve& a, const BezierCurve& b,
                                double adherent_tol = 0.18) {
    const std::vector<Vec2d> pa = a.sample(64);
    const std::vector<Vec2d> pb = b.sample(64);
    double worst = -1.0;
    for (std::size_t ai = 0; ai + 1 < pa.size(); ++ai) {
        const Vec2d a0 = pa[ai], a1 = pa[ai + 1];
        const Vec2d ad = a1 - a0;
        const double alen = ad.norm();
        if (alen < 1e-8) continue;
        const Vec2d au = ad / alen;
        for (std::size_t bi = 0; bi + 1 < pb.size(); ++bi) {
            const Vec2d b0 = pb[bi], b1 = pb[bi + 1];
            const Vec2d bd = b1 - b0;
            const double blen = bd.norm();
            if (blen < 1e-8) continue;
            const Vec2d bu = bd / blen;
            if (std::abs(au.dot(bu)) < 0.96) continue;
            if (std::abs(cross2dLocal(au, b0 - a0)) > adherent_tol ||
                std::abs(cross2dLocal(au, b1 - a0)) > adherent_tol)
                continue;
            const double b0s = (b0 - a0).dot(au);
            const double b1s = (b1 - a0).dot(au);
            const double lo = std::max(0.0, std::min(b0s, b1s));
            const double hi = std::min(alen, std::max(b0s, b1s));
            if (lo > hi + adherent_tol) continue;
            const Vec2d w = a0 + 0.5 * (lo + hi) * au;
            double d = std::min(std::min((w - a.startPt()).norm(),
                                         (w - a.endPt()).norm()),
                                std::min((w - b.startPt()).norm(),
                                         (w - b.endPt()).norm()));
            worst = std::max(worst, d);
        }
    }
    return worst;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 3) {
        fprintf(stderr,
                "用法: %s <data.json> <conn_id> [entry_bias exit_bias]\n",
                argv[0]);
        return 2;
    }
    IntersectionInput input = IntersectionIO::loadFromFile(argv[1]);
    input = InputNormalizer(input);
    ConnectivityDirectionNormalizer(input, ConnectivityDirectionConfig());
    const ConnId target(argv[2]);
    const double entry_bias = argc > 3 ? std::atof(argv[3]) : 0.0;
    const double exit_bias = argc > 4 ? std::atof(argv[4]) : 0.0;
    const double middle_handle_scale0 = envOverride(
        "ISG_DIAG_MIDDLE_HANDLE0", 1.0);
    const double middle_handle_scale1 = envOverride(
        "ISG_DIAG_MIDDLE_HANDLE1", 1.0);
    const bool scan = std::getenv("ISG_DIAG_SCAN") != nullptr;
    const bool scan_bias_only = std::getenv("ISG_DIAG_SCAN_BIAS_ONLY") != nullptr;
    // 可选：直接覆盖家族分档量，用来观察"若家族给本成员另一个档位会怎样"，
    // 绕开 buildSegmented 内部对反向偏置的槽位钳制。
    const bool override_stagger = argc > 6;
    const double stagger0_override = argc > 5 ? std::atof(argv[5]) : 0.0;
    const double stagger1_override = argc > 6 ? std::atof(argv[6]) : 0.0;

    const Connectivity* conn = nullptr;
    for (const auto& c : input.connectivities)
        if (c.id == target)
            conn = &c;
    if (conn == nullptr) {
        fprintf(stderr, "conn %s not found\n", target.c_str());
        return 1;
    }

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
    const UTurnFamilyInfo family =
        builder.build(*conn, input, scope, &solver, false);
    const std::pair<Vec2d, Vec2d> entry =
        input.entryPtDir(conn->entry_lane_id);
    const std::pair<Vec2d, Vec2d> exit =
        input.exitPtDir(conn->exit_lane_id);
    const double station_offset = envOverride(
        "ISG_DIAG_STATION_OFFSET", 0.0);
    double min_lead0 = family.lead0;
    double min_lead1 = family.lead1;
    const double target_station = std::isfinite(family.aligned_station)
        ? family.aligned_station + station_offset
        : family.aligned_station;
    builder.enforceAlignmentStation(
        entry.first, entry.second, exit.first, exit.second,
        target_station, min_lead0, min_lead1);
    const bool physical = !pureGeometricInput(input);
    const double step = physical && family.rank.family_size >= 3 ? 0.25 : 0.01;
    const UTurnFamilyLadder ladder = builder.familyLateralLadder(
        *conn, input, scope, step, &solver, false);

    const double station = family.rank.family_size > 1
        ? target_station : std::numeric_limits<double>::quiet_NaN();
    printf("==== %s ====\n", target.c_str());
    printf("  family size=%zu rank=%zu radius=%.4f station=%.4f\n",
           family.rank.family_size, family.rank.reverse_radius_rank,
           family.radius, station);
    const std::vector<const Connectivity*> family_members =
        builder.alignmentComponent(*conn, input, scope, &solver, false);
    printf("  family members:");
    for (const Connectivity* member : family_members)
        printf(" %s", member->id.c_str());
    printf("\n");
    printf("  lead floors: %.4f/%.4f  ladder: step=%.4f entry=%.4f exit=%.4f\n",
           family.lead0, family.lead1, ladder.step,
           ladder.entry_stagger, ladder.exit_stagger);
    printf("  bias: entry=%.4f exit=%.4f\n", entry_bias, exit_bias);

    // 同簇受约束的兄弟集合（结构性豁免对不参与）。
    std::vector<ConnId> partners;
    for (const auto& p : solver.pairs()) {
        if (p.exempt == CrossExemption::StructuralCross)
            continue;
        if (p.id_a == target) partners.push_back(p.id_b);
        else if (p.id_b == target) partners.push_back(p.id_a);
    }
    printf("  constrained partners: %zu\n", partners.size());

    const double max_curvature = segmentedUTurnMaxCurvatureLimit(
        entry.first, entry.second, exit.first, exit.second);
    Vec2d T0 = entry.second.normalized();
    Vec2d axis = (T0 - exit.second.normalized()).normalized();
    if (axis.dot(T0) < 0.0) axis = -axis;

    const std::vector<double> arc_alphas =
        {2.0 / 3.0, 0.75, 0.85, 1.0, 1.25, 1.50, 1.75, 2.0,
         2.5, 3.0, 4.0, 5.0, 0.58, 0.50, 0.38, 0.28, 0.16};
    const Vec2d lateral_ref{-axis.y(), axis.x()};
    const double side = (exit.first - entry.first).dot(lateral_ref) >= 0.0
        ? 1.0 : -1.0;
    const std::vector<double> scan_biases =
        {-2.0, -1.5, -1.0, -0.70, -0.35, 0.0, 0.35, 0.70, 1.0, 1.5, 2.0};
    const std::vector<double> scan_station_offsets =
        {-1.0, 0.0, 0.75, 1.5, 2.5, 3.5};
    const std::vector<double> scan_handle_scales =
        {0.50, 0.75, 1.0, 1.25, 1.50};
    const std::vector<double>& use_alphas = scan ? arc_alphas : arc_alphas;
    for (double alpha : use_alphas) {
      const std::vector<double>& entry_biases = scan
          ? scan_biases : std::vector<double>{entry_bias};
      const std::vector<double>& exit_biases = scan
          ? scan_biases : std::vector<double>{exit_bias};
      const std::vector<double> default_station_offsets = {station_offset};
      const std::vector<double>& station_offsets = scan && !scan_bias_only
          ? scan_station_offsets : default_station_offsets;
      const std::vector<double>& handle_scales = scan && !scan_bias_only
          ? scan_handle_scales : std::vector<double>{1.0};
      for (double use_entry_bias : entry_biases) {
       for (double use_exit_bias : exit_biases) {
        for (double use_station_offset : station_offsets) {
         for (double use_handle_scale0 : handle_scales) {
          for (double use_handle_scale1 : handle_scales) {
        double scan_min_lead0 = family.lead0;
        double scan_min_lead1 = family.lead1;
        const double scan_station = std::isfinite(family.aligned_station)
            ? family.aligned_station + use_station_offset
            : family.aligned_station;
        builder.enforceAlignmentStation(
            entry.first, entry.second, exit.first, exit.second,
            scan_station, scan_min_lead0, scan_min_lead1);
        const double use_entry_stagger =
            override_stagger ? stagger0_override : ladder.entry_stagger;
        const double use_exit_stagger =
            override_stagger ? stagger1_override : ladder.exit_stagger;
        BezierCurve candidate = UTurnCurveInitializer().buildSegmented(
            entry.first, entry.second, exit.first, exit.second,
            scan_min_lead0, scan_min_lead1, alpha, 0.0, 0.0, 0.0,
            use_entry_stagger, use_exit_stagger,
            use_entry_bias, use_exit_bias, ladder.step);
        if (candidate.numSegments() == 3) {
            BezierSegment& middle = candidate.segs[1];
            const Vec2d q0 = middle.ctrl[0];
            const Vec2d q1 = middle.ctrl[3];
            Vec2d d0 = middle.ctrl[1] - q0;
            Vec2d d1 = q1 - middle.ctrl[2];
            if (d0.norm() > 1e-8)
                middle.ctrl[1] = q0 + d0.normalized() *
                    d0.norm() * (scan ? use_handle_scale0 : middle_handle_scale0);
            if (d1.norm() > 1e-8)
                middle.ctrl[2] = q1 - d1.normalized() *
                    d1.norm() * (scan ? use_handle_scale1 : middle_handle_scale1);
        }
        printf("  alpha=%.3f bias=%.2f/%.2f station_offset=%.2f handle=%.2f/%.2f segs=%d",
               alpha, use_entry_bias, use_exit_bias, use_station_offset,
               use_handle_scale0, use_handle_scale1,
               (int)candidate.numSegments());
        if (candidate.numSegments() != 3) { printf(" (not segmented)\n"); continue; }
        const BoundarySafetyResult boundary_safety =
            curveBoundarySafetyForInput(candidate, input, 128);
        const std::vector<std::string> boundary_hits =
            rawBoundaryHits(candidate, input);
        const Vec2d q0 = candidate.segs.front().ctrl[3];
        const Vec2d q1 = candidate.segs.back().ctrl[0];
        const double kmax = candidate.maxCurvature(40);
        printf(" boundary=%d outside=%d raw=",
               (int)(boundary_safety.intersects || boundary_safety.outside_road_edge),
               (int)boundary_safety.outside_road_edge);
        for (const auto& boundary_id : boundary_hits)
            printf("%s,", boundary_id.c_str());
        printf(" lead=%.3f/%.3f gap=%.4f q0lat=%+.4f q1lat=%+.4f"
               " d_station=%.4f/%.4f kmax=%.3f(<%.3f) self=%d round=%d",
               (q0 - entry.first).norm(), (exit.first - q1).norm(),
               (q1 - q0).norm(),
               (q0 - exit.first).dot(lateral_ref) * side,
               (q1 - exit.first).dot(lateral_ref) * side,
               std::isfinite(scan_station) ? q0.dot(axis) - scan_station : 0.0,
               std::isfinite(scan_station) ? q1.dot(axis) - scan_station : 0.0,
               kmax, max_curvature,
               (int)curveSelfIntersectsBusiness(candidate, 1.0),
               (int)segmentedUTurnMiddleArcLooksRound(candidate, axis));
        if (scan || std::getenv("ISG_DIAG_POINTS"))
            printf(" q0=(%.3f,%.3f) q1=(%.3f,%.3f)",
                   q0.x(), q0.y(), q1.x(), q1.y());
        printf(" crosses:");
        int count = 0;
        std::vector<ConnId> hits;
        for (const ConnId& other : partners) {
            auto it = curves.find(other);
            if (it == curves.end()) continue;
            if (curvesIntersectBusiness(
                    candidate, *it->second, kClusterEndpointTol)) {
                printf(" %s", other.c_str());
                hits.push_back(other);
                ++count;
            }
        }
        if (count == 0) printf(" none");
        printf("\n");
        if (scan && boundary_hits.empty() && count == 0 &&
            kmax < max_curvature && !curveSelfIntersectsBusiness(candidate, 1.0) &&
            candidate.arcLength() / std::max(1e-6, (entry.first - exit.first).norm()) >= 1.35) {
            printf("    SCAN_SAFE alpha=%.3f bias=%.2f/%.2f station_offset=%.2f handle=%.2f/%.2f\n",
                   alpha, use_entry_bias, use_exit_bias, use_station_offset,
                   use_handle_scale0, use_handle_scale1);
        }
        {
            const Vec2d lat_dir{-axis.y(), axis.x()};
            const double sgn =
                (exit.first - entry.first).dot(lat_dir) >= 0.0 ? 1.0 : -1.0;
            for (const ConnId& other : hits) {
                const BezierCurve& ob = *curves.find(other)->second;
                printCrossings(candidate, ob,
                               exit.first, axis, lat_dir, sgn);
                printf("        %s: uturn(self)=%d leftright(other)=%d"
                       " uturn(other)=%d adherence_witness=%.3f -> %s\n",
                       other.c_str(),
                       (int)curveLooksUTurnForClusterExemption(candidate),
                       (int)looksLeftRightTurn(ob),
                       (int)curveLooksUTurnForClusterExemption(ob),
                       adherenceWitnessDistance(candidate, ob),
                       (curveLooksUTurnForClusterExemption(candidate) &&
                        looksLeftRightTurn(ob))
                           ? "结构性豁免候选(是否成立取决于贴行)"
                           : "无结构性豁免");
            }
        }
          }
         }
        }
       }
      }
    }

    // 同簇伙伴在"出口车道回退方向"上的有符号横向剖面：用于判断掉头的首尾直段
    // 应该整体偏向车道中心线的哪一侧才能与伙伴彻底分侧。
    // 参考系取本条掉头的 U 轴法线 lateral，正方向与 (p1-p0) 同号（与
    // buildSegmented 里的 side 同源）。
    printf("  frame: p0 back_axis=%+.4f lat=%+.4f | p1 back_axis=0 lat=0"
           " | station back_axis=%+.4f\n",
           (entry.first - exit.first).dot(axis),
           (entry.first - exit.first).dot(lateral_ref) * side,
           std::isfinite(station) ? station - exit.first.dot(axis) : 0.0);
    // 各伙伴曲线在家族平齐站位处的有符号横向坐标（可能多次穿越该站位）。
    if (std::isfinite(station)) {
        for (const ConnId& other : partners) {
            auto it = curves.find(other);
            if (it == curves.end()) continue;
            printf("    %-10s lat@station:", other.c_str());
            const std::vector<Vec2d> pts = it->second->sample(400);
            int printed = 0;
            for (std::size_t i = 0; i + 1 < pts.size(); ++i) {
                const double a = pts[i].dot(axis) - station;
                const double b = pts[i + 1].dot(axis) - station;
                if (a == 0.0 || (a < 0.0) != (b < 0.0)) {
                    const double t = std::abs(b - a) < 1e-12
                        ? 0.0 : -a / (b - a);
                    const Vec2d hit = pts[i] + (pts[i + 1] - pts[i]) * t;
                    printf(" %+.4f", (hit - exit.first).dot(lateral_ref) * side);
                    ++printed;
                }
            }
            if (printed == 0) printf(" (never reaches station)");
            printf("\n");
        }
    }
    printf("  exit-lane lateral profile (side=%+.0f, 正=远离入口直段一侧)\n",
           side);
    for (const ConnId& other : partners) {
        auto it = curves.find(other);
        if (it == curves.end()) continue;
        printf("    %-10s", other.c_str());
        for (double back : {1.5, 3.0, 5.0, 7.0, 9.24}) {
            // 找伙伴曲线上轴向站位最接近 p1 回退 back 米处的采样点。
            const Vec2d exit_back = -exit.second.normalized();
            const double want =
                exit.first.dot(axis) + back * exit_back.dot(axis);
            double best_lat = std::numeric_limits<double>::quiet_NaN();
            double best_err = 1e18;
            for (int i = 0; i <= 400; ++i) {
                const Vec2d pt = it->second->evaluate(
                    static_cast<double>(i) / 400.0);
                const double err = std::abs(pt.dot(axis) - want);
                if (err < best_err) {
                    best_err = err;
                    best_lat = (pt - exit.first).dot(lateral_ref) * side;
                }
            }
            if (best_err > 0.5)
                printf("  back%.1f=  --  ", back);
            else
                printf("  back%.1f=%+.3f", back, best_lat);
        }
        printf("\n");
    }
    return 0;
}
