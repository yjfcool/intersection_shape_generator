// 围栏（粗糙路口面）越界诊断：把 fence.containment / generator.fence_overflow 的
// 布尔结论拆成可定标的量。
//
// 动机：`fence.containment` 用零容差的 polygonContains 判断，而每条曲线的首尾连接点
// 就在车道端头，粗糙路口面的外环恰好从那里穿过——端点是"在内还是在外"变成浮点抛硬币。
// 审计报出的 generator.fence_overflow 里既有 `2.27m` 这种真实外溢，也有 `0.0017m`
// 这种毫米级噪声，必须先分清两类样本再决定容差，否则要么放过真实外溢、要么继续误报。
//
// 输出每条曲线：
//   out_pts     零容差判定为界外的采样点数 / 总采样数
//   max_dist    界外点到外环的最大距离
//   worst_s     最差点的弧长站位（占全长比例）
//   ep_dist     最差点到最近连接点的距离
//   near_ep     所有界外点是否都落在连接点附近（≤ 阈值）
//
// 用法: diag_fence [-v] [数据文件...]，缺省扫描 datas/ 全部 json。
//       -v 额外打印按段参数采样(sps=64,与审计一致)的越界点，即 curveInsideFence 的实际判据；
//       --sps N 指定该采样密度。
#include "constraints/fence_check.h"
#include "curve/curve_utils.h"
#include "intersection_shape_generator.h"
#include "io/iodata_json.h"
#include "utils.h"
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <dirent.h>
#include <string>
#include <vector>

using namespace isg;

namespace {

struct FenceProbe {
    int n_out = 0;
    int n_total = 0;
    double max_dist = 0.0;
    double worst_frac = 0.0;
    double worst_ep_dist = 0.0;
    double max_ep_dist_over_out = 0.0;  // 所有界外点里离端点最远的距离
};

// 按弧长均匀采样，逐点做零容差 polygonContains，与 curveInsideFence 的判据一致
// （后者按段参数采样，这里按弧长采样以便报告站位；两者对"是否存在界外点"等价）。
// 按段参数采样，复现 curveInsideFence 的取点方式（它对每段取 i/sps，
// 因此还会命中多段曲线的内部接缝，与按弧长采样命中的点集不同）。
// 布尔审计用参数采样、数值外溢用弧长采样，两者会互相矛盾——把差异打印出来。
void dumpParametricOffenders(const BezierCurve& curve, const Polygon2d& fence,
                             int sps, const char* id) {
    if (curve.empty() || fence.outer.empty())
        return;
    const Vec2d p0 = curve.startPt();
    const Vec2d p1 = curve.endPt();
    for (std::size_t s = 0; s < curve.segs.size(); ++s) {
        for (int i = 0; i <= sps; ++i) {
            const Vec2d pt = curve.segs[s].evaluate((double)i / sps);
            if (polygonContains(fence, pt))
                continue;
            const double d = pointToPolygonDist(pt, fence);
            const double ep = std::min((pt - p0).norm(), (pt - p1).norm());
            printf("      %-8s seg=%zu t=%.4f dist=%9.6fm ep_d=%9.6fm at (%.3f,%.3f)\n",
                   id, s, (double)i / sps, d, ep, pt.x(), pt.y());
        }
    }
}

// 连接两端点的直线弦相对围栏的最大外溢。
//
// 判据意义：若连 p0→p1 的直线弦都在粗糙面外，说明外环在这两个连接点之间是**内凹**的
// （典型情形：粗糙路口面把各进出口"路口喉部"用直线弦连起来，于是把转角切掉了），
// 那么任何贴着弦走的曲线都不可能待在面内——外溢不是生成器的形态缺陷，而是粗糙面
// 本身不覆盖该转向的可行走廊。这一栏把"生成器画歪了"和"面画小了"区分开。
double chordOverflow(const BezierCurve& curve, const Polygon2d& fence, int samples) {
    if (curve.empty() || fence.outer.empty())
        return 0.0;
    const Vec2d p0 = curve.startPt();
    const Vec2d p1 = curve.endPt();
    double worst = 0.0;
    for (int i = 0; i <= samples; ++i) {
        const Vec2d pt = p0 + (p1 - p0) * ((double)i / samples);
        if (polygonContains(fence, pt))
            continue;
        if (fenceOutsideIsConnectionPointRounding(fence, pt, p0, p1))
            continue;
        worst = std::max(worst, pointToPolygonDist(pt, fence));
    }
    return worst;
}

FenceProbe probeCurve(const BezierCurve& curve, const Polygon2d& fence, int samples) {
    FenceProbe probe;
    if (curve.empty() || fence.outer.empty())
        return probe;
    const std::vector<Vec2d> pts = curve.sampleByArcLength(samples);
    probe.n_total = (int)pts.size();
    if (pts.size() < 2)
        return probe;
    std::vector<double> station(pts.size(), 0.0);
    for (std::size_t i = 1; i < pts.size(); ++i)
        station[i] = station[i - 1] + (pts[i] - pts[i - 1]).norm();
    const double total = std::max(1e-9, station.back());
    const Vec2d p0 = curve.startPt();
    const Vec2d p1 = curve.endPt();
    for (std::size_t i = 0; i < pts.size(); ++i) {
        if (polygonContains(fence, pts[i]))
            continue;
        const double d = pointToPolygonDist(pts[i], fence);
        const double ep = std::min((pts[i] - p0).norm(), (pts[i] - p1).norm());
        ++probe.n_out;
        probe.max_ep_dist_over_out = std::max(probe.max_ep_dist_over_out, ep);
        if (d > probe.max_dist) {
            probe.max_dist = d;
            probe.worst_frac = station[i] / total;
            probe.worst_ep_dist = ep;
        }
    }
    return probe;
}

void runOne(const std::string& path, int verbose_sps) {
    IntersectionInput input = IntersectionIO::loadFromFile(path);
    if (input.area.geometry.outer.empty()) {
        printf("%-26s (no area polygon)\n", path.c_str());
        return;
    }
    IntersectionShapeGenerator gen;
    IntersectionOutput output;
    if (!gen.generate(input, output)) {
        printf("%-26s GENERATION FAILED\n", path.c_str());
        return;
    }
    const std::string base = path.substr(path.find_last_of('/') + 1);
    int n_curves = 0, n_violating = 0, n_endpoint_only = 0;
    double worst = 0.0;
    std::vector<std::string> lines;
    for (const auto& cc : output.connectivity_curves) {
        if (!cc.curve)
            continue;
        ++n_curves;
        const FenceProbe probe = probeCurve(*cc.curve, input.area.geometry, 160);
        if (probe.n_out == 0)
            continue;
        ++n_violating;
        worst = std::max(worst, probe.max_dist);
        const bool endpoint_only = probe.max_ep_dist_over_out <= 0.50;
        if (endpoint_only)
            ++n_endpoint_only;
        const double chord = chordOverflow(*cc.curve, input.area.geometry, 160);
        char buf[400];
        snprintf(buf, sizeof(buf),
                 "  %-10s turn=%d out=%3d/%3d max=%8.4fm worst_s=%.3f ep_d=%7.3fm "
                 "out_ep_max=%7.3fm chord_out=%8.4fm %s%s",
                 cc.id.c_str(), (int)cc.turn_type, probe.n_out, probe.n_total,
                 probe.max_dist, probe.worst_frac, probe.worst_ep_dist,
                 probe.max_ep_dist_over_out, chord,
                 endpoint_only ? "[ENDPOINT-ONLY]" : "[INTERIOR]",
                 chord > 1e-6 ? " [FENCE-CUTS-CORNER]" : "");
        lines.push_back(buf);
    }
    printf("%-26s rough=%d curves=%d violating=%d (endpoint-only=%d) worst=%.4fm\n",
           base.c_str(), (int)input.area.is_rough, n_curves, n_violating,
           n_endpoint_only, worst);
    for (const auto& l : lines)
        printf("%s\n", l.c_str());
    if (verbose_sps > 0) {
        // 审计判定落在这些参数采样点上，与上面按弧长的报告不是同一个点集：
        // 两者对同一条曲线可以给出相反结论。diag_all_violations 用 profile.samples=64，
        // 库内默认 25，越密越容易在端点附近多打出一个"不是连接点"的采样点。
        printf("    -- parametric samples (sps=%d, curveInsideFence 判据) --\n",
               verbose_sps);
        for (const auto& cc : output.connectivity_curves)
            if (cc.curve)
                dumpParametricOffenders(*cc.curve, input.area.geometry, verbose_sps,
                                        cc.id.c_str());
    }
}

}  // namespace

int main(int argc, char** argv) {
    std::vector<std::string> paths;
    int verbose_sps = 0;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "-v" || a == "--verbose")
            verbose_sps = 64;  // 与 diag_all_violations 的 profile.samples 一致
        else if (a == "--sps" && i + 1 < argc)
            verbose_sps = std::atoi(argv[++i]);
        else
            paths.push_back(a);
    }
    if (paths.empty()) {
        for (const char* dir : {"datas", "../datas"}) {
            DIR* d = opendir(dir);
            if (!d)
                continue;
            struct dirent* ent = nullptr;
            while ((ent = readdir(d)) != nullptr) {
                const std::string name = ent->d_name;
                if (name.size() > 5 && name.substr(name.size() - 5) == ".json")
                    paths.push_back(std::string(dir) + "/" + name);
            }
            closedir(d);
            break;
        }
        std::sort(paths.begin(), paths.end());
    }
    for (const auto& p : paths)
        runOne(p, verbose_sps);
    return 0;
}
