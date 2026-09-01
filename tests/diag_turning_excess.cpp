// 单向绕转跨度诊断：打印每条最终曲线的段数、弧长/弦长、最大曲率与
// curveTurningSpan（有符号累加的 max-min，单向绕转跨度）。
//
// 用途：`kMaxNonUTurnTurningSpan` 门禁是按"正常转向的绕转跨度上限"标定的。
// 一旦某个数据集的合法曲线跨度逼近阈值，门禁就会误伤，需要用实测分布重新定标。
//
// 用法：diag_turning_excess <数据文件> [--id 连接id]...
#include "curve/curve_utils.h"
#include "intersection_shape_generator.h"
#include "io/iodata_json.h"
#include "utils.h"

#include <cstdio>
#include <cstring>
#include <string>
#include <unordered_set>

using namespace isg;

int main(int argc, char** argv) {
    std::string path;
    std::unordered_set<std::string> only_ids;
    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "--id") == 0 && i + 1 < argc)
            only_ids.insert(argv[++i]);
        else if (path.empty())
            path = argv[i];
    }
    if (path.empty()) {
        std::fprintf(stderr, "用法: diag_turning_excess <数据文件> [--id 连接id]...\n");
        return 2;
    }

    IntersectionInput input = IntersectionIO::loadFromFile(path);
    IntersectionShapeGenerator gen;
    IntersectionOutput output;
    if (!gen.generate(input, output)) {
        std::fprintf(stderr, "generate() 返回 false\n");
        return 1;
    }

    double worst = 0.0;
    std::string worst_id;
    for (const auto& cc : output.connectivity_curves) {
        if (!cc.curve || cc.curve->empty())
            continue;
        if (!only_ids.empty() && only_ids.count(cc.id) == 0)
            continue;
        const BezierCurve& c = *cc.curve;
        const double chord = (c.endPt() - c.startPt()).norm();
        const double span_deg = curveTurningSpan(c) * RAD2DEG;
        if (span_deg > worst) {
            worst = span_deg;
            worst_id = cc.id;
        }
        std::printf("%-10s segs=%d arc/chord=%.3f maxk=%.4f span=%.1fdeg status=%d\n",
                    cc.id.c_str(), c.numSegments(),
                    chord > 1e-8 ? c.arcLength() / chord : 0.0,
                    c.maxCurvature(60), span_deg, (int)cc.status);
    }
    std::printf("\n最大绕转跨度: %s = %.1fdeg\n", worst_id.c_str(), worst);
    return 0;
}
