// 性能基准：逐份数据集测量单路口生成耗时。
//
// 输出两路：
//   1) stdout —— 与历史格式兼容的 `<path> ok=1 123.4 ms` 行（多次重复时取中位数）。
//   2) CSV（`--csv <路径>`）—— 与 diag_all_violations 同一套发现记录列：
//      source,category,check,dataset,subject,severity,metric,value,threshold,detail,location
//      由 tools/test_report.py 汇总进 Excel 的“性能”sheet。
//
// 判级：单路口生成耗时上限 15000ms（产品要求），超限记 error；
// 另外记录 3000ms 的交互体验参考线，超过记 warn，便于提前发现性能劣化趋势。
//
// 用法：
//   diag_gen_time [数据文件...] [--dir 数据目录] [--csv 输出.csv] [--repeat 次数]
//   不指定数据文件与目录时，默认扫描 PROJECT_ROOT_DIR/datas 下全部 .json。
#include "intersection_shape_generator.h"
#include "io/iodata_json.h"
#include "utils.h"

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <dirent.h>
#include <string>
#include <vector>

using namespace isg;

namespace {

const double kHardLimitMs = 15000.0;  ///< 产品硬上限
const double kSoftLimitMs = 3000.0;   ///< 交互体验参考线

/// 取中位数（重复测量抑制单次抖动）。
double median(std::vector<double> v) {
    if (v.empty()) return 0.0;
    std::sort(v.begin(), v.end());
    const size_t n = v.size();
    return n % 2 ? v[n / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]);
}

/// 取文件名（作为按文件统计的键）。
std::string baseName(const std::string& path) {
    const size_t p = path.find_last_of("/\\");
    return p == std::string::npos ? path : path.substr(p + 1);
}

std::string csvEscape(const std::string& s) {
    if (s.find_first_of(",\"\n\r") == std::string::npos) return s;
    std::string out = "\"";
    for (size_t i = 0; i < s.size(); ++i) {
        if (s[i] == '"') out += "\"\"";
        else if (s[i] == '\n' || s[i] == '\r') out += ' ';
        else out += s[i];
    }
    return out + "\"";
}

/// 单条发现记录（列顺序必须与 diag_all_violations 的 CSV 保持一致）。
struct Row {
    std::string check, dataset, severity, metric, detail;
    double value = 0.0, threshold = 0.0;
};

std::vector<Row> g_rows;

void emit(const std::string& check, const std::string& dataset,
          const std::string& severity, const std::string& metric, double value,
          double threshold, const std::string& detail) {
    Row r;
    r.check = check; r.dataset = dataset; r.severity = severity;
    r.metric = metric; r.value = value; r.threshold = threshold; r.detail = detail;
    g_rows.push_back(r);
}

bool writeCsv(const std::string& path) {
    FILE* f = std::fopen(path.c_str(), "w");
    if (!f) return false;
    std::fprintf(f, "source,category,check,dataset,subject,severity,metric,value,"
                    "threshold,detail,location\n");
    for (size_t i = 0; i < g_rows.size(); ++i) {
        const Row& r = g_rows[i];
        std::fprintf(f, "perf,perf,%s,%s,-,%s,%s,%.4f,%.4f,%s,%s\n",
                     csvEscape(r.check).c_str(), csvEscape(r.dataset).c_str(),
                     r.severity.c_str(), csvEscape(r.metric).c_str(), r.value,
                     r.threshold, csvEscape(r.detail).c_str(),
                     "tests/diag_gen_time.cpp");
    }
    std::fclose(f);
    return true;
}

/// 枚举目录下全部 .json（固定顺序，便于逐版本对比）。
std::vector<std::string> listJson(const std::string& dir) {
    std::vector<std::string> files;
    DIR* d = opendir(dir.c_str());
    if (!d) {
        std::fprintf(stderr, "Cannot open %s\n", dir.c_str());
        return files;
    }
    struct dirent* ent = nullptr;
    while ((ent = readdir(d)) != nullptr) {
        const std::string name = ent->d_name;
        if (name.size() < 6) continue;
        if (name.substr(name.size() - 5) != ".json") continue;
        files.push_back(dir + "/" + name);
    }
    closedir(d);
    std::sort(files.begin(), files.end());
    return files;
}

}  // namespace

int main(int argc, char** argv) {
    std::vector<std::string> paths;
    std::string csv, dir;
    int repeat = 1;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--csv" && i + 1 < argc) csv = argv[++i];
        else if (a == "--dir" && i + 1 < argc) dir = argv[++i];
        else if (a == "--repeat" && i + 1 < argc) repeat = std::atoi(argv[++i]);
        else if (a.size() > 2 && a.compare(0, 2, "--") == 0) {
            std::fprintf(stderr, "未知参数: %s\n", a.c_str());
            return 2;
        } else paths.push_back(a);
    }
    if (repeat < 1) repeat = 1;
    if (paths.empty()) {
        if (dir.empty()) dir = std::string(PROJECT_ROOT_DIR) + "/datas";
        paths = listJson(dir);
    }

    int over_hard = 0;
    for (size_t i = 0; i < paths.size(); ++i) {
        const std::string name = baseName(paths[i]);
        IntersectionInput input = IntersectionIO::loadFromFile(paths[i]);
        std::vector<double> samples;
        bool ok = true;
        for (int r = 0; r < repeat; ++r) {
            IntersectionShapeGenerator gen;
            IntersectionOutput output;
            IntersectionInput copy = input;  // generate 会规范化输入，逐次用干净副本
            const std::chrono::steady_clock::time_point t0 =
                std::chrono::steady_clock::now();
            ok = gen.generate(copy, output) && ok;
            samples.push_back(std::chrono::duration<double, std::milli>(
                std::chrono::steady_clock::now() - t0).count());
        }
        const double ms = median(samples);
        std::printf("%s ok=%d %.1f ms\n", paths[i].c_str(), ok ? 1 : 0, ms);

        const char* severity = ms > kHardLimitMs ? "error"
                             : (ms > kSoftLimitMs ? "warn" : "info");
        if (ms > kHardLimitMs) ++over_hard;
        char detail[160];
        std::snprintf(detail, sizeof(detail),
                      "单路口生成耗时中位数（%d 次测量，min=%.0fms max=%.0fms）",
                      repeat, *std::min_element(samples.begin(), samples.end()),
                      *std::max_element(samples.begin(), samples.end()));
        emit("perf.generation_ms", name, severity, "ms", ms, kHardLimitMs, detail);
        if (!ok)
            emit("perf.generation_failed", name, "error", "", 0.0, 0.0,
                 "generate() 返回 false，耗时数据仅供参考");
    }

    std::printf("\n==== 超过 %.0fms 硬上限: %d / %zu ====\n", kHardLimitMs, over_hard,
                paths.size());
    if (!csv.empty() && !writeCsv(csv)) {
        std::fprintf(stderr, "写入 CSV 失败: %s\n", csv.c_str());
        return 1;
    }
    return 0;
}
