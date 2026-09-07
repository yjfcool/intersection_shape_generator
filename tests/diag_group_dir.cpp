// 同组切向统一（ConnectivityDirectionNormalizer / GroupUnified）的强制转角定标。
//
// 动机：该工具的目的是抹掉**同一进出口内平行车道之间**的切向抖动，让同组曲线端点切向
// 一致。但它的实现是"给整组选一个方向，然后覆盖组内每条车道的端点切向"，一条车道自己
// 的朝向与被强加的方向差多少都不检查。一旦输入把**两个物理方向不同的路口臂**塞进同一
// 个 groupId（`110004764` 的 `43107602` 就是：五条车道朝 `-11.5°`、两条朝 `-83°~-88°`），
// 组内少数派车道的端点切向会被整整扭转几十度，之后所有以它为端点的曲线都在"错误的切向"
// 上做 G1，于是必然自交、越过路缘、与同簇兄弟成片相交——而审计用的是同一份规范化输入，
// 反而看不出 G1 异常。
//
// 本工具直接对比 `ConnectivityDirectionNormalizer` **前后**每条车道的端点切向，给出被
// 强制转过的角度，用来定标"合法抖动"与"跨臂误并"的分界。
//
// 用法: diag_group_dir [--limit 度] [--baseline 米] [数据文件...]，缺省扫描 datas/ 全部
// json。`--limit`/`--baseline` 直接透传到 `ConnectivityDirectionConfig`，可在不重编库的
// 前提下扫阈值；`--limit 1e9` 等价于关闭跨臂保护，此时"强制转角"即组方向与车道自身朝向
// 的原始偏差，便于重新定标。
#include "io/iodata_json.h"
#include "toolkits/toolkits.h"
#include "types.h"
#include "utils.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <dirent.h>
#include <string>
#include <unordered_map>
#include <vector>

using namespace isg;

namespace {

struct LaneDir {
    Vec2d dir = Vec2d(1, 0);
    Vec2d robust = Vec2d(0, 0);  // 长基线弦方向；退化时为零向量
    bool valid = false;
};

// 自路口端点回溯至累计长度达 baseline，取该点到端点的弦；与 toolkits.cpp 的判据一致。
Vec2d robustDir(const std::vector<Vec3d>& points, GroupRole role, double baseline) {
    if (baseline <= 0.0 || points.size() < 2)
        return Vec2d(0, 0);
    const bool from_back = role == GroupRole::Entry;
    const Vec2d endpoint = from_back ? Vec2d(points.back()) : Vec2d(points.front());
    Vec2d anchor = endpoint;
    double accumulated = 0.0;
    for (std::size_t step = 1; step < points.size(); ++step) {
        const std::size_t i = from_back ? points.size() - 1 - step : step;
        const Vec2d current(points[i]);
        const std::size_t prev = from_back ? i + 1 : i - 1;
        accumulated += (current - Vec2d(points[prev])).norm();
        anchor = current;
        if (accumulated >= baseline)
            break;
    }
    const Vec2d chord = from_back ? endpoint - anchor : anchor - endpoint;
    return chord.norm() > 1e-10 ? chord.normalized() : Vec2d(0, 0);
}

// 按车道在组内的角色取端点切向：Entry 取末点切向，Exit 取首点切向。
std::unordered_map<LaneId, LaneDir> snapshot(const IntersectionInput& input, double baseline) {
    std::unordered_map<LaneId, GroupRole> role;
    for (std::size_t g = 0; g < input.lane_groups.size(); ++g)
        for (std::size_t i = 0; i < input.lane_groups[g].lanes.size(); ++i)
            role[input.lane_groups[g].lanes[i]] = input.lane_groups[g].role;
    std::unordered_map<LaneId, LaneDir> out;
    for (std::size_t i = 0; i < input.lanes.size(); ++i) {
        const Lane& lane = input.lanes[i];
        const std::unordered_map<LaneId, GroupRole>::const_iterator it = role.find(lane.id);
        if (it == role.end() || lane.geometry.points.size() < 2)
            continue;
        const Vec2d d = it->second == GroupRole::Entry
            ? entryLineTangent(lane.geometry.points)
            : exitLineTangent(lane.geometry.points);
        LaneDir rec;
        rec.valid = d.norm() > 1e-10;
        if (rec.valid)
            rec.dir = d.normalized();
        rec.robust = robustDir(lane.geometry.points, it->second, baseline);
        out[lane.id] = rec;
    }
    return out;
}

double angleDeg(const Vec2d& a, const Vec2d& b) {
    double c = a.dot(b);
    c = std::max(-1.0, std::min(1.0, c));
    return std::acos(c) * 180.0 / M_PI;
}

struct Row {
    std::string dataset;
    LaneGroupId group;
    LaneId lane;
    double deg = 0.0;         // 两点式切向被强制转过的角度
    double robust_deg = 0.0;  // 长基线弦与统一后切向的夹角；-1 表示弦退化
};

void runOne(const std::string& path, std::vector<Row>& rows,
            const ConnectivityDirectionConfig& cfg) {
    const IntersectionInput raw = IntersectionIO::loadFromFile(path);
    IntersectionInput input = InputNormalizer(raw);
    const std::unordered_map<LaneId, LaneDir> before =
            snapshot(input, cfg.group_robust_baseline_m);
    ConnectivityDirectionNormalizer(input, cfg);
    const std::unordered_map<LaneId, LaneDir> after =
            snapshot(input, cfg.group_robust_baseline_m);

    std::unordered_map<LaneId, LaneGroupId> owner;
    for (std::size_t g = 0; g < input.lane_groups.size(); ++g)
        for (std::size_t i = 0; i < input.lane_groups[g].lanes.size(); ++i)
            owner[input.lane_groups[g].lanes[i]] = input.lane_groups[g].id;

    const std::string base = path.substr(path.find_last_of('/') + 1);
    std::vector<Row> local;
    for (std::unordered_map<LaneId, LaneDir>::const_iterator it = before.begin();
         it != before.end(); ++it) {
        const std::unordered_map<LaneId, LaneDir>::const_iterator jt = after.find(it->first);
        if (jt == after.end() || !it->second.valid || !jt->second.valid)
            continue;
        Row row;
        row.dataset = base;
        row.group = owner.count(it->first) ? owner[it->first] : LaneGroupId();
        row.lane = it->first;
        row.deg = angleDeg(it->second.dir, jt->second.dir);
        row.robust_deg = it->second.robust.norm() > 1e-10
                ? angleDeg(it->second.robust, jt->second.dir) : -1.0;
        local.push_back(row);
    }
    std::sort(local.begin(), local.end(),
              [](const Row& a, const Row& b) { return a.deg > b.deg; });
    double worst = local.empty() ? 0.0 : local.front().deg;
    int over5 = 0, over20 = 0;
    for (std::size_t i = 0; i < local.size(); ++i) {
        if (local[i].deg > 5.0) ++over5;
        if (local[i].deg > 20.0) ++over20;
    }
    std::printf("%-26s lanes=%2d worst=%7.2f°  >5°=%d  >20°=%d\n", base.c_str(),
                (int)local.size(), worst, over5, over20);
    for (std::size_t i = 0; i < local.size() && local[i].deg > 5.0; ++i)
        std::printf("    lane %-10s group %-10s forced %7.2f°  (长基线弦 %7.2f°)\n",
                    local[i].lane.c_str(), local[i].group.c_str(), local[i].deg,
                    local[i].robust_deg);
    rows.insert(rows.end(), local.begin(), local.end());
}

}  // namespace

int main(int argc, char** argv) {
    std::vector<std::string> paths;
    ConnectivityDirectionConfig cfg;  // 缺省与生成器默认档一致
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--limit" && i + 1 < argc)
            cfg.group_force_limit_deg = std::atof(argv[++i]);
        else if (arg == "--baseline" && i + 1 < argc)
            cfg.group_robust_baseline_m = std::atof(argv[++i]);
        else
            paths.push_back(arg);
    }
    std::printf("跨臂保护阈值 limit=%.3g°  抗抖动基线 baseline=%.3gm\n",
                cfg.group_force_limit_deg, cfg.group_robust_baseline_m);
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
    std::vector<Row> rows;
    for (std::size_t i = 0; i < paths.size(); ++i)
        runOne(paths[i], rows, cfg);

    std::sort(rows.begin(), rows.end(),
              [](const Row& a, const Row& b) { return a.deg > b.deg; });
    const double edges[] = {0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 45.0, 90.0, 181.0};
    int bins[9] = {0};
    for (std::size_t i = 0; i < rows.size(); ++i)
        for (int b = 0; b < 9; ++b)
            if (rows[i].deg < edges[b]) { ++bins[b]; break; }
    std::printf("\n==== 全 corpus 强制转角分布（车道数 %d）====\n", (int)rows.size());
    double lo = 0.0;
    for (int b = 0; b < 9; ++b) {
        std::printf("  [%6.1f°, %6.1f°)  %4d\n", lo, edges[b], bins[b]);
        lo = edges[b];
    }
    std::printf("\n==== 转角最大的 20 条车道（末列为长基线弦与统一后切向夹角）====\n");
    for (std::size_t i = 0; i < rows.size() && i < 20; ++i)
        std::printf("  %-26s lane %-10s group %-10s %7.2f°  弦 %7.2f°\n",
                    rows[i].dataset.c_str(), rows[i].lane.c_str(), rows[i].group.c_str(),
                    rows[i].deg, rows[i].robust_deg);
    return 0;
}
