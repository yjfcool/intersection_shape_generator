#include "intersection_shape_generator.h"
#include "utils.h"
#include "generator/connectivity_generator.h"
#include "constraints/infeasibility_detector.h"
#include "generator/polygon_builder.h"
#include "toolkits/toolkits.h"
#include "optimizer/sdf_field.h"
#include <chrono>
#include <cmath>

namespace isg {

IntersectionShapeGenerator::IntersectionShapeGenerator(): config_{} {
}

IntersectionShapeGenerator::IntersectionShapeGenerator(const Config& config): config_(config) {
}

/// 拓扑校验：检查车道组边界数量、连通关系转向类型与几何方位是否一致等
static ValidationReport validateTopology(const IntersectionInput& input) {
    auto turnTypeName = [&](ConnTurnType t) -> const char*  {
        switch (t) {
            case ConnTurnType::UTurnLeft: return "UTurnLeft";
            case ConnTurnType::TurnLeft: return "TurnLeft";
            case ConnTurnType::Straight: return "Straight";
            case ConnTurnType::TurnRight: return "TurnRight";
            case ConnTurnType::UTurnRight: return "UTurnRight";
            default: return "Unknown";
        }
    };
    ValidationReport r;
    for (auto& conn : input.connectivities) {
        ConnTurnType inf = TurnPreprocessor(conn, input);
        if (std::abs((int)conn.turn_type - (int)inf) > 1)
            r.warnings.push_back("Connectivity " + conn.id + ": declared=" + turnTypeName(conn.turn_type)
                + " geometry=" + turnTypeName(inf));
    }
    if (!input.area.geometry.outer.empty() && !isSimplePolygon(input.area.geometry))
        r.errors.push_back("Coarse intersection area is self-intersecting");
    for (auto& sl : input.stop_lines)
        if (!sl.associated_group_id.empty() && !input.laneGroupExists(sl.associated_group_id))
            r.warnings.push_back("StopLine " + sl.id + ": group '" + sl.associated_group_id + "' not found");
    return r;
}

/// 获取数据外轮廓矩形框
static BoundingBox2d getOuterBoundingBox(const IntersectionInput& norm_input) {
    BoundingBox2d roi;
    if (!norm_input.area.geometry.outer.empty())
        roi = norm_input.area.geometry.bbox();
    else {
        for (auto& l : norm_input.lanes) {
            for (auto& p : l.geometry.points)
                roi.expand(p);
        }
        roi.min_pt -= Vec2d(20, 20);
        roi.max_pt += Vec2d(20, 20);
    }
    return roi;
}

bool IntersectionShapeGenerator::generate(const IntersectionInput& input, IntersectionOutput& output) {
    try {
        /// 初始化原数据缺失信息
        const IntersectionInput norm_input = InputNormalizer(input);
        report_ = validateTopology(norm_input);
        if (!report_.is_valid())
            return false;

        /// 构建避让空间距离查询器
        auto t0 = std::chrono::steady_clock::now();
        SDFField sdf;
        if (!norm_input.obstacles.empty()) {
            double obstacle_clearance = obstacleAvoidanceClearanceForMode(norm_input.mode);
            if (norm_input.mode == 2)
                obstacle_clearance = std::max(obstacle_clearance, config_.obstacle_buffer);
            sdf.build(getOuterBoundingBox(norm_input), norm_input.obstacles, config_.sdf_cell_size, obstacle_clearance);
        }
        output.perf.sdf_build_ms =
            std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();

        /// 曲线生成
        double opt_ms = 0;
        ConnectivityGenerator cgen(config_.lbfgs, config_.connectivity_direction);
        output.connectivity_curves = cgen.generate(norm_input, sdf, &opt_ms);
        ElevationInterpolator(output.connectivity_curves, norm_input); //高程插值
        output.perf.optimize_ms = opt_ms;

        /// 车道边线生成：按兼容基线保持车道边线生成关闭，lane_edges 继续为空。
        auto te = std::chrono::steady_clock::now();
        // EdgeLineGenerator elgen;
        // output.lane_edges = elgen.generate(norm_input, output.connectivity_curves);
        output.perf.edge_gen_ms =
            std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - te).count();

        /// 路口面生成
        auto ta = std::chrono::steady_clock::now();
        IntersectionAreaBuilder areabuilder(intersectionAreaExtendDistance(norm_input.mode), 4.0);
        output.area = areabuilder.build(norm_input, output.connectivity_curves, output.lane_edges);
        output.perf.area_gen_ms =
            std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - ta).count();

    } catch (const std::exception& error) {
        report_.errors.push_back(std::string("generation failed: ") + error.what());
        output.connectivity_curves.clear();
        output.lane_edges.clear();
        output.area = IntersectionArea();
        return false;
    } catch (...) {
        report_.errors.push_back("generation failed: unknown exception");
        output.connectivity_curves.clear();
        output.lane_edges.clear();
        output.area = IntersectionArea();
        return false;
    }
    return true;
}

}
