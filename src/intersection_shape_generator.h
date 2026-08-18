#pragma once

#include "types.h"

namespace isg {

/// 拓扑校验报告：包含错误与警告列表
struct ValidationReport {
    std::vector<std::string> errors;
    std::vector<std::string> warnings;
    bool is_valid() const { return errors.empty(); }
};

/// 路口形态生成器：顶层入口，负责调度连通曲线生成、车道边线生成与路口面生成
class IntersectionShapeGenerator {
public:
    struct Config {
        double sdf_cell_size = 0.2;       ///< SDF网格单元尺寸（米）
        double obstacle_buffer = 1.0;     ///< mode=2障碍物绕行清距（米）；其他mode固定为0m
        double kappa_max = 0.25;          ///< 最大允许曲率
        LBFGSConfig lbfgs;                ///< LBFGS优化器配置
        ConnectivityDirectionConfig connectivity_direction;  ///< 连通方向配置
    };

    IntersectionShapeGenerator();
    explicit IntersectionShapeGenerator(const Config& config);

    /// 生成路口形态：输入路口数据，输出连通曲线、车道边线与精细路口面
    bool generate(const IntersectionInput& input, IntersectionOutput& output);

    /// 获取最近一次校验报告
    const ValidationReport& lastReport() const { return report_; }

private:
    Config config_;
    ValidationReport report_;
};

}  // 命名空间 isg
