//
// Created by yuanjinfa on 2026/8/14.
//

#ifndef ISG_TOOLKITS_H
#define ISG_TOOLKITS_H

#pragma once

#include "types.h"

namespace isg {

    /// 插值曲线z值：根据输入数据的z值
    void ElevationInterpolator(std::vector<ConnectivityCurve>& curves, const IntersectionInput& input);

    /// 根据端点切向计算转向元数据，同时保持输入转向类型兼容。
    ConnTurnType TurnPreprocessor(const Connectivity& connectivity, const IntersectionInput& input);

    /// 完善输入数据缺失信息：复制输入并补齐兼容的组引用(负责groupID补全和可配置方向统一)，不修改调用方原对象
    IntersectionInput InputNormalizer(const IntersectionInput& input);

    /// 按配置统一同组车道的路口端切向，不修改车道端点位置。
    void ConnectivityDirectionNormalizer(IntersectionInput& input, const ConnectivityDirectionConfig& config);

};


#endif //ISG_TOOLKITS_H
