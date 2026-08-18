#pragma once

#include "types.h"

namespace isg {

/// 将输入固有折线按原有规则转换为逐段三次 Bezier。
/// 该模块只负责表达转换，不判断障碍物、Boundary 或是否允许保留。
class FixedShapeInitializer {
public:
    BezierCurve build(const Connectivity& connectivity) const;
    bool hasGeometry(const Connectivity& connectivity) const;
};

}  // 命名空间 isg
