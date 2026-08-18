#pragma once

#include "domain/generation_types.h"

#include <cstddef>
#include <vector>

namespace isg {

/// 从已完成当前几何审计的候选中按业务键稳定选择。
/// 未审计或存在硬违规的候选不会进入比较；完全相同时保留枚举顺序。
class CandidateSelector {
public:
    static const std::size_t npos;

    std::size_t selectBest(const std::vector<CurveCandidate>& candidates) const;
};

}  // 命名空间 isg
