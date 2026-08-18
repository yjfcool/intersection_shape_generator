#pragma once

#include "types.h"

#include <functional>
#include <unordered_map>
#include <vector>

namespace isg {

struct CurvePatchEntry {
    ConnId id;
    BezierCurve curve;
};

typedef std::function<void(ConnectivityCurve&, const BezierCurve&)>
    CurvePatchPreparer;

/// 在所有目标可用且候选准备完成后一次性提交多条曲线更新。
class AtomicCurvePatch {
public:
    bool apply(const std::vector<CurvePatchEntry>& entries,
               std::vector<ConnectivityCurve>& results,
               const CurvePatchPreparer& prepare) const;
};

}  // 命名空间 isg
