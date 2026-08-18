#pragma once

#include "types.h"

namespace isg {

/**
 * Generates left/right connectivity lane edges, including shared edges for
 * adjacent centerlines. The application facade intentionally does not invoke
 * this component until the separate feature review enables it.
 */
class EdgeLineGenerator {
public:
    EdgeLineGenerator() = default;

    std::vector<ConnectivityLaneEdge> generate(
            const IntersectionInput& input, std::vector<ConnectivityCurve>& centerlines);
};

}  // 命名空间 isg
