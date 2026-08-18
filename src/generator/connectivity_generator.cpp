#include "connectivity_generator.h"
#include "generator/connectivity_generation_session.h"

namespace isg {

ConnectivityGenerator::ConnectivityGenerator(
    const LBFGSConfig& config, const ConnectivityDirectionConfig& direction_config)
    : config_(config), direction_config_(direction_config) {}

std::vector<ConnectivityCurve> ConnectivityGenerator::generate(
    const IntersectionInput& input, SDFField& sdf, double* out_ms) {
    ConnectivityGenerationSession session(config_, direction_config_);
    return session.run(input, sdf, out_ms);
}

}  // 命名空间 isg
