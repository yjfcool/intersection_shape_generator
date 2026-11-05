#include "nonmotor/nonmotor_connectivity_generator.h"
#include "preprocessing/crosswalk_clearance_calculator.h"
#include "utils.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace isg {

NonMotorConnectivityGenerator::NonMotorConnectivityGenerator(const NonmotorConifg& cfg)
    : config_(cfg)
{
}

std::vector<ConnectivityCurve> NonMotorConnectivityGenerator::generate(
    const IntersectionInput& input,
    const std::vector<ConnectivityCurve>& motor_curves,
    NonMotorBackbone& backbone) {
    std::vector<ConnectivityCurve> result;
    for (size_t i = 0; i < input.connectivities.size(); ++i) {
        ConnectivityCurve curve;
        // TODO:
        //curve.geometry = generateOne(backbone, ...);
        result.push_back(curve);
    }
    return result;
}

LineString2d NonMotorConnectivityGenerator::generateOne(
    const NonMotorBackbone& backbone,
    const Vec2d& entry_pt, const Vec2d& entry_dir,
    const Vec2d& exit_pt, const Vec2d& exit_dir,
    int path_index, double corner_radius) {
    LineString2d geom;
    return geom;
}

}  // namespace isg
