#include "nonmotor_backbone_extractor.h"
#include "io/iodata_shapefile.h"
#include "utils.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>


namespace isg {

NonMotorBackbone NonMotorBackboneExtractor::run(const IntersectionInput& input) const {
    NonMotorBackbone bb;
    // TODO:
    return bb;
}

bool NonMotorBackboneExtractor::output(const NonMotorBackbone& backbone, std::string out_dir) const {
    return save(backbone, out_dir, "");
}

}  // namespace isg
