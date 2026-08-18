#pragma once

#include "types.h"
#include "optimizer/penalty_cost.h"

namespace isg {

struct SampledCurve {
    std::vector<Vec2d> pts;
    BoundingBox2d bbox;
    Vec2d start{0, 0};
    Vec2d end{0, 0};
};

struct SampledSiblingCurve {
    ConnId id;
    SampledCurve sampled;
    BezierCurve curve;
    bool exempt_a1 = false;
    int expected_side = 0;
    Vec2d ref_perp{0, 0};
    bool shared_endpoint = false;
};

SampledCurve sampleCurveForIntersections(
    const BezierCurve& curve, int minimum_samples = 32);
std::vector<SampledSiblingCurve> sampleSiblingsForIntersections(
    const std::vector<SiblingCurve>& siblings, int minimum_samples = 32,
    bool sample_exempt = true);

}  // 命名空间 isg
