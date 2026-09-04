#pragma once

#include "types.h"
#include "optimizer/penalty_cost.h"

namespace isg {

struct SampledCurve {
    std::vector<Vec2d> pts;
    BoundingBox2d bbox;
    // 与采样折线段一一对应的包围盒和中点。共享族/闭包矩阵会反复检查
    // 同一批候选，预计算这些值可避免每次成对相交时重复构造。
    std::vector<BoundingBox2d> segment_boxes;
    std::vector<Vec2d> segment_midpoints;
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
