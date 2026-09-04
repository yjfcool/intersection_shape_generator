#pragma once

#include "types.h"
#include "optimizer/penalty_cost.h"

namespace isg {

// 采样折线的分块粒度。成对相交是 O(Na×Nb) 的双层循环，两条 240 点曲线
// 就有约 5.7 万次配对。把 b 侧每 kSampleBlockSize 段合并成一个粗包围盒，
// 可以整块跳过“分段包围盒必然不相交”的区间——这些配对在逐段判断里也会
// 被同一条包围盒条件拒掉，因此跳过不改变结果，只是不再逐个访问。
constexpr int kSampleBlockSize = 16;

struct SampledCurve {
    std::vector<Vec2d> pts;
    BoundingBox2d bbox;
    // 与采样折线段一一对应的包围盒和中点。共享族/闭包矩阵会反复检查
    // 同一批候选，预计算这些值可避免每次成对相交时重复构造。
    std::vector<BoundingBox2d> segment_boxes;
    std::vector<Vec2d> segment_midpoints;
    // 每 kSampleBlockSize 段合并的粗包围盒，用于成段剪枝。
    std::vector<BoundingBox2d> block_boxes;
    Vec2d start{0, 0};
    Vec2d end{0, 0};
};

// 由已填好的 pts 补齐 bbox、分段包围盒、分段中点与分块包围盒。
// 所有构造 SampledCurve 的位置都应走这里，保证剪枝索引始终一致。
void buildSampledCurveIndex(SampledCurve& sampled);

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
