#include "generation/curve_sampling.h"

#include <algorithm>
#include <cmath>

namespace isg {

SampledCurve sampleCurveForIntersections(
    const BezierCurve& curve, int minimum_samples) {
    SampledCurve sampled;
    sampled.start = curve.startPt();
    sampled.end = curve.endPt();
    int count = std::max(
        minimum_samples, static_cast<int>(std::ceil(curve.arcLength() / 0.35)) + 1);
    sampled.pts = curve.sampleByArcLength(std::min(count, 160));
    for (const auto& point : sampled.pts)
        sampled.bbox.expand(point);
    if (sampled.pts.size() >= 2) {
        sampled.segment_boxes.resize(sampled.pts.size() - 1);
        sampled.segment_midpoints.reserve(sampled.pts.size() - 1);
        for (std::size_t i = 0; i + 1 < sampled.pts.size(); ++i) {
            sampled.segment_boxes[i].expand(sampled.pts[i]);
            sampled.segment_boxes[i].expand(sampled.pts[i + 1]);
            sampled.segment_midpoints.push_back(
                0.5 * (sampled.pts[i] + sampled.pts[i + 1]));
        }
    }
    return sampled;
}

std::vector<SampledSiblingCurve> sampleSiblingsForIntersections(
    const std::vector<SiblingCurve>& siblings, int minimum_samples,
    bool sample_exempt) {
    std::vector<SampledSiblingCurve> sampled;
    sampled.reserve(siblings.size());
    for (const auto& sibling : siblings) {
        SampledSiblingCurve item;
        item.id = sibling.id;
        item.exempt_a1 = sibling.exempt_a1;
        item.expected_side = sibling.expected_side;
        item.ref_perp = sibling.ref_perp;
        item.shared_endpoint = sibling.shared_endpoint;
        item.curve = sibling.curve;
        if (sample_exempt || !sibling.exempt_a1)
            item.sampled = sampleCurveForIntersections(sibling.curve, minimum_samples);
        sampled.push_back(std::move(item));
    }
    return sampled;
}

}  // 命名空间 isg
