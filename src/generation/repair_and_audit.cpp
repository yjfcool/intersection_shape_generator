#include "generation/repair_and_audit.h"

#include <algorithm>
#include <unordered_map>
#include <unordered_set>

namespace isg {

// 物理修复按 Boundary、Obstacle 的既有优先级执行，并在每步后重做风险审计。
PhysicalRepairResult PhysicalRepairCoordinator::repair(
    const BezierCurve& initial, const BezierCurve& current,
    bool geometric_uturn, const PhysicalRiskAuditor& auditor,
    const PhysicalRepairAttempt& boundary_repair,
    const PhysicalRepairAttempt& obstacle_repair) const {
    PhysicalRepairResult result;
    result.curve = current;
    result.risk = auditor(result.curve);
    if (geometric_uturn && result.risk.boundary) {
        const PhysicalRiskSnapshot initial_risk = auditor(initial);
        if (!initial_risk.physical()) {
            result.curve = initial;
            result.risk = initial_risk;
        }
    }
    if (result.risk.boundary && !geometric_uturn) {
        BezierCurve repaired = result.curve;
        if (boundary_repair(result.curve, repaired)) {
            result.curve = repaired;
            result.boundary_repaired = true;
            result.risk = auditor(result.curve);
        }
    }
    if (result.risk.physical() && !geometric_uturn) {
        BezierCurve repaired = result.curve;
        if (obstacle_repair(result.curve, repaired)) {
            result.curve = repaired;
            result.risk = auditor(result.curve);
        }
    }
    return result;
}

// 影响闭包只纳入变更曲线关联的非结构性曲线对。
RepairImpactClosure RepairImpactClosureBuilder::build(
    const std::vector<CurvePair>& pairs,
    const std::unordered_set<ConnId>& changed_ids) const {
    RepairImpactClosure closure;
    closure.curve_ids = changed_ids;
    for (const auto& pair : pairs) {
        if (pair.exempt == CrossExemption::StructuralCross)
            continue;
        if (!changed_ids.count(pair.id_a) && !changed_ids.count(pair.id_b))
            continue;
        closure.pairs.push_back(pair);
        closure.curve_ids.insert(pair.id_a);
        closure.curve_ids.insert(pair.id_b);
        closure.neighbors[pair.id_a].push_back(pair.id_b);
        closure.neighbors[pair.id_b].push_back(pair.id_a);
    }
    return closure;
}

std::unordered_map<ConnId, std::vector<ConnId> >
RepairImpactClosureBuilder::neighbors(
    const std::vector<CurvePair>& pairs,
    const std::unordered_set<ConnId>& available_ids) const {
    std::unordered_map<ConnId, std::vector<ConnId> > result;
    for (const auto& pair : pairs) {
        if (pair.exempt == CrossExemption::StructuralCross ||
            !available_ids.count(pair.id_a) ||
            !available_ids.count(pair.id_b)) {
            continue;
        }
        result[pair.id_a].push_back(pair.id_b);
        result[pair.id_b].push_back(pair.id_a);
    }
    return result;
}

// 所有候选先完成准备，再统一写回，保证家族修复不会部分提交。
bool AtomicCurvePatch::apply(
    const std::vector<CurvePatchEntry>& entries,
    std::vector<ConnectivityCurve>& results,
    const CurvePatchPreparer& prepare) const {
    if (entries.empty() || !prepare)
        return false;
    std::unordered_map<ConnId, std::size_t> result_index;
    for (std::size_t i = 0; i < results.size(); ++i)
        result_index[results[i].id] = i;

    std::unordered_set<ConnId> seen;
    std::vector<std::size_t> indices;
    std::vector<ConnectivityCurve> prepared;
    indices.reserve(entries.size());
    prepared.reserve(entries.size());
    for (const auto& entry : entries) {
        const auto found = result_index.find(entry.id);
        if (found == result_index.end() || !seen.insert(entry.id).second)
            return false;
        indices.push_back(found->second);
    }
    for (std::size_t i = 0; i < entries.size(); ++i) {
        ConnectivityCurve candidate = results[indices[i]];
        prepare(candidate, entries[i].curve);
        prepared.push_back(candidate);
    }
    for (std::size_t i = 0; i < indices.size(); ++i)
        results[indices[i]] = prepared[i];
    return true;
}

// 最终单曲线审计统一回写状态和违规原因，避免各修复器各自标注。
void FinalCurveAuditor::apply(const FinalCurveAuditSnapshot& snapshot,
                              ConnectivityCurve& curve) const {
    curve.violation.reason.clear();
    curve.violation.max_obstacle_penetration = snapshot.max_obstacle_penetration;
    if (snapshot.obstacle_intersection) {
        curve.violation.max_obstacle_penetration = std::max(
            curve.violation.max_obstacle_penetration, 0.01);
        curve.violation.reason = "curve intersects obstacle";
    }
    if (snapshot.update_fence_overflow)
        curve.violation.max_fence_overflow = snapshot.max_fence_overflow;
    if (snapshot.self_intersection)
        curve.violation.reason = "curve self-intersects away from endpoints";
    if (snapshot.boundary_intersection) {
        curve.violation.reason = snapshot.uturn
            ? "U-turn curve intersects boundary away from endpoints"
            : "curve intersects boundary away from endpoints";
    }
    if (snapshot.road_edge_clearance_violation)
        curve.violation.reason = "curve violates RoadEdge clearance";
    if (snapshot.ordinary_single_axis_violation) {
        curve.violation.reason =
            "ordinary single cubic control points leave endpoint direction axes";
    }
    if (snapshot.self_intersection || snapshot.boundary_intersection ||
        snapshot.road_edge_clearance_violation ||
        snapshot.ordinary_single_axis_violation ||
        curve.violation.max_obstacle_penetration > 0.05) {
        curve.status = CurveStatus::Degraded;
    } else if (!curve.violation.exempt_crosses.empty()) {
        curve.status = CurveStatus::WarnA2;
    } else {
        curve.status = CurveStatus::OK;
    }
}

// 最终曲线对审计跳过结构性交叉，其余发现交由统一状态回写处理。
std::vector<FinalPairFinding> FinalPairAuditor::audit(
    const std::vector<CurvePair>& pairs,
    const std::vector<ConnectivityCurve>& results,
    const FinalPairViolationDetector& detector) const {
    std::vector<FinalPairFinding> findings;
    if (!detector)
        return findings;
    std::unordered_map<ConnId, std::size_t> index;
    for (std::size_t i = 0; i < results.size(); ++i)
        index[results[i].id] = i;
    for (const auto& pair : pairs) {
        if (pair.exempt == CrossExemption::StructuralCross)
            continue;
        const auto ia = index.find(pair.id_a);
        const auto ib = index.find(pair.id_b);
        if (ia == index.end() || ib == index.end())
            continue;
        const ConnectivityCurve& a = results[ia->second];
        const ConnectivityCurve& b = results[ib->second];
        if (!a.curve || !b.curve)
            continue;
        Vec2d location(0, 0);
        if (!detector(pair, *a.curve, *b.curve, location))
            continue;
        FinalPairFinding finding;
        finding.id_a = pair.id_a;
        finding.id_b = pair.id_b;
        finding.exemption = pair.exempt;
        finding.location = location;
        findings.push_back(finding);
    }
    return findings;
}

void FinalPairAuditor::apply(
    const std::vector<FinalPairFinding>& findings,
    std::vector<ConnectivityCurve>& results) const {
    std::unordered_map<ConnId, std::size_t> index;
    for (std::size_t i = 0; i < results.size(); ++i)
        index[results[i].id] = i;
    for (const auto& finding : findings) {
        const auto ia = index.find(finding.id_a);
        const auto ib = index.find(finding.id_b);
        if (ia == index.end() || ib == index.end())
            continue;
        ConnectivityCurve& a = results[ia->second];
        ConnectivityCurve& b = results[ib->second];
        a.violation.exempt_crosses.push_back(finding.location);
        b.violation.exempt_crosses.push_back(finding.location);
        const CurveStatus status =
            finding.exemption == CrossExemption::ObstacleCross
                ? CurveStatus::WarnA2 : CurveStatus::Degraded;
        if (a.status == CurveStatus::OK || status == CurveStatus::Degraded)
            a.status = status;
        if (b.status == CurveStatus::OK || status == CurveStatus::Degraded)
            b.status = status;
    }
}

}  // 命名空间 isg
