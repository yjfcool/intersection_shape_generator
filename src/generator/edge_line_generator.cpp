#include "generator/edge_line_generator.h"

#include <iostream>
#include <map>
#include <set>

#include "curve/bezier.h"
#include "utils.h"

namespace isg {

/**
 * 车道边线生成器(共享边优化版本)。
 *
 * 策略:
 *  1. 按(进入组, 退出组)聚合已生成中心线。
 *  2. 组内按进入车道顺序从内到外排序。
 *  3. 标记组内与其它车道中心线相交的特殊连通关系。
 *  4. 相邻非特殊车道若进入侧和退出侧都共享边界，则生成一条共享中线，
 *     避免左右两侧分别独立生成造成缝隙或重叠。
 *  5. 最外侧右边、最内侧左边和特殊车道边线仍独立生成。
 *
 * 共享中线特性:
 *  - 首尾端点G1连续。
 *  - 端点精确对齐共享边界连接点。
 *  - 使用alpha=0.38的三次Bezier拟合。
 */
class EdgeLineGeneratorImpl {
    double bezier_maxCurvature = 1.0;
    double edgeLine_defaultLaneWidth = 1.0;

public:
    std::vector<ConnectivityLaneEdge> generate(
        const IntersectionInput& inp, std::vector<ConnectivityCurve>& centerlines) {
        std::vector<ConnectivityLaneEdge> result;
        std::map<std::string, std::pair<std::string, std::string>> edgeIdMap;

        // 构建连通ID到生成中心线的映射。
        std::map<std::string, const ConnectivityCurve*> clMap;
        for (auto& gcl : centerlines) clMap[gcl.id] = &gcl;

        // 构建连通ID到输入连通关系的映射。
        std::map<std::string, const Connectivity*> connMap;
        for (auto& conn : inp.connectivities) connMap[conn.id] = &conn;

        // 按(进入组, 退出组)聚合中心线。
        using GroupKey = std::pair<std::string, std::string>;
        std::map<GroupKey, std::vector<const ConnectivityCurve*>> pairGroups;

        for (auto& gcl : centerlines) {
            auto cIt = connMap.find(gcl.id);
            if (cIt == connMap.end()) continue;
            const Connectivity* conn = cIt->second;
            GroupKey key = {conn->enterGroupId, conn->exitGroupId};
            pairGroups[key].push_back(&gcl);
        }

        // 逐组生成边线。
        for (auto& kv : pairGroups) {
            auto& key = kv.first;
            (void)key;
            auto& groupCls = kv.second;
            // 按进入车道顺序排序(内侧到外侧)。
            std::vector<const ConnectivityCurve*> sorted = groupCls;
            std::sort(sorted.begin(), sorted.end(),
                      [&](const ConnectivityCurve* a, const ConnectivityCurve* b) {
                          auto aIt = inp.findLane(a->entry_lane_id);
                          auto bIt = inp.findLane(b->entry_lane_id);
                          int aOrder = (aIt) ? aIt->laneOrder : 0;
                          int bOrder = (bIt) ? bIt->laneOrder : 0;
                          return aOrder < bOrder;
                      });

            // 标记组内与其它车道中心线相交的特殊连通关系。
            std::set<std::string> specialIds;
            if (sorted.size() > 1) {
                for (size_t i = 0; i < sorted.size(); ++i) {
                    for (size_t j = i + 1; j < sorted.size(); ++j) {
                        auto icurve = sorted[i]->curve;
                        auto jcurve = sorted[j]->curve;
                        if (!icurve || icurve->empty()) continue;
                        if (!jcurve || jcurve->empty()) continue;
                        auto ismps = icurve->sample(50);
                        auto jsmps = jcurve->sample(50);
                        if (ismps.size() < 2 || jsmps.size() < 2) continue;
                        if (polylinesIntersectExcludeEndpoints(ismps, jsmps)) {
                            specialIds.insert(sorted[i]->id);
                            specialIds.insert(sorted[j]->id);
                        }
                    }
                }
            }

            // 判定哪些相邻车道在进入侧和退出侧都共享边界。
            // sharedBetween[i] 为true表示第i条与第i+1条车道可共用中线。
            std::vector<bool> sharedBetween(sorted.size() > 0 ? sorted.size() - 1 : 0, false);
            // 记录每对相邻车道的共享边界ID。
            struct SharedEdgeInfo {
                std::string enterEdgeId;
                std::string exitEdgeId;
            };
            std::vector<SharedEdgeInfo> sharedEdgeInfos(sharedBetween.size());

            for (size_t i = 0; i + 1 < sorted.size(); ++i) {
                if (specialIds.count(sorted[i]->id) || specialIds.count(sorted[i + 1]->id))
                    continue;

                const std::string& enterLineI = sorted[i]->entry_lane_id;
                const std::string& enterLineI1 = sorted[i + 1]->entry_lane_id;
                const std::string& exitLineI = sorted[i]->exit_lane_id;
                const std::string& exitLineI1 = sorted[i + 1]->exit_lane_id;

                // 进入侧: 第i条车道右边界与第i+1条车道左边界一致。
                // ID不同时再比较几何连接点，兼容同一物理边界被不同ID表达的输入。
                auto enterI = inp.findLane(enterLineI);
                auto enterI1 = inp.findLane(enterLineI1);
                bool enterShared = false;
                std::string sharedEnterEdgeId;
                if (enterI && enterI1) {
                    const std::string& rightOfI = enterI->right_edge_id;
                    const std::string& leftOfI1 = enterI1->left_edge_id;
                    if (!rightOfI.empty() && !leftOfI1.empty()) {
                        // 优先按ID直接匹配。
                        if (rightOfI == leftOfI1) {
                            enterShared = true;
                            sharedEnterEdgeId = rightOfI;
                        } else {
                            // ID不同但连接点重合时也认为共享。
                            auto rightEdgeIt = inp.findEdge(rightOfI);
                            auto leftEdgeIt = inp.findEdge(leftOfI1);
                            if (rightEdgeIt && leftEdgeIt) {
                                auto rEdgePt = getConnPoint(rightEdgeIt->geometry.points, true);
                                auto lEdgePt = getConnPoint(leftEdgeIt->geometry.points, true);
                                double ptDist = dist(rEdgePt, lEdgePt);
                                if (ptDist < 0.01) {
                                    enterShared = true;
                                    sharedEnterEdgeId = rightOfI;
                                }
                            }
                        }
                    }
                }

                // 退出侧采用同样规则: 相邻车道之间的共享边界应重合。
                auto exitI = inp.findLane(exitLineI);
                auto exitI1 = inp.findLane(exitLineI1);
                bool exitShared = false;
                std::string sharedExitEdgeId;
                if (exitI && exitI1) {
                    const std::string& rightOfI = exitI->right_edge_id;
                    const std::string& leftOfI1 = exitI1->left_edge_id;
                    if (!rightOfI.empty() && !leftOfI1.empty()) {
                        // 优先按ID直接匹配。
                        if (rightOfI == leftOfI1) {
                            exitShared = true;
                            sharedExitEdgeId = rightOfI;
                        } else {
                            // ID不同但连接点重合时也认为共享。
                            auto rightEdgeIt = inp.findEdge(rightOfI);
                            auto leftEdgeIt = inp.findEdge(leftOfI1);
                            if (rightEdgeIt && leftEdgeIt) {
                                auto rEdgePt = getConnPoint(rightEdgeIt->geometry.points, false);
                                auto lEdgePt = getConnPoint(leftEdgeIt->geometry.points, false);
                                double ptDist = dist(rEdgePt, lEdgePt);
                                if (ptDist < 0.01) {
                                    exitShared = true;
                                    sharedExitEdgeId = rightOfI;
                                }
                            }
                        }
                    }
                }

                if (enterShared && exitShared) {
                    sharedBetween[i] = true;
                    sharedEdgeInfos[i] = {sharedEnterEdgeId, sharedExitEdgeId};
                }
            }

            // 为组内每条车道生成左右边线。
            for (size_t i = 0; i < sorted.size(); ++i) {
                const ConnectivityCurve& gcl = *sorted[i];
                if (!gcl.curve) continue;

                auto connIt = connMap.find(gcl.id);
                if (connIt == connMap.end()) continue;
                const Connectivity* conn = connIt->second;

                bool isSpecial = specialIds.count(gcl.id) > 0;

                // 左边界是否由上一条车道的共享右边界提供。
                bool leftShared = (!isSpecial && i > 0 && sharedBetween[i - 1]);
                // 右边界是否与下一条车道的左边界共享。
                bool rightShared = (!isSpecial && i + 1 < sorted.size() && sharedBetween[i]);

                // 生成左边界。
                std::string leftElId;
                if (leftShared) {
                    // 左边界已由(i-1, i)共享中线生成。
                    leftElId = "gen_el_shared_" + sorted[i - 1]->id + "_" + gcl.id;
                } else {
                    // 生成独立左边界。
                    leftElId = "gen_el_left_" + conn->id;
                    std::vector<Vec2d> leftPts = generateIndependentEdge(gcl, conn, true, inp);
                    ConnectivityLaneEdge el;
                    el.id = leftElId;
                    el.geometry.points = toVec3dArray(leftPts);
                    result.push_back(el);
                }

                // 生成右边界。
                std::string rightElId;
                if (rightShared) {
                    // 生成第i条与第i+1条车道之间的共享中线。
                    rightElId = "gen_el_shared_" + gcl.id + "_" + sorted[i + 1]->id;
                    std::vector<Vec2d> sharedPts = generateSharedMidline(
                        gcl, *sorted[i + 1],
                        sharedEdgeInfos[i].enterEdgeId,
                        sharedEdgeInfos[i].exitEdgeId,
                        inp);
                    ConnectivityLaneEdge el;
                    el.id = rightElId;
                    el.geometry.points = toVec3dArray(sharedPts);
                    // 共享边线ID记录在内侧车道右边界，外侧车道稍后回填到左边界。
                    result.push_back(el);
                } else {
                    // 生成独立右边界。
                    rightElId = "gen_el_right_" + conn->id;
                    std::vector<Vec2d> rightPts = generateIndependentEdge(gcl, conn, false, inp);
                    ConnectivityLaneEdge el;
                    el.id = rightElId;
                    el.geometry.points = toVec3dArray(rightPts);
                    result.push_back(el);
                }

                edgeIdMap[gcl.id] = {leftElId, rightElId};
            }
        }

        // 回填中心线的左右边线ID。
        for (auto& gcl : centerlines) {
            auto it = edgeIdMap.find(gcl.id);
            if (it != edgeIdMap.end()) {
                gcl.left_edge_id = it->second.first;
                gcl.right_edge_id = it->second.second;
            }
        }

        return result;
    }

private:
    // 生成相邻两条车道之间的共享中线。
    std::vector<Vec2d> generateSharedMidline(
        const ConnectivityCurve& gclI,
        const ConnectivityCurve& gclI1,
        const std::string& sharedEnterEdgeId,
        const std::string& sharedExitEdgeId,
        const IntersectionInput& inp) const {
        // 起点为进入侧共享边界连接点。
        Vec2d startPt{0, 0};
        auto enterEdgeIt = inp.findEdge(sharedEnterEdgeId);
        if (enterEdgeIt) {
            startPt = getConnPoint(enterEdgeIt->geometry.points, true);
        }

        // 终点为退出侧共享边界连接点。
        Vec2d endPt{0, 0};
        auto exitEdgeIt = inp.findEdge(sharedExitEdgeId);
        if (exitEdgeIt) {
            endPt = getConnPoint(exitEdgeIt->geometry.points, false);
        }

        // 起点切向取两条进入车道切向的平均方向。
        auto enterClI = inp.findLane(gclI.entry_lane_id);
        auto enterClI1 = inp.findLane(gclI1.entry_lane_id);
        Vec2d startTangI = (enterClI) ? getConnTangent(enterClI->geometry.points, true) : Vec2d{0, 1};
        Vec2d startTangI1 = (enterClI1) ? getConnTangent(enterClI1->geometry.points, true) : Vec2d{0, 1};
        Vec2d startTangSum = startTangI + startTangI1;
        Vec2d startTang = (startTangSum.norm() > EPS) ? startTangSum.normalized() : startTangI;

        // 终点切向取两条退出车道切向的平均方向。
        auto exitClI = inp.findLane(gclI.exit_lane_id);
        auto exitClI1 = inp.findLane(gclI1.exit_lane_id);
        Vec2d endTangI = (exitClI) ? getConnTangent(exitClI->geometry.points, false) : Vec2d{0, 1};
        Vec2d endTangI1 = (exitClI1) ? getConnTangent(exitClI1->geometry.points, false) : Vec2d{0, 1};
        Vec2d endTangSum = endTangI + endTangI1;
        Vec2d endTang = (endTangSum.norm() > EPS) ? endTangSum.normalized() : endTangI;

        // 使用alpha=0.38构造三次Bezier。
        double d = dist(startPt, endPt);
        if (d < EPS) return {startPt}; // 退化时返回单点折线，避免零长度线段。
        BezierSegment cb = makeCubicG1(startPt, startTang, endPt, endTang, 0.38);

        // 自适应采样以匹配既有边线质量。
        std::vector<Vec2d> pts = cb.sampleAdaptive(3.0, 2.0, 0.1);
        if (pts.empty()) {
            pts = cb.sampleCount(30);
        }

        // 强制端点精确对齐。
        if (!pts.empty()) {
            pts.front() = startPt;
            pts.back() = endPt;
        }

        return pts;
    }

    // 生成独立边界：先按中心线偏移，再平滑中段，并用Bezier重新贴合首尾端。
    std::vector<Vec2d> generateIndependentEdge(
        const ConnectivityCurve& gcl, const Connectivity* conn, bool isLeft, const IntersectionInput& inp) const {
        // 查找进入组。
        auto grpIt = inp.findGroup(conn->enterGroupId);
        if (!grpIt) {
            // 缺少组信息时退化为简单偏移。
            double hw = edgeLine_defaultLaneWidth * 0.5;
            return offsetPolyline(gcl.curve->sample(50), isLeft ? hw : -hw);
        }
        const LaneGroup& grp = *grpIt;

        // 估计半车道宽。
        double hw = estimateHalfWidth(conn->entry_lane_id, isLeft, grp, inp, true);

        // 生成偏移折线。
        std::vector<Vec2d> pts = offsetPolyline(gcl.curve->sample(50), isLeft ? hw : -hw);

        // 删除近重复点。
        removeDuplicates(pts, 0.02);

        if (pts.size() < 2) return pts;

        // 查找边界在进入/退出侧的端点与切向。
        std::string enterGrpId = conn->enterGroupId;
        std::string exitGrpId = conn->exitGroupId;

        Vec2d startPt, endPt;
        if (isLeft) {
            startPt = findEdgePtInGroup(conn->entry_lane_id, true, enterGrpId, inp, hw, true);
            endPt = findEdgePtInGroup(conn->exit_lane_id, true, exitGrpId, inp, hw, false);
        } else {
            startPt = findEdgePtInGroup(conn->entry_lane_id, false, enterGrpId, inp, hw, true);
            endPt = findEdgePtInGroup(conn->exit_lane_id, false, exitGrpId, inp, hw, false);
        }
        Vec2d enterTang = getEdgeTangent(conn->entry_lane_id, inp, true);
        Vec2d exitTang = getEdgeTangent(conn->exit_lane_id, inp, false);

        // 平滑中段。
        smoothMiddle(pts, 5);

        // 首尾用Bezier贴合到指定端点和切向。
        int smoothPts = std::max(6, (int)pts.size() / 4);
        bezierAlignEnds(pts, startPt, enterTang, endPt, exitTang, smoothPts);

        return pts;
    }

    // 删除连续近重复点。
    void removeDuplicates(std::vector<Vec2d>& pts, double minDist) const {
        if (pts.size() < 3) return;
        std::vector<Vec2d> out;
        out.push_back(pts.front());
        for (size_t i = 1; i < pts.size() - 1; ++i) {
            if (dist(pts[i], out.back()) >= minDist)
                out.push_back(pts[i]);
        }
        out.push_back(pts.back());
        pts = out;
    }

    // 使用角平分线法向偏移折线；offset为正表示左侧，负表示右侧。
    std::vector<Vec2d> offsetPolyline(const std::vector<Vec2d>& pts, double offset) const {
        if (pts.size() < 2) return pts;
        int n = (int)pts.size();
        std::vector<Vec2d> out;
        out.reserve(n);

        for (int i = 0; i < n; ++i) {
            Vec2d normal;
            if (i == 0) {
                Vec2d dir = (pts[1] - pts[0]).normalized();
                normal = rotLeft(dir);
            } else if (i == n - 1) {
                Vec2d dir = (pts[n - 1] - pts[n - 2]).normalized();
                normal = rotLeft(dir);
            } else {
                Vec2d d1 = (pts[i] - pts[i - 1]).normalized();
                Vec2d d2 = (pts[i + 1] - pts[i]).normalized();
                Vec2d n1 = rotLeft(d1);
                Vec2d n2 = rotLeft(d2);
                Vec2d avg = n1 + n2;
                if (avg.norm() < EPS) {
                    normal = n1;
                } else {
                    normal = avg.normalized();
                    double cosHalf = n1.dot(normal);
                    if (cosHalf > 0.3) {
                        double miter = std::min(1.0 / cosHalf, 3.0);
                        out.push_back(pts[i] + normal * (offset * miter));
                        continue;
                    }
                }
            }
            out.push_back(pts[i] + normal * offset);
        }
        return out;
    }

    // 移动平均平滑中段，保留首尾一定长度不动。
    void smoothMiddle(std::vector<Vec2d>& pts, int passes) const {
        int n = (int)pts.size();
        if (n < 5) return;
        int margin = std::max(2, n / 8);

        for (int pass = 0; pass < passes; ++pass) {
            std::vector<Vec2d> tmp = pts;
            for (int i = margin; i < n - margin; ++i) {
                Vec2d sum = pts[i] * 4.0;
                int cnt = 4;
                if (i - 1 >= 0) {
                    sum += pts[i - 1] * 2.0;
                    cnt += 2;
                }
                if (i + 1 < n) {
                    sum += pts[i + 1] * 2.0;
                    cnt += 2;
                }
                if (i - 2 >= 0) {
                    sum += pts[i - 2] * 1.0;
                    cnt += 1;
                }
                if (i + 2 < n) {
                    sum += pts[i + 2] * 1.0;
                    cnt += 1;
                }
                tmp[i] = sum / (double)cnt;
            }
            pts = tmp;
        }
    }

    // 首尾G1 Bezier重拟合，并强制端点对齐。
    void bezierAlignEnds(
        std::vector<Vec2d>& pts, const Vec2d& startPt, const Vec2d& startTang,
        const Vec2d& endPt, const Vec2d& endTang, int K) const {
        int n = (int)pts.size();
        if (n < 6) return;

        K = std::max(4, std::min(K, n / 3));

        pts.front() = startPt;
        pts.back() = endPt;

        // 头部区间[0..K]重拟合。
        {
            Vec2d p0 = startPt;
            Vec2d p3 = pts[K];
            double d = dist(p0, p3);
            if (d > EPS) {
                Vec2d t3;
                if (K + 1 < n) {
                    t3 = (pts[K + 1] - pts[K - 1]).normalized();
                } else {
                    t3 = (pts[K] - pts[K - 1]).normalized();
                }
                double alpha = 0.38;
                Vec2d p1 = p0 + startTang * (alpha * d);
                Vec2d p2 = p3 - t3 * (alpha * d);
                BezierSegment cb;
                cb.ctrl = {p0, p1, p2, p3};
                for (int i = 1; i < K; ++i) {
                    double t = (double)i / K;
                    pts[i] = cb.evaluate(t);
                }
            }
        }

        // 尾部区间[n-1-K..n-1]重拟合。
        {
            int startIdx = n - 1 - K;
            if (startIdx < K) startIdx = K;
            Vec2d p0 = pts[startIdx];
            Vec2d p3 = endPt;
            double d = dist(p0, p3);
            if (d > EPS) {
                Vec2d t0;
                if (startIdx > 0 && startIdx + 1 < n) {
                    t0 = (pts[startIdx] - pts[startIdx - 1]).normalized();
                } else {
                    t0 = (p3 - p0).normalized();
                }
                double alpha = 0.38;
                Vec2d p1 = p0 + t0 * (alpha * d);
                Vec2d p2 = p3 - endTang * (alpha * d);
                BezierSegment cb;
                cb.ctrl = {p0, p1, p2, p3};
                int count = n - 1 - startIdx;
                for (int i = 1; i < count; ++i) {
                    double t = (double)i / count;
                    pts[startIdx + i] = cb.evaluate(t);
                }
            }
        }

        pts.front() = startPt;
        pts.back() = endPt;
    }

    // 根据组边界端点估计半车道宽。
    double estimateHalfWidth(const std::string& entry_lane_id, bool isLeft,
                             const LaneGroup& grp, const IntersectionInput& inp, bool is_entryline) const {
        auto clit = inp.findLane(entry_lane_id);
        if (!clit) return edgeLine_defaultLaneWidth * 0.5;

        const Vec2d& clPt = getConnPoint(clit->geometry.points, is_entryline);
        Vec2d tangent = getConnTangent(clit->geometry.points, is_entryline);
        Vec2d normal = rotLeft(tangent);

        double bestDist = -1;
        for (auto& eid : grp.boundaries) {
            auto elit = inp.findEdge(eid);
            if (!elit) continue;
            if (elit->geometry.points.empty()) continue;
            const Vec2d& ep = getConnPoint(elit->geometry.points, is_entryline);
            double lateral = (ep - clPt).dot(normal);
            if (isLeft && lateral > 0.01 && lateral < 10.0) {
                if (bestDist < 0 || lateral < bestDist) bestDist = lateral;
            } else if (!isLeft && lateral < -0.01 && lateral > -10.0) {
                if (bestDist < 0 || (-lateral) < bestDist) bestDist = -lateral;
            }
        }

        if (bestDist > 0) return bestDist;
        return edgeLine_defaultLaneWidth * 0.5;
    }

    // 在指定组内查找边界端点。
    Vec2d findEdgePtInGroup(
        const std::string& lineId, bool isLeft, const std::string& groupId,
        const IntersectionInput& inp, double hw, bool is_entryline) const {
        auto clit = inp.findLane(lineId);
        if (!clit) return {0, 0};

        const Vec2d& clPt = getConnPoint(clit->geometry.points, is_entryline);
        const Vec2d& tang = getConnTangent(clit->geometry.points, is_entryline);
        Vec2d normal = rotLeft(tang);

        auto git = inp.findGroup(groupId);
        if (git) {
            for (auto& eid : git->boundaries) {
                auto elit = inp.findEdge(eid);
                if (!elit) continue;
                if (elit->geometry.points.empty()) continue;
                const Vec2d& ep = getConnPoint(elit->geometry.points, is_entryline);
                double lat = (ep - clPt).dot(normal);
                if (isLeft && lat > 0.01 && std::abs(lat - hw) < hw * 0.8) return ep;
                if (!isLeft && lat < -0.01 && std::abs(-lat - hw) < hw * 0.8) return ep;
            }
        }
        return clPt + normal * (isLeft ? hw : -hw);
    }

    // 获取车道中心线端点切向。
    Vec2d getEdgeTangent(const std::string& lineId, const IntersectionInput& inp, bool is_entryline) const {
        auto clit = inp.findLane(lineId);
        if (!clit) return {0, 1};
        return getConnTangent(clit->geometry.points, is_entryline);
    }
};

std::vector<ConnectivityLaneEdge> EdgeLineGenerator::generate(
    const IntersectionInput& input,
    std::vector<ConnectivityCurve>& centerlines) {
    return EdgeLineGeneratorImpl().generate(input, centerlines);
}

}
