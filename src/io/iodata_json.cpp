#include "io/iodata_json.h"

#include "utils/json.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <memory>
#include <sstream>

namespace isg {

// 适配自定义的轻量级JSON库的序列化/反序列化功能

// 为std::vector<std::string>、std::map等类型添加JSON序列化助手
namespace detail {
    inline nlohmann::json vectorToStringJson(const std::vector<std::string>& vec) {
        nlohmann::json j = nlohmann::json::array();
        for (const auto& item : vec) {
            j.push_back(item);
        }
        return j;
    }

    inline std::vector<std::string> vectorFromStringJson(const nlohmann::json& j) {
        std::vector<std::string> result;
        if (j.is_array()) {
            for (size_t i = 0; i < j.size(); ++i) {
                result.push_back(j[i].get<std::string>());
            }
        }
        return result;
    }

    inline nlohmann::json mapToAttrMapJson(const std::map<std::string, std::string>& m) {
        nlohmann::json j = nlohmann::json::object();
        for (auto it = m.begin(); it != m.end(); ++it) {
            j[it->first] = it->second;
        }
        return j;
    }

    inline std::map<std::string, std::string> attrMapFromJson(const nlohmann::json& j) {
        std::map<std::string, std::string> result;
        if (j.is_object()) {
            for (auto it = j.items().begin(); it != j.items().end(); ++it) {
                result[it.key()] = it.value().get<std::string>();
            }
        }
        return result;
    }
}

// 序列化和反序列化 Vec2d（方向向量仍为二维）
inline nlohmann::json vec2dToJson(const Vec2d& v) {
    nlohmann::json j = nlohmann::json::object();
    j["x"] = v.x();
    j["y"] = v.y();
    return j;
}

inline Vec2d vec2dFromJson(const nlohmann::json& j) {
    Vec2d v;
    if (j.contains("x") && j.contains("y")) {
        v << j["x"].get<double>(), j["y"].get<double>();
    } else {
        v << 0.0, 0.0;  // 默认值
    }
    return v;
}

// 序列化和反序列化 Vec3d（坐标点包含高程，旧二维JSON缺省z=0）
nlohmann::json vec3dToJson(const Vec3d& v) {
    nlohmann::json j = nlohmann::json::object();
    j["x"] = v.x();
    j["y"] = v.y();
    j["z"] = v.z();
    return j;
}

Vec3d vec3dFromJson(const nlohmann::json& j) {
    if (j.contains("x") && j.contains("y")) {
        double z = j.contains("z") ? j["z"].get<double>() : 0.0;
        return Vec3d(j["x"].get<double>(), j["y"].get<double>(), z);
    }
    return Vec3d(0.0, 0.0, 0.0);
}

// 序列化和反序列化 BoundingBox2d
inline nlohmann::json boundingBox2dToJson(const BoundingBox2d& bbox) {
    nlohmann::json j = nlohmann::json::object();
    j["min_pt"] = vec2dToJson(bbox.min_pt);
    j["max_pt"] = vec2dToJson(bbox.max_pt);
    return j;
}

inline BoundingBox2d boundingBox2dFromJson(const nlohmann::json& j) {
    BoundingBox2d bbox;
    if (j.contains("min_pt")) {
        bbox.min_pt = vec2dFromJson(j["min_pt"]);
    }
    if (j.contains("max_pt")) {
        bbox.max_pt = vec2dFromJson(j["max_pt"]);
    }
    return bbox;
}

// 序列化和反序列化 LineString2d
nlohmann::json lineString2dToJson(const LineString2d& ls) {
    nlohmann::json j = nlohmann::json::object();
    nlohmann::json points = nlohmann::json::array();
    for (const auto& pt : ls.points) {
        points.push_back(vec3dToJson(pt));
    }
    j["points"] = points;
    return j;
}

inline LineString2d lineString2dFromJson(const nlohmann::json& j) {
    LineString2d ls;
    if (j.contains("points")) {
        const nlohmann::json& points = j["points"];
        for (size_t i = 0; i < points.size(); ++i) {
            ls.points.push_back(vec3dFromJson(points[i]));
        }
    }
    return ls;
}

// 序列化和反序列化 Polygon2d
inline nlohmann::json polygon2dToJson(const Polygon2d& poly) {
    nlohmann::json j = nlohmann::json::object();
    nlohmann::json outer = nlohmann::json::array();
    for (const auto& pt : poly.outer) {
        outer.push_back(vec3dToJson(pt));
    }
    j["outer"] = outer;

    nlohmann::json holes = nlohmann::json::array();
    for (const auto& hole : poly.holes) {
        nlohmann::json hole_array = nlohmann::json::array();
        for (const auto& pt : hole) {
            hole_array.push_back(vec3dToJson(pt));
        }
        holes.push_back(hole_array);
    }
    j["holes"] = holes;
    return j;
}

inline Polygon2d polygon2dFromJson(const nlohmann::json& j) {
    Polygon2d poly;

    if (j.contains("outer")) {
        const nlohmann::json& outer = j["outer"];
        for (size_t i = 0; i < outer.size(); ++i) {
            poly.outer.push_back(vec3dFromJson(outer[i]));
        }
    }

    if (j.contains("holes")) {
        const nlohmann::json& holes = j["holes"];
        for (size_t h_idx = 0; h_idx < holes.size(); ++h_idx) {
            std::vector<Vec3d> hole;
            const nlohmann::json& hole_json = holes[h_idx];
            for (size_t i = 0; i < hole_json.size(); ++i) {
                hole.push_back(vec3dFromJson(hole_json[i]));
            }
            poly.holes.push_back(hole);
        }
    }
    return poly;
}

// 序列化和反序列化 LaneEdge
inline nlohmann::json laneEdgeToJson(const LaneEdge& le) {
    nlohmann::json j = nlohmann::json::object();
    j["id"] = le.id;
    j["geometry"] = lineString2dToJson(le.geometry);
    j["is_shared"] = le.is_shared;
    j["attrs"] = detail::mapToAttrMapJson(le.attrs);
    j["groupId"] = le.groupId;
    j["lineOrder"] = le.lineOrder;

    if (le.shared_by) {
        nlohmann::json shared_by_obj = nlohmann::json::object();
        shared_by_obj["first"] = le.shared_by->first;
        shared_by_obj["second"] = le.shared_by->second;
        j["shared_by"] = shared_by_obj;
    }

    return j;
}

inline LaneEdge laneEdgeFromJson(const nlohmann::json& j) {
    LaneEdge le;
    if (j.contains("id")) {
        le.id = j["id"].get<std::string>();
    }
    if (j.contains("geometry")) {
        le.geometry = lineString2dFromJson(j["geometry"]);
    }
    if (j.contains("is_shared")) {
        le.is_shared = j["is_shared"].get<bool>();
    }
    if (j.contains("attrs")) {
        le.attrs = detail::attrMapFromJson(j["attrs"]);
    }
    if (j.contains("groupId")) {
        le.groupId = j["groupId"].get<std::string>();
    }
    if (j.contains("lineOrder")) {
        le.lineOrder = j["lineOrder"].get<int>();
    }

    if (j.contains("shared_by")) {
        const nlohmann::json& shared_data = j["shared_by"];
        auto shared_pair = std::make_shared<std::pair<LaneId, LaneId>>();
        if (shared_data.contains("first")) {
            shared_pair->first = shared_data["first"].get<std::string>();
        }
        if (shared_data.contains("second")) {
            shared_pair->second = shared_data["second"].get<std::string>();
        }
        le.shared_by = shared_pair;
    }

    return le;
}

// 序列化和反序列化 Lane
inline nlohmann::json laneToJson(const Lane& l) {
    nlohmann::json j = nlohmann::json::object();
    j["id"] = l.id;
    j["left_edge_id"] = l.left_edge_id;
    j["right_edge_id"] = l.right_edge_id;
    j["width"] = l.width;
    j["geometry"] = lineString2dToJson(l.geometry);
    j["attrs"] = detail::mapToAttrMapJson(l.attrs);
    j["groupId"] = l.groupId;
    j["laneOrder"] = l.laneOrder;

    return j;
}

inline Lane laneFromJson(const nlohmann::json& j) {
    Lane l;
    if (j.contains("id")) {
        l.id = j["id"].get<std::string>();
    }
    if (j.contains("left_edge_id")) {
        l.left_edge_id = j["left_edge_id"].get<std::string>();
    }
    if (j.contains("right_edge_id")) {
        l.right_edge_id = j["right_edge_id"].get<std::string>();
    }
    if (j.contains("width")) {
        l.width = j["width"].get<double>();
    }
    if (j.contains("geometry")) {
        l.geometry = lineString2dFromJson(j["geometry"]);
    }
    if (j.contains("attrs")) {
        l.attrs = detail::attrMapFromJson(j["attrs"]);
    }
    if (j.contains("groupId")) {
        l.groupId = j["groupId"].get<std::string>();
    }
    if (j.contains("laneOrder")) {
        l.laneOrder = j["laneOrder"].get<int>();
    }

    return l;
}

// 序列化和反序列化 LaneGroup
inline nlohmann::json laneGroupToJson(const LaneGroup& lg) {
    nlohmann::json j = nlohmann::json::object();
    j["id"] = lg.id;
    j["role"] = static_cast<int>(lg.role);
    j["lanes"] = detail::vectorToStringJson(lg.lanes);
    j["boundaries"] = detail::vectorToStringJson(lg.boundaries);
    j["attrs"] = detail::mapToAttrMapJson(lg.attrs);

    return j;
}

inline LaneGroup laneGroupFromJson(const nlohmann::json& j) {
    LaneGroup lg;
    if (j.contains("id")) {
        lg.id = j["id"].get<std::string>();
    }
    if (j.contains("role")) {
        lg.role = static_cast<GroupRole>(j["role"].get<int>());
    }
    if (j.contains("lanes")) {
        lg.lanes = detail::vectorFromStringJson(j["lanes"]);
    }
    if (j.contains("boundaries")) {
        lg.boundaries = detail::vectorFromStringJson(j["boundaries"]);
    }
    if (j.contains("attrs")) {
        lg.attrs = detail::attrMapFromJson(j["attrs"]);
    }

    return lg;
}

// 序列化和反序列化 Connectivity
inline nlohmann::json connectivityToJson(const Connectivity& conn) {
    nlohmann::json j = nlohmann::json::object();
    j["id"] = conn.id;
    j["entry_lane_id"] = conn.entry_lane_id;
    j["exit_lane_id"] = conn.exit_lane_id;
    j["turn_type"] = static_cast<int>(conn.turn_type);
    j["enterGroupId"] = conn.enterGroupId;
    j["exitGroupId"] = conn.exitGroupId;
    j["lane_type"] = static_cast<int>(conn.lane_type);
    j["fixed_shape"] = (conn.fixed_shape ? 1 : 0);
    if (!conn.geometry.points.empty()) {
        j["geometry"] = lineString2dToJson(conn.geometry);
    }

    return j;
}

inline Connectivity connectivityFromJson(const nlohmann::json& j) {
    Connectivity conn;
    if (j.contains("id")) {
        conn.id = j["id"].get<std::string>();
    }
    if (j.contains("entry_lane_id")) {
        conn.entry_lane_id = j["entry_lane_id"].get<std::string>();
    }
    if (j.contains("exit_lane_id")) {
        conn.exit_lane_id = j["exit_lane_id"].get<std::string>();
    }
    if (j.contains("turn_type")) {
        conn.turn_type = static_cast<ConnTurnType>(j["turn_type"].get<int>());
    }
    if (j.contains("enterGroupId")) {
        conn.enterGroupId = j["enterGroupId"].get<std::string>();
    }
    if (j.contains("exitGroupId")) {
        conn.exitGroupId = j["exitGroupId"].get<std::string>();
    }
    if (j.contains("lane_type")) {
        conn.lane_type = static_cast<ConnLaneType>(j["lane_type"].get<int>());
    }
    if (j.contains("fixed_shape")) {
        conn.fixed_shape = j["fixed_shape"].get<int>();
    }
    if (j.contains("geometry")) {
        conn.geometry = lineString2dFromJson(j["geometry"]);
    }

    return conn;
}

// 序列化和反序列化 Obstacle
inline nlohmann::json obstacleToJson(const Obstacle& obs) {
    nlohmann::json j = nlohmann::json::object();
    j["id"] = obs.id;
    j["geometry"] = polygon2dToJson(obs.geometry);
    j["buffered_geometry"] = polygon2dToJson(obs.buffered_geometry);

    return j;
}

inline Obstacle obstacleFromJson(const nlohmann::json& j) {
    Obstacle obs;
    if (j.contains("id")) {
        obs.id = j["id"].get<std::string>();
    }
    if (j.contains("geometry")) {
        obs.geometry = polygon2dFromJson(j["geometry"]);
    }
    if (j.contains("buffered_geometry")) {
        obs.buffered_geometry = polygon2dFromJson(j["buffered_geometry"]);
    }

    return obs;
}

// 序列化和反序列化 Boundary
inline nlohmann::json boundaryToJson(const Boundary& bd) {
    nlohmann::json j = nlohmann::json::object();
    j["id"] = bd.id;
    j["type"] = static_cast<int>(bd.type);
    j["geometry"] = lineString2dToJson(bd.geometry);

    return j;
}

inline Boundary boundaryFromJson(const nlohmann::json& j) {
    Boundary bd;
    if (j.contains("id")) {
        bd.id = j["id"].get<std::string>();
    }
    if (j.contains("type")) {
        bd.type = static_cast<Boundary::Type>(j["type"].get<int>());
    }
    if (j.contains("geometry")) {
        bd.geometry = lineString2dFromJson(j["geometry"]);
    }

    return bd;
}

// 序列化和反序列化 StopLine
inline nlohmann::json stopLineToJson(const StopLine& sl) {
    nlohmann::json j = nlohmann::json::object();
    j["id"] = sl.id;
    j["geometry"] = lineString2dToJson(sl.geometry);
    j["associated_group_id"] = sl.associated_group_id;
    j["normal_direction"] = vec2dToJson(sl.normal_direction);

    return j;
}

inline StopLine stopLineFromJson(const nlohmann::json& j) {
    StopLine sl;
    if (j.contains("id")) {
        sl.id = j["id"].get<std::string>();
    }
    if (j.contains("geometry")) {
        sl.geometry = lineString2dFromJson(j["geometry"]);
    }
    if (j.contains("associated_group_id")) {
        sl.associated_group_id = j["associated_group_id"].get<std::string>();
    }
    if (j.contains("normal_direction")) {
        sl.normal_direction = vec2dFromJson(j["normal_direction"]);
    }

    return sl;
}

// 序列化和反序列化 Crosswalk
inline nlohmann::json crosswalkToJson(const Crosswalk& cw) {
    nlohmann::json j = nlohmann::json::object();
    j["id"] = cw.id;
    j["geometry"] = polygon2dToJson(cw.geometry);
    j["crossing_direction"] = vec2dToJson(cw.crossing_direction);

    return j;
}

inline Crosswalk crosswalkFromJson(const nlohmann::json& j) {
    Crosswalk cw;
    if (j.contains("id")) {
        cw.id = j["id"].get<std::string>();
    }
    if (j.contains("geometry")) {
        cw.geometry = polygon2dFromJson(j["geometry"]);
    }
    if (j.contains("crossing_direction")) {
        cw.crossing_direction = vec2dFromJson(j["crossing_direction"]);
    }

    return cw;
}

// 序列化和反序列化 IntersectionArea
inline nlohmann::json intersectionAreaToJson(const IntersectionArea& ia) {
    nlohmann::json j = nlohmann::json::object();
    j["id"] = ia.id;
    j["geometry"] = polygon2dToJson(ia.geometry);
    j["is_rough"] = ia.is_rough;

    return j;
}

inline IntersectionArea intersectionAreaFromJson(const nlohmann::json& j) {
    IntersectionArea ia;
    if (j.contains("id")) {
        ia.id = j["id"].get<std::string>();
    }
    if (j.contains("geometry")) {
        ia.geometry = polygon2dFromJson(j["geometry"]);
    }
    if (j.contains("is_rough")) {
        ia.is_rough = j["is_rough"].get<bool>();
    }

    return ia;
}

// 序列化和反序列化 ViolationInfo
inline nlohmann::json violationInfoToJson(const ViolationInfo& vi) {
    nlohmann::json j = nlohmann::json::object();
    j["type"] = static_cast<int>(vi.type);
    j["max_obstacle_penetration"] = vi.max_obstacle_penetration;
    j["max_fence_overflow"] = vi.max_fence_overflow;
    j["fence_expansion_applied"] = vi.fence_expansion_applied;

    nlohmann::json exempt_crosses = nlohmann::json::array();
    for (const auto& ec : vi.exempt_crosses) {
        exempt_crosses.push_back(vec2dToJson(ec));
    }
    j["exempt_crosses"] = exempt_crosses;

    j["reason"] = vi.reason;

    return j;
}

inline ViolationInfo violationInfoFromJson(const nlohmann::json& j) {
    ViolationInfo vi;
    if (j.contains("type")) {
        vi.type = static_cast<ViolationInfo::InfeasibilityType>(j["type"].get<int>());
    }
    if (j.contains("max_obstacle_penetration")) {
        vi.max_obstacle_penetration = j["max_obstacle_penetration"].get<double>();
    }
    if (j.contains("max_fence_overflow")) {
        vi.max_fence_overflow = j["max_fence_overflow"].get<double>();
    }
    if (j.contains("fence_expansion_applied")) {
        vi.fence_expansion_applied = j["fence_expansion_applied"].get<double>();
    }

    if (j.contains("exempt_crosses")) {
        const nlohmann::json& crosses = j["exempt_crosses"];
        for (size_t i = 0; i < crosses.size(); ++i) {
            vi.exempt_crosses.push_back(vec2dFromJson(crosses[i]));
        }
    }

    if (j.contains("reason")) {
        vi.reason = j["reason"].get<std::string>();
    }

    return vi;
}

// 序列化和反序列化 ConnectivityCurve
inline nlohmann::json connectivityCurveToJson(const ConnectivityCurve& cc) {
    nlohmann::json j = nlohmann::json::object();
    j["id"] = cc.id;
    j["entry_lane_id"] = cc.entry_lane_id;
    j["exit_lane_id"] = cc.exit_lane_id;
    j["turn_type"] = static_cast<int>(cc.turn_type);
    j["lane_type"] = static_cast<int>(cc.lane_type);
    j["fixed_shape"] = (cc.fixed_shape ? 1 : 0);
    LineString2d geometry = cc.geometry;
    if (geometry.points.empty() && cc.curve) {
        int n = std::max(2, std::min(240, (int)std::ceil(cc.curve->arcLength() / 0.3) + 1));
        geometry.points = toVec3dArray(cc.curve->sampleByArcLength(n));
    }
    if (!geometry.points.empty()) {
        j["geometry"] = lineString2dToJson(geometry);
    }
    j["status"] = static_cast<int>(cc.status);
    j["violation"] = violationInfoToJson(cc.violation);
    j["left_edge_id"] = cc.left_edge_id;
    j["right_edge_id"] = cc.right_edge_id;

    return j;
}

inline ConnectivityCurve connectivityCurveFromJson(const nlohmann::json& j) {
    ConnectivityCurve cc;
    if (j.contains("id")) {
        cc.id = j["id"].get<std::string>();
    }
    if (j.contains("entry_lane_id")) {
        cc.entry_lane_id = j["entry_lane_id"].get<std::string>();
    }
    if (j.contains("exit_lane_id")) {
        cc.exit_lane_id = j["exit_lane_id"].get<std::string>();
    }
    if (j.contains("turn_type")) {
        cc.turn_type = static_cast<ConnTurnType>(j["turn_type"].get<int>());
    }
    if (j.contains("geometry")) {
        cc.geometry = lineString2dFromJson(j["geometry"]);
    }
    if (j.contains("fixed_shape")) {
        cc.fixed_shape = j["fixed_shape"].get<int>();
    }
    if (j.contains("status")) {
        cc.status = static_cast<CurveStatus>(j["status"].get<int>());
    }
    if (j.contains("violation")) {
        cc.violation = violationInfoFromJson(j["violation"]);
    }
    if (j.contains("left_edge_id")) {
        cc.left_edge_id = j["left_edge_id"].get<std::string>();
    }
    if (j.contains("right_edge_id")) {
        cc.right_edge_id = j["right_edge_id"].get<std::string>();
    }

    return cc;
}

// 序列化和反序列化 ConnectivityLaneEdge
inline nlohmann::json connectivityLaneEdgeToJson(const ConnectivityLaneEdge& cle) {
    nlohmann::json j = nlohmann::json::object();
    j["id"] = cle.id;
    j["geometry"] = lineString2dToJson(cle.geometry);
    j["is_shared"] = cle.is_shared;
    j["attrs"] = detail::mapToAttrMapJson(cle.attrs);
    j["groupId"] = cle.groupId;
    j["lineOrder"] = cle.lineOrder;

    if (cle.shared_by) {
        nlohmann::json shared_by_obj = nlohmann::json::object();
        shared_by_obj["first"] = cle.shared_by->first;
        shared_by_obj["second"] = cle.shared_by->second;
        j["shared_by"] = shared_by_obj;
    }

    return j;
}

inline ConnectivityLaneEdge connectivityLaneEdgeFromJson(const nlohmann::json& j) {
    ConnectivityLaneEdge cle;
    if (j.contains("id")) {
        cle.id = j["id"].get<std::string>();
    }
    if (j.contains("geometry")) {
        cle.geometry = lineString2dFromJson(j["geometry"]);
    }
    if (j.contains("is_shared")) {
        cle.is_shared = j["is_shared"].get<bool>();
    }
    if (j.contains("attrs")) {
        cle.attrs = detail::attrMapFromJson(j["attrs"]);
    }
    if (j.contains("groupId")) {
        cle.groupId = j["groupId"].get<std::string>();
    }
    if (j.contains("lineOrder")) {
        cle.lineOrder = j["lineOrder"].get<int>();
    }

    if (j.contains("shared_by")) {
        const nlohmann::json& shared_data = j["shared_by"];
        auto shared_pair = std::make_shared<std::pair<LaneId, LaneId>>();
        if (shared_data.contains("first")) {
            shared_pair->first = shared_data["first"].get<std::string>();
        }
        if (shared_data.contains("second")) {
            shared_pair->second = shared_data["second"].get<std::string>();
        }
        cle.shared_by = shared_pair;
    }

    return cle;
}

// 序列化和反序列化 IntersectionOutput::PerfStats
inline nlohmann::json perfStatsToJson(const IntersectionOutput::PerfStats& ps) {
    nlohmann::json j = nlohmann::json::object();
    j["sdf_build_ms"] = ps.sdf_build_ms;
    j["precheck_ms"] = ps.precheck_ms;
    j["optimize_ms"] = ps.optimize_ms;
    j["smooth_ms"] = ps.smooth_ms;
    j["edge_gen_ms"] = ps.edge_gen_ms;
    j["area_gen_ms"] = ps.area_gen_ms;

    return j;
}

inline IntersectionOutput::PerfStats perfStatsFromJson(const nlohmann::json& j) {
    IntersectionOutput::PerfStats ps;
    if (j.contains("sdf_build_ms")) {
        ps.sdf_build_ms = j["sdf_build_ms"].get<double>();
    }
    if (j.contains("precheck_ms")) {
        ps.precheck_ms = j["precheck_ms"].get<double>();
    }
    if (j.contains("optimize_ms")) {
        ps.optimize_ms = j["optimize_ms"].get<double>();
    }
    if (j.contains("smooth_ms")) {
        ps.smooth_ms = j["smooth_ms"].get<double>();
    }
    if (j.contains("edge_gen_ms")) {
        ps.edge_gen_ms = j["edge_gen_ms"].get<double>();
    }
    if (j.contains("area_gen_ms")) {
        ps.area_gen_ms = j["area_gen_ms"].get<double>();
    }

    return ps;
}

// 序列化和反序列化 IntersectionOutput
inline nlohmann::json intersectionOutputToJson(const IntersectionOutput& io) {
    nlohmann::json j = nlohmann::json::object();

    nlohmann::json curves = nlohmann::json::array();
    for (const auto& curve : io.connectivity_curves) {
        curves.push_back(connectivityCurveToJson(curve));
    }
    j["connectivity_curves"] = curves;

    nlohmann::json edges = nlohmann::json::array();
    for (const auto& edge : io.lane_edges) {
        edges.push_back(connectivityLaneEdgeToJson(edge));
    }
    j["lane_edges"] = edges;

    j["area"] = intersectionAreaToJson(io.area);
    j["perf"] = perfStatsToJson(io.perf);

    return j;
}

inline IntersectionOutput intersectionOutputFromJson(const nlohmann::json& j) {
    IntersectionOutput io;

    if (j.contains("connectivity_curves")) {
        const nlohmann::json& curves = j["connectivity_curves"];
        for (size_t i = 0; i < curves.size(); ++i) {
            io.connectivity_curves.push_back(connectivityCurveFromJson(curves[i]));
        }
    }

    if (j.contains("lane_edges")) {
        const nlohmann::json& edges = j["lane_edges"];
        for (size_t i = 0; i < edges.size(); ++i) {
            io.lane_edges.push_back(connectivityLaneEdgeFromJson(edges[i]));
        }
    }

    if (j.contains("area")) {
        io.area = intersectionAreaFromJson(j["area"]);
    }
    if (j.contains("perf")) {
        io.perf = perfStatsFromJson(j["perf"]);
    }

    return io;
}

// 主要的IntersectionInput序列化函数
inline nlohmann::json intersectionInputToJson(const IntersectionInput& input) {
    nlohmann::json j = nlohmann::json::object();

    // 序列化lane_groups
    nlohmann::json lane_groups = nlohmann::json::array();
    for (const auto& lg : input.lane_groups) {
        lane_groups.push_back(laneGroupToJson(lg));
    }
    j["lane_groups"] = lane_groups;

    // 序列化lanes
    nlohmann::json lanes = nlohmann::json::array();
    for (const auto& lane : input.lanes) {
        lanes.push_back(laneToJson(lane));
    }
    j["lanes"] = lanes;

    // 序列化lane_edges
    nlohmann::json lane_edges = nlohmann::json::array();
    for (const auto& edge : input.lane_edges) {
        lane_edges.push_back(laneEdgeToJson(edge));
    }
    j["lane_edges"] = lane_edges;

    // 序列化connectivities
    nlohmann::json connectivities = nlohmann::json::array();
    for (const auto& conn : input.connectivities) {
        connectivities.push_back(connectivityToJson(conn));
    }
    j["connectivities"] = connectivities;

    // 序列化obstacles
    nlohmann::json obstacles = nlohmann::json::array();
    for (const auto& obs : input.obstacles) {
        obstacles.push_back(obstacleToJson(obs));
    }
    j["obstacles"] = obstacles;

    // 序列化boundaries
    nlohmann::json boundaries = nlohmann::json::array();
    for (const auto& bd : input.boundaries) {
        boundaries.push_back(boundaryToJson(bd));
    }
    j["boundaries"] = boundaries;

    // 序列化stop_lines
    nlohmann::json stop_lines = nlohmann::json::array();
    for (const auto& sl : input.stop_lines) {
        stop_lines.push_back(stopLineToJson(sl));
    }
    j["stop_lines"] = stop_lines;

    // 序列化crosswalks
    nlohmann::json crosswalks = nlohmann::json::array();
    for (const auto& cw : input.crosswalks) {
        crosswalks.push_back(crosswalkToJson(cw));
    }
    j["crosswalks"] = crosswalks;

    // 序列化area
    j["area"] = intersectionAreaToJson(input.area);

    // 序列化id
    j["id"] = input.id;

    // 序列化mode
    j["mode"] = input.mode;

    return j;
}

// 主要的IntersectionInput反序列化函数
inline IntersectionInput intersectionInputFromJson(const nlohmann::json& j) {
    IntersectionInput input;

    // 反序列化lane_groups
    if (j.contains("lane_groups")) {
        const nlohmann::json& lane_groups = j["lane_groups"];
        for (size_t i = 0; i < lane_groups.size(); ++i) {
            auto lanegroup = laneGroupFromJson(lane_groups[i]);
            auto it_lanegroup = std::find_if(input.lane_groups.begin(), input.lane_groups.end(), [lanegroup](isg::LaneGroup& l){
                return l.id == lanegroup.id;
            });
            if (it_lanegroup == input.lane_groups.end()) {
                input.lane_groups.push_back(lanegroup);
            }
        }
    }

    // 反序列化lanes
    if (j.contains("lanes")) {
        const nlohmann::json& lanes = j["lanes"];
        for (size_t i = 0; i < lanes.size(); ++i) {
            auto lane = laneFromJson(lanes[i]);
            auto it_lane = std::find_if(input.lanes.begin(), input.lanes.end(), [lane](isg::Lane& l){
                return l.id == lane.id;
            });
            if (it_lane == input.lanes.end()) {
                input.lanes.push_back(lane);
            }
        }
    }

    // 反序列化lane_edges
    if (j.contains("lane_edges")) {
        const nlohmann::json& lane_edges = j["lane_edges"];
        for (size_t i = 0; i < lane_edges.size(); ++i) {
            auto lane_edge = laneEdgeFromJson(lane_edges[i]);
            auto it_lane_edge = std::find_if(input.lane_edges.begin(), input.lane_edges.end(), [lane_edge](isg::LaneEdge& l){
                return l.id == lane_edge.id;
            });
            if (it_lane_edge == input.lane_edges.end()) {
                input.lane_edges.push_back(lane_edge);
            }
        }
    }

    // 反序列化connectivities
    if (j.contains("connectivities")) {
        const nlohmann::json& connectivities = j["connectivities"];
        for (size_t i = 0; i < connectivities.size(); ++i) {
            auto conn = connectivityFromJson(connectivities[i]);
            auto it_conn = std::find_if(input.connectivities.begin(), input.connectivities.end(), [conn](isg::Connectivity& l){
                return l.id == conn.id;
            });
            if (it_conn == input.connectivities.end()) {
                input.connectivities.push_back(conn);
            }
        }
    }

    // 反序列化obstacles
    if (j.contains("obstacles")) {
        const nlohmann::json& obstacles = j["obstacles"];
        for (size_t i = 0; i < obstacles.size(); ++i) {
            auto obstacle = obstacleFromJson(obstacles[i]);
            auto it_obstacle = std::find_if(input.obstacles.begin(), input.obstacles.end(), [obstacle](isg::Obstacle& l){
                return l.id == obstacle.id;
            });
            if (it_obstacle == input.obstacles.end()) {
                input.obstacles.push_back(obstacle);
            }
        }
    }

    // 反序列化boundaries
    if (j.contains("boundaries")) {
        const nlohmann::json& boundaries = j["boundaries"];
        for (size_t i = 0; i < boundaries.size(); ++i) {
            auto boundary = boundaryFromJson(boundaries[i]);
            auto it_boundary = std::find_if(input.boundaries.begin(), input.boundaries.end(), [boundary](isg::Boundary& l){
                return l.id == boundary.id;
            });
            if (it_boundary == input.boundaries.end()) {
                input.boundaries.push_back(boundary);
            }
        }
    }

    // 反序列化stop_lines
    if (j.contains("stop_lines")) {
        const nlohmann::json& stop_lines = j["stop_lines"];
        for (size_t i = 0; i < stop_lines.size(); ++i) {
            auto stopline = stopLineFromJson(stop_lines[i]);
            auto it_stopline = std::find_if(input.stop_lines.begin(), input.stop_lines.end(), [stopline](isg::StopLine& l){
                return l.id == stopline.id;
            });
            if (it_stopline == input.stop_lines.end()) {
                input.stop_lines.push_back(stopline);
            }
        }
    }

    // 反序列化crosswalks
    if (j.contains("crosswalks")) {
        const nlohmann::json& crosswalks = j["crosswalks"];
        for (size_t i = 0; i < crosswalks.size(); ++i) {
            auto crosswalk = crosswalkFromJson(crosswalks[i]);
            auto it_crosswalk = std::find_if(input.crosswalks.begin(), input.crosswalks.end(), [crosswalk](isg::Crosswalk& l){
                return l.id == crosswalk.id;
            });
            if (it_crosswalk == input.crosswalks.end()) {
                input.crosswalks.push_back(crosswalk);
            }
        }
    }

    // 反序列化area
    if (j.contains("area")) {
        input.area = intersectionAreaFromJson(j["area"]);
    }

    // 反序列化id
    if (j.contains("id")) {
        input.id = j["id"].get<std::string>();
    }

    // 反序列化mode
    if (j.contains("mode")) {
        input.mode = j["mode"].get<int>();
    }
    return input;
}

// 文件和字符串门面统一复用同一套字段编解码，避免两条路径产生格式差异。
void IntersectionIO::saveToFile(const IntersectionInput& input,
                                const std::string& filepath) {
    const nlohmann::json json = intersectionInputToJson(input);
    std::ofstream file(filepath);
    if (!file.is_open())
        throw std::runtime_error("Cannot open file for writing: " + filepath);
    file << json.dump(4);  // 使用四空格缩进，保持既有输出格式。
}

IntersectionInput IntersectionIO::loadFromFile(const std::string& filepath) {
    std::ifstream file(filepath);
    if (!file.is_open())
        throw std::runtime_error("Cannot open file for reading: " + filepath);
    const std::string content((std::istreambuf_iterator<char>(file)),
                              std::istreambuf_iterator<char>());
    return intersectionInputFromJson(nlohmann::json::parse(content));
}

std::string IntersectionIO::toJsonString(const IntersectionInput& input) {
    return intersectionInputToJson(input).dump(4);
}

IntersectionInput IntersectionIO::fromJsonString(
    const std::string& json_string) {
    return intersectionInputFromJson(nlohmann::json::parse(json_string));
}

}  // 命名空间 isg
