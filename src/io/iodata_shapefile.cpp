#include "iodata_shapefile.h"
#include "utils/filesystem.hpp"
#include "utils/shapefile.hpp"
#include <algorithm>
#include <vector>

namespace isg {
namespace {

namespace fs = ghc::filesystem;

const std::vector<IsgDbfField> kLaneFields = {
        {"ID", 'C', 64, 0},
        {"LANE_ORDER", 'C', 64, 0},
        {"LEDGE", 'C', 64, 0},
        {"REDGE", 'C', 64, 0},
        {"GROUP_ID", 'C', 64, 0},
        {"GROUP_TYPE", 'C', 64, 0},
};
const std::vector<IsgDbfField> kConnectivityFields = {
        {"ID", 'C', 64, 0},
        {"TURN_TYPE", 'C', 64, 0},
        {"FLANE", 'C', 64, 0},
        {"TLANE", 'C', 64, 0},
        {"LANE_TYPE", 'C', 64, 0},
        {"FIXED_SHAPE", 'C', 64, 0},
};
const std::vector<IsgDbfField> kLaneEdgeFields = {
        {"ID", 'C', 64, 0},
        {"GROUP_ID", 'C', 64, 0},
        {"GROUP_TYPE", 'C', 64, 0},
        {"LEFT_CLINE_ID", 'C', 64, 0},
        {"RIGHT_CLINE_ID", 'C', 64, 0},
};
const std::vector<IsgDbfField> kIdTypeFields = {
        {"ID", 'C', 64, 0},
        {"TYPE", 'C', 64, 0},
};
const std::vector<IsgDbfField> kStopLineFields = {
        {"ID", 'C', 64, 0},
        {"ENTRY_GROUP_ID", 'C', 64, 0},
};
const std::vector<IsgDbfField> kControlPointFields = {
        {"ID", 'C', 64, 0},
        {"NUM_SEGS", 'C', 64, 0},
        {"TURN_TYPE", 'C', 64, 0},
        {"TLANE", 'C', 64, 0},
        {"LANE_TYPE", 'C', 64, 0},
        {"FIXED_SHAPE", 'C', 64, 0},
};

std::vector<ShapePoint> toShapePoints(const std::vector<Vec2d>& points) {
    std::vector<ShapePoint> result;
    result.reserve(points.size());
    for (const auto& point : points)
        result.push_back({point[0], point[1], 0.0});
    return result;
}

std::vector<ShapePoint> toShapePoints(const std::vector<Vec3d>& points) {
    std::vector<ShapePoint> result;
    result.reserve(points.size());
    for (const auto& point : points)
        result.push_back({point[0], point[1], point[2]});
    return result;
}

template <typename T, typename AttributeBuilder, typename GeometryBuilder>
void writeShapes(const std::string& dir,
                 const std::string& name,
                 int shape_type,
                 const std::vector<IsgDbfField>& fields,
                 const std::vector<T>& items,
                 AttributeBuilder attributes,
                 GeometryBuilder geometry) {
    fs::create_directories(dir);
    std::vector<ShapeRecord> records;
    records.reserve(items.size());
    for (std::size_t i = 0; i < items.size(); ++i) {
        std::vector<ShapePoint> points;
        if (!geometry(items[i], points))
            continue;
        records.emplace_back(static_cast<int>(records.size()), shape_type,
                             attributes(items[i], i), points, std::vector<int32_t>{0});
    }
    ShapefileEngine::write(dir, name, shape_type, fields, records);
}

template <typename T>
bool lineGeometry(const T& value, std::vector<ShapePoint>& points) {
    points = toShapePoints(value.geometry.points);
    return true;
}

template <typename T>
bool polygonGeometry(const T& value, std::vector<ShapePoint>& points) {
    points = toShapePoints(value.geometry.outer);
    return true;
}

void writeArea(const std::string& dir, const std::string& name, const IntersectionArea& area) {
    const std::vector<Polygon2d> areas{area.geometry};
    writeShapes(dir, name, SHP_POLYGONZ, kIdTypeFields, areas,
        [](const Polygon2d&, std::size_t i) {
            return std::vector<std::string>{std::to_string(i), "-1"};
        },
        [](const Polygon2d& polygon, std::vector<ShapePoint>& points) {
            points = toShapePoints(polygon.outer);
            return true;
        });
}

}  // namespace

bool save(IntersectionInput& input, std::string out_dir, std::string prefix) {
    writeShapes(out_dir, prefix + "_lanes", SHP_POLYLINEZ, kLaneFields, input.lanes,
        [](const Lane& lane, std::size_t) {
            return std::vector<std::string>{lane.id, std::to_string(lane.laneOrder),
                lane.left_edge_id, lane.right_edge_id, lane.groupId, ""};
        }, lineGeometry<Lane>);

    writeShapes(out_dir, prefix + "_connectivities", SHP_POLYLINEZ, kConnectivityFields,
        input.connectivities,
        [](const Connectivity& connectivity, std::size_t) {
            return std::vector<std::string>{connectivity.id,
                std::to_string(static_cast<int>(connectivity.turn_type)), connectivity.entry_lane_id,
                connectivity.exit_lane_id, std::to_string(static_cast<int>(connectivity.lane_type)),
                connectivity.fixed_shape ? "1" : "0"};
        },
        [](const Connectivity& connectivity, std::vector<ShapePoint>& points) {
            if (connectivity.fixed_shape && !connectivity.geometry.points.empty())
                points = toShapePoints(connectivity.geometry.points);
            return true;
        });

    writeShapes(out_dir, prefix + "_laneedges", SHP_POLYLINEZ, kLaneEdgeFields, input.lane_edges,
        [](const LaneEdge& edge, std::size_t) {
            return std::vector<std::string>{edge.id, "-1", edge.groupId, "", ""};
        }, lineGeometry<LaneEdge>);

    writeShapes(out_dir, prefix + "_stoplines", SHP_POLYLINEZ, kStopLineFields, input.stop_lines,
        [](const StopLine& line, std::size_t) {
            return std::vector<std::string>{line.id, line.associated_group_id};
        }, lineGeometry<StopLine>);

    writeShapes(out_dir, prefix + "_roadedges", SHP_POLYLINEZ, kIdTypeFields, input.boundaries,
        [](const Boundary& boundary, std::size_t) {
            return std::vector<std::string>{boundary.id,
                std::to_string(static_cast<int>(boundary.type))};
        }, lineGeometry<Boundary>);

    writeShapes(out_dir, prefix + "_obstacles", SHP_POLYGONZ, kIdTypeFields, input.obstacles,
        [](const Obstacle& obstacle, std::size_t) {
            return std::vector<std::string>{obstacle.id, "-1"};
        }, polygonGeometry<Obstacle>);

    writeShapes(out_dir, prefix + "_crosswalks", SHP_POLYGONZ, kIdTypeFields, input.crosswalks,
        [](const Crosswalk& crosswalk, std::size_t) {
            return std::vector<std::string>{crosswalk.id, "-1"};
        }, polygonGeometry<Crosswalk>);

    writeArea(out_dir, prefix + "_areas", input.area);
    return true;
}

bool save(IntersectionOutput& output, std::string out_dir, std::string prefix) {
#ifdef PROJECT_ROOT_DIR
    std::string project_dir = PROJECT_ROOT_DIR;
    const std::string qgis_name = "intersection_gen.qgs";
    std::string qgis_template = project_dir + "/datas/" + qgis_name;
    if (!fs::exists(qgis_template))
        qgis_template = project_dir + "/" + qgis_name;
    if (fs::exists(qgis_template)) {
        const std::string destination =
            fs::path(out_dir).parent_path().parent_path().string() + "/" + qgis_name;
        fs::copy_file(qgis_template, destination, fs::copy_options::skip_existing);
    }
#endif

    writeShapes(out_dir, prefix + "_ctrlpts", SHP_POLYLINEZ, kControlPointFields,
        output.connectivity_curves,
        [](const ConnectivityCurve& curve, std::size_t) {
            const int segment_count = curve.curve ? curve.curve->numSegments() : 0;
            return std::vector<std::string>{curve.id, std::to_string(segment_count),
                std::to_string(static_cast<int>(curve.turn_type)), curve.entry_lane_id,
                curve.exit_lane_id, std::to_string(static_cast<int>(curve.lane_type)),
                curve.fixed_shape ? "1" : "0"};
        },
        [](const ConnectivityCurve& curve, std::vector<ShapePoint>& points) {
            if (curve.curve && curve.curve->numSegments() > 0) {
                for (const auto& segment : curve.curve->segs) {
                    for (const auto& point : segment.ctrl) {
                        const ShapePoint shape_point{point[0], point[1], 0.0};
                        if (std::find(points.begin(), points.end(), shape_point) == points.end())
                            points.push_back(shape_point);
                    }
                }
            }
            if (points.empty() && !curve.geometry.points.empty())
                points = toShapePoints(curve.geometry.points);
            return true;
        });

    writeShapes(out_dir, prefix + "_lanes", SHP_POLYLINEZ, kConnectivityFields,
        output.connectivity_curves,
        [](const ConnectivityCurve& curve, std::size_t) {
            return std::vector<std::string>{curve.id,
                std::to_string(static_cast<int>(curve.turn_type)), curve.entry_lane_id,
                curve.exit_lane_id, std::to_string(static_cast<int>(curve.lane_type)),
                curve.fixed_shape ? "1" : "0"};
        },
        [](const ConnectivityCurve& curve, std::vector<ShapePoint>& points) {
            if (!curve.geometry.points.empty()) {
                points = toShapePoints(curve.geometry.points);
                return true;
            }
            if (!curve.curve || curve.curve->numSegments() == 0)
                return false;
            points = toShapePoints(curve.curve->sampleByShape());
            return true;
        });

    writeShapes(out_dir, prefix + "_laneedges", SHP_POLYLINEZ, kLaneEdgeFields,
        output.lane_edges,
        [](const ConnectivityLaneEdge& edge, std::size_t) {
            return std::vector<std::string>{edge.id, "-1", "", "", ""};
        }, lineGeometry<ConnectivityLaneEdge>);

    writeArea(out_dir, prefix + "_areas", output.area);
    return true;
}

}  // namespace isg
