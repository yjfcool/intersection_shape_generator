#pragma once
// ─────────────────────────────────────────────────────────────────────────────
// ScenarioBuilder: construct IntersectionInput for various intersection types.
//
// Lane geometry convention (from types.h):
//   Entry lane: outside → junction edge   (last pt = junction entry point)
//   Exit  lane: junction edge → outside   (first pt = junction exit point)
//
// Right-hand-traffic positioning:
//   Entry lanes are on the RIGHT of entry_dir (right_perp = (dir.y, -dir.x))
//   Exit  lanes are on the LEFT  of entry_dir
// ─────────────────────────────────────────────────────────────────────────────
#include "types.h"
#include <sstream>
#include <cmath>

using namespace isg;

// ── helpers ──────────────────────────────────────────────────────────────────
inline Vec2d rightPerp(const Vec2d& dir) {
    return Vec2d(dir[1], -dir.x()); // rotate -90° (right side)
}
inline Polygon2d makeRectObstacle(Vec2d centre, double hw, double hh) {
    Polygon2d p;
    p.outer = { {centre[0]-hw, centre[1]-hh}, {centre[0]+hw, centre[1]-hh},
                {centre[0]+hw, centre[1]+hh}, {centre[0]-hw, centre[1]+hh} };
    return p;
}
inline Polygon2d makeCircleObstacle(Vec2d centre, double r, int segs=16) {
    Polygon2d p;
    for (int i=0;i<segs;++i) {
        double a = 2*M_PI*i/segs;
        p.outer.push_back(Vec3d(centre + Vec2d(r*std::cos(a), r*std::sin(a))));
    }
    return p;
}

// ─────────────────────────────────────────────────────────────────────────────
// ArmDef: one road arm attached to the junction
// ─────────────────────────────────────────────────────────────────────────────
struct ArmDef {
    std::string name;       // e.g. "N", "S", "E", "W"
    Vec2d  junction_pt;     // arm–junction connection point
    Vec2d  entry_dir;       // unit vector pointing INTO junction
    int    n_entry = 1;     // number of entry (inbound) lanes
    int    n_exit  = 1;     // number of exit  (outbound) lanes
    double lane_w  = 3.5;   // lane width (m)
    double arm_len = 15.0;  // arm length beyond junction edge (m)
};

// Build one arm: populate lane_groups, lanes, lane_edges of `inp`
// Returns: {all_entry_lane_ids, all_exit_lane_ids}
inline std::pair<std::vector<LaneId>, std::vector<LaneId>>
buildArm(const ArmDef& arm, IntersectionInput& inp)
{
    Vec2d rp = rightPerp(arm.entry_dir);   // right perpendicular = entry side
    Vec2d exit_dir = -arm.entry_dir;       // direction out of junction for exit lanes

    auto mkLaneId = [&](bool entry, int i) {
        return arm.name + (entry?"_E_":"_X_") + std::to_string(i);
    };
    auto mkEdgeId = [&](bool entry, int i) {
        return arm.name + (entry?"_EedgeB":"_XedgeB") + std::to_string(i);
    };

    // ── Entry group ──────────────────────────────────────────────────────────
    LaneGroup eg;
    eg.id   = arm.name + "_ENTRY";
    eg.role = GroupRole::Entry;
    //eg.direction  = arm.entry_dir;                 // retained but unused by algorithm
    //eg.ref_point  = arm.junction_pt;

    std::vector<LaneId> entry_ids, exit_ids;

    for (int i = 0; i < arm.n_entry; ++i) {
        double offset = arm.lane_w * (i + 0.5);   // RHT: positive rp = right side
        Vec2d centre_jct = arm.junction_pt + offset * rp;
        Vec2d centre_out = centre_jct - arm.entry_dir * arm.arm_len; // outside end

        Lane l;
        l.id    = mkLaneId(true, i);
        l.width = arm.lane_w;
        l.groupId = eg.id;
        l.geometry.points = { centre_out, centre_jct };   // outside → junction
        eg.lanes.push_back(l.id);
        inp.lanes.push_back(l);
        entry_ids.push_back(l.id);
    }
    // N+1 boundary IDs
    for (int i=0; i<=arm.n_entry; ++i) eg.boundaries.push_back(mkEdgeId(true,i));
    inp.lane_groups.push_back(eg);

    // ── Exit group ───────────────────────────────────────────────────────────
    LaneGroup xg;
    xg.id   = arm.name + "_EXIT";
    xg.role = GroupRole::Exit;
    //xg.direction = exit_dir;
    //xg.ref_point = arm.junction_pt;

    for (int i = 0; i < arm.n_exit; ++i) {
        double offset = arm.lane_w * (i + 0.5);   // left side = -rp
        Vec2d centre_jct = arm.junction_pt - offset * rp;
        Vec2d centre_out = centre_jct + exit_dir * arm.arm_len;

        Lane l;
        l.id    = mkLaneId(false, i);
        l.width = arm.lane_w;
        l.groupId = xg.id;
        l.geometry.points = { centre_jct, centre_out };   // junction → outside
        xg.lanes.push_back(l.id);
        inp.lanes.push_back(l);
        exit_ids.push_back(l.id);
    }
    for (int i=0; i<=arm.n_exit; ++i) xg.boundaries.push_back(mkEdgeId(false,i));
    inp.lane_groups.push_back(xg);

    return {entry_ids, exit_ids};
}

// ─────────────────────────────────────────────────────────────────────────────
// addConn: add a connectivity record
// ─────────────────────────────────────────────────────────────────────────────
static int gConnIdx = 0;
inline void addConn(IntersectionInput& inp,
                    const LaneId& entry, const LaneId& exit,
                    ConnTurnType tt)
{
    auto findGroupForLane = [&](const LaneId& lane_id, GroupRole role) {
        for (const auto& g : inp.lane_groups) {
            if (g.role != role) continue;
            for (const auto& lid : g.lanes)
                if (lid == lane_id)
                    return g.id;
        }
        return LaneGroupId{};
    };

    Connectivity c;
    c.id            = "C" + std::to_string(gConnIdx++);
    c.entry_lane_id = entry;
    c.exit_lane_id  = exit;
    c.turn_type     = tt;
    c.enterGroupId  = findGroupForLane(entry, GroupRole::Entry);
    c.exitGroupId   = findGroupForLane(exit, GroupRole::Exit);
    inp.connectivities.push_back(c);
}

// ─────────────────────────────────────────────────────────────────────────────
// Standard square coarse area and road-edge boundaries
// ─────────────────────────────────────────────────────────────────────────────
inline void addStandardArea(IntersectionInput& inp, double half = 12.0) {
    Polygon2d fence;
    fence.outer = { {-half,-half},{half,-half},{half,half},{-half,half} };
    inp.area.geometry = fence;
}
inline void addRoadEdgeBoundary(IntersectionInput& inp,
                                 const std::string& id,
                                 Vec2d a, Vec2d b) {
    Boundary bnd;
    bnd.id   = id;
    bnd.type = Boundary::Type::RoadEdge;
    bnd.geometry.points = {a, b};
    inp.boundaries.push_back(bnd);
}

// ─────────────────────────────────────────────────────────────────────────────
// addObstacle helper
// ─────────────────────────────────────────────────────────────────────────────
inline void addObstacle(IntersectionInput& inp,
                         const std::string& id,
                         const Polygon2d& poly) {
    Obstacle obs;
    obs.id       = id;
    obs.geometry = poly;
    inp.obstacles.push_back(obs);
}
