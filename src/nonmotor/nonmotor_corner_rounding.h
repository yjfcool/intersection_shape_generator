#pragma once

#include "types.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <vector>

namespace isg {
namespace nonmotor_detail {

/// 圆角的三维控制点；平面弧长用于采样，高程随同一个参数插值。
struct CornerBezier {
    std::array<Vec3d, 4> ctrl;
    Vec3d evaluate(double t) const {
        const double u = 1.0 - t;
        return ctrl[0] * (u*u*u) + ctrl[1] * (3*u*u*t) +
               ctrl[2] * (3*u*t*t) + ctrl[3] * (t*t*t);
    }
    Vec3d tangent(double t) const {
        const double u = 1.0 - t;
        return (ctrl[1] - ctrl[0]) * (3*u*u) +
               (ctrl[2] - ctrl[1]) * (6*u*t) + (ctrl[3] - ctrl[2]) * (3*t*t);
    }
    Vec3d acceleration(double t) const {
        return (ctrl[2]-ctrl[1]*2+ctrl[0])*(6*(1-t)) +
               (ctrl[3]-ctrl[2]*2+ctrl[1])*(6*t);
    }
};

/// 对解析导数作有界细分审计；拒绝驻点和额外反向小圈，不以输出采样密度代替曲率检查。
inline double attachmentBezierCurvature(const CornerBezier& curve) {
    double maximum = 0, turn = 0;
    Vec2d previous = curve.tangent(0).xy();
    for (int i=0; i<=128; ++i) {
        const double t=i/128.0;
        const Vec2d d=curve.tangent(t).xy(), dd=curve.acceleration(t).xy();
        if (d.norm()<1e-6) return std::numeric_limits<double>::infinity();
        maximum=std::max(maximum,std::abs(d.x()*dd.y()-d.y()*dd.x())/std::pow(d.norm(),3));
        if (i) turn+=std::abs(std::atan2(previous.x()*d.y()-previous.y()*d.x(),previous.dot(d)));
        previous=d;
    }
    return turn > 1.5*PI ? std::numeric_limits<double>::infinity() : maximum;
}

/// 只删除连续重复点，圆角拼接和不足一个采样间距的末段均保留精确端点。
inline void appendCornerPoint(std::vector<Vec3d>& points, const Vec3d& point) {
    if (points.empty() || (points.back() - point).norm() > 1e-9) points.push_back(point);
}

/// 相邻偏移段以两条支撑线的交点拼接，保持两侧直段平行；原先端点均值会把折角拉向原骨干。
/// 近平行或过长尖角采用均值回退，最大拼接移动为原始法向偏移的4倍，避免无界尖刺。
inline Vec3d joinOffsetSegments(const Vec3d& before, const Vec3d& incoming_end,
    const Vec3d& outgoing_start, const Vec3d& after, const Vec3d& vertex,
    bool require_forward = false) {
    const Vec2d first = incoming_end.xy()-before.xy(), last = after.xy()-outgoing_start.xy();
    const Vec3d middle = (incoming_end+outgoing_start)*0.5;
    const double determinant = first.x()*last.y()-first.y()*last.x();
    if (first.norm() < 1e-9 || last.norm() < 1e-9 || std::abs(determinant) < 1e-8*first.norm()*last.norm())
        return middle;
    const Vec2d delta = outgoing_start.xy()-incoming_end.xy();
    const double station = (delta.x()*last.y()-delta.y()*last.x())/determinant;
    // 前方以完整相邻段的起点/终点为基准：交点可从incoming_end回退，
    // 但不能越过before；从交点出发也不能越过after再折返。
    // 旧门禁station>=0误拒绝正常内侧miter，均值回退会破坏平行轨位。
    const double outgoing_station = (delta.x()*first.y()-delta.y()*first.x())/determinant;
    const Vec2d point = incoming_end.xy()+first*station;
    const double offset = std::max((incoming_end.xy()-vertex.xy()).norm(), (outgoing_start.xy()-vertex.xy()).norm());
    if (!std::isfinite(station) || !std::isfinite(outgoing_station) ||
        (require_forward && (station <= -1.0+1e-7 || outgoing_station >= 1.0-1e-7)) ||
        (point-vertex.xy()).norm() > 4*std::max(offset, 1e-7)) return middle;
    return Vec3d(point, middle.z());
}

/// 用两侧截取点及前驱/后继方向拟合三次贝塞尔。
/// 等截距时采用圆弧的标准把手比例；不等截距仍保持两端G1且控制点位于折角三角形内。
inline CornerBezier fitCornerBezier(const Vec3d& before, const Vec3d& vertex, const Vec3d& after) {
    const Vec2d incoming = vertex.xy() - before.xy(), outgoing = after.xy() - vertex.xy();
    const double cosine = std::max(-1.0, std::min(1.0, incoming.normalized().dot(outgoing.normalized())));
    const double half_cosine = std::sqrt(0.5 * (1.0 + cosine));
    const double ratio = (4.0 / 3.0) * half_cosine / (1.0 + half_cosine);
    CornerBezier curve;
    curve.ctrl = {{before, before + (vertex - before) * ratio,
                   after + (vertex - after) * ratio, after}};
    return curve;
}

/// 平面弧长近似等距采样：按控制多边形长度/0.01选择等参数查找折线段数（64至16384），
/// 再按查找折线弧长取0.2米站距；该分辨率不代表每段弦长严格为1厘米。
/// 每个圆角独立起算，尾部不足0.2米时补精确终点，Z沿查找折线同步插值。
inline std::vector<Vec3d> sampleCornerBezier(const CornerBezier& curve) {
    double polygon_length = 0.0;
    for (size_t i = 1; i != 4; ++i) polygon_length += (curve.ctrl[i].xy() - curve.ctrl[i-1].xy()).norm();
    const int count = static_cast<int>(std::max(64.0, std::min(16384.0, std::ceil(polygon_length / 0.01))));
    std::vector<Vec3d> points;
    appendCornerPoint(points, curve.ctrl[0]);
    Vec3d previous = curve.ctrl[0];
    double traveled = 0.0, next_sample = 0.2;
    for (int i = 1; i <= count; ++i) {
        const Vec3d current = curve.evaluate(static_cast<double>(i) / count);
        const double length = (current.xy() - previous.xy()).norm();
        if (length > 1e-12) {
            while (next_sample <= traveled + length) {
                const double f = (next_sample - traveled) / length;
                appendCornerPoint(points, previous + (current - previous) * f);
                next_sample += 0.2;
            }
        }
        traveled += length; previous = current;
    }
    appendCornerPoint(points, curve.ctrl[3]);
    return points;
}

/// 密集采样路径的局部尖折修复。只替换尖折前后radius米弧长内的点，其他点按原序完整保留。
/// 截取点按弧长在线段内部插值，不把跨越截取位置的整条长段当作圆角范围；
/// 相邻尖折的修复窗口在两折点弧长中点处分界。两端控制柄沿保留路径的切向，保持拼接G1。
/// 无可用截距或精确反向的退化折点保持原样，交给最终自交和空间约束审计。
inline LineString2d smoothSharpCorners(const LineString2d& source, double radius) {
    if (radius <= 0.0 || source.points.size() < 3) return source;
    std::vector<Vec3d> points;
    for (const auto& point : source.points) appendCornerPoint(points, point);
    if (points.size() < 3) return source;
    std::vector<double> stations(points.size(), 0.0);
    std::vector<size_t> corners;
    for (size_t i = 1; i < points.size(); ++i)
        stations[i] = stations[i-1] + (points[i].xy()-points[i-1].xy()).norm();
    for (size_t i = 1; i+1 < points.size(); ++i) {
        const Vec2d before = points[i].xy()-points[i-1].xy();
        const Vec2d after = points[i+1].xy()-points[i].xy();
        if (before.norm() > 1e-7 && after.norm() > 1e-7 &&
            before.normalized().dot(after.normalized()) < 0.70) corners.push_back(i);
    }
    if (corners.empty()) return source;
    // 截取点恰好落在已有顶点时，左端采用保留前缀的入切向，右端采用保留后缀的出切向。
    const auto sample = [&](double station, bool incoming, Vec2d& tangent) {
        size_t next = static_cast<size_t>(std::lower_bound(stations.begin(), stations.end(), station)-stations.begin());
        if (!incoming && next < stations.size() && std::abs(stations[next]-station) < 1e-9) ++next;
        next = std::max(size_t(1), std::min(next, points.size()-1));
        while (next+1 < points.size() && stations[next]-stations[next-1] < 1e-9) ++next;
        const double length = stations[next]-stations[next-1];
        const double t = length > 1e-9 ? std::max(0.0, std::min(1.0, (station-stations[next-1])/length)) : 0.0;
        tangent = length > 1e-9 ? Vec2d((points[next].xy()-points[next-1].xy())/length) : Vec2d(0, 0);
        return Vec3d(points[next-1]+(points[next]-points[next-1])*t);
    };
    LineString2d result;
    size_t cursor = 0;
    for (size_t k = 0; k < corners.size(); ++k) {
        const double center = stations[corners[k]];
        const double left = std::max(center-radius, k == 0 ? 0.0 :
            (stations[corners[k-1]]+center)*0.5);
        const double right = std::min(center+radius, k+1 == corners.size() ? stations.back() :
            (center+stations[corners[k+1]])*0.5);
        Vec2d incoming, outgoing;
        const Vec3d before = sample(left, true, incoming);
        const Vec3d after = sample(right, false, outgoing);
        const double cosine = std::max(-1.0, std::min(1.0, incoming.dot(outgoing)));
        if (center-left < 1e-7 || right-center < 1e-7 || incoming.norm() < 1e-7 ||
            outgoing.norm() < 1e-7 || cosine < -1.0+1e-10) continue;
        const double half_cosine = std::sqrt(0.5*(1.0+cosine));
        const double ratio = (4.0/3.0)*half_cosine/(1.0+half_cosine);
        const Vec3d first(before.xy()+incoming*(ratio*(center-left)),
            before.z()+(after.z()-before.z())/3.0);
        const Vec3d last(after.xy()-outgoing*(ratio*(right-center)),
            after.z()-(after.z()-before.z())/3.0);
        const CornerBezier curve{{{before, first, last, after}}};
        // 扫描尖折与写出前缀使用不同游标：不能把未进入修复窗口的骨干点漏写成一条直连线。
        while (cursor < points.size() && stations[cursor] < left-1e-9)
            appendCornerPoint(result.points, points[cursor++]);
        for (const auto& point : sampleCornerBezier(curve)) appendCornerPoint(result.points, point);
        while (cursor < points.size() && stations[cursor] <= right+1e-9) ++cursor;
    }
    for (; cursor < points.size(); ++cursor) appendCornerPoint(result.points, points[cursor]);
    if (!result.points.empty()) {
        result.points.front() = source.points.front();
        result.points.back() = source.points.back();
    }
    return result;
}

/// 给定两端坐标与方向建立切向过渡；可连接固定车道端点和偏移后的骨干挂接点。
inline CornerBezier fitAttachmentTransition(const Vec3d& a, const Vec2d& start_direction,
                                            const Vec3d& b, const Vec2d& end_direction, double radius,
                                            const Vec2d* monotone_support = nullptr) {
    const Vec2d start = start_direction.norm() > 1e-9 ? start_direction.normalized() : Vec2d(1, 0);
    const Vec2d end = end_direction.norm() > 1e-9 ? end_direction.normalized() : start;
    double handle = std::min((b.xy() - a.xy()).norm() / 3.0, radius);
    // 端部过渡的第二控制点必须留在起始支撑线两端之间。否则控制多边形会先
    // 向骨干内侧移动，再返回终点，形成用户看到的“向骨干鼓肚子”。
    // 控制点仍位于终点切线，故只缩短把手，不改变端点G1方向。
    const Vec2d support = monotone_support && monotone_support->norm() > 1e-9
        ? monotone_support->normalized() : start;
    const Vec2d normal = rotLeft(support);
    const double y0 = a.xy().dot(normal), y1 = b.xy().dot(normal);
    const double lateral = end.dot(normal);
    if (std::abs(lateral) > 1e-9) {
        const double lower = std::min(y0, y1), upper = std::max(y0, y1);
        const double max_handle = lateral > 0 ? (y1 - lower) / lateral : (upper - y1) / (-lateral);
        handle = std::min(handle, std::max(0.0, max_handle));
    }
    // 极端相向切向下仍保留非零把手，确保端点解析切向存在；其位移小于微米级
    // 容差，不会重新形成可见的法向鼓肚。
    if (handle < 1e-6 && (b.xy()-a.xy()).norm() > 1e-6) handle = 1e-6;
    CornerBezier curve;
    const Vec3d first(a.xy() + start * handle, a.z() + (b.z()-a.z()) / 3.0);
    const Vec3d last(b.xy() - end * handle, b.z() - (b.z()-a.z()) / 3.0);
    curve.ctrl = {{a, first, last, b}};
    return curve;
}

/// 以公共点前后各corner_radius为目标截距替换折角，并按0.2米平面弧长采样贝塞尔。
/// 短段截距不超过本段长度；同段两端圆角相遇时同比缩短截距，避免反向或重叠。
/// 首尾保留调用者传入的坐标（可以是偏移后的骨干交点）；可选方向用于首尾切向过渡。
/// 首尾切向过渡仅占用不超过radius的邻段长度，中间直线保持原轨位，不将整条长直段拉成弓形。
/// 进出挂接点使用受限切向过渡，控制多边形在支撑线法向上不回退；生产布局会在
/// 真实折点前后屏蔽corner_radius保护带内的占用事件短段，但不会删除有效骨干段。
inline LineString2d roundPolylineCorners(const std::vector<Vec3d>& source, double radius,
                                         const Vec2d* start_direction = nullptr,
                                         const Vec2d* end_direction = nullptr,
                                         const Vec2d* start_support_direction = nullptr,
                                         const Vec2d* end_support_direction = nullptr) {
    LineString2d result;
    if (radius <= 0.0 || source.size() <= 1) { result.points = source; return result; }
    std::vector<Vec3d> points;
    for (const auto& p : source) {
        appendCornerPoint(points, p);
        while (points.size() >= 3) {
            const size_t n = points.size();
            const Vec3d a = points[n-3], b = points[n-2], c = points[n-1];
            const Vec2d ab = b.xy()-a.xy(), bc = c.xy()-b.xy();
            const double first = ab.norm(), last = bc.norm();
            if (first < 1e-9 || last < 1e-9 || ab.dot(bc) <= 0 ||
                std::abs(ab.x()*bc.y()-ab.y()*bc.x()) > 1e-10*first*last ||
                std::abs(b.z() - (a.z() + (c.z()-a.z()) * first/(first+last))) > 1e-6) break;
            points[n-2] = c; points.pop_back();
        }
    }
    const size_t n = points.size();
    if (n < 2) { result.points = points; return result; }
    std::vector<double> lengths(n-1), incoming(n, 0.0), outgoing(n, 0.0);
    for (size_t i = 0; i+1 < n; ++i) lengths[i] = (points[i+1].xy()-points[i].xy()).norm();
    for (size_t i = 1; i+1 < n; ++i) {
        if (lengths[i-1] < 1e-9 || lengths[i] < 1e-9) continue;
        const double dot = (points[i].xy()-points[i-1].xy()).dot(points[i+1].xy()-points[i].xy()) / (lengths[i-1]*lengths[i]);
        // 直线无需圆角；精确180度折返无法在零横向空间内拟合，保留交给最终自交审计。
        if (dot > 1.0-1e-10 || dot < -1.0+1e-10) continue;
        incoming[i] = std::min(radius, lengths[i-1]);
        outgoing[i] = std::min(radius, lengths[i]);
    }
    const auto changes_direction = [](const Vec2d* direction, const Vec2d& segment) {
        return direction && direction->norm() > 1e-9 && segment.norm() > 1e-9 &&
               direction->normalized().dot(segment.normalized()) < 1.0-1e-10;
    };
    const bool start_transition = changes_direction(start_direction, points[1].xy()-points[0].xy());
    const bool end_transition = changes_direction(end_direction, points[n-1].xy()-points[n-2].xy());
    if (n > 2 && start_transition) incoming[1] = std::min(incoming[1], lengths[0]*0.5);
    if (n > 2 && end_transition) outgoing[n-2] = std::min(outgoing[n-2], lengths[n-2]*0.5);
    for (size_t i = 0; i+1 < n; ++i) {
        const double sum = outgoing[i] + incoming[i+1];
        if (sum > lengths[i] && sum > 0) {
            const double scale = lengths[i]/sum;
            outgoing[i] *= scale; incoming[i+1] *= scale;
        }
    }
    const auto connect = [&](const Vec3d& a, const Vec3d& b, const Vec2d& first,
                             const Vec2d& last, bool first_transition, bool last_transition) {
        const double length = (b.xy()-a.xy()).norm();
        if (length < 1e-9) {
            appendCornerPoint(result.points, b); return;
        }
        if (!first_transition && !last_transition) {
            appendCornerPoint(result.points, b); return;
        }
        double start_length = first_transition ? std::min(radius, length) : 0;
        double end_length = last_transition ? std::min(radius, length) : 0;
        if (start_length+end_length > length) {
            const double scale = length/(start_length+end_length);
            start_length *= scale; end_length *= scale;
        }
        const Vec2d direction = (b.xy()-a.xy())/length;
        const Vec3d start = a+(b-a)*(start_length/length);
        const Vec3d end = b+(a-b)*(end_length/length);
        if (start_length > 1e-9)
            for (const auto& p : sampleCornerBezier(fitAttachmentTransition(a, first, start, direction, radius)))
                appendCornerPoint(result.points, p);
        // 端部过渡之间按0.2米补齐直线点；保持轨位和高程，不引入整段横向拱起。
        const double straight = (end.xy()-start.xy()).norm();
        for (double station = 0.2; station < straight; station += 0.2)
            appendCornerPoint(result.points, start+(end-start)*(station/straight));
        appendCornerPoint(result.points, end);
        if (end_length > 1e-9)
            for (const auto& p : sampleCornerBezier(fitAttachmentTransition(end, direction, b, last, radius)))
                appendCornerPoint(result.points, p);
        else appendCornerPoint(result.points, b);
    };
    appendCornerPoint(result.points, points.front());
    for (size_t i = 1; i+1 < n; ++i) {
        const Vec3d before = lengths[i-1] > 1e-9 ? points[i] + (points[i-1]-points[i]) * (incoming[i]/lengths[i-1]) : points[i];
        const Vec3d after = lengths[i] > 1e-9 ? points[i] + (points[i+1]-points[i]) * (outgoing[i]/lengths[i]) : points[i];
        const Vec2d first = start_direction ? *start_direction : points[1].xy()-points[0].xy();
        connect(result.points.back(), before, first, points[i].xy()-points[i-1].xy(), i == 1 && start_transition, false);
        if (incoming[i] > 1e-9 && outgoing[i] > 1e-9) {
            // 进出挂接点已经由车道射线与偏移骨干支撑线求得；若继续把该点作为
            // 普通折点的控制顶点，偏移量变化会把曲线拉回骨干后再离开。端部只
            // 在截取点两侧之间做一次切向过渡，内部真实骨干折点仍使用标准圆角。
            const bool endpoint_corner = (i == 1 && start_direction && n > 2) ||
                (i + 2 == n && end_direction && n > 2);
            const auto corner = endpoint_corner
                ? fitAttachmentTransition(before,
                    (i == 1 && start_direction) ? *start_direction : points[i].xy()-points[i-1].xy(),
                    after,
                    (i + 2 == n && end_direction) ? *end_direction : points[i+1].xy()-points[i].xy(),
                    radius,
                    (i == 1 ? start_support_direction : end_support_direction))
                : fitCornerBezier(before, points[i], after);
            for (const auto& p : sampleCornerBezier(corner))
                appendCornerPoint(result.points, p);
        } else appendCornerPoint(result.points, points[i]);
    }
    const Vec2d last = end_direction ? *end_direction : points[n-1].xy()-points[n-2].xy();
    const Vec2d first = n == 2 && start_direction ? *start_direction : points[n-1].xy()-points[n-2].xy();
    connect(result.points.back(), points.back(), first, last, n == 2 && start_transition, end_transition);
    // 弧长采样中的最后一次线性插值可能产生数值舍入，公共接口仍须返回精确端点。
    if (!result.points.empty()) {
        result.points.front() = points.front();
        result.points.back() = points.back();
    }
    return result;
}

/// 以车道延长线和偏移骨干端段的交点作为公共转角，一次性圆角整个挂接路径。
/// 骨干端点不再强制车道切向；仅车道端点固定。近乎平行、反向或交点越过邻点时保留原挂接点，
/// 由车道端点的局部切向过渡处理。求交移动限制在4倍截距内，避免近乎平行的远处交点。
inline LineString2d roundAttachedPolyline(const std::vector<Vec3d>& backbone,
    const Vec3d& entry, const Vec2d& entry_direction, const Vec3d& exit,
    const Vec2d& exit_direction, double radius) {
    std::vector<Vec3d> knots = backbone;
    if (radius > 0 && knots.size() >= 2) {
        const auto intersection = [&](bool at_entry) {
            const size_t index = at_entry ? 0 : knots.size()-1;
            const size_t neighbor = at_entry ? 1 : knots.size()-2;
            const Vec3d a = backbone[index], b = backbone[neighbor];
            const Vec2d tangent = b.xy()-a.xy();
            const Vec2d lane_direction = at_entry ? entry_direction : -exit_direction;
            const Vec3d origin = at_entry ? entry : exit;
            const double determinant = lane_direction.norm() > 1e-9 && tangent.norm() > 1e-9
                ? (lane_direction.x()*tangent.y()-lane_direction.y()*tangent.x()) : 0;
            if (std::abs(determinant) < 1e-8*lane_direction.norm()*tangent.norm() || std::abs(determinant) < 1e-12) return;
            const Vec2d delta = a.xy()-origin.xy();
            const double fraction = (delta.x()*lane_direction.y()-delta.y()*lane_direction.x())/determinant;
            const Vec3d candidate = a+(b-a)*fraction;
            if (!std::isfinite(fraction) || fraction >= 1-1e-7 ||
                (candidate.xy()-a.xy()).norm() > 4*radius ||
                (candidate.xy()-origin.xy()).dot(lane_direction) < -1e-7) return;
            knots[index] = candidate;
        };
        intersection(true);
        intersection(false);
        // 双端求交不能使短骨干段反向；退回原点列，由圆角统一限制两端截距。
        if (knots.size() == 2 && (knots[1].xy()-knots[0].xy()).dot(backbone[1].xy()-backbone[0].xy()) <= 0)
            knots = backbone;

    }
    std::vector<Vec3d> full;
    appendCornerPoint(full, entry);
    for (const auto& point : knots) appendCornerPoint(full, point);
    appendCornerPoint(full, exit);
    const Vec2d entry_support = knots.size() >= 2 ? knots[1].xy()-knots[0].xy() : entry_direction;
    const Vec2d exit_support = knots.size() >= 2 ? knots[knots.size()-1].xy()-knots[knots.size()-2].xy() : exit_direction;
    auto result = roundPolylineCorners(full, radius, &entry_direction, &exit_direction,
        &entry_support, &exit_support);
    // 末个采样点可能与精确端点仅差浮点舍入量，仍应精确保留调用者给定的车道端点XYZ。
    if (result.points.size() >= 2) { result.points.front() = entry; result.points.back() = exit; }
    return result;
}

}  // namespace nonmotor_detail
}  // namespace isg
