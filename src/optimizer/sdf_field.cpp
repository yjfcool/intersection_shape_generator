#include "sdf_field.h"
#include "constraints/fence_check.h"
#include "utils/clipper.hpp"
#include <cmath>
#include <algorithm>
#include <limits>
#include <functional>

namespace isg {

static double ptSegDist(const Vec2d& p, const Vec2d& a, const Vec2d& b) {
    Vec2d ab = b - a, ap = p - a;
    double t = ab.dot(ap);
    double l2 = ab.squaredNorm();
    if (l2 < 1e-20) return ap.norm();
    t = std::max(0.0, std::min(1.0, t / l2));
    return (p - (a + t * ab)).norm();
}

double SDFField::distToPolygons(const Vec2d& pt, const std::vector<Polygon2d>& polys) const {
    double m = std::numeric_limits<double>::max();
    for (auto& poly : polys) {
        auto& ring = poly.outer;
        int n = (int)ring.size();
        for (int i = 0; i < n; ++i)
            m = std::min(m, ptSegDist(pt, ring[i], ring[(i + 1) % n]));
    }
    return m;
}

bool SDFField::insideAny(const Vec2d& pt, const std::vector<Polygon2d>& polys) const {
    for (auto& p : polys)
        if (polygonContains(p, pt))
            return true;
    return false;
}

void SDFField::build(const BoundingBox2d& roi, const std::vector<Obstacle>& obs, double cs, double buf) {
    CacheKey key{roi, obs, cs, buf};

    // 尝试从缓存获取
    auto it = cache_map_.find(key);
    if (it != cache_map_.end()) {
        grid_ = it->second;
        cs_ = cs;
        roi_ = roi;
        buffer_radius_ = buf;
        // 重建原始障碍物列表；缓冲在网格值中以signed_distance - buffer_radius_表达。
        std::vector<Polygon2d> buffered;
        for (auto& o : obs) {
            if (o.geometry.outer.empty()) continue;
            buffered.push_back(o.geometry);
        }
        if (buffered.size() > 1) {
            std::vector<std::vector<std::array<double,2>>> subs;
            for (auto& p : buffered) if (!p.outer.empty()) subs.push_back(toArray(p.outer));
            auto sol = ClipperUtil::UnionPaths(subs, ClipperLib::pftNonZero);
            if (!sol.empty()) {
                buffered_.clear();
                for (auto& path : sol) {
                    Polygon2d pg;
                    pg.outer = toArray(path);
                    buffered_.push_back(pg);
                }
            } else {
                buffered_ = buffered;
            }
        } else {
            buffered_ = buffered;
        }
        // 根据ROI和网格尺寸计算行数列数
        double mg = cs_;
        BoundingBox2d ext;
        ext.min_pt = roi_.min_pt - Vec2d(mg, mg);
        ext.max_pt = roi_.max_pt + Vec2d(mg, mg);
        roi_ = ext;
        cols_ = std::max(2, (int)std::ceil(roi_.width() / cs_));
        rows_ = std::max(2, (int)std::ceil(roi_.height() / cs_));
        return;
    }

    // 缓存未命中,正常构建并存入缓存
    buildInternal(roi, obs, cs, buf);

    // 存入缓存
    cache_map_[key] = grid_;
}

void SDFField::buildFromPolygons(const BoundingBox2d& roi, const std::vector<Polygon2d>& polys, double cs) {
    cs_ = cs;
    buffer_radius_ = 0.0;
    roi_ = roi;
    buffered_ = polys;
    rebuildGrid(polys);
}

void SDFField::rebuildGrid(const std::vector<Polygon2d>& ops) {
    double mg = cs_;
    BoundingBox2d ext;
    ext.min_pt = roi_.min_pt - Vec2d(mg, mg);
    ext.max_pt = roi_.max_pt + Vec2d(mg, mg);
    roi_ = ext;
    cols_ = std::max(2, (int)std::ceil(roi_.width() / cs_));
    rows_ = std::max(2, (int)std::ceil(roi_.height() / cs_));
    grid_.assign(rows_ * cols_, 0.0);
    for (int r = 0; r < rows_; ++r) {
        for (int c = 0; c < cols_; ++c) {
            Vec2d pt = cellToWorld(r, c);
            double d = distToPolygons(pt, ops);
            bool ins = insideAny(pt, ops);
            double signed_d = ins ? -d : d;
            grid_[idx(r, c)] = signed_d - buffer_radius_;
        }
    }
}

void SDFField::updateRegion(const BoundingBox2d& dirty, const std::vector<Obstacle>& obs, double buf) {
    std::vector<Polygon2d> lp;
    for (auto& o : obs)
        lp.push_back(o.geometry);
    std::pair<int, int> _mc = worldToCell(dirty.min_pt);
    int rm = _mc.first;
    int cm = _mc.second;
    std::pair<int, int> _xc = worldToCell(dirty.max_pt);
    int rx = _xc.first;
    int cx = _xc.second;
    rm = std::max(0, rm - 1);
    cm = std::max(0, cm - 1);
    rx = std::min(rows_ - 1, rx + 1);
    cx = std::min(cols_ - 1, cx + 1);
    for (int r = rm; r <= rx; ++r)
        for (int c = cm; c <= cx; ++c) {
            Vec2d pt = cellToWorld(r, c);
            double d = distToPolygons(pt, lp);
            bool ins = insideAny(pt, lp);
            double signed_d = ins ? -d : d;
            grid_[idx(r, c)] = signed_d - buf;
        }
}

double SDFField::rawAt(int r, int c) const {
    r = std::max(0, std::min(rows_ - 1, r));
    c = std::max(0, std::min(cols_ - 1, c));
    return grid_[idx(r, c)];
}

std::pair<int, int> SDFField::worldToCell(const Vec2d& p) const {
    return {(int)((p[1] - roi_.min_pt.y()) / cs_), (int)((p[0] - roi_.min_pt.x()) / cs_)};
}

Vec2d SDFField::cellToWorld(int r, int c) const {
    return Vec2d(roi_.min_pt[0] + (c + 0.5) * cs_, roi_.min_pt[1] + (r + 0.5) * cs_);
}

std::pair<double, Vec2d> SDFField::queryWithGrad(const Vec2d& pt) const {
    if (grid_.empty()) return {1e18, Vec2d(0, 0)};
    double fx = (pt[0] - roi_.min_pt.x()) / cs_ - 0.5, fy = (pt[1] - roi_.min_pt.y()) / cs_ - 0.5;
    int x0 = (int)std::floor(fx), y0 = (int)std::floor(fy), x1 = x0 + 1, y1 = y0 + 1;
    double tx = fx - x0, ty = fy - y0;
    x0 = std::max(0, std::min(cols_ - 1, x0));
    x1 = std::max(0, std::min(cols_ - 1, x1));
    y0 = std::max(0, std::min(rows_ - 1, y0));
    y1 = std::max(0, std::min(rows_ - 1, y1));
    double v00 = rawAt(y0, x0), v10 = rawAt(y1, x0), v01 = rawAt(y0, x1), v11 = rawAt(y1, x1);
    double val = (1 - tx) * (1 - ty) * v00 + tx * (1 - ty) * v01 + (1 - tx) * ty * v10 + tx * ty * v11;
    double gx = ((1 - ty) * (v01 - v00) + ty * (v11 - v10)) / cs_;
    double gy = ((1 - tx) * (v10 - v00) + tx * (v11 - v01)) / cs_;
    return {val, Vec2d(gx, gy)};
}

bool SDFField::isSafe(const Vec2d& pt, double cl) const {
    std::pair<double,Vec2d> _q = queryWithGrad(pt);
    return _q.first >= cl;
}

double SDFField::obstaclePenalty(const Vec2d& pt, double cl) const {
    std::pair<double,Vec2d> _q2 = queryWithGrad(pt);
    double s = _q2.first - cl;
    return s >= 0 ? 0.0 : s * s;
}

Vec2d SDFField::obstaclePenaltyGrad(const Vec2d& pt, double cl) const {
    std::pair<double, Vec2d> _q3 = queryWithGrad(pt);
    double d = _q3.first;
    Vec2d gd = _q3.second;
    double s = d - cl;
    return s >= 0 ? Vec2d(0, 0) : 2.0 * s * gd;
}

// ─── 静态缓存实现 ──────────────────────────────────────────────────────

// 定义静态成员
std::unordered_map<SDFField::CacheKey, std::vector<double>, SDFField::CacheKeyHash> SDFField::cache_map_;

bool SDFField::CacheKey::operator==(const CacheKey& other) const {
    // 比较ROI
    if (roi.min_pt.x() != other.roi.min_pt.x() || roi.min_pt.y() != other.roi.min_pt.y() ||
        roi.max_pt.x() != other.roi.max_pt.x() || roi.max_pt.y() != other.roi.max_pt.y()) {
        return false;
    }

    // 比较网格尺寸与缓冲
    if (cell_size != other.cell_size || buffer != other.buffer) {
        return false;
    }

    // 比较障碍物
    if (obstacles.size() != other.obstacles.size()) {
        return false;
    }

    // 逐顶点比较。早期实现只比对障碍物个数与首个障碍物的首个顶点，
    // 两组"个数相同、首顶点相同"的不同障碍物会命中同一份栅格，
    // 从而拿到错误的距离场并直接破坏障碍物避让判定。
    for (size_t i = 0; i < obstacles.size(); ++i) {
        const auto& ring1 = obstacles[i].geometry.outer;
        const auto& ring2 = other.obstacles[i].geometry.outer;
        if (ring1.size() != ring2.size()) {
            return false;
        }
        for (size_t v = 0; v < ring1.size(); ++v) {
            if (std::abs(ring1[v].x() - ring2[v].x()) > 1e-9 ||
                std::abs(ring1[v].y() - ring2[v].y()) > 1e-9) {
                return false;
            }
        }
    }

    return true;
}

std::size_t SDFField::CacheKeyHash::operator()(const SDFField::CacheKey& k) const {
    std::size_t h1 = std::hash<double>{}(k.roi.min_pt.x());
    std::size_t h2 = std::hash<double>{}(k.roi.min_pt.y());
    std::size_t h3 = std::hash<double>{}(k.roi.max_pt.x());
    std::size_t h4 = std::hash<double>{}(k.roi.max_pt.y());
    std::size_t h5 = std::hash<double>{}(k.cell_size);
    std::size_t h6 = std::hash<double>{}(k.buffer);
    std::size_t h7 = std::hash<size_t>{}(k.obstacles.size());

    // 全顶点折叠，与 operator== 的比较范围保持一致。
    std::size_t h8 = 0;
    for (const auto& obs : k.obstacles) {
        h8 = h8 * 1000003u + std::hash<size_t>{}(obs.geometry.outer.size());
        for (const auto& v : obs.geometry.outer) {
            h8 = h8 * 1000003u + std::hash<double>{}(v.x());
            h8 = h8 * 1000003u + std::hash<double>{}(v.y());
        }
    }

    // 组合哈希值
    return h1 ^ (h2 << 1) ^ (h3 << 2) ^ (h4 << 3) ^ (h5 << 4) ^ (h6 << 5) ^ (h7 << 6) ^ (h8 << 7);
}

void SDFField::buildInternal(const BoundingBox2d& roi, const std::vector<Obstacle>& obs, double cs, double buf) {
    cs_ = cs;
    buffer_radius_ = buf;
    roi_ = roi;
    std::vector<Polygon2d> buffered;
    for (auto& o : obs) {
        if (o.geometry.outer.empty()) continue;
        buffered.push_back(o.geometry);
    }
    if (buffered.size() > 1) {
        std::vector<std::vector<std::array<double,2>>> subs;
        for (auto& p : buffered) if (!p.outer.empty()) subs.push_back(toArray(p.outer));
        auto sol = ClipperUtil::UnionPaths(subs, ClipperLib::pftNonZero);
        if (!sol.empty()) {
            buffered_.clear();
            for (auto& path : sol) {
                Polygon2d pg;
                pg.outer = toArray(path);
                buffered_.push_back(pg);
            }
        } else {
            buffered_ = buffered;
        }
    } else {
        buffered_ = buffered;
    }
    rebuildGrid(buffered_);
}

// ─── 文件末尾补充 ────────────────────────────────────────────────────────────

}
