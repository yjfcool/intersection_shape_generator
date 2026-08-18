#pragma once
#include "types.h"
#include <unordered_map>

namespace isg {

/// 符号距离场(SDF)：用于快速查询任意点到障碍物最近表面的有符号距离
/// 网格化存储，支持带梯度的查询与静态缓存复用
class SDFField {
public:
    /// 由障碍物多边形构建SDF
    /// cs=网格单元尺寸(米), buf=障碍物缓冲距离(米)
    void build(const BoundingBox2d&, const std::vector<Obstacle>&, double cs = 0.2, double buf = 0.0);

    /// 由通用多边形构建SDF
    void buildFromPolygons(const BoundingBox2d&, const std::vector<Polygon2d>&, double cs = 0.2);

    /// 查询: 返回(SDF值, 梯度向量)
    std::pair<double, Vec2d> queryWithGrad(const Vec2d&) const;

    /// 点是否安全（SDF >= cl）
    bool isSafe(const Vec2d&, double cl = 0.1) const;

    /// 障碍物惩罚值：穿透深度平方
    double obstaclePenalty(const Vec2d&, double cl = 0.1) const;

    /// 障碍物惩罚梯度
    Vec2d obstaclePenaltyGrad(const Vec2d&, double cl = 0.1) const;

    /// 局部更新SDF区域
    void updateRegion(const BoundingBox2d&, const std::vector<Obstacle>&, double buf = 0.0);

    double cellSize() const { return cs_; }
    int rows() const { return rows_; }
    int cols() const { return cols_; }
    std::pair<int, int> worldToCell(const Vec2d&) const;
    Vec2d cellToWorld(int r, int c) const;
    bool valid() const { return !grid_.empty(); }

private:
    void rebuildGrid(const std::vector<Polygon2d>&);
    double distToPolygons(const Vec2d&, const std::vector<Polygon2d>&) const;
    bool insideAny(const Vec2d&, const std::vector<Polygon2d>&) const;
    inline int idx(int r, int c) const { return r * cols_ + c; }
    double rawAt(int r, int c) const;

    /// 缓存相关方法与数据
    struct CacheKey {
        BoundingBox2d roi;
        std::vector<Obstacle> obstacles;
        double cell_size;
        double buffer;
        bool operator==(const CacheKey& other) const;
    };
    struct CacheKeyHash {
        std::size_t operator()(const SDFField::CacheKey& k) const;
    };
    void buildInternal(const BoundingBox2d& roi, const std::vector<Obstacle>& obs, double cs, double buf);

private:
    double cs_ = 0.2;
    double buffer_radius_ = 0.0;
    BoundingBox2d roi_;
    int rows_ = 0, cols_ = 0;
    std::vector<double> grid_;
    std::vector<Polygon2d> buffered_;
    static std::unordered_map<CacheKey, std::vector<double>, CacheKeyHash> cache_map_;  ///< 静态缓存存储
};

}
