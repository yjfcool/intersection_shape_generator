#pragma once

#include "types.h"

namespace isg {

/// U-turn 候选相对参考曲线的长度、空间半径、拱高和曲率包络检测。
class UTurnEnvelopeConstraint {
public:
    bool exceeds(const BezierCurve& curve,
                 const BezierCurve& reference,
                 const Vec2d& entry,
                 const Vec2d& exit) const;

    bool collapses(const BezierCurve& curve,
                   const BezierCurve& reference,
                   const Vec2d& entry,
                   const Vec2d& exit) const;

private:
    double lateralBulge(const BezierCurve& curve,
                        const Vec2d& entry,
                        const Vec2d& exit) const;
};

}  // 命名空间 isg
