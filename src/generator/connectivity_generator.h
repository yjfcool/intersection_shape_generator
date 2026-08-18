#pragma once
#include "types.h"

namespace isg {

/// 无状态门面；单路口可变状态保存在独立生成会话中。
class ConnectivityGenerator {
public:
    explicit ConnectivityGenerator(
        const LBFGSConfig& cfg = {},
        const ConnectivityDirectionConfig& direction_cfg = {});

    /// 生成全部连通曲线
    /// @param out_ms 输出优化耗时(毫秒)
    std::vector<ConnectivityCurve> generate(const IntersectionInput&, SDFField&, double* out_ms = nullptr);

private:
    LBFGSConfig config_;
    ConnectivityDirectionConfig direction_config_;
};

}
