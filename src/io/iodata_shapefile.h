// 路口输入输出的 Shapefile 调查数据导出接口。

#pragma once

#include "types.h"
#include <string>

namespace isg {

bool save(IntersectionInput& input, std::string out_dir, std::string prefix);
bool save(IntersectionOutput& output, std::string out_dir, std::string prefix);

}
