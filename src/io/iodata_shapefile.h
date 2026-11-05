// 路口输入输出的 Shapefile 调查数据导出接口。

#pragma once

#include <string>
#include "types.h"
#include "nonmotor/nonmotor_type.h"

namespace isg {

bool save(const IntersectionInput& input, std::string out_dir, std::string prefix);
bool save(const IntersectionOutput& output, std::string out_dir, std::string prefix);
bool save(const NonMotorBackbone& nonmotor_backbone, std::string dir, std::string prefix);
}
