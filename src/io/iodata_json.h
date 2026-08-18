#pragma once

#include "types.h"
#include "utils/json.hpp"

#include <string>

namespace isg {

// 三维点 JSON 转换保留为兼容公共 API，二维旧数据缺省高程为零。
Vec3d vec3dFromJson(const nlohmann::json& json);
nlohmann::json vec3dToJson(const Vec3d& point);
nlohmann::json lineString2dToJson(const LineString2d& line);

// 路口输入的文件和字符串 JSON 编解码门面。
class IntersectionIO {
public:
    static void saveToFile(const IntersectionInput& input, const std::string& filepath);
    static IntersectionInput loadFromFile(const std::string& filepath);
    static std::string toJsonString(const IntersectionInput& input);
    static IntersectionInput fromJsonString(const std::string& json_string);
};

}  // 命名空间 isg
