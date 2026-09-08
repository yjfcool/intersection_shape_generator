#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <cmath>
#include <algorithm>
#include <cstdint>
#include <cstring>
#include <iomanip>

// =============================================================================
// 数据结构定义
// =============================================================================

typedef std::array<double, 3> ShapePoint;

// Shapefile 几何类型
static constexpr int SHP_NULL = 0;
static constexpr int SHP_POINT = 1;
static constexpr int SHP_POLYLINE = 3;
static constexpr int SHP_POLYGON = 5;
static constexpr int SHP_POINTZ = 11;
static constexpr int SHP_POLYLINEZ = 13;
static constexpr int SHP_POLYGONZ = 15;

// DBF 字段描述符
struct IsgDbfField {
    std::string name; // 字段名（最大 11 个字节）
    char type = 'C'; // 'C'=字符串(Character), 'N'=数值(Numeric), 'F'=浮点数(Float)
    uint8_t length = 254; // 字段总长度
    uint8_t decimals = 0; // 小数位数（仅在 N 或 F 类型时有效）

    IsgDbfField() {};

    IsgDbfField(std::string _name, char _type, uint8_t _length, uint8_t _decimals)
            : name(_name), type(_type), length(_length), decimals(_decimals) {}
};

struct ShapeRecord {
    int32_t id = 0; // 1-based 记录编号
    int32_t shapeType = 0; // 0:None, 1:Point, 3:PolyLine, 5:Polygon, 11:PointZ, 13:PolyLineZ, 15:PolygonZ
    std::vector<std::string> attributes; // 属性值列表: 与DbfField顺序一一对应的属性值（全存为字符串，方便统一管理）
    std::vector<ShapePoint> points; // 顶点坐标集合
    std::vector<int32_t> parts = {0}; // 每个 Part 的起始点索引
    double box[4] = {0, 0, 0, 0}; // 局部二维边界框 [Xmin, Ymin, Xmax, Ymax]

    ShapeRecord() {};

    ShapeRecord(int32_t _id, int32_t _shapeType, const std::vector<std::string> &_attributes,
                const std::vector<ShapePoint> &_points, std::vector<int32_t> _parts = {0})
            : id(_id), shapeType(_shapeType), attributes(_attributes), points(_points), parts(_parts) {}
};

class ShapefileEngine {
private:
    // =============================================================================
    // 基础字节序与工具辅助函数
    // =============================================================================
    static inline int32_t readBE32(std::ifstream &f) {
        uint8_t b[4];
        f.read(reinterpret_cast<char *>(b), 4);
        return (b[0] << 24) | (b[1] << 16) | (b[2] << 8) | b[3];
    }

    static inline int32_t readLE32(std::ifstream &f) {
        uint8_t b[4];
        f.read(reinterpret_cast<char *>(b), 4);
        return b[0] | (b[1] << 8) | (b[2] << 16) | (b[3] << 24);
    }

    static inline int16_t readLE16(std::ifstream &f) {
        uint8_t b[2];
        f.read(reinterpret_cast<char *>(b), 2);
        return b[0] | (b[1] << 8);
    }

    static inline double readLEDouble(std::ifstream &f) {
        double val;
        f.read(reinterpret_cast<char *>(&val), 8);
        return val;
    }

    static inline void writeBE32(std::ofstream &f, int32_t val) {
        uint8_t b[4] = {
                static_cast<uint8_t>((val >> 24) & 0xFF),
                static_cast<uint8_t>((val >> 16) & 0xFF),
                static_cast<uint8_t>((val >> 8) & 0xFF),
                static_cast<uint8_t>(val & 0xFF)
        };
        f.write(reinterpret_cast<char *>(b), 4);
    }

    static inline void writeLE32(std::ofstream &f, int32_t val) {
        uint8_t b[4] = {
                static_cast<uint8_t>(val & 0xFF),
                static_cast<uint8_t>((val >> 8) & 0xFF),
                static_cast<uint8_t>((val >> 16) & 0xFF),
                static_cast<uint8_t>((val >> 24) & 0xFF)
        };
        f.write(reinterpret_cast<char *>(b), 4);
    }

    static inline void writeLE16(std::ofstream &f, int16_t val) {
        uint8_t b[2] = {static_cast<uint8_t>(val & 0xFF), static_cast<uint8_t>((val >> 8) & 0xFF)};
        f.write(reinterpret_cast<char *>(b), 2);
    }

    static inline void writeLEDouble(std::ofstream &f, double val) {
        f.write(reinterpret_cast<char *>(&val), 8);
    }

    static inline std::string trim(const std::string &str) {
        size_t first = str.find_first_not_of(" \t\r\n");
        if (first == std::string::npos)
            return "";
        size_t last = str.find_last_not_of(" \t\r\n");
        return str.substr(first, (last - first + 1));
    }

    static inline std::string changeExtension(const std::string &path, const std::string &newExt) {
        size_t dot = path.find_last_of('.');
        if (dot == std::string::npos)
            return path + newExt;
        return path.substr(0, dot) + newExt;
    }

    static inline bool fileExists(const std::string &path) {
        std::ifstream f(path);
        return f.good();
    }

    static void calculateRecordBox(ShapeRecord &rec) {
        if (rec.points.empty())
            return;
        rec.box[0] = rec.box[2] = rec.points[0][0];
        rec.box[1] = rec.box[3] = rec.points[0][1];
        for (const auto &p : rec.points) {
            if (p[0] < rec.box[0])
                rec.box[0] = p[0];
            if (p[1] < rec.box[1])
                rec.box[1] = p[1];
            if (p[0] > rec.box[2])
                rec.box[2] = p[0];
            if (p[1] > rec.box[3])
                rec.box[3] = p[1];
        }
    }

    static bool writePrj(const std::string &path) {
        std::ofstream f(path);
        if (!f) {
            std::cerr << "[ERROR] Cannot open prj: " << path << "\n";
            return false;
        }
        // 默认使用WGS84地理坐标系（投影坐标系时此处需定制）
        f << "GEOGCS[\"GCS_WGS_1984\","
          << "DATUM[\"D_WGS_1984\","
          << "SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],"
          << "PRIMEM[\"Greenwich\",0.0],"
          << "UNIT[\"Degree\",0.017453292519943295]]";
        return true;
    }

public:
    // =============================================================================
    // Reader（自适应同步读取 .shp 与 .dbf）
    // =============================================================================
    static bool read(const std::string &dir, const std::string &fname,
                     std::vector<IsgDbfField> &outFields, std::vector<ShapeRecord> &outRecords,
                     int32_t &outGlobalShapeType, double outGlobalBoxXYZ[6] = 0) {
        std::string shpPath = dir + "/" + fname + ".shp";
        std::string shxPath = changeExtension(shpPath, ".shx");
        std::string dbfPath = changeExtension(shpPath, ".dbf");

        if (fileExists(shpPath)) {
            // 1. 读取几何数据 (.shp)
            std::ifstream fShp(shpPath, std::ios::binary);
            if (!fShp.is_open())
                return false;

            // 解析 100 字节的 Main File Header
            int32_t fileCode = readBE32(fShp);
            if (fileCode != 9994)
                return false; // Shapefile 固定的 Magic Number

            //f.seekg(24, std::ios::beg); // 跳过未使用的 20 字节
            //int32_t fileLengthInWords = readBE32(f); // 16位字(Word)单位的文件长度
            fShp.seekg(32, std::ios::beg);

            int32_t version = readLE32(fShp);
            if (version != 1000)
                return false;

            outGlobalShapeType = readLE32(fShp);
            outGlobalBoxXYZ[0] = readLEDouble(fShp); // Xmin
            outGlobalBoxXYZ[1] = readLEDouble(fShp); // Ymin
            outGlobalBoxXYZ[3] = readLEDouble(fShp); // Xmax
            outGlobalBoxXYZ[4] = readLEDouble(fShp); // Ymax
            outGlobalBoxXYZ[2] = readLEDouble(fShp); // Zmin
            outGlobalBoxXYZ[5] = readLEDouble(fShp); // Zmax

            fShp.seekg(100, std::ios::beg); // 彻底移至数据区起点
            outRecords.clear();

            // 2. 循环读取记录体
            while (fShp.peek() != EOF) {
                ShapeRecord rec;
                rec.id = readBE32(fShp);
                int32_t contentLengthInWords = readBE32(fShp);
                std::streampos recordContentStart = fShp.tellg();

                rec.shapeType = readLE32(fShp);

                // 处理 Point (2D: 1 | 3D: 11)
                if (rec.shapeType == 1 || rec.shapeType == 11) {
                    ShapePoint p;
                    p[0] = readLEDouble(fShp);
                    p[1] = readLEDouble(fShp);
                    p[2] = (rec.shapeType == 11) ? readLEDouble(fShp) : 0.0;
                    rec.points.push_back(p);
                }
                    // 处理 PolyLine / Polygon (2D: 3, 5 | 3D: 13, 15)
                else if (rec.shapeType == 3 || rec.shapeType == 5 || rec.shapeType == 13 || rec.shapeType == 15) {
                    rec.box[0] = readLEDouble(fShp);
                    rec.box[1] = readLEDouble(fShp);
                    rec.box[2] = readLEDouble(fShp);
                    rec.box[3] = readLEDouble(fShp);

                    int32_t numParts = readLE32(fShp);
                    int32_t numPoints = readLE32(fShp);

                    rec.parts.resize(numParts);
                    for (int i = 0; i < numParts; ++i)
                        rec.parts[i] = readLE32(fShp);

                    rec.points.resize(numPoints);
                    for (int i = 0; i < numPoints; ++i) {
                        rec.points[i][0] = readLEDouble(fShp);
                        rec.points[i][1] = readLEDouble(fShp);
                        rec.points[i][2] = 0.0; // 先初始化为默认值
                    }
                    // 如果是 3D 几何，继续回读 Z 轴信息
                    if (rec.shapeType == 13 || rec.shapeType == 15) {
                        double recZMin = readLEDouble(fShp); // 跳过/读取单要素 Zmin
                        double recZMax = readLEDouble(fShp); // 跳过/读取单要素 Zmax
                        for (int i = 0; i < numPoints; ++i) {
                            rec.points[i][2] = readLEDouble(fShp); // 填充实际Z值
                        }
                    }
                }
                outRecords.push_back(rec);

                // 移动到下一条记录的边界位置（规避可能存在的填充字节）
                fShp.seekg(recordContentStart + static_cast<std::streamoff>(contentLengthInWords * 2), std::ios::beg);
            }
            fShp.close();
        }

        if (fileExists(dbfPath)) {
            // 2. 读取属性数据 (.dbf)
            std::string dbfPath = changeExtension(shpPath, ".dbf");
            std::ifstream fDbf(dbfPath, std::ios::binary);
            if (!fDbf.is_open())
                return true; // DBF 丢失时允许只保留几何结构

            fDbf.seekg(4, std::ios::beg);
            int32_t dbfRecordsCount = readLE32(fDbf);
            int16_t headerLength = readLE16(fDbf);
            int16_t recordLength = readLE16(fDbf);

            // 设置有效记录空间
            bool hasRecords = outRecords.size() > 0;
            if (!hasRecords) {
                outRecords.resize(dbfRecordsCount);
            }

            // 计算字段数量
            int32_t fieldCount = (headerLength - 32 - 1) / 32;
            outFields.resize(fieldCount);
            fDbf.seekg(32, std::ios::beg);

            // 解析描述符定义
            for (int32_t i = 0; i < fieldCount; ++i) {
                char nameBuf[12] = {0};
                fDbf.read(nameBuf, 11);
                outFields[i].name = trim(std::string(nameBuf));
                fDbf.read(&outFields[i].type, 1);
                fDbf.seekg(4, std::ios::cur); // 跳过数据地址
                fDbf.read(reinterpret_cast<char *>(&outFields[i].length), 1);
                fDbf.read(reinterpret_cast<char *>(&outFields[i].decimals), 1);
                fDbf.seekg(14, std::ios::cur); // 跳过保留字段
            }

            // 定位到数据记录起点
            fDbf.seekg(headerLength, std::ios::beg);
            for (int32_t i = 0; i < dbfRecordsCount && i < static_cast<int32_t>(outRecords.size()); ++i) {
                char deleteFlag;
                fDbf.read(&deleteFlag, 1); // 0x20 有效，0x2A 已删除

                if (!hasRecords)
                    outRecords[i].id = i;
                outRecords[i].attributes.resize(fieldCount);
                for (int32_t j = 0; j < fieldCount; ++j) {
                    std::vector<char> valBuf(outFields[j].length + 1, 0);
                    fDbf.read(valBuf.data(), outFields[j].length);
                    outRecords[i].attributes[j] = trim(std::string(valBuf.data()));
                }
            }
            fDbf.close();
        }

        return true;
    }

    // =============================================================================
    // Writer（同步生成 .shp, .shx .dbf）
    // =============================================================================
    static bool write(const std::string &dir, const std::string &fname, int32_t globalShapeType,
                      const std::vector<IsgDbfField> &fields, std::vector<ShapeRecord> &records) {
        std::string shpPath = dir + "/" + fname + ".shp";
        std::string shxPath = changeExtension(shpPath, ".shx");
        std::string dbfPath = changeExtension(shpPath, ".dbf");
        std::string prjPath = changeExtension(shpPath, ".prj");

        // ---------------------------------------------------------
        // 核心一：生成几何数据 (.shp & .shx)
        // ---------------------------------------------------------
        if (globalShapeType != 0) {
            std::ofstream fShp(shpPath, std::ios::binary);
            std::ofstream fShx(shxPath, std::ios::binary);
            if (!fShp.is_open() || !fShx.is_open())
                return false;

            double globalBox[4] = {0, 0, 0, 0};
            double globalZBox[2] = {0, 0};

            if (!records.empty()) {
                // 计算局部包围盒
                for (auto &rec : records)
                    calculateRecordBox(rec);

                // 初始化全局全空间外包络
                bool inited = false;
                for (const auto &rec : records) {
                    for (const auto &p : rec.points) {
                        if (!inited) {
                            globalBox[0] = globalBox[2] = p[0];
                            globalBox[1] = globalBox[3] = p[1];
                            globalZBox[0] = globalZBox[1] = p[2];
                            inited = true;
                        }
                        if (p[0] < globalBox[0])
                            globalBox[0] = p[0];
                        if (p[1] < globalBox[1])
                            globalBox[1] = p[1];
                        if (p[0] > globalBox[2])
                            globalBox[2] = p[0];
                        if (p[1] > globalBox[3])
                            globalBox[3] = p[1];
                        if (p[2] < globalZBox[0])
                            globalZBox[0] = p[2];
                        if (p[2] > globalZBox[1])
                            globalZBox[1] = p[2];
                    }
                }
            }

            // 2. 预留 100 字节的文件头空间
            char zeroHeader[100] = {0};
            fShp.write(zeroHeader, 100);
            fShx.write(zeroHeader, 100);

            std::vector<std::pair<int32_t, int32_t>> indexMap; // 存储 shx 的 [Offset, Length] 对

            // 3. 循环写入每条记录
            int32_t currentRecordId = 1;
            for (auto &rec : records) {
                if (rec.points.empty()) {
                    continue;
                }
                // 记录相对 shp 文件头的起始位置 (Word 单位)
                int32_t shpOffsetInWords = static_cast<int32_t>(fShp.tellp()) / 2;

                // 计算该记录的内容长度 (不含记录头的 8 字节)
                int32_t contentBytes = 4; // Shape Type (int32) 占用 4 字节

                // 计算精准字节长度
                if (globalShapeType == 1) {
                    contentBytes += 16; // X(8) + Y(8)
                } else if (globalShapeType == 11) {
                    // PointZ
                    contentBytes += 24; // X(8) + Y(8) + Z(8)
                } else if (globalShapeType == 3 || globalShapeType == 5) {
                    // PolyLine, Polygon
                    contentBytes += 32; // Box (4 * 8)
                    contentBytes += 4; // NumParts
                    contentBytes += 4; // NumPoints
                    contentBytes += static_cast<int32_t>(rec.parts.size()) * 4;
                    contentBytes += static_cast<int32_t>(rec.points.size()) * 16;
                } else if (globalShapeType == 13 || globalShapeType == 15) {
                    // PolyLineZ, PolygonZ
                    contentBytes += 32; // Box (4 * 8)
                    contentBytes += 4; // NumParts
                    contentBytes += 4; // NumPoints
                    contentBytes += static_cast<int32_t>(rec.parts.size()) * 4;
                    contentBytes += static_cast<int32_t>(rec.points.size()) * 16;
                    contentBytes += 16; // Record Zmin(8) + Zmax(8)
                    contentBytes += static_cast<int32_t>(rec.points.size()) * 8; // Z Array
                }
                int32_t contentLengthInWords = contentBytes / 2;

                // 写入记录头 (Big Endian)
                writeBE32(fShp, currentRecordId++);
                writeBE32(fShp, contentLengthInWords);

                // 写入记录内容 (Little Endian)
                writeLE32(fShp, globalShapeType);
                if (globalShapeType == 1 || globalShapeType == 11) {
                    writeLEDouble(fShp, rec.points[0][0]);
                    writeLEDouble(fShp, rec.points[0][1]);
                    if (globalShapeType == 11) {
                        writeLEDouble(fShp, rec.points[0][2]);
                    }
                } else {
                    // 写入共通的 2D 描述块
                    writeLEDouble(fShp, rec.box[0]);
                    writeLEDouble(fShp, rec.box[1]);
                    writeLEDouble(fShp, rec.box[2]);
                    writeLEDouble(fShp, rec.box[3]);

                    writeLE32(fShp, static_cast<int32_t>(rec.parts.size()));
                    writeLE32(fShp, static_cast<int32_t>(rec.points.size()));

                    for (int32_t pIdx : rec.parts)
                        writeLE32(fShp, pIdx);
                    for (const auto &p : rec.points) {
                        writeLEDouble(fShp, p[0]);
                        writeLEDouble(fShp, p[1]);
                    }

                    // 如果是 3D 类型，继续追加写入 Z 轴特有描述块
                    if (globalShapeType == 13 || globalShapeType == 15) {
                        double recZMin = rec.points[0][2];
                        double recZMax = rec.points[0][2];
                        for (const auto &p : rec.points) {
                            if (p[2] < recZMin)
                                recZMin = p[2];
                            if (p[2] > recZMax)
                                recZMax = p[2];
                        }
                        writeLEDouble(fShp, recZMin);
                        writeLEDouble(fShp, recZMax);
                        for (const auto &p : rec.points)
                            writeLEDouble(fShp, p[2]);
                    }
                }

                // 保存索引信息
                indexMap.push_back({shpOffsetInWords, contentLengthInWords});
            }

            // 4. 写入 shx 索引文件
            for (const auto &idx : indexMap) {
                writeBE32(fShx, idx.first);
                writeBE32(fShx, idx.second);
            }

            // 5. 回溯并填写完整的 SHP 与 SHX 文件头
            int32_t totalShpLengthInWords = static_cast<int32_t>(fShp.tellp()) / 2;
            int32_t totalShxLengthInWords = static_cast<int32_t>(fShx.tellp()) / 2;

            // 重写完整的 100 字节全局文件头
            auto writeHeaderLambda = [&](std::ofstream &stream, int32_t totalLengthWords) {
                stream.seekp(0, std::ios::beg);
                writeBE32(stream, 9994); // File Code
                for (int i = 0; i < 5; ++i)
                    writeBE32(stream, 0); // Unused
                writeBE32(stream, totalLengthWords); // File Length
                writeLE32(stream, 1000); // Version
                writeLE32(stream, globalShapeType);
                writeLEDouble(stream, globalBox[0]); //Xmin
                writeLEDouble(stream, globalBox[1]); //Ymin
                writeLEDouble(stream, globalBox[2]); //Xmax
                writeLEDouble(stream, globalBox[3]); //Ymax
                // Z 和 M 边界框 (共计 4 double)，规范允许在不使用时将其填为 0.0
                writeLEDouble(stream, globalZBox[0]); // Zmin
                writeLEDouble(stream, globalZBox[1]); // Zmax
                writeLEDouble(stream, 0.0); // Mmin
                writeLEDouble(stream, 0.0); // Mmax
            };
            writeHeaderLambda(fShp, totalShpLengthInWords);
            writeHeaderLambda(fShx, totalShxLengthInWords);
            fShp.close();
            fShx.close();

            // 添加投影文件
            writePrj(prjPath);
        }

        // ---------------------------------------------------------
        // 核心二：构建并写入属性表数据 (.dbf)
        // ---------------------------------------------------------
        std::ofstream fDbf(dbfPath, std::ios::binary);
        if (!fDbf.is_open())
            return false;
        int16_t dbfHeaderLength = static_cast<int16_t>(32 + fields.size() * 32 + 1);
        int16_t dbfRecordLength = 1; // 1 字节的删除标记位
        for (const auto &field : fields)
            dbfRecordLength += field.length;

        // 1. 写入 DBF 基础头信息 (32 字节)
        fDbf.put(0x03); // Version: dBase III w/o memo
        fDbf.put(126); // Last update Year (2026 - 1900 = 126)
        fDbf.put(5); // Month
        fDbf.put(23); // Day
        writeLE32(fDbf, static_cast<int32_t>(records.size()));
        writeLE16(fDbf, dbfHeaderLength);
        writeLE16(fDbf, dbfRecordLength);
        char dbfReserved[20] = {0};
        fDbf.write(dbfReserved, 20);

        // 2. 写入字段描述列表 (每条 32 字节)
        for (const auto &field : fields) {
            char nameBuf[11] = {0};
            std::strncpy(nameBuf, field.name.c_str(), 11);
            fDbf.write(nameBuf, 11);
            fDbf.put(field.type);
            writeLE32(fDbf, 0); // Data Address offset
            fDbf.put(field.length);
            fDbf.put(field.decimals);
            char fieldReserved[14] = {0};
            fDbf.write(fieldReserved, 14);
        }
        fDbf.put(0x0D); // 字段区域终止标记符

        // 3. 写入各个要素的属性记录体
        for (const auto &rec : records) {
            fDbf.put(0x20); // 正常活跃数据标记 (空格)
            for (size_t j = 0; j < fields.size(); ++j) {
                std::string s = (j < rec.attributes.size()) ? rec.attributes[j] : "";
                if (s.length() > fields[j].length) {
                    s = s.substr(0, fields[j].length);
                } else {
                    // 对齐填充规范：数字右对齐，文本左对齐
                    if (fields[j].type == 'N' || fields[j].type == 'F') {
                        s = std::string(fields[j].length - s.length(), ' ') + s;
                    } else {
                        s = s + std::string(fields[j].length - s.length(), ' ');
                    }
                }
                fDbf.write(s.c_str(), fields[j].length);
            }
        }
        fDbf.put(0x1A); // DBF 文件终止标记符 (Ctrl+Z)
        fDbf.close();

        return true;
    }
};

inline std::string toWKT(std::vector<std::array<double, 2>> points, std::string id = "") {
    std::stringstream ss;
    ss << std::setprecision(12) << (id.length() > 0 ? id + ", " : "") << "\"LINESTRING (";
    for (int i = 0; i < points.size(); ++i)
        ss << (i == 0 ? "" : ", ") << points[i][0] << " " << points[i][1];
    ss << ")\"";
    return ss.str();
};

inline std::string toWKT(std::vector<std::array<double, 3>> points, std::string id = "") {
    std::stringstream ss;
    ss << std::setprecision(12) << (id.length() > 0 ? id + ", " : "") << "\"LINESTRING (";
    for (int i = 0; i < points.size(); ++i)
        ss << (i == 0 ? "" : ", ") << points[i][0] << " " << points[i][1];
    ss << ")\"";
    return ss.str();
};

inline bool WKT2Shapefile(std::string wkt_file, std::string shp_dir, std::string prefix = "") {
    // csvline to vector: support "," parse
    auto parse_csvline = [&](const std::string &line) -> std::vector<std::string> {
        std::vector<std::string> result;
        std::string current_token = "";
        bool in_quotes = false;
        for (size_t i = 0; i < line.length(); ++i) {
            char c = line[i];
            if (c == '"') {
                in_quotes = !in_quotes; // change status of ""
            } else if (c == ',' && !in_quotes) {
                result.push_back(current_token);
                current_token = "";
            } else {
                current_token += c;
            }
        }
        result.push_back(current_token);
        return result;
    };
    // WKT to coords
    auto parse_wkt = [&](const std::string &wkt, std::string &type) -> std::vector<std::array<double, 3>> {
        std::vector<std::array<double, 3>> coords;
        if (wkt.empty())
            return coords;
        // get GEOTYPE: POINT, LINESTRING
        size_t first_bracket = wkt.find('(');
        if (first_bracket == std::string::npos)
            return coords;
        type = wkt.substr(0, first_bracket);
        type.erase(std::remove_if(type.begin(), type.end(), [](unsigned char c) { return std::isspace(c); }),
                   type.end());
        std::transform(type.begin(), type.end(), type.begin(), ::toupper);
        // filter string
        std::string clean_wkt = wkt.substr(first_bracket);
        std::string numbers_only = "";
        for (char c : clean_wkt) {
            if (std::isdigit(c) || c == '.' || c == '-' || c == ' ') {
                numbers_only += c;
            } else if (c == ',') {
                numbers_only += " , ";
            }
        }
        // get X Y
        std::stringstream ss(numbers_only);
        double x, y;
        std::string separator;
        while (ss >> x >> y) {
            coords.push_back({x, y});
            if (!(ss >> separator) || separator != ",") {
                // ',' of LINESTRING | POLYGON
                if (ss.eof())
                    break;
                ss.clear();
                size_t pos = ss.tellg();
                ss.seekg(pos - separator.length());
            }
        }
        return coords;
    };

    // convert csv to shapefile
    std::vector<std::string> def_fields;
    std::map<std::string, std::vector<int32_t>> dt_ftrs;
    std::map<int32_t, std::vector<std::string>> ftr_attrs;
    std::map<int32_t, std::pair<std::string, std::vector<std::array<double, 3>>>> ftr_geoms;
    std::ifstream file(wkt_file);
    if (!file.is_open()) {
        std::cerr << "[ERROR]：open csv failed ！" << std::endl;
        return false;
    } else {
        std::string line;
        bool is_header = true;
        int32_t row = -1;
        while (std::getline(file, line)) {
            if (is_header) {
                //pass header
                is_header = false;
                continue;
            }
            if (line.empty())
                continue;
            std::vector<std::string> values = parse_csvline(line);
            if (values.empty())
                continue;
            ++row;
            for (const auto &val : values) {
                if (val.find("POINT") != std::string::npos
                    || val.find("LINE") != std::string::npos
                    || val.find("POLYGON") != std::string::npos) {
                    std::string wkt_type;
                    auto coords = parse_wkt(val, wkt_type);
                    ftr_geoms[row] = {wkt_type, coords};
                } else {
                    ftr_attrs[row].emplace_back(val);
                }
            }
        }
        file.close();
    };
    return true;
}

/**
 * @brief 根据起点和方向向量生成一条向量线
 * @param pos 起点坐标
 * @param v 方向向量
 * @param length 向量线长度
 * @return 向量线坐标串
 */
inline std::vector<std::array<double, 3>> genVectorline(
        const std::array<double, 3> &pos, const std::array<double, 3> &v, double length) {
    double vx = v[0], vy = v[1], vz = v[2];
    double vlen = std::sqrt(vx * vx + vy * vy + vz * vz);
    if (vlen < 1e-9) { return {pos, pos}; }
    // end =（pos + v[i]/vlen * length） //v[i]/vlen: unit vector
    std::array<double, 3> end;
    end[0] = pos[0] + vx / vlen * length;
    end[1] = pos[1] + vy / vlen * length;
    end[2] = pos[2] + vz / vlen * length;
    return {pos, end};
}
/**
 * @brief 根据起点和方向角生成一条向量线
 * @param pos 起点坐标
 * @param theta 方向角
 * @param length 向量线长度
 * @return 向量线坐标串
 */
inline std::vector<std::array<double, 3>> genVectorline(
        const std::array<double, 3> &pos, double theta, double length) {
    std::array<double, 3> end = {pos[0] + length * std::cos(theta), pos[1] + length * std::sin(theta), pos[2]};
    return {pos, end};
}
