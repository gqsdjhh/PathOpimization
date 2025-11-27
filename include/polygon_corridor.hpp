#pragma once
#include <vector>
#include "grid.hpp"

// 2D 浮点数点
struct P2 {
    double x, y;
};

using Polygon = std::vector<P2>;

// 主接口：计算路径段的凸安全走廊
Polygon computeConvexCorridor(
    const std::vector<Point>& path,
    int L_idx,
    int R_idx,
    const std::vector<std::vector<int>>& grid
);

// -----------------------------------------------------------
// 数学辅助函数
// -----------------------------------------------------------

// 向量点积
inline double dot(P2 a, P2 b) {
    return a.x * b.x + a.y * b.y;
}

// 向量减法
inline P2 sub(P2 a, P2 b) {
    return {a.x - b.x, a.y - b.y};
}

// 向量加法
inline P2 add(P2 a, P2 b) {
    return {a.x + b.x, a.y + b.y};
}

// 向量乘标量
inline P2 mul(P2 a, double s) {
    return {a.x * s, a.y * s};
}