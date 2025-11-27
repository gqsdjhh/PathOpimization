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