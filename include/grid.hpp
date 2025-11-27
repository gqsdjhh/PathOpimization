#pragma once
#include <vector>

struct Point {
    int x, y;
    bool operator==(Point const& o) const noexcept { return x == o.x && y == o.y; }
};

struct PointHash {
    size_t operator()(Point const& p) const noexcept {
        return (std::hash<int>()(p.x) * 73856093u) ^ 
               (std::hash<int>()(p.y) * 19349663u);
    }
};

// 将 RAW_MAP 转为 vector
std::vector<std::vector<int>> loadHandWrittenMap(int choice);