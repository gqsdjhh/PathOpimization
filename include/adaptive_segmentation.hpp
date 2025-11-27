#pragma once
#include "grid.hpp"

// 表示路径点
struct Segment {
    Point start;
    Point end;
    int start_idx;
    int end_idx;
};

// 计算两点之间的距离
double dist(const Point& a, const Point& b);
// 计算路径长度
double calculatePathLength(const std::vector<Point>& pts, int start, int end);
// 自适应路径分割函数
std::vector<Segment> adaptiveSegmentation(const std::vector<Point>& path, const std::vector<std::vector<int>>& obstacle_map, double max_segment_length = 5, double complexity_threshold = 3.0);
