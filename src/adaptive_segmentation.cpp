#include "adaptive_segmentation.hpp"
#include <cmath>
#include <set>
#include <vector>
#include <iostream>

double dist(const Point& a, const Point& b) {
    return std::hypot(a.x - b.x, a.y - b.y);
}

double calculatePathLength(const std::vector<Point>& pts, int start, int end) {
    double length = 0.0;
    for (int i = start; i < end; ++i) {
        length += dist(pts[i], pts[i + 1]);
    }
    return length;
}

double calculateEnvironmentComplexity(
        const std::vector<Point>& pts,
        int start, int end,
        const std::vector<std::vector<int>>& map,
        double radius)
{
    int H = map.size();
    int W = map[0].size();
    
    std::set<std::pair<int, int>> obstacles;

    for (int i = start; i <= end; ++i) {
        int cx = pts[i].x;
        int cy = pts[i].y;

        for (int dx = -radius; dx <= radius; ++dx) {
            for (int dy = -radius; dy <= radius; ++dy) {

                int nx = cx + dx;
                int ny = cy + dy;

                if (nx >= 0 && nx < W && ny >= 0 && ny < H) {
                    if (map[ny][nx] == 1) {
                        obstacles.insert({nx, ny});
                    }
                }
            }
        }
    }

    int obstacle_count = obstacles.size();
    int point_count = end - start + 1;

    return static_cast<double>(obstacle_count) / (radius * point_count);
}


std::vector<Segment> adaptiveSegmentation(
        const std::vector<Point>& path,
        const std::vector<std::vector<int>>& map,
        double max_segment_length,
        double complexity_threshold)
{
    std::vector<Segment> segments;
    int N = path.size();

    if (N < 2) return segments;

    int current_start = 0;

    while (current_start < N - 1) {

        int j = current_start + 1;
        int valid_end = j;   // 最近一次满足要求的 j

        while (j < N) {

            double len = calculatePathLength(path, current_start, j);
            double comp = calculateEnvironmentComplexity(path, current_start, j, map, 1.0);

            if (len > max_segment_length || comp > complexity_threshold)
                break;

            valid_end = j;
            j++;
        }

        // 确保至少两个点
        if (valid_end == current_start) valid_end = current_start + 1;

        // 生成段
        segments.push_back(
            { path[current_start], path[valid_end], current_start, valid_end }
        );

        // 更新下一段起点
        current_start = valid_end;
    }

    return segments;
}
