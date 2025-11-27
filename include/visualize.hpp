#pragma once
#include <vector>
#include <string>

#include "grid.hpp"
#include "adaptive_segmentation.hpp"
#include "polygon_corridor.hpp"

void print_grid_with_path(
    const std::vector<std::vector<int>>& grid,
    Point start, Point goal,
    const std::vector<Point>& path
);

void write_ppm_scaled(
    const std::string &filename,
    const std::vector<std::vector<int>>& grid,
    Point start, Point goal,
    const std::vector<Point>& path,
    int scale = 20
);

void write_ppm_segments(
    const std::string& filename,
    const std::vector<std::vector<int>>& grid,
    const std::vector<Point>& path,
    const std::vector<Segment>& segments,
    int scale = 20
);

void draw_polygon_ppm(
    const std::string& file,
    const Polygon& poly,
    const std::vector<std::vector<int>>& grid,
    int scale = 20
);
