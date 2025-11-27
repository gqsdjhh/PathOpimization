#pragma once
#include <vector>
#include <unordered_map>
#include "grid.hpp"

struct Node {
    Point p;
    double f;
};

struct NodeCmp {
    bool operator()(Node const& a, Node const& b) const {
        return a.f > b.f;
    }
};

double heuristic(Point a, Point b);

std::pair<
    std::unordered_map<Point, Point, PointHash>,
    std::unordered_map<Point, double, PointHash>
>
astar_grid(const std::vector<std::vector<int>>& grid, Point start, Point goal);

std::vector<Point> reconstruct_path(
    std::unordered_map<Point, Point, PointHash> &came_from,
    Point start, Point goal
);
