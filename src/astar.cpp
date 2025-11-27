#include "astar.hpp"
#include <queue>
#include <cmath>
#include <algorithm>

double heuristic(Point a, Point b) {
    return std::abs(a.x - b.x) + std::abs(a.y - b.y);
}

std::vector<Point> reconstruct_path(
    std::unordered_map<Point, Point, PointHash> &came_from,
    Point start, Point goal
){
    std::vector<Point> path;
    Point cur = goal;

    if (came_from.find(cur) == came_from.end()) return path;

    while (!(cur == start)) {
        path.push_back(cur);
        cur = came_from[cur];
    }

    path.push_back(start);
    std::reverse(path.begin(), path.end());
    return path;
}

std::pair<
    std::unordered_map<Point, Point, PointHash>,
    std::unordered_map<Point, double, PointHash>
>
astar_grid(const std::vector<std::vector<int>>& grid, Point start, Point goal)
{
    int H = grid.size();
    int W = grid[0].size();

    auto in_bounds = [&](Point p){ return p.x >= 0 && p.x < W && p.y >= 0 && p.y < H; };
    auto passable = [&](Point p){ return in_bounds(p) && grid[p.y][p.x] == 0; };

    std::priority_queue<Node, std::vector<Node>, NodeCmp> frontier;
    frontier.push({start, 0.0});

    std::unordered_map<Point, Point, PointHash> came_from;
    std::unordered_map<Point, double, PointHash> cost_so_far;

    came_from[start] = start;
    cost_so_far[start] = 0.0;

    std::vector<Point> dirs = {{1,0},{-1,0},{0,1},{0,-1}};

    while (!frontier.empty()) {
        Node curNode = frontier.top();
        frontier.pop();

        Point current = curNode.p;
        if (current == goal) break;

        for (auto d : dirs) {
            Point next{ current.x + d.x, current.y + d.y };
            if (!passable(next)) continue;

            double new_cost = cost_so_far[current] + 1.0;

            if (!cost_so_far.count(next) || new_cost < cost_so_far[next]) {
                cost_so_far[next] = new_cost;
                double priority = new_cost + heuristic(goal, next);
                frontier.push({ next, priority });
                came_from[next] = current;
            }
        }
    }
    return {came_from, cost_so_far};
}
