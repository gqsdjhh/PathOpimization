#include "polygon_corridor.hpp"
#include "grid.hpp"
#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <unordered_set>

// 计算点到线段的最短距离对应的线段上的点
// 线段: p1 -> p2, 目标点: p
P2 closestPointOnSegment(P2 p, P2 p1, P2 p2) {
    P2 segment = sub(p2, p1);
    double len2 = dot(segment, segment);
    
    if (len2 < 1e-6) return p1; // 线段退化为一个点

    P2 to_p = sub(p, p1);
    double t = dot(to_p, segment) / len2;

    if (t < 0.0) return p1;
    if (t > 1.0) return p2;
    return add(p1, mul(segment, t));
}

// -----------------------------------------------------------
// 凸多边形裁剪 (Sutherland-Hodgman 算法)
// -----------------------------------------------------------

// 定义一个半平面: normal.x * x + normal.y * y + constant >= 0 为 "内部"
struct HyperPlane {
    P2 normal;
    double c; // constant offset
    
    // 判断点是否在半平面内（或边缘）
    bool isInside(P2 p) const {
        return (normal.x * p.x + normal.y * p.y + c) >= -1e-9; 
    }

    // 计算直线与线段(p1-p2)的交点
    P2 intersect(P2 p1, P2 p2) const {
        double d1 = normal.x * p1.x + normal.y * p1.y + c;
        double d2 = normal.x * p2.x + normal.y * p2.y + c;
        
        // 避免除以零（虽然逻辑上只有跨越平面时才调用，不会平行）
        double t = d1 / (d1 - d2); 
        return add(p1, mul(sub(p2, p1), t));
    }
};

// 裁剪函数：输入多边形，用平面裁剪，返回新多边形
Polygon clipPolygon(const Polygon& subject, const HyperPlane& plane) {
    Polygon output;
    if (subject.empty()) return output;

    for (size_t i = 0; i < subject.size(); ++i) {
        P2 curr = subject[i];
        P2 prev = subject[(i > 0) ? (i - 1) : (subject.size() - 1)];

        bool currIn = plane.isInside(curr);
        bool prevIn = plane.isInside(prev);

        if (currIn) {
            if (!prevIn) {
                // 从外进内：添加交点，再添加当前点
                output.push_back(plane.intersect(prev, curr));
            }
            output.push_back(curr);
        } else {
            if (prevIn) {
                // 从内出外：添加交点
                output.push_back(plane.intersect(prev, curr));
            }
            // 都在外：什么都不做
        }
    }
    return output;
}

// -----------------------------------------------------------
// 主函数：生成凸安全走廊
// -----------------------------------------------------------
Polygon computeConvexCorridor(
    const std::vector<Point>& path,
    int L_idx,
    int R_idx,
    const std::vector<std::vector<int>>& grid
) {
    int H = grid.size();
    int W = grid[0].size();
    
    // 1. 初始化多边形为整个地图的大矩形
    // 为了防止边界问题，稍微扩一点点
    Polygon poly;
    poly.push_back({0.0, 0.0});
    poly.push_back({(double)W, 0.0});
    poly.push_back({(double)W, (double)H});
    poly.push_back({0.0, (double)H});

    // 确定搜索范围 (Bounding Box of path segment + Buffer)
    int min_x = W, max_x = 0, min_y = H, max_y = 0;
    for (int i = L_idx; i <= R_idx; ++i) {
        min_x = std::min(min_x, path[i].x);
        max_x = std::max(max_x, path[i].x);
        min_y = std::min(min_y, path[i].y);
        max_y = std::max(max_y, path[i].y);
    }
    
    int search_dist = 8; // 搜索半径，决定走廊最大可能的宽度
    int s_min_x = std::max(0, min_x - search_dist);
    int s_max_x = std::min(W - 1, max_x + search_dist);
    int s_min_y = std::max(0, min_y - search_dist);
    int s_max_y = std::min(H - 1, max_y + search_dist);

    // 2. 遍历范围内的所有障碍物
    // 对于每个障碍物，生成一个分割平面，切掉多边形的一部分
    for (int y = s_min_y; y <= s_max_y; ++y) {
        for (int x = s_min_x; x <= s_max_x; ++x) {
            if (grid[y][x] == 0) continue; // 也可以跳过非障碍物

            P2 obs_p = {x + 0.5, y + 0.5}; // 障碍物中心

            // 2.1 找到路径段上离该障碍物最近的点
            // 这种方法比单纯找每个路径点更鲁棒，防止把弯曲路径切断
            double min_dist_sq = 1e9;
            P2 closest_path_p = {(double)path[L_idx].x, (double)path[L_idx].y};

            // 简化：在离散路径点中找最近点 (对于A*路径通常足够密集)
            // 如果需要极高精度，可以遍历每两个点组成的线段调用 closestPointOnSegment
            for (int i = L_idx; i <= R_idx; ++i) {
                P2 p = {(double)path[i].x, (double)path[i].y};
                double d2 = (p.x - obs_p.x)*(p.x - obs_p.x) + (p.y - obs_p.y)*(p.y - obs_p.y);
                if (d2 < min_dist_sq) {
                    min_dist_sq = d2;
                    closest_path_p = p;
                }
            }
            
            // 也可以选择开启“线段插值最近点”以获得更平滑的结果：
            /*
            for (int i = L_idx; i < R_idx; ++i) {
                P2 p1 = {(double)path[i].x, (double)path[i].y};
                P2 p2 = {(double)path[i+1].x, (double)path[i+1].y};
                P2 cp = closestPointOnSegment(obs_p, p1, p2);
                double d2 = (cp.x - obs_p.x)*(cp.x - obs_p.x) + (cp.y - obs_p.y)*(cp.y - obs_p.y);
                if (d2 < min_dist_sq) {
                    min_dist_sq = d2;
                    closest_path_p = cp;
                }
            }
            */

            // 2.2 构建分割向量 (从障碍物指向路径)
            P2 diff = sub(closest_path_p, obs_p);
            double dist = std::sqrt(dot(diff, diff));

            if (dist < 1e-6) continue; // 路径点和障碍物重合（异常情况）

            // 法向量归一化
            P2 normal = {diff.x / dist, diff.y / dist};

            // 2.3 定义分割点
            // 我们希望平面贴近障碍物，给路径留出最大空间
            // 分割点位置 = 障碍物中心 + 法向量 * (安全边距)
            // 安全边距可以设为 0.5 (障碍物半径) + epsilon
            double margin = 0.7; 
            P2 plane_pt = add(obs_p, mul(normal, margin));

            // 2.4 构建平面参数: n_x * x + n_y * y + c >= 0
            // 代入点 plane_pt，使得 plane_pt 刚好在边界上 (result = 0)
            // c = -(n_x * px + n_y * py)
            double c = -(normal.x * plane_pt.x + normal.y * plane_pt.y);

            HyperPlane plane = {normal, c};

            // 2.5 执行裁剪
            poly = clipPolygon(poly, plane);
            
            if (poly.empty()) break; // 多边形已消失（路径被完全堵死）
        }
    }

    return poly;
}