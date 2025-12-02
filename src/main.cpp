#include <iostream>
#include <iomanip> // 用于控制输出格式
#include "grid.hpp"
#include "astar.hpp"
#include "visualize.hpp"
#include "adaptive_segmentation.hpp"
#include "polygon_corridor.hpp"
#include "optimizer.hpp"

int main()
{
    auto grid = loadHandWrittenMap(2);
    int W = grid[0].size();
    int H = grid.size();

    Point start{0, 0};
    Point goal{W - 1, H - 1};

    grid[start.y][start.x] = 0;
    grid[goal.y][goal.x] = 0;

    // -----------------------------
    // 1. A* 寻路
    // -----------------------------
    auto result = astar_grid(grid, start, goal);
    auto path = reconstruct_path(result.first, start, goal);

    if (path.empty())
    {
        std::cout << "未能找到通向目标的路径。\n";
        return 0;
    }
    else
    {
        std::cout << "找到路径！路径长度 = " << path.size() - 1 << " 步\n";
    }

    write_ppm_scaled("output.ppm", grid, start, goal, path, 20);
    std::cout << "已生成 output.ppm\n";

    // -----------------------------
    // 2. 自适应分段
    // -----------------------------
    double max_len = 4.0;
    double threshold = 1;

    auto segments = adaptiveSegmentation(path, grid, max_len, threshold);
    write_ppm_segments("segments.ppm", grid, path, segments, 20);
    std::cout << "已生成 segments.ppm(分段可视化)\n";

    // -----------------------------
    // 3. 为每一段生成凸多边形安全走廊
    // -----------------------------
    std::vector<Polygon> corridors;

    for (auto &s : segments)
    {
        Polygon corridor = computeConvexCorridor(
            path,
            s.start_idx,
            s.end_idx,
            grid);

        corridors.push_back(corridor);

        // std::cout << "计算走廊段: ("
        //           << s.start_idx << ", " << s.end_idx << ") ,"
        //           << " polygon size = " << corridor.size() << "\n";
    }

    // -----------------------------
    // 4. 保存每段的 corridor 图像
    // -----------------------------
    for (int i = 0; i < corridors.size(); ++i)
    {
        if (corridors[i].empty())
        {
            std::cout << "Corridor " << i << " is empty, skip.\n";
            continue;
        }
        draw_polygon_ppm("corridor_segment" + std::to_string(i) + ".ppm",
                         corridors[i], grid, 20);
        // std::cout << "已生成 corridor_segment" << i << ".ppm\n";
    }

    // -----------------------------
    // 5. 使用 OSQP 进行多项式优化
    // -----------------------------
    TrajectoryOptimizer optimizer;
    OptimizationConfig opt_config;
    opt_config.avg_vel = 2.0; 
    opt_config.max_vel = 3.0; 
    opt_config.max_acc = 2.0; 

    Trajectory trajectory = optimizer.solvePolynomial(path, corridors, segments, opt_config);

    // -----------------------------
    // 6. 结果处理与可视化
    // -----------------------------
    if (trajectory.pieces.empty()) {
        std::cerr << "优化失败：未生成轨迹 (Optimization Failed)\n";
        return -1;
    }

    std::cout << "\n================ [优化成功] 轨迹详细信息 ================\n";
    std::cout << "多项式基底: p(t) = c0 + c1*t + c2*t^2 + c3*t^3\n";

    double total_duration = 0.0;
    
    // 打印每一段的详细数学表达式
    for (size_t k = 0; k < trajectory.pieces.size(); ++k) {
        const auto& piece = trajectory.pieces[k];
        total_duration += piece.duration_;

        std::cout << "\n--- 第 " << k << " 段 (时长 T=" << std::fixed << std::setprecision(4) << piece.duration_ << "s) ---\n";
        
        // 打印 X 轴公式
        std::cout << "  x(t) = ";
        for(size_t i = 0; i < piece.polynomial_x_coefficients_.size(); ++i) {
            double c = piece.polynomial_x_coefficients_[i];
            if(i > 0 && c >= 0) std::cout << "+";
            
            // 打印系数
            std::cout << std::scientific << std::setprecision(2) << c;
            
            if(i > 0) {
                std::cout << " * (t/" << std::fixed << std::setprecision(2) << piece.duration_ << ")^" << i;
            }
            std::cout << " ";
        }
        std::cout << "\n";

        // Y 轴
        std::cout << "  y(t) = ";
        for(size_t i = 0; i < piece.polynomial_y_coefficients_.size(); ++i) {
            double c = piece.polynomial_y_coefficients_[i];
            if(i > 0 && c >= 0) std::cout << "+";
            std::cout << std::scientific << std::setprecision(2) << c;
            if(i > 0) {
                std::cout << " * (t/" << std::fixed << std::setprecision(2) << piece.duration_ << ")^" << i;
            }
            std::cout << " ";
        }
        std::cout << "\n";
    }
    std::cout << "========================================================\n";
    std::cout << "轨迹总时长: " << std::fixed << std::setprecision(2) << total_duration << " s\n";

    // 6.2 轨迹采样 (Sampling)
    std::vector<Point> final_path_points;
    double dt = 0.05; // 采样步长 50ms

    for (const auto& piece : trajectory.pieces) {
        for (double t = 0; t < piece.duration_; t += dt) {
            std::pair<double, double> pos = piece.evaluate(t);
            final_path_points.push_back({
                static_cast<int>(std::round(pos.first)), 
                static_cast<int>(std::round(pos.second))
            });
        }
        std::pair<double, double> end_pos = piece.evaluate(piece.duration_);
        final_path_points.push_back({
            static_cast<int>(std::round(end_pos.first)), 
            static_cast<int>(std::round(end_pos.second))
        });
    }

    std::cout << "采样点数: " << final_path_points.size() << "\n";

    // 6.3 生成最终结果图像
    write_ppm_scaled("final_trajectory.ppm", grid, start, goal, final_path_points, 20);
    std::cout << "已生成 final_trajectory.ppm\n";

    return 0;
}