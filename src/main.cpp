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
    auto grid = loadHandWrittenMap(1);
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

    print_grid_with_path(grid, start, goal, path);
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

        std::cout << "计算走廊段: ("
                  << s.start_idx << ", " << s.end_idx << ") ,"
                  << " polygon size = " << corridor.size() << "\n";
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
        std::cout << "已生成 corridor_segment" << i << ".ppm\n";
    }

    // -----------------------------
    // 5. 使用 OSQP 进行多项式优化 (最小化 Vel^2 + Acc^2)
    // -----------------------------
    TrajectoryOptimizer optimizer;
    OptimizationConfig opt_config;
    opt_config.avg_vel = 2.0; 
    opt_config.max_vel = 3.0; 
    opt_config.max_acc = 2.0; 
    opt_config.w_vel = 1.0;   // 速度权重
    opt_config.w_acc = 3.0;   // 加速度权重

    Trajectory traj = optimizer.solvePolynomial(path, corridors, segments, opt_config);

    if (!traj.pieces.empty()) {
        std::cout << "优化成功！总时长: " << traj.getTotalDuration() << "s\n";
        
        // ==========================================
        //  新增：打印每一段的数学表达式
        // ==========================================
        std::cout << "\n================ 轨迹数学表达式 ================\n";
        std::cout << "公式通式 (5阶贝塞尔): P(t) = Σ [ C_i * (1-u)^(5-i) * u^i ]\n";
        std::cout << "其中 u = t / Duration (归一化时间), C_i = binomial(5,i) * ControlPoint_i\n";
        
        // 5阶二项式系数
        static const int Binom[6] = {1, 5, 10, 10, 5, 1};

        for (size_t k = 0; k < traj.pieces.size(); ++k) {
            const auto& seg = traj.pieces[k];
            double T = seg.duration;

            std::cout << "\n--- 第 " << k << " 段 (时长 T = " << std::fixed << std::setprecision(2) << T << " s) ---\n";
            
            // 打印控制点信息
            std::cout << "  控制点 P0-P5: ";
            for(const auto& p : seg.control_points) {
                std::cout << "(" << p.x << "," << p.y << ") ";
            }
            std::cout << "\n\n";

            // Lambda 函数：打印单个维度的方程
            auto print_poly_equation = [&](const char* axis, auto get_val) {
                std::cout << "  " << axis << "(t) = ";
                for (int i = 0; i <= 5; ++i) {
                    double val = get_val(seg.control_points[i]);
                    double coef = Binom[i] * val;
                    
                    if (i > 0) {
                        if (coef >= 0) std::cout << " + ";
                        else std::cout << " "; // 负号会随数字打印出来
                    }
                    
                    // 打印项：Coef * (1-t/T)^(5-i) * (t/T)^i
                    std::cout << std::defaultfloat << coef;
                    
                    if (5 - i > 0) {
                        if (5 - i == 1) std::cout << "*(1 - t/" << std::fixed << std::setprecision(2) << T << ")";
                        else            std::cout << "*(1 - t/" << std::fixed << std::setprecision(2) << T << ")^" << (5 - i);
                    }
                    
                    if (i > 0) {
                        if (i == 1) std::cout << "*(t/" << std::fixed << std::setprecision(2) << T << ")";
                        else        std::cout << "*(t/" << std::fixed << std::setprecision(2) << T << ")^" << i;
                    }
                }
                std::cout << "\n";
            };

            print_poly_equation("x", [](const P2& p){ return p.x; });
            print_poly_equation("y", [](const P2& p){ return p.y; });
        }
        std::cout << "===============================================\n\n";

        // 采样并绘图
        std::vector<Point> display_path;
        double total_T = traj.getTotalDuration();
        for (double t = 0; t <= total_T; t += 0.05) {
            auto p = traj.at(t);
            display_path.push_back({(int)p.first, (int)p.second});
        }
        write_ppm_scaled("optimized_poly_path.ppm", grid, start, goal, display_path, 20);
        std::cout << "已保存 optimized_poly_path.ppm\n";
    } else {
        std::cout << "优化失败。\n";
    }

    return 0;
}