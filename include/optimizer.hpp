#pragma once
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "polygon_corridor.hpp" 
#include "adaptive_segmentation.hpp" 
#include "grid.hpp"
#include "polynomial.hpp"

struct OptimizationConfig {
    double avg_vel = 1.5;     // 期望平均速度 (m/s)
    double w_vel = 1.0;       // 速度代价权重
    double w_acc = 3.0;       // 加速度代价权重
};

class TrajectoryOptimizer {
public:
    TrajectoryOptimizer() = default;

    // 多项式轨迹优化入口
    Trajectory solvePolynomial(
        const std::vector<Point>& raw_path,
        const std::vector<Polygon>& corridors,
        const std::vector<Segment>& segments,
        OptimizationConfig config = OptimizationConfig()
    );

private:
    std::pair<Eigen::MatrixXd, Eigen::VectorXd> convertPolygonToConstraints(const Polygon& poly);
};