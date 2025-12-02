#include <iostream>
#include <OsqpEigen/OsqpEigen.h>
#include <cmath>

#include "optimizer.hpp"
#include "polygon_corridor.hpp"

// ---------------------------------------------------------
// 辅助函数
// ---------------------------------------------------------
inline double getDist(Point a, Point b) {
    return std::hypot(a.x - b.x, a.y - b.y);
}

// 辅助：将多边形转换为半平面 Ax + By <= C
std::pair<Eigen::MatrixXd, Eigen::VectorXd> 
TrajectoryOptimizer::convertPolygonToConstraints(const Polygon& poly) {
    int n = poly.size();
    Eigen::MatrixXd AB(n, 2);
    Eigen::VectorXd C(n);
    double center_x = 0.0, center_y = 0.0;
    for (const auto& p : poly) {
        center_x += p.x;
        center_y += p.y;
    }
    center_x /= n;
    center_y /= n;
    for (int i = 0; i < n; ++ i) {
        P2 p1 = poly[i];
        P2 p2 = poly[(i + 1) % n];
        double Kx = p1.x - p2.x;
        double Ky = p1.y - p2.y;
        double A = -Ky;
        double B = Kx;
        AB(i, 0) = A;
        AB(i, 1) = B;
        C(i) = (A * p1.x + B * p1.y);
        if (AB(i, 0) * center_x + AB(i, 1) * center_y > C(i)) {
            AB(i, 0) = - AB(i, 0);
            AB(i, 1) = - AB(i, 1);
            C(i) = - C(i);
        }
    }
    return std::make_pair(AB, C);
}

// 辅助：计算单段多项式曲线需要的放缩系数
double computeTimeScalingFactor(const PolynomialCurve& curve, double v_max) {
    double max_vel = 0.0;
    double t = curve.duration_;
    double dt = 0.05;               // 采样时间间隔

    for (double i = 0.0; i < t; i += dt) {
        std::pair<double, double> p1 = curve.evaluate(i);
        std::pair<double, double> p2 = curve.evaluate(std::min(i + dt, t)); 

        double vx = (p2.first - p1.first) / dt;
        double vy = (p2.second - p1.second) / dt;
        double vel = std::sqrt(vx*vx + vy*vy);
        if (vel > max_vel) max_vel = vel;
    }
    return max_vel / v_max;
}

// ---------------------------------------------------------
// 核心求解逻辑
// ---------------------------------------------------------
Trajectory TrajectoryOptimizer::solvePolynomial(
    const std::vector<Point>& raw_path,
    const std::vector<Polygon>& corridors,
    const std::vector<Segment>& segments,
    OptimizationConfig config
) {
    if (segments.empty()) return Trajectory();

    int segment_num = segments.size();
    int polynomial_order = 4; // 三次多项式 
    int total_vars = segment_num * polynomial_order * 2; 

    // 2. 时间分配
    std::vector<double> T(segment_num);
    for (int i = 0; i < segment_num; ++i) {
        double dist = getDist(segments[i].start, segments[i].end);

        T[i] = std::max(dist / config.avg_vel, 0.5); 
    }

    OsqpEigen::Solver solver;
    solver.settings()->setVerbosity(false);
    solver.settings()->setWarmStart(true);

    // 3. 构建目标函数矩阵 P 
    Eigen::SparseMatrix<double> P(total_vars, total_vars);
    std::vector<Eigen::Triplet<double>> p_triplets;

    // minimum jerk
    for (size_t i = 0; i < segment_num; ++i) {
        double t5 = std::pow(T[i], 5);

        double value = 72.0 / t5; 
        for (int axis = 0; axis < 2; ++axis) {  
            int row = i * polynomial_order * 2 + axis * polynomial_order + 3;

            p_triplets.emplace_back(row, row, value);
        }
    }
    P.setFromTriplets(p_triplets.begin(), p_triplets.end());

    Eigen::VectorXd q = Eigen::VectorXd::Zero(total_vars);

    // 4. 构建约束
    std::vector<Eigen::Triplet<double>> A_triplets;
    std::vector<double> l_vec, u_vec;
    int constraint_idx = 0;

    auto add_constraint = [&](int row, int col, double value) {
        A_triplets.emplace_back(row, col, value);
    };

    auto set_rng = [&](double lb, double ub) {
        l_vec.push_back(lb); 
        u_vec.push_back(ub); 
        constraint_idx++;
    };


    // 起点状态
    for (size_t axis = 0; axis < 2; ++axis) {
        int column = axis * polynomial_order;

        // Pos
        add_constraint(constraint_idx, column, 1.0);
        double pos = (axis == 0) ? segments[0].start.x : segments[0].start.y;
        set_rng(pos, pos);

        // Vel
        add_constraint(constraint_idx, column + 1, 1.0);
        set_rng(0.0, 0.0);  
    }


    // 终点状态
    for (size_t axis = 0; axis < 2; ++axis) {
        int column = (segment_num - 1) * polynomial_order * 2 + axis * polynomial_order;

        // Pos: c0 + c1*t + c2*t^2 + c3*t^3
        add_constraint(constraint_idx, column + 0, 1.0);
        add_constraint(constraint_idx, column + 1, 1.0);
        add_constraint(constraint_idx, column + 2, 1.0);
        add_constraint(constraint_idx, column + 3, 1.0);
        double pos = (axis == 0) ? segments.back().end.x : segments.back().end.y;
        set_rng(pos, pos);

        // Vel: c1 + 2*c2*t + 3*c3*t^2
        add_constraint(constraint_idx, column + 1, 1.0);
        add_constraint(constraint_idx, column + 2, 2.0);
        add_constraint(constraint_idx, column + 3, 3.0);
        set_rng(0.0, 0.0);  
    }

    
    // 连续性
    for (int i = 0; i < segment_num - 1; ++i) {
        double T_cur = T[i];
        double T_next = T[i+1];

        for (int axis = 0; axis < 2; ++axis) {
            int index_current = i * polynomial_order * 2 + axis * polynomial_order;
            int index_next = (i + 1) * polynomial_order * 2 + axis * polynomial_order;

            // 位置连续
            // c_{i,0} + c_{i,1} + c_{i,2} + c_{i,3}) - c_{i+1,0} = 0
            add_constraint(constraint_idx, index_current + 0, 1.0); 
            add_constraint(constraint_idx, index_current + 1, 1.0);   
            add_constraint(constraint_idx, index_current + 2, 1.0);  
            add_constraint(constraint_idx, index_current + 3, 1.0);  
            add_constraint(constraint_idx, index_next + 0, -1.0);   
            set_rng(0.0, 0.0);

            // 速度连续
            // T_{i+1} * ( c_{i,1} + 2c_{i,2} + 3c_{i,3} ) - T_i * c_{i+1,1} = 0
            add_constraint(constraint_idx, index_current + 1, 1.0 * T_next);
            add_constraint(constraint_idx, index_current + 2, 2.0 * T_next);
            add_constraint(constraint_idx, index_current + 3, 3.0 * T_next);
            add_constraint(constraint_idx, index_next + 1, -1.0 * T_cur);
            set_rng(0.0, 0.0);

            // 加速度连续
            // T_{i+1}^2 * ( 2c_{i,2} + 6c_{i,3}) - T_i^2 * 2c_{i + 1, c2} = 0
            add_constraint(constraint_idx, index_current + 2, 2.0 * T_next * T_next);
            add_constraint(constraint_idx, index_current + 3, 6.0 * T_next * T_next);
            add_constraint(constraint_idx, index_next + 2, -2.0 * T_cur * T_cur);
            set_rng(0.0, 0.0);
        }
    }


    // 安全走廊约束
    constexpr int samples_per_segment = 3;  

    for (int k = 0; k < segment_num; ++k) {
        auto [Mat, Ub] = convertPolygonToConstraints(corridors[k]); 
        int num_planes = Mat.rows();
        if (num_planes == 0) continue;

        // 在每个线段上均匀采样
        double dt_sample = T[k] / (samples_per_segment + 1);

        for (int s = 0; s <= samples_per_segment; ++s) {
            double t = s * dt_sample;
            double s_t = t / T[k];
            double s_pow1 = s_t;
            double s_pow2 = s_t * s_t;
            double s_pow3 = s_pow2 * s_t;

            for (int r = 0; r < num_planes; ++r) {
                double A_val = Mat(r, 0);
                double B_val = Mat(r, 1);
                double C_val = Ub(r);

                int idx_x = k * polynomial_order * 2;
                int idx_y = k * polynomial_order * 2 + polynomial_order;

                // X 轴
                add_constraint(constraint_idx, idx_x + 0, A_val * 1.0);
                add_constraint(constraint_idx, idx_x + 1, A_val * s_pow1);
                add_constraint(constraint_idx, idx_x + 2, A_val * s_pow2);
                add_constraint(constraint_idx, idx_x + 3, A_val * s_pow3);

                // Y 轴
                add_constraint(constraint_idx, idx_y + 0, B_val * 1.0);
                add_constraint(constraint_idx, idx_y + 1, B_val * s_pow1);
                add_constraint(constraint_idx, idx_y + 2, B_val * s_pow2);
                add_constraint(constraint_idx, idx_y + 3, B_val * s_pow3);

                set_rng(-OsqpEigen::INFTY, C_val);
            }
        }
    }

    // 5. 求解
    Eigen::SparseMatrix<double> A(constraint_idx, total_vars);
    A.setFromTriplets(A_triplets.begin(), A_triplets.end());
    
    Eigen::VectorXd l_e(l_vec.size()), u_e(u_vec.size());
    for(size_t i=0; i<l_vec.size(); ++i) {
        l_e(i)=l_vec[i]; u_e(i)=u_vec[i]; 
    }

    solver.data()->setNumberOfVariables(total_vars);
    solver.data()->setNumberOfConstraints(constraint_idx);
    if(!solver.data()->setHessianMatrix(P)) return {};
    if(!solver.data()->setGradient(q)) return {};
    if(!solver.data()->setLinearConstraintsMatrix(A)) return {};
    if(!solver.data()->setLowerBound(l_e)) return {};
    if(!solver.data()->setUpperBound(u_e)) return {};

    if(!solver.initSolver() || solver.solveProblem() != OsqpEigen::ErrorExitFlag::NoError) {
        std::cerr << "OSQP Solve Error!" << std::endl;
        return {};
    }

    // 6. 提取结果
    Eigen::VectorXd sol = solver.getSolution();
    Trajectory trajectory;
    for (int k = 0; k < segment_num; ++k) {
        PolynomialCurve curve;
        curve.duration_ = T[k];
        
        for (int axis = 0; axis < 2; ++axis) {
            int index = k * polynomial_order * 2 + axis * polynomial_order;
            for (int i = 0; i < polynomial_order; ++i) {
                double val = sol(index + i);

                if (axis == 0) curve.polynomial_x_coefficients_.push_back(val);
                else           curve.polynomial_y_coefficients_.push_back(val);
            }
        }
        trajectory.pieces.push_back(curve);
    }

    // 7. 执行时间放缩 (Time Scaling)
    double global_ratio = 1.0;
    for(const auto& piece : trajectory.pieces) {
        double r = computeTimeScalingFactor(piece, config.max_vel);
        if(r > global_ratio) global_ratio = r;
    }

    if(global_ratio > 1.01) {
        std::cout << "[Optimizer] Time Scaling Triggered: ratio = " << global_ratio << std::endl;
        for(auto& piece : trajectory.pieces) {
            piece.duration_ *= global_ratio;
        }
    }
    else {
        std::cout << "[Optimizer] No Time Scaling Needed." << std::endl;
    }

    return trajectory;
}

