#include <iostream>
#include <OsqpEigen/OsqpEigen.h>
#include <cmath>

#include "optimizer.hpp"
#include "polygon_corridor.hpp"

// ---------------------------------------------------------
// 辅助函数
// ---------------------------------------------------------
double getDist(Point a, Point b) {
    return std::hypot(a.x - b.x, a.y - b.y);
}

// 计算组合数 C(n, k)
double binomial(int n, int k) {
    if (k < 0 || k > n) return 0;
    double res = 1;
    for (int i = 1; i <= k; ++i) res = res * (n - i + 1) / i;
    return res;
}

// 计算 Bernstein 基函数 B_{i,n}(u) 在 u 处的 k 阶导数值
double bernstein_poly_deriv(int primitive_order, int i, double u, int target_derivative_order) {
    if (target_derivative_order == 0) {
        return binomial(primitive_order, i) * std::pow(u, i) * std::pow(1 - u, primitive_order - i);
    }
    // 递归公式: B'_{i,n}(u) = n * (B_{i-1, n-1}(u) - B_{i, n-1}(u))
    return primitive_order * (bernstein_poly_deriv(primitive_order - 1, i - 1, u, target_derivative_order - 1)
            - bernstein_poly_deriv(primitive_order - 1, i, u, target_derivative_order - 1));
}

// 通用：计算 k 阶导数的 Hessian 矩阵 (Degree N)
// H(i, j) = Integral(0~1) [ B^(k)_i(u) * B^(k)_j(u) ] du
Eigen::MatrixXd computeBernsteinHessian(int n_degree, int target_derivative_order) {
    //+1是因为N阶曲线有N+1个控制点
    Eigen::MatrixXd H(n_degree + 1, n_degree + 1);
    int steps = 100; // 积分步数
    double dt = 1.0 / steps;
    
    for (int i = 0; i < n_degree + 1; ++i) {
        for (int j = i; j < n_degree + 1; ++j) {
            double sum = 0.0;

            for(int t = 0; t < steps; ++t) {
                // 数值积分时中点法一般比左端点法精度更高，这里使用左端点法会导致误差过大，H矩阵失去正定，OSQP求解报错
                double u = (t + 0.5) * dt;
                double bi = bernstein_poly_deriv(n_degree, i, u, target_derivative_order);
                double bj = bernstein_poly_deriv(n_degree, j, u, target_derivative_order);

                sum += bi * bj * dt;
            }

            H(i, j) = sum;
            H(j, i) = sum; // 对称矩阵
        }
    }

    return H;
}

// 辅助：将多边形转换为半平面 Ax + By <= C
std::pair<Eigen::MatrixXd, Eigen::VectorXd> 
TrajectoryOptimizer::convertPolygonToConstraints(const Polygon& poly) {
    int n = poly.size();

    Eigen::MatrixXd AB(n, 2);
    Eigen::VectorXd C(n);

    // 计算多边形中心点
    double center_x = 0.0, center_y = 0.0;
    for (const auto& p : poly) {
        center_x += p.x;
        center_y += p.y;
    }
    center_x /= n;
    center_y /= n;

    // 循环访问每一个顶点p1和下一个顶点p2这两个点构成多边形的一条边。
    for (int i = 0; i < n; ++ i) {
        P2 p1 = poly[i];
        P2 p2 = poly[(i + 1) % n];

        // 求方向向量
        double Kx = p1.x - p2.x;
        double Ky = p1.y - p2.y;

        // 利用方向向量 * 法向量 = 0 求法向量
        double A = -Ky;
        double B = Kx;

        AB(i, 0) = A;
        AB(i, 1) = B;

        C(i) = (A * p1.x + B * p1.y);
        
        // 确保中心点在半平面内
        if (AB(i, 0) * center_x + AB(i, 1) * center_y > C(i)) {
            AB(i, 0) = - AB(i, 0);
            AB(i, 1) = - AB(i, 1);
            C(i) = - C(i);
        }
    }

    return std::make_pair(AB, C);
}

// 辅助：计算单段贝塞尔曲线需要的放缩系数
double computeTimeScalingFactor(const BezierCurve& curve, double v_max, double a_max) {
    // 1. 提取控制点
    std::vector<P2> P;
    for(auto& p : curve.control_points) P.push_back({(double)p.x, (double)p.y});
    
    int n = 5; // 5阶
    double T = curve.duration;

    // 2. 计算速度控制点 (Order 4)
    // V_i = n * (P_{i+1} - P_i)
    // 真实速度 V(t) = Bezier(V_cp, u) / T
    double max_vel_norm = 0.0;
    std::vector<P2> V_cp; 
    for(int i=0; i<n; ++i) {
        P2 v = sub(P[i+1], P[i]); // sub函数来自你的 polygon_corridor.cpp，或手动写
        // 注意：这里先不除以 T，最后统一算
        double val = std::sqrt(dot(v, v)) * n; 
        V_cp.push_back({v.x * n, v.y * n}); // 存下来用于算加速度
        
        if(val > max_vel_norm) max_vel_norm = val;
    }

    // 3. 计算加速度控制点 (Order 3)
    // A_i = (n-1) * (V_{i+1} - V_i)
    // 真实加速度 A(t) = Bezier(A_cp, u) / T^2
    double max_acc_norm = 0.0;
    for(int i=0; i<n-1; ++i) {
        P2 v_diff = sub(V_cp[i+1], V_cp[i]); // V_cp 已经是 P的差分乘以n了
        // 这里需要再乘以 (n-1) = 4
        // 实际上 A_control = n*(n-1)*(P_{i+2} - 2P_{i+1} + P_i)
        // 但我们要用 V_cp 算出的值
        // V_cp 存储的是 (P_{i+1}-P_i)
        // 所以我们算的 v_diff 实际上没有乘以 n，上面 V_cp push的时候我仅仅push了向量
        // 修正一下逻辑：
        
        // 更简单的写法：直接用原始点算
        // A_k = n * (n-1) * (P_{k+2} - 2P_{k+1} + P_k)
        P2 acc_vec = add(sub(P[i+2], mul(P[i+1], 2.0)), P[i]);
        double val = std::sqrt(dot(acc_vec, acc_vec)) * n * (n - 1);
        if(val > max_acc_norm) max_acc_norm = val;
    }

    // 4. 计算需要的比例 Ratio
    // 我们需要: max_vel_norm / (T * ratio) <= v_max
    //           max_acc_norm / (T * ratio)^2 <= a_max
    
    double ratio_v = max_vel_norm / (T * v_max);
    double ratio_a = std::sqrt(max_acc_norm / (T * T * a_max));

    return std::max(1.0, std::max(ratio_v, ratio_a));
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
    if (segments.empty()) return {};

    int K = segments.size();
    int N_poly = 5;         // 5阶贝塞尔
    int N_cp = N_poly + 1;  // 每段 6 个控制点
    int n_vars_per_axis = K * N_cp;
    int n_vars = 2 * n_vars_per_axis;

    // 1. 预计算 Hessian 矩阵 (Velocity & Acceleration)
    static Eigen::MatrixXd H_vel = computeBernsteinHessian(N_poly, 1);
    static Eigen::MatrixXd H_acc = computeBernsteinHessian(N_poly, 2);

    // 2. 时间分配
    std::vector<double> T(K);
    for (int i = 0; i < K; ++i) {
        double dist = getDist(segments[i].start, segments[i].end);
        T[i] = std::max(dist / config.avg_vel, 0.5); 
    }

    OsqpEigen::Solver solver;
    solver.settings()->setVerbosity(false);
    solver.settings()->setWarmStart(true);

    // 3. 构建目标函数矩阵 P
    Eigen::SparseMatrix<double> P(n_vars, n_vars);
    std::vector<Eigen::Triplet<double>> p_triplets;

    for (int k = 0; k < K; ++k) {
        // 积分换元系数:
        // Vel term: integral (p'(t))^2 dt = (1/T) * integral (p'(u))^2 du
        // Acc term: integral (p''(t))^2 dt = (1/T^3) * integral (p''(u))^2 du
        double c_vel = config.w_vel / T[k];
        double c_acc = config.w_acc / std::pow(T[k], 3); // 注意这里是T的三次方

        for (int r = 0; r < N_cp; ++r) {
            for (int c = 0; c < N_cp; ++c) {
                // 组合 Velocity 和 Acceleration 代价
                double val = c_vel * H_vel(r, c) + c_acc * H_acc(r, c);
                
                // X 轴
                int idx_x = k * N_cp + r;
                int idy_x = k * N_cp + c;
                p_triplets.emplace_back(idx_x, idy_x, val);
                
                // Y 轴
                p_triplets.emplace_back(idx_x + n_vars_per_axis, idy_x + n_vars_per_axis, val);
            }
        }
    }
    P.setFromTriplets(p_triplets.begin(), p_triplets.end());

    // 目标函数线性项 q
    // 设置为零向量
    Eigen::VectorXd q = Eigen::VectorXd::Zero(n_vars);

    // 4. 构建约束
    std::vector<Eigen::Triplet<double>> A_triplets;
    std::vector<double> l_vec, u_vec;
    int constraint_idx = 0;

    // 辅助函数 用来添加约束
    auto add_constraint = [&](int row, int col, double value) {
        A_triplets.emplace_back(row, col, value);
    };
    // 辅助函数 用来设置约束范围
    auto set_rng = [&](double lb, double ub) {
        l_vec.push_back(lb); u_vec.push_back(ub); constraint_idx++;
    };

    // --- (A) 起点终点约束 ---
    // Start X/Y
    add_constraint(constraint_idx, 0, 1.0); set_rng(raw_path.front().x, raw_path.front().x); // x0
    add_constraint(constraint_idx, n_vars_per_axis, 1.0); set_rng(raw_path.front().y, raw_path.front().y); // y0
    
    // End X/Y
    int last_idx = (K - 1) * N_cp + 5;
    add_constraint(constraint_idx, last_idx, 1.0); set_rng(raw_path.back().x, raw_path.back().x); // x_end
    add_constraint(constraint_idx, last_idx + n_vars_per_axis, 1.0); set_rng(raw_path.back().y, raw_path.back().y); // y_end

    // --- (B) 连续性约束 (C0, C1, C2) ---
    for (int k = 0; k < K - 1; ++k) {
        int i_end = k * N_cp + 5;
        int i_prev1 = k * N_cp + 4;
        int i_prev2 = k * N_cp + 3;
        
        int j_start = (k + 1) * N_cp + 0;
        int j_next1 = (k + 1) * N_cp + 1;
        int j_next2 = (k + 1) * N_cp + 2;

        double Tk = T[k], Tk1 = T[k+1];

        for (int axis = 0; axis < 2; ++axis) {
            int off = axis * n_vars_per_axis;

            // C0: Pos continuous
            add_constraint(constraint_idx, i_end + off, 1.0);
            add_constraint(constraint_idx, j_start + off, -1.0);
            set_rng(0, 0);

            // C1: Vel continuous -> (P_n - P_{n-1})/Tk = (Q_1 - Q_0)/Tk1
            double s1 = 1.0 / Tk;
            double s2 = 1.0 / Tk1;
            add_constraint(constraint_idx, i_end + off, s1);
            add_constraint(constraint_idx, i_prev1 + off, -s1);
            add_constraint(constraint_idx, j_next1 + off, -s2);
            add_constraint(constraint_idx, j_start + off, s2);
            set_rng(0, 0);

            // C2: Acc continuous -> (P_n - 2P_{n-1} + P_{n-2})/Tk^2 = ...
            double a1 = 1.0 / (Tk*Tk);
            double a2 = 1.0 / (Tk1*Tk1);
            add_constraint(constraint_idx, i_end + off, a1);
            add_constraint(constraint_idx, i_prev1 + off, -2*a1);
            add_constraint(constraint_idx, i_prev2 + off, a1);
            
            add_constraint(constraint_idx, j_next2 + off, -a2);
            add_constraint(constraint_idx, j_next1 + off, 2*a2);
            add_constraint(constraint_idx, j_start + off, -a2);
            set_rng(0, 0);
        }
    }

    // --- (C) 安全走廊约束 (Convex Hull Property) ---
    for (int k = 0; k < K; ++k) {

        // 获取 Corridor 对应的线性约束 Ax + By <= C
        auto [Mat, Ub] = convertPolygonToConstraints(corridors[k]);
        
        for (int i = 0; i < N_cp; ++i) {
            int idx = k * N_cp + i;
            // 获取每个系数 A, B, C
            for (int r = 0; r < Mat.rows(); ++r) {
                double A = Mat(r, 0), B = Mat(r, 1), C = Ub(r);
                add_constraint(constraint_idx, idx, A);                   // A*x
                add_constraint(constraint_idx, idx + n_vars_per_axis, B); // B*y
                set_rng(-OsqpEigen::INFTY, C);
            }
        }
    }

    // 5. 求解
    Eigen::SparseMatrix<double> A(constraint_idx, n_vars);
    A.setFromTriplets(A_triplets.begin(), A_triplets.end());
    
    Eigen::VectorXd l_e(l_vec.size()), u_e(u_vec.size());
    for(size_t i=0; i<l_vec.size(); ++i) { l_e(i)=l_vec[i]; u_e(i)=u_vec[i]; }

    solver.data()->setNumberOfVariables(n_vars);
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
    Trajectory traj;
    for (int k = 0; k < K; ++k) {
        BezierCurve curve;
        curve.duration = T[k];
        for (int i = 0; i < N_cp; ++i) {
            int idx = k * N_cp + i;
            curve.control_points.push_back({
                (int)std::round(sol(idx)), 
                (int)std::round(sol(idx + n_vars_per_axis))
            });
        }
        traj.pieces.push_back(curve);
    }

    // 7. 执行时间放缩 (Time Scaling)
    double global_ratio = 1.0;
    
    // 策略A：全局统一放缩 (保持各段相对时间比例不变，最平滑)
    for(const auto& piece : traj.pieces) {
        double r = computeTimeScalingFactor(piece, config.max_vel, config.max_acc);
        if(r > global_ratio) global_ratio = r;
    }

    // 应用放缩
    if(global_ratio > 1.001) {
        std::cout << "[Optimizer] Time Scaling Triggered: ratio = " << global_ratio << std::endl;
        for(auto& piece : traj.pieces) {
            piece.duration *= global_ratio;
        }
    }

    return traj;
}