#pragma once
#include <vector>
#include <cmath>
#include <iostream>
#include "grid.hpp"

// 5阶贝塞尔曲线段
// P(t) = (1-t)^5*P0 + 5(1-t)^4*t*P1 + 10(1-t)^3*t^2*P2 + 10(1-t)^2*t^3*P3 + 5(1-t)*t^4*P4 + t^5*P5, t in [0,1]
class BezierCurve {
public:
    std::vector<P2> control_points; // 6个控制点
    double duration;                   // 该段时长 T

    BezierCurve() = default;

    // 计算 t 时刻的位置 (t in [0, duration])
    std::pair<double, double> evaluate(double t) const {
        if (t <= 0) t = 0;
        if (t >= duration) t = duration;
        double u = t / duration;
        int n = 5;
        // 5阶系数: 1, 5, 10, 10, 5, 1
        static const int C[6] = {1, 5, 10, 10, 5, 1};
        double x = 0.0, y = 0.0;
        for (int i = 0; i <= n; ++i) {
            double b = C[i] * std::pow(1 - u, n - i) * std::pow(u, i);
            x += b * control_points[i].x;
            y += b * control_points[i].y;
        }
        return {x, y};
    }
};

struct Trajectory {
    std::vector<BezierCurve> pieces;

    // 计算 t 时刻的位置
    std::pair<double, double> at(double t) const {
        if (pieces.empty()) return {0,0};
        for (const auto& curve : pieces) {
            if (t <= curve.duration) return curve.evaluate(t);
            t -= curve.duration;
        }
        return pieces.back().evaluate(pieces.back().duration);
    }
    double getTotalDuration() const {
        double sum = 0;
        for(auto& p : pieces) sum += p.duration;
        return sum;
    }
};