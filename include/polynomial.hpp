#pragma once
#include <vector>
#include <cmath>
#include <iostream>
#include "grid.hpp"


class PolynomialCurve {
public:
    int degree_ = 3;
    std::vector<double> polynomial_x_coefficients_;
    std::vector<double> polynomial_y_coefficients_;
    double duration_;

    PolynomialCurve(int degree = 3) 
        : degree_(degree), 
          duration_(0.0) 
        {
            polynomial_x_coefficients_.reserve(degree_ + 1);
            polynomial_y_coefficients_.reserve(degree_ + 1);
        }

    std::pair<double, double> evaluate(double t) const {
        double x = 0.0;
        double y = 0.0;

        // 归一化处理
        if (duration_ < 1e-4) return {0,0}; 
        if (t < 0) t = 0;
        if (t > duration_) t = duration_;
        
        double s = t / duration_; // s in [0, 1]

        // x = c0 + s*(c1 + s*(c2 + s*c3))
        for (int i = polynomial_x_coefficients_.size() - 1; i >= 0; --i) {
            x = polynomial_x_coefficients_[i] + s * x;
            y = polynomial_y_coefficients_[i] + s * y;
        }
        
        return {x, y};
    }
};

struct Trajectory {
    std::vector<PolynomialCurve> pieces;

    // 计算 t 时刻的位置
    std::pair<double, double> at(double t) const {
        if (pieces.empty()) return {0,0};
        for (const auto& curve : pieces) {
            if (t <= curve.duration_) {
                return curve.evaluate(t);
            }
            t -= curve.duration_;
        }
        return pieces.back().evaluate(pieces.back().duration_);
    }

    double getTotalDuration() const {
        double sum = 0;
        for(auto& p : pieces) sum += p.duration_;
        return sum;
    }
};