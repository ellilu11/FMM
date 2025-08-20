#pragma once

#include <Eigen/Dense>

// using vec3d = std::array<double, 3>; 
using vec3d = Eigen::Vector3d;

using vec3cd = Eigen::Vector3cd;

using vecXcd = Eigen::VectorXcd;

using mat3d = Eigen::Matrix3d;

using matXcd = Eigen::MatrixXcd;

const vec3d zeroVec = vec3d::Zero();

std::ostream& operator<< (std::ostream& os, const vec3d& X) {
    os << X[0] << " " << X[1] << " " << X[2];
    return os;
}

std::istream& operator>>(std::istream& is, vec3d& X) {
    double x, y, z;
    if (is >> x >> y >> z)
        X = vec3d{ x, y, z };
    return is;
}

/*
vec3d operator+ (const vec3d& X0, const vec3d& X1) {
    auto [x0, y0, z0] = X0;
    auto [x1, y1, z1] = X1;
    return vec3d{ x0+x1, y0+y1, z0+z1 };
}

vec3d operator- (const vec3d& X0, const vec3d& X1) {
    auto [x0, y0, z0] = X0;
    auto [x1, y1, z1] = X1;
    return vec3d{ x0-x1, y0-y1, z0-z1 };
}

vec3d operator* (const double a, const vec3d& X) {
    auto [x, y, z] = X;
    return vec3d{ a*x, a*y, a*z };
}
*/

//std::vector<vec3d> operator+ (const std::vector<vec3d>& Xs, const std::vector<vec3d>& Ys) {
//    assert(Xs.size() == Ys.size());
//    std::vector<vec3d> sum;
//    for (size_t i = 0; i < Xs.size(); ++i)
//        sum.push_back(Xs[i] + Ys[i]);
//    return sum;
//}

std::ostream& operator<<(std::ostream& out, const std::vector<vec3d>& vec) {
    for (const auto& ele : vec)
        out << ele << ' ';
    out << '\n';
    return out;
}


