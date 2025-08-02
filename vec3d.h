#pragma once

#include <Eigen/Dense>

// using vec3d = std::array<double, 3>; 
using vec3d = Eigen::Vector3d;
using vec3dVec = std::vector<vec3d>;

using vec3cd = Eigen::Vector3cd;
using vec3cdVec = std::vector<vec3cd>;

using vecXcd = Eigen::VectorXcd;
using vecXcdVec = std::vector<vecXcd>;

using matXcd = Eigen::MatrixXcd;
using matXcdVec = std::vector<matXcd>;

const vec3d zeroVec = vec3d::Zero();
// constexpr vec3d zeroVec{ 0,0,0 };

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

vec3dVec operator+ (const vec3dVec& Xs, const vec3dVec& Ys) {
    assert(Xs.size() == Ys.size());
    vec3dVec sum;
    for (size_t i = 0; i < Xs.size(); ++i)
        sum.push_back(Xs[i] + Ys[i]);
    return sum;
}

std::ostream& operator<<(std::ostream& out, const vec3dVec& vec) {
    for (const auto& ele : vec)
        out << ele << ' ';
    out << '\n';
    return out;
}


