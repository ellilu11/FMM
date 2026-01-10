#pragma once

#include <Eigen/Dense>

using realVec = std::vector<double>;
using cmplx = std::complex<double>;
using cmplxVec = std::vector<cmplx>;

using pair2i = std::pair<int, int>;
using pair2d = std::pair<double, double>;

using vec3d = Eigen::Vector3d;
using vec3cd = Eigen::Vector3cd;
using vecXcd = Eigen::VectorXcd;

using mat3d = Eigen::Matrix3d;
using matXcd = Eigen::MatrixXcd;

constexpr cmplx iu(0, 1);
const double PI = std::acos(-1.0);
const vec3d zeroVec = vec3d::Zero();

template <typename T>
std::vector<T> operator+ (const std::vector<T>& zs, const std::vector<T>& ws) {
    std::vector<T> sum;
    for (size_t i = 0; i < zs.size(); ++i)
        sum.push_back(zs[i] + ws[i]);
    return sum;
}

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

std::array<bool, 3> operator> (const vec3d& x, const vec3d& y) {
    std::array<bool, 3> bools{ x[0] > y[0], x[1] > y[1], x[2] > y[2] };
    return bools;
}