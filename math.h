#pragma once

#include <cmath>
// #include <Eigen/Dense>
#include "fmm.h"

#define _USE_MATH_DEFINES

using vec3d = std::array<double, 3>;
using vec3dVec = std::vector<vec3d>;

std::array<bool,3> operator> (vec3d& x, vec3d& y) {
    std::array<bool,3> bools(x[0] > y[0], x[1] > y[1], x[2] > y[2]);
    return bools;
}

size_t bools2Idx(std::array<bool,3>& x) {
    return x[0] + 2 * x[1] + 4 * x[2];
}

std::ostream& operator<< (std::ostream& os, const vec3d x) {
    os << x[0] << " " << x[1] << " " << x[2];
    return os;
}

std::istream& operator>>(std::istream& is, vec3d& x) {
    double x0, x1, x2;
    if (is >> x0 >> x1 >> x2)
        x = vec3d(x0, x1, x2);
    return is;
}

//vec3d cartToSph(vec3d& x) {
//    auto ssq = x[0]*x[0] + x[1]*x[1];
//    auto r = std::sqrt(x[0]*x[0] + x[1]*x[1]+ x[2]*x[2] );
//    auto th = std::acos(x[2] / r );
//    auto ph = std::atan(x[0]/ );
//
//    return vec3d(r, th, ph);
//}

const uint64_t fallingFactorial(int n, int k) {
    return n == k ? 1 : n * fallingFactorial(n - 1, k);
}

const uint64_t binom(int n, int k) {
    return fallingFactorial(n, k) / fallingFactorial(n - k, 0); 
}

std::ostream& operator<<(std::ostream& out, const vec3dVec& vec) {
    for (const auto& ele : vec)
        out << ele << ' ';
    out << '\n';
    return out;
}

template <typename T>
bool contains(std::vector<T>& vec, T val) {
    return std::find(vec.begin(), vec.end(), val) != vec.end();
}

// return \sum_i (coeffs[i] * z^i)
// understand why passing coeffs by ref yields larger error
template <typename T>
const T evaluatePoly(std::vector<T> coeffs, const T z) {
    for (ptrdiff_t i = coeffs.size()-2; i >= 0; --i) 
        coeffs[i] += coeffs[i+1] * z;
    return coeffs[0];
}


