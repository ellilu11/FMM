#pragma once

#include <algorithm>
#include <cmath>
// #include <Eigen/Dense>
// #include <valarray>
#include "fmm.h"

#define _USE_MATH_DEFINES

using realVec = std::vector<double>;
using cmplx = std::complex<double>;
using cmplxVec = std::vector<cmplx>;
using vec3d = std::array<double, 3>; // = Eigen::Vector3d;
using vec3dVec = std::vector<vec3d>;

constexpr cmplx iu(0, 1);
constexpr vec3d zeroVec{ 0,0,0 }; // Eigen::Vector3d::Zero();

std::array<bool,3> operator> (const vec3d& x, const vec3d& y) {
    std::array<bool, 3> bools{ x[0] > y[0], x[1] > y[1], x[2] > y[2] };
    return bools;
}

size_t bools2Idx(const std::array<bool,3>& x) {
    return x[0] + 2 * x[1] + 4 * x[2];
}

std::ostream& operator<< (std::ostream& os, const vec3d& x) {
    os << x[0] << " " << x[1] << " " << x[2];
    return os;
}

std::istream& operator>>(std::istream& is, vec3d& X) {
    double x, y, z;
    if (is >> x >> y >> z)
        X = vec3d{ x, y, z };
    return is;
}

vec3d operator+ (const vec3d& X0, const vec3d& X1) {
    auto [x0, y0, z0] = X0;
    auto [x1, y1, z1] = X1;
    return vec3d{ x0+x1, y0+y1, z0+z1 };
}

vec3d operator* (const double a, const vec3d& X) {
    auto [x, y, z] = X;
    return vec3d{ a*x, a*y, a*z };
}

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

vec3d idx2pm(const int x) {
    assert(x < 8);
    auto xmod4 = x%4;
    vec3d res{ xmod4%2, xmod4/2, x/4 };

    //for_each(res.begin(), res.end(), [](double& x) {
    //    return pow(-1.0, x+1);
    //    } 
    //);
    for (auto& x : res) x = pow(-1.0, x+1);
    return res;
}

double abs(const vec3d& X) {
    auto [x, y, z] = X;
    return sqrt(x*x + y*y + z*z);
}

vec3d sph2Cart(const vec3d& X) {
    auto [r, th, ph] = X;
    return vec3d{ r* sin(th)* cos(ph), r* sin(th)* sin(ph), r* cos(th) };
}

vec3d cart2Sph(const vec3d& X) {
    auto [x, y, z] = X;
    auto r = abs(X);
    return vec3d{ r, acos(z/r), atan(y/x) }; // care with atan
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

const uint64_t fallingFactorial(int n, int k) {
    return n == k ? 1 : n * fallingFactorial(n - 1, k);
}

const uint64_t factorial(int n) {
    return fallingFactorial(n, n);
}

const double sphHarmonicCoeff(int l, int m) {
    auto absm = abs(m);
    assert(absm <= l);
    return sqrt(factorial(l-absm) / factorial(l+absm));
}

const double legendreLM(double x, int l, int m) {
    return pow(-1, m) * pow(1-x*x, static_cast<double>(m)/2.0); // *dmLegendreL
}
