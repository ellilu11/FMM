#pragma once

#include <algorithm>
#include <cmath>

#include "fmm.h"
#include "vec3d.h"

#define _USE_MATH_DEFINES

using realVec = std::vector<double>;
using cmplx = std::complex<double>;
using cmplxVec = std::vector<cmplx>;
const double PI = std::acos(-1.0);
constexpr cmplx iu(0, 1);

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

double pm(const int k) {
    return k % 2 ? -1 : 1;
}

double abs(const vec3d& X) {
    auto [x, y, z] = X;
    return sqrt(x*x + y*y + z*z);
}

vec3d sph2Cart(const vec3d& R) {
    auto [r, th, ph] = R;
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

// return factorial(n) / factorial(n-k)
const int fallingFactorial(int n, int k) {
    return k == 0 ? 1 : n * fallingFactorial(n - 1, k - 1);
}

const double fallingFactorial(double x, int k) {
    return k == 0 ? 1 : x * fallingFactorial(x - 1, k - 1);
}

const int factorial(int n) {
    return fallingFactorial(n, n);
}

const double sphHarmonicCoeff(int l, int absm) {
    assert(absm <= l);
    return sqrt(factorial(l-absm) / factorial(l+absm)) * 
        pm(absm) * pow(2.0, l);
}