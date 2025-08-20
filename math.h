#pragma once

#include <algorithm>
#include <cmath>

#include "fmm.h"
#include "vec3d.h"

#define _USE_MATH_DEFINES

using realVec = std::vector<double>;
using cmplx = std::complex<double>;
using cmplxVec = std::vector<cmplx>;

using pair2i = std::pair<int, int>;
using pair2d = std::pair<double, double>;
using pairSol = std::pair<double, vec3d>;
using solVec = std::vector<pairSol>;

const double PI = std::acos(-1.0);
constexpr cmplx iu(0,1);

pairSol operator+ (const pairSol& sol0, const pairSol& sol1) {
    return pairSol(sol0.first + sol1.first,
        sol0.second + sol1.second);
}

// TODO : Use a concept to enforce that the template type is summable
template <typename T>
std::vector<T> operator+ (const std::vector<T>& zs, const std::vector<T>& ws) {
    std::vector<T> sum;
    for (size_t i = 0; i < zs.size(); ++i)
        sum.push_back(zs[i] + ws[i]);
    return sum;
}

std::array<bool,3> operator> (const vec3d& x, const vec3d& y) {
    std::array<bool,3> bools{ x[0] > y[0], x[1] > y[1], x[2] > y[2] };
    return bools;
}

size_t bools2Idx(const std::array<bool, 3>& x) {
    return x[0] + 2 * x[1] + 4 * x[2];
}

// convert x to reverse binary, then replace each bit as 0 -> -1, 1 -> 1
vec3d idx2pm(const int x) {
    assert(x < 8);
    auto xmod4 = x%4;
    vec3d bits(xmod4%2, xmod4/2, x/4);

    for (auto& bit : bits) 
        bit = (bit == 0 ? -1 : 1);
    
    return bits;
}

// return pow(-1,k)
double pm(const int k) {
    return k % 2 ? -1.0 : 1.0;
}

vec3d fromSph(const vec3d& R) {
    auto r = R[0], th = R[1], ph = R[2];
    return vec3d(
        r * std::sin(th) * std::cos(ph), 
        r * std::sin(th) * std::sin(ph), 
        r * std::cos(th) );
}

vec3d toSph(const vec3d& X) {
    auto x = X[0], y = X[1], z = X[2], r = X.norm();
    assert(r != 0);

    auto toPhi = [](double x, double y) {
        if (x == 0 && y == 0) return 0.0; // pick phi = 0.0
        return std::atan2(y,x);
    };

    return vec3d( r, std::acos(z/r), toPhi(x,y) );
}

//size_t lm2Idx(const int l, const int m) {
//    return l*l + l + m;
//}

//constexpr int constPow(int base, int exp) {
//    return (exp == 0) ? 1 : base * constPow(base,exp-1);
//}

// return \sum_i (coeffs[i] * z^i)
template <typename T>
const T evaluatePoly(std::vector<T> coeffs, const T z) {
    for (ptrdiff_t i = coeffs.size()-2; i >= 0; --i)
        coeffs[i] += coeffs[i+1] * z;
    return coeffs[0];
}

const double fallingFactorial(double x, int k) {
    return k == 0 ? 1 : x * fallingFactorial(x - 1, k - 1);
}

//const uint64_t factorial(int n) {
//    return n == 0 ? 1 : n * factorial(n-1);
//}

const double factorial(double n) {
    return n == 0 ? 1 : n * factorial(n-1);
}

const double coeffYlm(int l, int abs_m) {
    assert(abs_m <= l);
    return
        std::sqrt(factorial(l-abs_m) / static_cast<double>(factorial(l+abs_m))) * // Ylm coeffs
        pm(abs_m) * std::pow(2.0, l); // legendreLM coeffs
}

const cmplx expI(const double arg) {
    return std::exp(iu*arg);
}

const cmplx powI(const uint32_t m) {
    switch (m % 4) {
        case 0: return 1;
        case 1: return iu;
        case 2: return -1.0;
        case 3: return -iu;
    }
}

mat3d rotationR(const pair2d angles) {
    auto [th, ph] = angles;
    return mat3d {
        {  cos(th)*cos(ph),  cos(th)*sin(ph), -sin(th) },
        { -sin(ph),          cos(ph),          0       },
        {  sin(th)*cos(ph),  sin(th)*sin(ph),  cos(th) }
    };
}

matXcd wignerD_l(const pair2d angles, const int l) {
    using namespace std;
    const auto [th, ph] = angles;

    auto sumCoeff = [th, l](int m, int n, int s) {
        int a0 = l+m-s, a1 = n-m+s, a2 = s, a3 = l-n-s;
        return pow(-1.0,n-m+s) * 
            pow(cos(th/2.0), a0+a3) * pow(sin(th/2.0), a1+a2) /
            ( factorial(a0) * factorial(a1) * factorial(a2) * factorial(a3) );
    };

    matXcd mat = matXcd::Zero(2*l+1, 2*l+1);
    for (int n = -l; n <= l; ++n) {
        int n_ = n+l;
        double pm_n = (n < 0) ? pm(n) : 1.0;
        cmplx exp_n = expI(static_cast<double>(-n)*ph);
        for (int m = -l; m <= l; ++m) {
            int m_ = m+l;
            double pm_m = (m < 0) ? pm(m) : 1.0;
            for (int s = max(m-n, 0); s <= min(l+m, l-n); ++s)
                mat(n_, m_) += sumCoeff(m, n, s);

            mat(n_, m_) *= exp_n 
                * pm_n / pm_m
                * sqrt(factorial(l+n)*factorial(l-n)*factorial(l+m)*factorial(l-m));
        }
    }

    return mat;
}