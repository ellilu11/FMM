#pragma once

#include <algorithm>
#include <cmath>

#include "fmm.h"
#include "vec3d.h"

#define _USE_MATH_DEFINES

using realVec = std::vector<double>;
using cmplx = std::complex<double>;
using cmplxVec = std::vector<cmplx>;
using pair2d = std::pair<double, double>;

const double PI = std::acos(-1.0);
const double TAU = 2.0 * PI;
constexpr cmplx iu(0, 1);

std::array<bool, 3> operator> (const vec3d& x, const vec3d& y) {
    std::array<bool, 3> bools{ x[0] > y[0], x[1] > y[1], x[2] > y[2] };
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
            // throw std::runtime_error("Azimuthal angle undefined");
        return std::atan2(y,x);
    };

    return vec3d{ r, std::acos(z/r), toPhi(x,y) };
}

//size_t lm2Idx(const int l, const int m) {
//    return l*l + l + m;
//}

constexpr int constPow(int base, int exp) {
    return (exp == 0) ? 1 : base * constPow(base,exp-1);
}

template <typename T>
bool contains(std::vector<T>& vec, T val) {
    return std::find(vec.begin(), vec.end(), val) != vec.end();
}

// return \sum_i (coeffs[i] * z^i)
template <typename T>
const T evaluatePoly(std::vector<T> coeffs, const T z) {
    for (ptrdiff_t i = coeffs.size()-2; i >= 0; --i)
        coeffs[i] += coeffs[i+1] * z;
    return coeffs[0];
}

//const uint64_t fallingFactorial(int n, int k) {
//    return k == 0 ? 1 : n * fallingFactorial(n - 1, k - 1);
//}

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
    return // std::sqrt(2*l+1) * 
        std::sqrt(factorial(l-abs_m) / static_cast<double>(factorial(l+abs_m))) * // Ylm coeffs
        pm(abs_m) * std::pow(2.0, l) ; // legendreLM coeffs
}

const cmplx expI(const double arg) {
    return std::exp(iu*arg);
}

//const double sumCoeff(const double th, const int l, int m, int n, int s) {
//    int a0 = l+m-s, a1 = n-m+s, a2 = s, a3 = l-n-s;
//    // std::cout << factorial(a0) << ' ' << factorial(a1) << ' ' << factorial(a2) << ' ' << factorial(a3) << '\n';
//    return pow(-1.0, n-m+s) *
//        pow(cos(th/2.0), a0+a3) * pow(sin(th/2.0), a1+a2) /
//        static_cast<double>(factorial(a0) * factorial(a1) * factorial(a2) * factorial(a3));
//}

Eigen::MatrixXcd wignerD_l(const pair2d& angles, const int l) {
    using namespace std;
    auto [th, ph] = angles;

    auto sumCoeff = [th, l](int m, int n, int s) {
        int a0 = l+m-s, a1 = n-m+s, a2 = s, a3 = l-n-s;
        return pow(-1.0,n-m+s) * 
            pow(cos(th/2.0), a0+a3) * pow(sin(th/2.0), a1+a2) /
            ( factorial(a0) * factorial(a1) * factorial(a2) * factorial(a3) );
    };

    /* double a0 = cos(th/2.0), a1 = sin(th/2.0)*sin(ph), a2 = -sin(th/2.0)*cos(ph);
    pair<double, cmplx> xi(cos(th/2.0), sin(th/2.0) * expI(-ph));
    auto sumCoeff = [xi, l](int m, int n, int s) {
        double alpha0 = l+m-s, alpha1 = n-m+s, alpha2 = s, alpha3 = l-n-s;
        return pow(xi.first, alpha0) * pow(xi.second, alpha1) *
            pow(-conj(xi.second), alpha2) * pow(xi.first, alpha3) /
            ( factorial(alpha0) * factorial(alpha1) * factorial(alpha2) * factorial(alpha3) );
        //return pow(a0, alpha0) * pow(-iu*a1-a2, alpha1) *
        //    pow(-iu*a1+a2, alpha2) * pow(a0, alpha3) /
        //    ( factorial(alpha0) * factorial(alpha1) * factorial(alpha2) * factorial(alpha3) );
    };*/

    Eigen::MatrixXcd mat = Eigen::MatrixXcd::Zero(2*l+1, 2*l+1);
    for (int n = -l; n <= l; ++n) {
        int _n = n+l;
        cmplx exp_n = expI(static_cast<double>(-n)*ph);
        auto pm_n = (n < 0) ? pm(n) : 1.0;
        for (int m = -l; m <= l; ++m) {
            int _m = m+l;
            auto pm_m = (m < 0) ? pm(m) : 1.0;
            for (int s = max(m-n, 0); s <= min(l+m, l-n); ++s)
                mat(_n, _m) += sumCoeff(m, n, s);

            mat(_n, _m) *= exp_n * pm_n / pm_m *
                sqrt(factorial(l+n)*factorial(l-n)*factorial(l+m)*factorial(l-m));
        }
    }

    return mat;
}

/*Eigen::MatrixXcd expI_l(const double ph, const int l) {
    Eigen::MatrixXcd mat = Eigen::MatrixXcd::Zero(2*l+1, 2*l+1);

    for (int m = -l; m <= l; ++m)
        mat(m+l,m+l) = expI(static_cast<double>(-m)*ph);

    return mat;
}*/