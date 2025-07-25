#pragma once

#include <cmath>
#include "fmm.h"

#define _USE_MATH_DEFINES

using cmplx = std::complex<double>;
using realVec = std::vector<double>;
using cmplxVec = std::vector<cmplx>;
constexpr cmplx iu(0, 1);

std::complex<bool> operator> (cmplx z, cmplx w) {
    std::complex<bool> bools(z.real() > w.real(), z.imag() > w.imag());
    return bools;
}

template <typename T>
size_t cmplx2Idx(std::complex<T> z) {
    return z.real() + 2*z.imag();
}

std::ostream& operator<< (std::ostream& os, const cmplx z) {
    os << z.real() << " " << z.imag();
    return os;
}

std::istream& operator>>(std::istream& is, cmplx& z) {
    double real, imag;
    if (is >> real >> imag)
        z = cmplx(real, imag);
    return is;
}

cmplxVec operator+= (cmplxVec& zs, const cmplxVec& ws) {
    for (size_t i = 0; i < zs.size(); ++i)
       zs[i] += ws[i];
    return zs;
}

cmplxVec operator+ (const cmplxVec& zs, const cmplxVec& ws) {
    cmplxVec sum;
    for (size_t i = 0; i < zs.size(); ++i)
        sum.emplace_back( zs[i] + ws[i] );
    return sum;
}

const uint64_t fallingFactorial(int n, int k) {
    return n == k ? 1 : n * fallingFactorial(n - 1, k);
}

const uint64_t binom(int n, int k) {
    // return fallingFactorial(n,k) / factorial(n - k);
    return fallingFactorial(n, k) / fallingFactorial(n - k, 0); 
}

std::ostream& operator<<(std::ostream& out, const cmplxVec& vec) {
    for (const auto& ele : vec)
        out << ele << ' ';
    out << '\n';
    return out;
}

template <typename T>
bool contains(std::vector<T>& vec, T val) {
    return std::find(vec.begin(), vec.end(), val) != vec.end();
}

// returns \sum_i coeffs[i] * z^i
template <typename T>
const T evaluatePoly(std::vector<T> coeffs, const T z) {
    for (ptrdiff_t i = coeffs.size()-2; i >= 0; --i) 
        coeffs[i] += coeffs[i+1] * z;
    return coeffs[0];
}


