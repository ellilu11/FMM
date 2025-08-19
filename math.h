#pragma once

#include <cmath>
#include "fmm.h"

#define _USE_MATH_DEFINES

using cmplx = std::complex<double>;
using cmplxVec = std::vector<cmplx>;

std::complex<bool> operator> (cmplx z, cmplx w) {
    std::complex<bool> bools(z.real() > w.real(), z.imag() > w.imag());
    return bools;
}

size_t bools2Idx(std::complex<bool> z) {
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
        sum.push_back( zs[i] + ws[i] );
    return sum;
}

// return factorial(n) / factorial(n-k)
const uint64_t fallingFactorial(int n, int k) {
    return k == 0 ? 1 : n * fallingFactorial(n - 1, k - 1);
}

const uint64_t binom(int n, int k) {
    return fallingFactorial(n, k) / fallingFactorial(k,k); 
}

std::ostream& operator<<(std::ostream& out, const cmplxVec& vec) {
    for (const auto& ele : vec)
        out << ele << ' ';
    out << '\n';
    return out;
}

// return \sum_i (coeffs[i] * z^i)
// understand why passing coeffs by ref yields larger error
template <typename T>
const T evaluatePoly(std::vector<T> coeffs, const T z) {
    for (ptrdiff_t i = coeffs.size()-2; i >= 0; --i) 
        coeffs[i] += coeffs[i+1] * z;
    return coeffs[0];
}


