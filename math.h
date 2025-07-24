#pragma once

#include <cmath>
#include <complex>

using namespace std;

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

std::ostream& operator<< (std::ostream& out, const cmplx z) {
    out << z.real() << " " << z.imag();
    return out;
}

std::istream& operator>>(std::istream& in, cmplx& z) {
    double real, imag;
    if (in >> real >> imag)
        z = cmplx(real, imag);
    return in;
}

cmplxVec operator+= (cmplxVec& zs, const cmplxVec& ws) {
    for (size_t i = 0; i < zs.size(); ++i)
       zs[i] += ws[i];
    return zs;
}

cmplxVec operator+ (cmplxVec& zs, cmplxVec& ws) {
    return zs += ws;
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





