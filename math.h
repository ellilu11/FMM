#pragma once

#include <cmath>
#include <complex>

using cmplx = std::complex<double>;
using cmplxVec = std::vector<cmplx>;
constexpr cmplx iu(0,1);

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

cmplxVec operator+= (cmplxVec& zs, const cmplxVec& ws) {
    for (size_t i = 0; i < zs.size(); ++i)
       zs[i] += ws[i];
    return zs;
}

cmplxVec operator+ (cmplxVec& zs, cmplxVec& ws) {
    return zs += ws;
}

//const double factorial(int n) {
//	return n == 0 ? 1 : n * factorial(n-1);
//}

const double fallingFactorial(int n, int k) {
    return n == k ? 1 : n * fallingFactorial(n - 1, k);
}

const double binom(int n, int k) {
    // return fallingFactorial(n,k) / factorial(n - k);
    return fallingFactorial(n, k) / fallingFactorial(n - k, 0);
}

template <typename T>
bool vecContains(std::vector<T>& vec, T val) {
    return std::find(vec.begin(), vec.end(), val) != vec.end();
}

//template <typename T>
//bool ptrContains(std::vector<std::shared_ptr<T>>& vec, std::shared_ptr<T>& ele) {
//    auto it = std::find_if(vec.begin(), vec.end(),
//        [&val_to_find](ele) {
//            return ele && *ele == val_to_find;
//        });
//    return it != vec.end();
//}



