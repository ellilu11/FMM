#ifndef MATH_H
#define MATH_H

#include <cmath>
#include <complex>

using cmplx = std::complex<double>;
constexpr cmplx iu(0,1);

std::complex<bool> operator> (cmplx z, cmplx w) {
	std::complex<bool> out(z.real() > w.real(), z.imag() > w.imag());
	return out;
}

size_t cmplx2Idx(std::complex<bool> z) {
	return z.real() + 2*z.imag();
}

double factorial(int n) {
	return n == 0 ? 1 : n * factorial(n-1);
}

double fallingFactorial(int n, int k) {
	return n == k ? 1 : n * fallingFactorial(n - 1, k);
}

double binom(int n, int k) {
	return fallingFactorial(n,k) / factorial(n - k);
}

#endif