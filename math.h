#ifndef MATH_H
#define MATH_H

#include <cmath>
#include <complex>

typedef std::complex<double> cmplx;
constexpr cmplx iu(0,1);

double factorial(int n) {
	return n == 0 ? 1 : n * factorial(n - 1);
}

//double fallingFactorial(n, k) {
//
//}

double binom(int n, int k) {
	return factorial(n) / (factorial(k) * factorial(n - k));
}

#endif