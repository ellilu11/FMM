#include <cmath>

double factorial(int n) {
	return n == 0 ? 1 : n * factorial(n - 1);
}

//double fallingFactorial(n, k) {
//
//}

double binom(int n, int k) {
	return factorial(n) / (factorial(k) * factorial(n - k));
}