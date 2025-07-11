#include <complex>
#include <random>
#include <iostream>
#include "node.h"

using namespace std;

int main()
{
	constexpr int N = 50;
	constexpr double L = 10.0;
	constexpr double Q = 1.0;
	constexpr double eps = 1.0E-6;

	const int Nlvl = ceil(log(N) / log(4.0));
	const int P = ceil(-log(eps) / log(2));

	random_device rd;
	mt19937 gen(rd());

	uniform_real_distribution<double> real(-L/2,L/2);
	uniform_real_distribution<double> imag(-L/2,L/2);

	vector<complex<double>> &pos;
	vector<double> qs;

	for (const auto& ele : pos) {
		complex<double> z(real(gen), imag(gen));
		pos.push_back(z);
		qs.push_back(Q);
	}

	shared_ptr<Node> masterNode = make_shared<Node>(pos,complex<double>(0,0),L,Nlvl);
	masterNode->buildCoeffs();

	return 0;
}