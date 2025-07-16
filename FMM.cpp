#define _USE_MATH_DEFINES
#include <cmath>
#include <fstream>
#include <random>
#include "node.h"
#include "node.cpp"
#include "stem.cpp"
#include "leaf.cpp"


using namespace std;

//namespace PARAM {
//
//}

const cmplx evaluateFfieldAnl(const vector<cmplx>& pos, vector<double>& qs, const cmplx z) {
	assert(pos.size() == qs.size());

	cmplx phi;
	for (size_t n = 0; n < pos.size(); ++n)
		phi += qs[n] * std::log(z - pos[n]);
	return phi;
}

//void multipoleTest(shared_ptr<Node> root, const int Nobs, const double L, const int P) {
//
//}

int main()
{
	// Define computational domain
	constexpr double L = 10.0;
	constexpr double eps = 1.0E-6;
	const int P = ceil(-log(eps) / log(2)); // # terms in multipole expansion

	random_device rd;
	mt19937 gen(rd());

	uniform_real_distribution<double> real(-L / 2, L / 2);
	uniform_real_distribution<double> imag(-L / 2, L / 2);

	//normal_distribution<double> u(0, 0.2 * L);
	//uniform_real_distribution<double> th(0, 2.0 * M_PI);

	// Populate domain with particles
	constexpr int N = 1000;
	constexpr double Q = 1.0;
	const int Nlvl = ceil(log(N) / log(4.0)); // # lvls of refinement

	vector<cmplx> pos;
	vector<double> qs;

	for (int n = 0; n < N; ++n) {
		cmplx z(real(gen), imag(gen));
		//const auto R = std::abs(u(gen));
		//cmplx z(R * std::cos(th(gen)), R * std::sin(th(gen)));
		pos.push_back(z);
		qs.push_back(Q);
	}

	// Refine domain
	shared_ptr<Node> root;
	if (N > 1)
		root = make_shared<Stem>(pos, qs, 0, L, Nlvl, 0, nullptr);
	else
		root = make_shared<Leaf>(pos, qs, 0, L, Nlvl, 0, nullptr);

	// Upward pass: Aggregate multipole expansions
	root->buildCoeffs(P);

	// Tests (move later into test cpp files)
	root->iListTest(); // interaction list finding test

	// Farfield test
	constexpr int Nobs = 1000;
	// multipoleTest(root, Nobs, L, P);

	const double c = 10.0;
	const double R = c * L;

	ofstream obsFile, phiFile;
	obsFile.open("observers.txt");
	phiFile.open("farfield.txt");

	for (int n = 0; n < Nobs; ++n) {
		const double th = 2 * M_PI * static_cast<double>(n) / static_cast<double>(Nobs);
		const double x = R * cos(th);
		const double y = R * sin(th);
		const cmplx z(x, y);

		const auto phi = root->evaluateFfield(z, P);
		const auto phiAnl = evaluateFfieldAnl(pos, qs, z);
		obsFile << x << " " << y << std::endl;
		phiFile << phi.real() << " " << phi.imag() << " " << phiAnl.real() << " " << phiAnl.imag() << std::endl;
	}

	return 0;
}