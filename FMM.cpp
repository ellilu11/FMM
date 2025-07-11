// FMM.cpp : Defines the entry point for the application.

#include <complex>
#include <fstream>
#include <memory>
#include <random>
#include <vector>
#include "node.h"
#include "stem.h"
#include "stem.cpp"
#include "leaf.h"
#include "leaf.cpp"

using namespace std;

int main()
{
	constexpr int N = 1000;
	constexpr double Q = 1.0;
	constexpr double eps = 1.0E-6;
	const double L = 10.0;

	const int Nlvl = ceil(log(N) / log(4.0));
	const int P = ceil(-log(eps) / log(2));

	random_device rd;
	mt19937 gen(rd());

	uniform_real_distribution<double> real(-L / 2, L / 2);
	uniform_real_distribution<double> imag(-L / 2, L / 2);

	vector<cmplx> pos;
	vector<double> qs;

	for (int n = 0; n < N; ++n) {
		cmplx z(real(gen), imag(gen));
		pos.push_back(z);
		qs.push_back(Q);
	}
	shared_ptr<Node> master;
	if (N > 1 && Nlvl)
		master = make_shared<Stem>(pos, qs, 0, L, Nlvl, P);
	else
		master = make_shared<Leaf>(pos, qs, 0, L, Nlvl, P);

	// master->buildCoeffs();

	ofstream posFile, nodeFile;
	posFile.open("positions.txt");
	nodeFile.open("nodes.txt");
	
	master->printPos(posFile);
	master->printNode(nodeFile);

	return 0;
}