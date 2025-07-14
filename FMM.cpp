// FMM.cpp : Defines the entry point for the application.

#include <complex>
#include <fstream>
#include <random>
#include <vector>
#include "node.h"
#include "node.cpp"
#include "stem.cpp"
#include "leaf.cpp"

using namespace std;

int main()
{
	// Define computational domain
	const double L = 10.0;
	constexpr double eps = 1.0E-6;
	const int P = ceil(-log(eps) / log(2)); // # terms in multipole expansion

	random_device rd;
	mt19937 gen(rd());

	uniform_real_distribution<double> real(-L / 2, L / 2);
	uniform_real_distribution<double> imag(-L / 2, L / 2);

	// Populate particles
	constexpr int N = 1000;
	constexpr double Q = 1.0;
	const int Nlvl = ceil(log(N) / log(4.0)); // # lvls of refinement

	vector<cmplx> pos;
	vector<double> qs;

	for (int n = 0; n < N; ++n) {
		cmplx z(real(gen), imag(gen));
		pos.push_back(z);
		qs.push_back(Q);
	}

	// Refine computational domain
	shared_ptr<Node> root;
	if (N > 1 && Nlvl)
		root = make_shared<Stem>(pos, qs, 0, L, Nlvl, 0, nullptr);
	else
		root = make_shared<Leaf>(pos, qs, 0, L, Nlvl, 0, nullptr);

	root->buildCoeffs(P);

	// Near neighbor test
	uniform_int_distribution<> branchIdx(0, 3);
	shared_ptr<Node> node = root;

	// while (typeid(*node) == typeid(Stem))
	while (node->getLvl() > 1)
		node = node->getBranches(branchIdx(gen));
	node->setNborFlag(2);

	auto nbors = node->getNearNeighbors();
	for (const auto& nbor : nbors)
		nbor->setNborFlag(1);

	ofstream posFile, nodeFile;
	posFile.open("positions.txt");
	nodeFile.open("nodes.txt");

	root->printPos(posFile);
	root->printNode(nodeFile);

	return 0;
}