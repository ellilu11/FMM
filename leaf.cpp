#include "leaf.h"
#include <iostream>

Leaf::Leaf(std::vector<cmplx>& pos,
	std::vector<double>& qs,
	cmplx z0,
	const double L,
	const int lvl,
	const int P)
	: Node(pos, qs, z0, L, lvl, P)
{
}

void Leaf::buildCoeffs() {
	for (int i = 0; i < pos.size(); ++i) {
		coeffs[0] += qs[i];
		for (int k = 1; k < P; ++k)
			coeffs[k] -= qs[i] * pow(pos[i] - z0, k) / static_cast<double>(k); // pos[i]-z0 or pos[i] ?
	}
}

void Leaf::printNode(std::ofstream& f) {
	f << z0.real() << " " << z0.imag() << " " << L << std::endl;
}