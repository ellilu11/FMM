#include "leaf.h"
#include <iostream>

Leaf::Leaf(std::vector<cmplx>& pos,
	std::vector<double>& qs,
	cmplx z0,
	const double L,
	const int lvl,
	const int branchIdx,
	Node* const root)
	: Node(pos, qs, z0, L, lvl, branchIdx, root)
{
}

void Leaf::buildCoeffs(const int P) {

	for (int k = 0; k < P; ++k) {
		cmplx a_k;
		for (size_t i = 0; i < pos.size(); ++i)
			a_k += qs[i] * 
				k == 0 ? 1 : -pow(pos[i]-z0, k) / static_cast<double>(k); // pos[i]-z0 or pos[i] ?
		coeffs.push_back(a_k);
	}

}

void Leaf::buildLocalCoeffs(const int P) {

}

void Leaf::printNode(std::ofstream& f) {
	f << z0.real() << " " << z0.imag() << " " << L << " " << nborFlag << std::endl;
}