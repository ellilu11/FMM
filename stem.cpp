#include "stem.h"
#include "leaf.h"
#include <cassert>
#include <iostream>
#include <memory>

Stem::Stem(std::vector<cmplx>& pos,
	std::vector<double>& qs,
	cmplx z0,
	const double L,
	const int lvl,
	const int branchIdx,
	Node* const root)
	: Node(pos, qs, z0, L, lvl, branchIdx, root)
{
	std::vector<std::vector<cmplx>> branchPos(4);
	std::vector<std::vector<double>> branchQs(4);
	// for (const auto& ele : pos) {
	for (size_t n = 0; n < pos.size(); ++n) {
		size_t k = (pos[n].real() > z0.real()) + 2 * (pos[n].imag() > z0.imag());
		assert(k < 4);
		branchPos[k].push_back(pos[n]);
		branchQs[k].push_back(qs[n]);
	}
	for (size_t k = 0; k < 4; ++k) {
		cmplx dz0 = (pow(-1, k % 2 + 1) + iu * pow(-1, k / 2 + 1)) * L / 4.0;
		std::shared_ptr<Node> branch;

		// if branch has at least two particles and is not lvl 0, then further subdivide it
		if ((branchPos[k]).size() > 1 && lvl-1)
			branch = std::make_shared<Stem>(branchPos[k], branchQs[k], z0+dz0, L/2, lvl-1, k, this);
		else
			branch = std::make_shared<Leaf>(branchPos[k], branchQs[k], z0+dz0, L/2, lvl-1, k, this);

		branches.push_back(branch);
	}
}

void Stem::buildCoeffs(const int P) {
	cmplx b_0, b_k;

	// recursively combine coeffs from child nodes into coeffs for parent node
	// only combine coeffs of nodes on the same lvl!
	for (const auto& branch : branches) {
		branch->buildCoeffs(P);
		auto branchCoeffs = branch->getCoeffs();
		b_0 += branchCoeffs[0];
	}
	coeffs.push_back(b_0);

	for (size_t k = 1; k < P; ++k) {
		for (const auto& branch : branches) {
			auto branchCoeffs = branch->getCoeffs();
			for (size_t l = 0; l < k; ++l)
				b_k += branchCoeffs[l] * pow(branch->getCenter(), k - l) * binom(k - 1, l - 1)
				- branchCoeffs[0] * pow(branch->getCenter(), k) / static_cast<double>(k);
		}
		coeffs.push_back(b_k);
	}
}

void Stem::buildLocalCoeffs(const int P) {

}

void Stem::printNode(std::ofstream& f) {
	f << z0.real() << " " << z0.imag() << " " << L << " " << nborFlag << std::endl;
	for (const auto& branch : branches)
		branch->printNode(f);
}