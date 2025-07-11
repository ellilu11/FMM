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
	const int P)
	: Node(pos, qs, z0, L, lvl, P)
{
	std::vector<std::vector<cmplx>> branchPos(4);
	std::vector<std::vector<double>> branchQs(4);
	// for (const auto& ele : pos) {
	for (int n = 0; n < pos.size(); ++n) {
		int k = (pos[n].real() > z0.real()) + 2 * (pos[n].imag() > z0.imag());
		assert(k < 4);
		branchPos[k].push_back(pos[n]);
		branchQs[k].push_back(qs[n]);
	}
	for (int k = 0; k < 4; ++k) {
		cmplx dz0 = (pow(-1, k % 2 + 1) + iu * pow(-1, k / 2 + 1)) * L / 4.0;
		std::shared_ptr<Node> branch;

		// if child has at least two particles and is not lvl 0, then subdivide it
		if (pos.size() > 1 && lvl)
			branch = std::make_shared<Stem>(branchPos[k], branchQs[k], z0 + dz0, L / 2, lvl - 1, P);
		else
			branch = std::make_shared<Leaf>(branchPos[k], branchQs[k], z0 + dz0, L / 2, lvl - 1, P);

		branches.push_back(branch);
	}
}

void Stem::buildCoeffs() {
	// recursively combine coeffs from child nodes into coeffs for parent node
	// only combine coeffs of nodes on the same lvl!
	for (const auto& branch : branches) {
		branch->buildCoeffs();
		auto branchCoeffs = branch->getCoeffs();
		coeffs[0] += branchCoeffs[0];
		for (int k = 1; k < P; ++k)
			for (int l = 0; l < k; ++l)
				coeffs[k] += branchCoeffs[l] * pow(branch->getCenter(), k - l) * binom(k - 1, l - 1)
				- branchCoeffs[0] * pow(branch->getCenter(), k) / static_cast<double>(k);
	}
}

void Stem::printNode(std::ofstream& f) {
	f << z0.real() << " " << z0.imag() << " " << L << std::endl;
	for (const auto& branch : branches)
		branch->printNode(f);
}