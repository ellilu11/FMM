#include <bitset>
#include <cassert>
#include <iostream>
#include <memory>
#include "leaf.h"
#include "math.h"
#include "stem.h"

Stem::Stem(std::vector<cmplx>& pos,
	std::vector<double>& qs,
	const cmplx z0,
	const double L,
	const int lvl,
	const int branchIdx,
	Node* const base)
	: Node(pos, qs, z0, L, lvl, branchIdx, base)
{
	std::vector<std::vector<cmplx>> branchPos(4);
	std::vector<std::vector<double>> branchQs(4);
	// for (const auto& ele : pos) {
	for (size_t n = 0; n < pos.size(); ++n) {
		size_t k = cmplx2Idx(pos[n] > z0);
		assert(k < 4);

		branchPos[k].push_back(pos[n]);
		branchQs[k].push_back(qs[n]);
	}
	for (size_t k = 0; k < branchPos.size(); ++k) {
		cmplx dz0( pow(-1,k%2+1), pow(-1,k/2+1) );
		dz0 *= L / 4.0; // center of branch is shifted by dz0 from center of node

		std::shared_ptr<Node> branch;
		// if branch has at least two particles and is not lvl 0, then further subdivide it
		if (branchPos[k].size() > 1 && lvl-1)
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
	setInteractionList(); // place in constructor instead if possible

	cmplx b_0, b_k;

	for (const auto& iNode : iList) {
		cmplx z0 = iNode->getCenter();
		b_0 += iNode->getCoeffs()[0] * std::log(-z0);
		for (size_t l = 1; l < P; ++l)
			b_0 += iNode->getCoeffs()[l] / pow(-z0, l);
	}
	localCoeffs.push_back(b_0);

	for (size_t k = 1; k < P; ++k) {
		for (const auto& iNode : iList) {
			cmplx z0 = iNode->getCenter();
			b_k -= iNode->getCoeffs()[0] / (static_cast<double>(k) * pow(z0, k));
			for (size_t l = 1; l < P; ++l)
				b_k += iNode->getCoeffs()[l] / pow(-z0, l) / pow (z0, k) * binom(k+l-1,l-1);
		}
		localCoeffs.push_back(b_k);
	}

}

void Stem::printNode(std::ofstream& f) {
	f << z0.real() << " " << z0.imag() << " " << L << " " << nodeStat << std::endl;
	for (const auto& branch : branches)
		branch->printNode(f);
}