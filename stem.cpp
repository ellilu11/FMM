#include <bitset>
#include <cassert>
#include <iostream>
#include <memory>
#include <ranges>
#include "leaf.h"
#include "math.h"
#include "stem.h"

Stem::Stem(std::vector<cmplx>& pos,
	std::vector<double>& qs,
	const cmplx zk,
	const double L,
	const int lvl,
	const int branchIdx,
	Node* const base)
	: Node(pos, qs, zk, L, lvl, branchIdx, base)
{
	std::vector<std::vector<cmplx>> branchPos(4);
	std::vector<std::vector<double>> branchQs(4);
	// for (auto [ele, elq] : std::views::zip(pos, qs)) {
	for (size_t n = 0; n < pos.size(); ++n) {
		size_t k = cmplx2Idx<bool>(pos[n] > zk);
		assert(k < branchPos.size());

		branchPos[k].push_back(pos[n]);
		branchQs[k].push_back(qs[n]);
	}
	for (size_t k = 0; k < branchPos.size(); ++k) {
		cmplx dzk( pow(-1,k%2+1), pow(-1,k/2+1) );
		dzk *= L / 4.0;

		std::shared_ptr<Node> branch;
		// if branch has at least two particles and is not lvl 0, then further subdivide it
		if (branchPos[k].size() > 1 && lvl-1)
			branch = std::make_shared<Stem>(branchPos[k], branchQs[k], zk+dzk, L/2, lvl-1, k, this);
		else
			branch = std::make_shared<Leaf>(branchPos[k], branchQs[k], zk+dzk, L/2, lvl-1, k, this);

		branches.push_back(branch);
	}
}

void Stem::buildCoeffs(const int P) {
	cmplx b_0, b_k;

	// recursively combine coeffs from branch nodes
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
			auto z0 = branch->getCenter();
			for (size_t l = 0; l < k; ++l)
				b_k += branchCoeffs[l] * pow(z0-zk, k - l) * binom(k - 1, l - 1)
				- branchCoeffs[0] * pow(z0-zk, k) / static_cast<double>(k);
		}
		coeffs.push_back(b_k);
	}
}

std::vector<cmplx> Stem::shiftLocalCoeffs(const std::vector<cmplx>& coeffs, const int P) {
	return localCoeffs;
}

void Stem::buildLocalCoeffs(const int P) {
	localCoeffs = shiftLocalCoeffs( base->getLocalCoeffs(), P ); // add local coeffs to those from base (after shifting to z0)
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

	for (const auto& branch : branches)
		branch->buildLocalCoeffs(P);

}

void Stem::printNode(std::ofstream& f) {
	f << zk.real() << " " << zk.imag() << " " << L << " " << nodeStat << std::endl;
	for (const auto& branch : branches)
		branch->printNode(f);
}

void Stem::iListTest() {
	std::random_device rd;
	std::mt19937 gen(rd());

	std::uniform_int_distribution<> branchIdx(0,3);
	std::shared_ptr<Node> node = std::make_shared<Stem>(*this);

	while (node->isNodeType<Stem>())
	//	// while (node->getLvl() > 3)
		node = node->getBranches(branchIdx(gen));
	node->setNodeStat(3);

	node->setInteractionList();
	auto nbors = node->getInteractionList();

	for (const auto& nbor : nbors)
		nbor->setNodeStat(2);

	std::ofstream posFile, nodeFile;
	posFile.open("positions.txt");
	nodeFile.open("nodes.txt");

	printPos(posFile);
	printNode(nodeFile);
}