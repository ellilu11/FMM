#include <bitset>
#include <cassert>
#include <iostream>
#include <memory>
#include <ranges>
#include "leaf.h"
#include "math.h"
#include "stem.h"

Stem::Stem(cmplxVec& psn,
	std::vector<double>& qs,
	const cmplx zk,
	const double L,
	const int lvl,
	const int branchIdx,
	Stem* const base)
	: Node(psn, qs, zk, L, lvl, branchIdx, base)
{
	std::vector<cmplxVec> branchpsn(4);
	std::vector<std::vector<double>> branchQs(4);
	// for (auto [ele, elq] : std::views::zip(psn, qs)) {
	for (size_t n = 0; n < psn.size(); ++n) {
		size_t k = cmplx2Idx<bool>(psn[n] > zk);
		assert(k < branchpsn.size());

		branchpsn[k].push_back(psn[n]);
		branchQs[k].push_back(qs[n]);
	}
	for (size_t k = 0; k < branchpsn.size(); ++k) {
		cmplx dzk( pow(-1,k%2+1), pow(-1,k/2+1) );
		dzk *= L / 4.0;

		std::shared_ptr<Node> branch;
		// if branch has at least two particles and is not lvl 0, then further subdivide it
		if (branchpsn[k].size() > 1 && lvl-1)
			branch = std::make_shared<Stem>(branchpsn[k], branchQs[k], zk+dzk, L/2, lvl-1, k, this);
		else
			branch = std::make_shared<Leaf>(branchpsn[k], branchQs[k], zk+dzk, L/2, lvl-1, k, this);

		branches.push_back(branch);
	}
}

void Stem::buildCoeffs(const int P) {
	cmplx b_0, b_k;

    // std::vector<cmplxVec> branchCoeffss;

	// recursively combine coeffs from branch nodes
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

void Stem::buildLocalCoeffs(const int P) {
	
	buildInteractionList();
	cmplx b_0, b_k;

	for (const auto& iNode : iList) {
		cmplx z0 = iNode->getCenter();
		b_0 += iNode->getCoeffs()[0] * std::log(-z0-zk);
		for (size_t k = 1; k < P; ++k)
			b_0 += iNode->getCoeffs()[k] * pow(-1.0, k) / pow(z0-zk, k);
	}
	localCoeffs.push_back(b_0);

	for (size_t k = 1; k < P; ++k) {
		for (const auto& iNode : iList) {
			cmplx z0 = iNode->getCenter();
			b_k -= iNode->getCoeffs()[0] / (static_cast<double>(k) * pow(z0-zk, k));
			for (size_t l = 1; l < P; ++l)
				b_k += iNode->getCoeffs()[l] * pow(-1.0,l) / pow(z0-zk, k+l) * binom(k+l-1,l-1);
		}
		localCoeffs.push_back(b_k);
	}

    localCoeffs += shiftBaseLocalCoeffs(P);

	for (const auto& branch : branches)
		branch->buildLocalCoeffs(P);

    // nbors.clear();
    iList.clear();
}

void Stem::printNode(std::ofstream& f) {
	f << zk << " " << L << " " << nodeStat << std::endl;
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

    node->buildInteractionList();
    auto nbors = node->getInteractionList();

    for (const auto& nbor : nbors)
	    nbor->setNodeStat(2);

    std::ofstream psnFile, nodeFile;
    psnFile.open("out/psnitions.txt");
    nodeFile.open("out/nodes.txt");

    printpsn(psnFile);
    printNode(nodeFile);
}