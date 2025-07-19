#include "leaf.h"
#include <iostream>

Leaf::Leaf(cmplxVec& psn,
	std::vector<double>& qs,
	const cmplx zk,
	const double L,
	const int lvl,
	const int branchIdx,
	Stem* const base)
	: Node(psn, qs, zk, L, lvl, branchIdx, base)
{
    phis.reserve(psn.size());
}

void Leaf::buildCoeffs(const int P) {
	//auto sumOverpsn =
	//	[this](const double q, const cmplx psn, const int idx) {
	//		return q *
	//			idx == 0 ? 1 : -pow(psn-z0, idx) / static_cast<double>(idx); 
	//	};
	cmplx a_k;
	for (int k = 0; k < P; ++k) {
		for (size_t i = 0; i < psn.size(); ++i)
			a_k += qs[i] * 
				k == 0 ? 1 : -pow(psn[i]-zk, k) / static_cast<double>(k);
		coeffs.push_back(a_k);
	}
}

void Leaf::buildLocalCoeffs(const int P) {
    buildInteractionList();
    cmplx b_0, b_k;

    for (const auto& iNode : iList) {
        cmplx z0 = iNode->getCenter();
        b_0 += iNode->getCoeffs()[0] * std::log(-z0 - zk);
        for (size_t k = 1; k < P; ++k)
            b_0 += iNode->getCoeffs()[k] * pow(-1.0, k) / pow(z0 - zk, k);
    }
    localCoeffs.push_back(b_0);

    for (size_t k = 1; k < P; ++k) {
        for (const auto& iNode : iList) {
            cmplx z0 = iNode->getCenter();
            b_k -= iNode->getCoeffs()[0] / (static_cast<double>(k) * pow(z0 - zk, k));
            for (size_t l = 1; l < P; ++l)
                b_k += iNode->getCoeffs()[l] * pow(-1.0, l) / pow(z0 - zk, k + l) * binom(k + l - 1, l - 1);
        }
        localCoeffs.push_back(b_k);
    }

    localCoeffs += shiftBaseLocalCoeffs(P);
    iList.clear();

    evaluatePhi(P);
    // evaluateFld(P);
}

void Leaf::evaluatePhiLocalExp(const int P) {
    cmplx phi;
    cmplxVec phisLocalExp;
    for (const auto& ele : psn) {
        for (size_t k = 0; k < P; ++k)
            phi += localCoeffs[k] * pow(ele, k);
        phisLocalExp.push_back(phi);
    }
    phis += phisLocalExp;
}

void Leaf::evaluatePhiDirect(const int P) {
    cmplx phi;
    cmplxVec phisDirect;
    // cmplxVec 

    for (const auto& nbor : nbors) {

    }

    for (const auto& obs : psn) {
        for (const auto& nbor : nbors) {
            cmplxVec srcs = nbor->getPsn();
            for (size_t n = 0; n < srcs.size(); ++n)
                phi -= (nbor->getQs())[n] * std::log(obs - srcs[n]);
        }
        phisDirect.push_back(phi);
    }
    phis += phisDirect;
}

void Leaf::evaluatePhi(const int P) {
    evaluatePhiLocalExp(P);
    evaluatePhiDirect(P);
}

void Leaf::printNode(std::ofstream& f) {
	f << zk << " " << L << " " << nodeStat << std::endl;
}

void Leaf::iListTest() {
	setNodeStat(3);
	buildInteractionList();
	auto nbors = getInteractionList();

	for (const auto& nbor : nbors)
		nbor->setNodeStat(2);

	std::ofstream psnFile, nodeFile;
	psnFile.open("out/psnitions.txt");
	nodeFile.open("out/nodes.txt");

	printpsn(psnFile);
	printNode(nodeFile);
}