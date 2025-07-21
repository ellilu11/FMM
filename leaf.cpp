#include "leaf.h"
#include "node.h"
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
    phis.resize(psn.size());
}

void Leaf::buildMpoleCoeffs() {
    for (int k = 0; k <= P_; ++k) {
        cmplx a_k;
        for (size_t n = 0; n < psn.size(); ++n)
            a_k += qs[n] *
                k == 0 ? 1 : -pow(psn[n]-zk, k) / static_cast<double>(k);
        coeffs.push_back(a_k);
    }
}

void Leaf::buildLocalCoeffs() {
    buildNearNeighbors();

    if (!isRoot()) {
        buildInteractionList();
        cmplx b_0;

        for (const auto& iNode : iList) {
            auto z0 = iNode->getCenter();
            auto mpoleCoeffs = iNode->getMpoleCoeffs();
            b_0 += mpoleCoeffs[0] * std::log(-(z0-zk));
            for (size_t k = 1; k <= P_; ++k)
                b_0 += mpoleCoeffs[k] * pow(-1.0, k) / pow(z0 - zk, k);
        }
        localCoeffs.push_back(b_0);

        for (size_t k = 1; k <= P_; ++k) {
            cmplx b_k;
            for (const auto& iNode : iList) {
                auto z0 = iNode->getCenter();
                auto mpoleCoeffs = iNode->getMpoleCoeffs();
                b_k -= mpoleCoeffs[0] / (static_cast<double>(k) * pow(z0 - zk, k));
                for (size_t l = 1; l <= P_; ++l)
                    b_k += mpoleCoeffs[l] * pow(-1.0, l) / pow(z0 - zk, k + l) *
                        static_cast<double>(binomTable[k+l-1][l-1]);
            }
            localCoeffs.push_back(b_k);
        }

        if (!base->isRoot()) localCoeffs += shiftBaseLocalCoeffs();
        iList.clear();
    }

    evaluatePhi();
    // evaluateFld();
}

void Leaf::evaluatePhiLocalExp() {
    cmplx phi;
    cmplxVec phisLocalExp;
    for (const auto& ele : psn) {
        for (size_t k = 0; k < P_; ++k)
            phi += localCoeffs[k] * pow(ele, k);
        phisLocalExp.push_back(phi);
    }
    phis += phisLocalExp;
}

void Leaf::evaluatePhiDirect() {
    cmplx phi;
    cmplxVec phisDirect;
    // cmplxVec 

    // phi due to other particles in this node (apply reciprocity later)
    for (size_t m = 0; m < psn.size(); ++m)
        for (size_t n = 0; n < psn.size(); ++n)
            if (n != m) phi -= qs[n] * std::log(psn[m] - psn[n]);

    // phi due to particles in neighboring nodes
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

void Leaf::evaluatePhi() {
    evaluatePhiLocalExp();
    evaluatePhiDirect();
}

void Leaf::iListTest() {
    setNodeStat(3);
    buildInteractionList();
    auto nbors = getInteractionList();

    for (const auto& nbor : nbors)
        nbor->setNodeStat(2);

    std::ofstream psnFile, nodeFile;
    psnFile.open("out/srcs.txt");
    nodeFile.open("out/nodes.txt");

    printPsn(psnFile);
    printNode(nodeFile);
}