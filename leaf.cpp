#include "leaf.h"
#include "math.h"
#include "node.h"
#include <iostream>

Leaf::Leaf(cmplxVec& psn,
    realVec& qs,
    const cmplx zk,
    const double L,
    const int lvl,
    const int branchIdx,
    Stem* const base)
    : Node(psn, qs, zk, L, lvl, branchIdx, base)
{
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
            b_0 += mpoleCoeffs[0] * std::log(zk-z0);
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

        // comment out for mpoleToLocalTest()
        if (!base->isRoot()) localCoeffs += base->getShiftedLocalCoeffs(zk); 
        // iList.clear();
    }

    evaluatePhi();
    // evaluateFld();
}

cmplxVec Leaf::getPhiFarSrc() {
    cmplxVec phis;
    for (const auto& obs : psn) {
        cmplx phi;
        for (size_t k = 0; k <= P_; ++k)
            phi -= localCoeffs[k] * pow(obs-zk, k);
        phis.push_back(phi);
    }
    return phis;
}

cmplxVec Leaf::getPhiNearSrc() {
    cmplxVec phis;

    for (size_t obs = 0; obs < psn.size(); ++obs) {
        cmplx phi;

        // phi due to other particles in this node (apply reciprocity later)
        for (size_t src = 0; src < psn.size(); ++src)
            if (src != obs) phi -= qs[src] * std::log(psn[obs] - psn[src]);

        // phi due to particles in neighboring nodes
        for (const auto& nbor : nbors) {
            cmplxVec psnNbor = nbor->getPsn();
            for (size_t src = 0; src < psnNbor.size(); ++src)
                phi -= (nbor->getQs())[src] * std::log(psn[obs] - psnNbor[src]);
        }
        phis.push_back(phi);
    }
    return phis;
}

void Leaf::evaluatePhi() {
    phis.resize(psn.size());
    for (size_t n = 0; n < psn.size(); ++n)
        phis[n] = getPhiFarSrc()[n] + getPhiNearSrc()[n];
    // phis = getPhiFarSrc() + getPhiNearSrc();
}

