#include <cassert>
#include <iostream>
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
        dzk *= L_ / 4.0;

        std::shared_ptr<Node> branch;
        // if branch has at least two particles and is not lvl 0, then further subdivide it
        if (branchpsn[k].size() > 1 && lvl-1)
            branch = std::make_shared<Stem>(branchpsn[k], branchQs[k], zk+dzk, L_/2, lvl-1, k, this);
        else
            branch = std::make_shared<Leaf>(branchpsn[k], branchQs[k], zk+dzk, L_/2, lvl-1, k, this);

        branches.push_back(branch);
    }
}

void Stem::buildMpoleCoeffs() {
    cmplx b_0;
    for (const auto& branch : branches) {
        branch->buildMpoleCoeffs();
        auto branchCoeffs = branch->getMpoleCoeffs();
        b_0 += branchCoeffs[0];
    }
    coeffs.push_back(b_0);

    for (size_t k = 1; k <= P_; ++k) {
        cmplx b_k;
        for (const auto& branch : branches) {
            auto branchCoeffs = branch->getMpoleCoeffs();
            auto z0 = branch->getCenter();
            b_k -= branchCoeffs[0] * pow(z0-zk, k) / static_cast<double>(k);
            for (size_t l = 1; l <= k; ++l)
                b_k += branchCoeffs[l] * pow(z0-zk, k - l) * static_cast<double>(binomTable[k-1][l-1]);
        }
        coeffs.push_back(b_k);
    }
}

void Stem::buildLocalCoeffs() {
    buildNearNeighbors();

    if (!isRoot()) {
        buildInteractionList();
        cmplx b_0;

        for (const auto& iNode : iList) {
            auto z0 = iNode->getCenter();
            auto mpoleCoeffs = iNode->getMpoleCoeffs();
            b_0 += mpoleCoeffs[0] * std::log(zk-z0);
            for (size_t k = 1; k <= P_; ++k)
                b_0 += mpoleCoeffs[k] * pow(-1.0, k) / pow(z0-zk, k);
        }
        localCoeffs.push_back(b_0);

        for (size_t k = 1; k <= P_; ++k) {
            cmplx b_k;
            for (const auto& iNode : iList) {
                auto z0 = iNode->getCenter();
                auto mpoleCoeffs = iNode->getMpoleCoeffs();
                b_k -= mpoleCoeffs[0] / (static_cast<double>(k) * pow(z0-zk, k));
                for (size_t l = 1; l <= P_; ++l)
                    b_k += mpoleCoeffs[l] * pow(-1.0, l) / pow(z0-zk, k + l) * 
                        static_cast<double>(binomTable[k+l-1][l-1]);
            }
            localCoeffs.push_back(b_k);
        }
        // comment out for mpoleToLocalTest()
        // if (!base->isRoot()) localCoeffs += shiftBaseLocalCoeffs(); 
        // iList.clear();
    }

    for (const auto& branch : branches)
        branch->buildLocalCoeffs();
}

