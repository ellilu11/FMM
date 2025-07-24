#include <cassert>
#include <iostream>
#include "leaf.h"
#include "math.h"
#include "stem.h"

Stem::Stem(
    ParticleVec& particles,
    const cmplx zk,
    const double L,
    const int lvl,
    const int branchIdx,
    Stem* const base)
    : Node(particles, zk, L, lvl, branchIdx, base)
{
    // Assign every particle in node to a branch node
    std::vector<ParticleVec> branchParts(4);
    for (const auto& particle : particles) {
        size_t k = cmplx2Idx<bool>(particle->getPos() > zk);
        assert(k < branchParts.size());
        branchParts[k].push_back(particle);
    }

    // Construct branch nodes
    for (size_t k = 0; k < branchParts.size(); ++k) {
        cmplx dzk( pow(-1,k%2+1), pow(-1,k/2+1) );
        dzk *= L_ / 4.0;

        std::shared_ptr<Node> branch;
        // if branch has at least two particles and is not lvl 0, then further subdivide it
        if (branchParts[k].size() > 1 && lvl-1)
            branch = std::make_shared<Stem>(branchParts[k], zk+dzk, L_/2, lvl-1, k, this);
        else
            branch = std::make_shared<Leaf>(branchParts[k], zk+dzk, L_/2, lvl-1, k, this);

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

        if (!base->isRoot()) localCoeffs += base->getShiftedLocalCoeffs(zk); 
        iList.clear();
    }

    for (const auto& branch : branches)
        branch->buildLocalCoeffs();
}