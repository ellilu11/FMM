#include "stem.h"

Stem::Stem(
    const ParticleVec& particles,
    const int branchIdx,
    Stem* const base)
    : Node(particles, branchIdx, base)
{
    // Assign every particle in node to a branch based on its position relative to center
    std::array<ParticleVec,8> branchParts;
    for (const auto& p : particles)
        branchParts[bools2Idx(p->getPos() > center)].push_back(p);
 
    // Construct branch nodes
    for (size_t k = 0; k < branchParts.size(); ++k) {
        std::shared_ptr<Node> branch;

        if (branchParts[k].size() > maxNodeParts)
            branch = std::make_shared<Stem>(branchParts[k], k, this);
        else
            branch = std::make_shared<Leaf>(branchParts[k], k, this);

        branches.push_back(branch);
    }
}

void Stem::buildNbors() {
    assert(!isRoot());

    for (int i = 0; i < numDir; ++i) {
        Dir dir = static_cast<Dir>(i);
        auto nbor = getNeighborGeqSize(dir);
        if (nbor != nullptr)
            nbors.push_back(nbor);
    }
    assert(nbors.size() <= numDir);
}

void Stem::buildLists() {
    if (!isRoot()) {
        buildNbors();
        buildInteractionList();
        buildOuterInteractionList();
    }

    for (const auto& branch : branches)
        branch->buildLists();
}

void Stem::buildMpoleCoeffs() {
    for (int l = 0; l <= order; ++l)
        coeffs.emplace_back(vecXcd::Zero(2*l+1));

    for (const auto& branch : branches) {
        branch->buildMpoleCoeffs();
        auto branchCoeffs = branch->getMpoleCoeffs();
        auto branchIdx = branch->getBranchIdx();
        double r = (branch->getCenter() - center).norm();

        for (int j = 0; j <= order; ++j) {
            branchCoeffs[j] = wignerD[branchIdx][j] * branchCoeffs[j];

            vecXcd shiftedBranchCoeffs_j = vecXcd::Zero(2*j+1);
            for (int k = -j; k <= j; ++k) {
                int k_ = k + j;
                double r2n = 1.0;
                for (int n = 0; n <= min(j+k, j-k); ++n) {
                    // if ( max(k+n-j, -n) <= 0 && 0 <= min(k+j-n, n) )
                    shiftedBranchCoeffs_j[k_] += branchCoeffs[j-n][k_-n] *
                        tables.A_[n][n] * tables.A_[j-n][k_-n] / tables.A_[j][k_] *
                        r2n; // legendreCos(0.0, n, 0) = 1 for all n;
                    r2n *= r;
                }
            }

            coeffs[j] += wignerDInv[branchIdx][j] * shiftedBranchCoeffs_j;
        }
    }
}

void Stem::propagateExpCoeffs() {
    if (!isRoot()) {
        for (int dir = 0; dir < 6; ++dir) {
            auto start = chrono::high_resolution_clock::now();

            // build exp coeffs of branches, merge them, and propagate to outer dirlist
            auto mergedExpCoeffs = getMergedExpCoeffs(dir);

            
            t_M2X += chrono::high_resolution_clock::now() - start;

            start = chrono::high_resolution_clock::now();

            for (const auto& iNode : outerDirList[dir])
                iNode->addShiftedExpCoeffs(mergedExpCoeffs, center, dir);

            // for lvl > 1, propagate own exp coeffs to inner dirlist
            if (!base->isRoot())
                for (const auto& iNode : dirList[dir])
                    iNode->addShiftedExpCoeffs(expCoeffsOut[dir], center, dir);

            t_X2X += chrono::high_resolution_clock::now() - start;

        }
        expCoeffsOut = {};
    }

    for (const auto& branch : branches)
        branch->propagateExpCoeffs();
}

void Stem::buildLocalCoeffs() {
    if (!isRoot()) {
        auto start = std::chrono::high_resolution_clock::now();

        buildLocalCoeffsFromLeafIlist();

        t_X2L_l4 += std::chrono::high_resolution_clock::now() - start;

        start = std::chrono::high_resolution_clock::now();

        buildLocalCoeffsFromDirList();

        t_X2L += std::chrono::high_resolution_clock::now() - start;

        start = std::chrono::high_resolution_clock::now();

        if (!base->isRoot()) {
            auto shiftedLocalCoeffs = base->getShiftedLocalCoeffs(branchIdx);
            for (int l = 0; l <= order; ++l)
                localCoeffs[l] += shiftedLocalCoeffs[l];
        }

        t_L2L += std::chrono::high_resolution_clock::now() - start;
    }

    for (const auto& branch : branches)
        branch->buildLocalCoeffs();
}

