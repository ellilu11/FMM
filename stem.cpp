#include "stem.h"

Stem::Stem(
    const ParticleVec& particles,
    const int branchIdx,
    Stem* const base)
    : Node(particles, branchIdx, base)
{
    // Assign every particle in node to a branch based on its position relative to center
    constexpr int nbranch = 8; //  std::pow(2, DIM);
    std::vector<ParticleVec> branchParts(nbranch);
    for (const auto& p : particles) {
        auto k = bools2Idx(p->getPos() > center);
        branchParts[k].push_back(p);
    }
    // Construct branch nodes
    for (size_t k = 0; k < branchParts.size(); ++k) {
        std::shared_ptr<Node> branch;

        if (branchParts[k].size() >= maxNodeParts)
            branch = std::make_shared<Stem>(branchParts[k], k, this);
        else
            branch = std::make_shared<Leaf>(branchParts[k], k, this);

        branches.push_back(branch);
    }
}

void Stem::buildMpoleCoeffs() {
    coeffs.resize(order+1);

    //for (const auto& branch : branches) {
    //    branch->buildMpoleCoeffs();
    //    auto branchCoeffs = branch->getMpoleCoeffs();
    //    coeffs[0] += branchCoeffs[0];

    //    for (size_t l = 1; l <= order; ++l) {
    //        auto branchCoeffs = branch->getMpoleCoeffs();
    //        auto dz = branch->getCenter() - center;

    //        coeffs[l] -= branchCoeffs[0] * pow(dz, l) / static_cast<double>(l);

    //        for (size_t k = 1; k <= l; ++k)
    //            coeffs[l] += branchCoeffs[k] * pow(dz, l - k) * 
    //                            static_cast<double>(binomTable[l-1][k-1]);
    //    }
    //}
}

void Stem::buildLocalCoeffs() {
    buildMpoleToLocalCoeffs();

    for (const auto& branch : branches)
        branch->buildLocalCoeffs();
}