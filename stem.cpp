#include "stem.h"

Stem::Stem(
    ParticleVec& particles,
    const cmplx center,
    const int branchIdx,
    Stem* const base)
    : Node(particles, center, branchIdx, base)
{
    // Assign every particle in node to a branch based on its position relative to center
    std::vector<ParticleVec> branchParts(4);
    for (const auto& p : particles)
        branchParts[bools2Idx(p->getPos() > center)].push_back(p);

    // Construct branch nodes
    for (size_t k = 0; k < branchParts.size(); ++k) {
        cmplx dz( pow(-1,k%2+1), pow(-1,k/2+1) );
        dz *= nodeLeng / 4.0;

        std::shared_ptr<Node> branch;
        // if branch has at least two particles and is not max lvl, then further subdivide it
        if (branchParts[k].size() > 1 && lvl < maxLvl)
            branch = std::make_shared<Stem>(branchParts[k], center+dz, k, this);
        else
            branch = std::make_shared<Leaf>(branchParts[k], center+dz, k, this);

        branches.push_back(branch);
    }
}

void Stem::buildMpoleCoeffs() {
    coeffs.resize(order+1);

    for (const auto& branch : branches) {
        branch->buildMpoleCoeffs();
        auto branchCoeffs = branch->getMpoleCoeffs();
        coeffs[0] += branchCoeffs[0];

        for (size_t l = 1; l <= order; ++l) {
            auto branchCoeffs = branch->getMpoleCoeffs();
            auto dz = branch->getCenter() - center;

            coeffs[l] -= branchCoeffs[0] * pow(dz, l) / static_cast<double>(l);

            //cmplxVec innerCoeffs;
            //for (ptrdiff_t n = l-1; n >=0; --n)
            //    innerCoeffs.push_back( branchCoeffs[l-n] 
            //        * static_cast<double>(binomTable[l-1][l-n-1]) );
            //coeffs[l] += evaluatePoly<cmplx>(innerCoeffs, dz);
            for (size_t k = 1; k <= l; ++k)
                coeffs[l] += branchCoeffs[k] * pow(dz, l - k) * 
                                static_cast<double>(binomTable[l-1][k-1]);
        }
    }
}

void Stem::buildLocalCoeffs() {
    buildMpoleToLocalCoeffs();

    for (const auto& branch : branches)
        branch->buildLocalCoeffs();
}