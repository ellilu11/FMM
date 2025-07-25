#include "stem.h"

Stem::Stem(
    ParticleVec& particles,
    const cmplx center,
    const double L,
    const int lvl,
    const int branchIdx,
    Stem* const base)
    : Node(particles, center, L, lvl, branchIdx, base)
{
    // Assign every particle in node to a branch node
    std::vector<ParticleVec> branchParts(4);
    for (const auto& particle : particles) {
        size_t k = cmplx2Idx<bool>(particle->getPos() > center);
        assert(k < branchParts.size());
        branchParts[k].push_back(particle);
    }

    // Construct branch nodes
    for (size_t k = 0; k < branchParts.size(); ++k) {
        cmplx dcenter( pow(-1,k%2+1), pow(-1,k/2+1) );
        dcenter *= L_ / 4.0;

        std::shared_ptr<Node> branch;
        // if branch has at least two particles and is not lvl 0, then further subdivide it
        if (branchParts[k].size() > 1 && lvl-1)
            branch = std::make_shared<Stem>(branchParts[k], center+dcenter, L_/2, lvl-1, k, this);
        else
            branch = std::make_shared<Leaf>(branchParts[k], center+dcenter, L_/2, lvl-1, k, this);

        branches.push_back(branch);
    }
}

void Stem::buildMpoleCoeffs() {
    coeffs.resize(P_+1);

    for (const auto& branch : branches) {
        branch->buildMpoleCoeffs();
        auto branchCoeffs = branch->getMpoleCoeffs();
        coeffs[0] += branchCoeffs[0];

        for (size_t l = 1; l <= P_; ++l) {
            auto branchCoeffs = branch->getMpoleCoeffs();
            auto dz = branch->getCenter() - center;

            coeffs[l] -= branchCoeffs[0] * pow(dz, l) / static_cast<double>(l);

            //cmplxVec innerCoeffs;
            //for (ptrdiff_t n = l-1; n >=0; --n)
            //    innerCoeffs.push_back( branchCoeffs[l-n] 
            //        * static_cast<double>(binomTable[l-1][l-n-1]) );
            //coeffs[l] += evaluatePoly<cmplx>(innerCoeffs, dz);
            for (size_t k = 1; k <= l; ++k)
                coeffs[l] += branchCoeffs[k] * pow(dz, l - k) * static_cast<double>(binomTable[l-1][k-1]);
        }
    }
}

void Stem::buildLocalCoeffs() {
    buildLocalCoeffsFromIList();

    for (const auto& branch : branches)
        branch->buildLocalCoeffs();
}