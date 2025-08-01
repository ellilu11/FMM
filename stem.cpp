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
    for (const auto& p : particles)
        branchParts[bools2Idx(p->getPos() > center)].push_back(p);
 
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
    coeffs.resize(pow(order+1, 2));

    for (const auto& branch : branches) {
        branch->buildMpoleCoeffs();
        auto branchCoeffs = branch->getMpoleCoeffs();
        auto [r, th, ph] = cart2Sph(branch->getCenter() - center);
        // auto Rmatrix = rotationMatrix(branch->getBranchIdx());
        // branchCoeffs *= Rmatrix;
        int idx = 0;

        /*for (int j = 0; j <= order; ++j) {
            for (int k = -j; k <= j; ++k) {

                for (int n = 0; n < j; ++n) {
                    for (int m = -n; m < n; ++m) {
                        auto idx2 = lm2Idx(j-n, k-m);
                        coeffs[idx] += 
                            branchCoeffs[idx2] * pow(iu, )
                    }
                }
                idx++;
            }
        }*/
    }
}

void Stem::buildLocalCoeffs() {
    buildMpoleToLocalCoeffs();

    for (const auto& branch : branches)
        branch->buildLocalCoeffs();
}

