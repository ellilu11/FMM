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

        if (branchParts[k].size() > maxNodeParts)
            branch = std::make_shared<Stem>(branchParts[k], k, this);
        else
            branch = std::make_shared<Leaf>(branchParts[k], k, this);

        branches.push_back(branch);
    }
}

void Stem::buildMpoleCoeffs() {
    for (int l = 0; l <= order; ++l)
        coeffs.emplace_back(vecXcd::Zero(2*l+1)) ;

    // precompute as LUT
    realVec legendreLMCoeffs;
    for (int n = 0; n <= order; ++n)
        legendreLMCoeffs.push_back(legendreLM(0, n, 0));

    for (const auto& branch : branches) {
        branch->buildMpoleCoeffs();
        auto branchCoeffs = branch->getMpoleCoeffs();
        auto branchIdx = branch->getBranchIdx();

        // apply rotation
        for (int l = 0; l < branchCoeffs.size(); ++l) {
            //auto mat = rotationMat[branchIdx][l];
            //std::cout << branchIdx << ' ' << l << ' '
            //    << '(' << mat.rows() << ',' << mat.cols() << ") "
            //    << '(' << branchCoeffs[l].rows() << ',' << branchCoeffs[l].cols() << ")\n";    
            branchCoeffs[l] = rotationMat[branchIdx][l] * branchCoeffs[l];
        }

        //auto dR = toSph(branch->getCenter() - center);
        //auto r = dR[0];
        double r = (branch->getCenter() - center).norm();

        for (size_t j = 0; j <= order; ++j) {
            vecXcd rotatedBranchCoeffs_j(2*j+1);

            for (int k = -j; k <= j; ++k) {
                size_t k_ = k+j;

                for (size_t n = 0; n <= j; ++n)
                    rotatedBranchCoeffs_j[k_] +=
                        branchCoeffs[j-n][k_] * A[n][0] * A[j-n][k_] / A[j][k_] *
                        pow(r, n) * legendreLMCoeffs[n];
            }
            // apply inverse rotation
            // assume rotationMat is unitary, so its adjoint is its inverse
            coeffs[j] += rotationMat[branchIdx][j].adjoint() * rotatedBranchCoeffs_j;
        }
    }
}

void Stem::buildLocalCoeffs() {
    buildMpoleToLocalCoeffs();

    for (const auto& branch : branches)
        branch->buildLocalCoeffs();
}

