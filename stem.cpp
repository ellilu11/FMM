#include "stem.h"

Stem::Stem(
    const ParticleVec& particles,
    const int branchIdx,
    Stem* const base)
    : Node(particles, branchIdx, base)
{
    // std::cout << "Stem has " << particles.size() << " particles\n";

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

void Stem::buildLists() {
    if (!isRoot()) {
        buildNearNeighbors();
        buildInteractionList();
        buildDirectedIList();
    }

    for (const auto& branch : branches)
        branch->buildLists();
}

// rotation matrices
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
                        r2n; // legendreLM(0.0, n, 0) = 1 for all n;
                    r2n *= r;
                }
            }

            coeffs[j] += wignerDInv[branchIdx][j] * shiftedBranchCoeffs_j;
        }
    }
}

// no rotation matrices
/*void Stem::buildMpoleCoeffs() {
    using namespace std;
    for (int l = 0; l <= order; ++l)
        coeffs.emplace_back(vecXcd::Zero(2*l+1));

    for (const auto& branch : branches) {
        branch->buildMpoleCoeffs();
        auto branchCoeffs = branch->getMpoleCoeffs();
        auto branchIdx = branch->getBranchIdx();

        auto dR = toSph(branch->getCenter() - center);
        double r = dR[0], th = dR[1], ph = dR[2];

        for (int j = 0; j <= order; ++j) {
            for (int k = -j; k <= j; ++k) {
                int k_ = k + j;

                for (int n = 0; n <= j; ++n) {
                    for (int m = max(k+n-j,-n); m <= min(k+j-n,n); ++m) {
                        int m_ = m + n;

                        coeffs[j][k_] +=
                            branchCoeffs[j-n][k_-m_] * pow(iu, abs(k)-abs(m)-abs(k-m)) *
                            tables.A[n][m_] * tables.A[j-n][k_-m_] / tables.A[j][k_] *
                            pow(r, n) * legendreLM(th, n, abs(-m)) * expI(static_cast<double>(-m)*ph);
                    }
                }
            }
        }
    }
}*/

void Stem::propagateExpCoeffs() {
    if (!isRoot()) {
        for (int dir = 0; dir < 6; ++dir) {
            auto expCoeffs = getMpoleToExpCoeffs(dir);
            auto iList = dirList[dir];
            for (const auto& iNode : iList)
                iNode->addShiftedExpCoeffs(expCoeffs, center, dir);
        }
    }

    for (const auto& branch : branches)
        branch->propagateExpCoeffs();
}

void Stem::buildLocalCoeffs() {
    if (!isRoot()) {
        buildLocalCoeffsFromLeafIlist();
        buildLocalCoeffsFromDirList();

        if (!base->isRoot())
            for (int l = 0; l <= order; ++l)
                localCoeffs[l] += (base->getShiftedLocalCoeffs(branchIdx))[l];
    }

    for (const auto& branch : branches)
        branch->buildLocalCoeffs();
}

