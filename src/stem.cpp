#include "stem.h"

Stem::Stem(
    const ParticleVec& particles,
    const int branchIdx,
    Stem* const base)
    : Node(particles, branchIdx, base)
{
    // Assign every particle to a branch from its position relative to center
    std::array<ParticleVec,8> branchParts;

    for (const auto& p : particles)
        branchParts[Math::bools2Idx(p->getPos() > center)].push_back(p);
 
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

/* buildNeighbors()
 * Find all neighbor nodes of equal or greater size
 */
void Stem::buildNeighbors() {
    assert(!isRoot());

    for (int i = 0; i < numDir; ++i) {
        Dir dir = static_cast<Dir>(i);
        auto nbor = getNeighborGeqSize(dir);
        if (nbor != nullptr)
            nbors.push_back(nbor);
    }
    assert(nbors.size() <= numDir);
}

/* buildLists()
 * Find neighbor and interaction lists
 * Add self as near non-neighbor (list 3 node) of any list 4 nodes
 */
void Stem::buildLists() {
    if (!isRoot()) {
        buildNeighbors();

        buildInteractionList();

        pushSelfToNearNonNbors();
    }

    for (const auto& branch : branches)
        branch->buildLists();
}

/* buildMpoleCoeffs()
 * (M2M) Build mpole coeffs by merging branch mpole coeffs 
 */
void Stem::buildMpoleCoeffs() {
    for (int l = 0; l <= order; ++l)
        coeffs.push_back(vecXcd::Zero(2*l+1));

    for (const auto& branch : branches) {
        branch->buildMpoleCoeffs();
        auto branchCoeffs = branch->getMpoleCoeffs();
        const auto branchIdx = branch->getBranchIdx();

        const double r = (branch->getCenter() - center).norm();

        for (int j = 0; j <= order; ++j) {
            branchCoeffs[j] = wignerD[branchIdx][j] * branchCoeffs[j];

            vecXcd shiftedBranchCoeffs_j = vecXcd::Zero(2*j+1);
            for (int k = -j; k <= j; ++k) {
                int k_j = k + j;
                double r2n = 1.0;

                for (int n = 0; n <= min(j+k, j-k); ++n) {
                    shiftedBranchCoeffs_j[k_j] += branchCoeffs[j-n][k_j-n] *
                        tables.A_[n][n] * tables.A_[j-n][k_j-n] / tables.A_[j][k_j] *
                        r2n; // legendreCos(0.0, n, 0) = 1 for all n;

                    r2n *= r;
                }
            }

            coeffs[j] += wignerDInv[branchIdx][j] * shiftedBranchCoeffs_j;
        }
    }
}

/* propagateExpCoeffs()
 * (M2X) Convert mpole coeffs into outgoing exp coeffs
 * (X2X) Translate outgoing exp coeffs to nodes in all dirlists
 */ 
void Stem::propagateExpCoeffs() {
    if (!isRoot()) {
        for (int dir = 0; dir < 6; ++dir) {
            auto start = Clock::now();

            auto expCoeffs = getMpoleToExpCoeffs(dir);
            
            t.M2X += Clock::now() - start;

            auto iList = dirList[dir];
            
            start = Clock::now();
            
            for (const auto& node : iList)
                node->addShiftedExpCoeffs(expCoeffs, center, dir);
            
            t.X2X += Clock::now() - start;
        }
    }

    for (const auto& branch : branches)
        branch->propagateExpCoeffs();
}

/* buildLocalCoeffs() 
 * (X2L) Receive incoming exp coeffs and add to local coeffs
 * (P2L) Add contribution from list 4 nodes t to local coeffs
 * (L2L) Shift base local coeffs to center and add to local coeffs
 */
void Stem::buildLocalCoeffs() {
    if (!isRoot()) {
        auto start = Clock::now();

        evalExpToLocalCoeffs();

        t.X2L += Clock::now() - start;

        start = Clock::now();

        evalLeafIlistSols();

        t.P2L += Clock::now() - start;

        start = Clock::now();
        
        if (!base->isRoot()) {
            auto shiftedLocalCoeffs = base->getShiftedLocalCoeffs(branchIdx);

            for (int l = 0; l <= order; ++l)
                localCoeffs[l] += shiftedLocalCoeffs[l];
        }
        
        t.L2L += Clock::now() - start;
    }

    for (const auto& branch : branches)
        branch->buildLocalCoeffs();
}

