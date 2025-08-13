#include "node.h"

using enum Dir;

int Node::order;
int Node::orderExp;
int Node::maxNodeParts;
double Node::rootLeng;
Tables Node::tables;
std::array<std::vector<matXcd>,14> Node::wignerD;
std::array<std::vector<matXcd>,14> Node::wignerDInv;
std::array<mat3d,6> Node::rotMatR;

void Node::setNodeParams(const Config& config) {
    order = ceil(-std::log(config.EPS) / std::log(2));
    maxNodeParts = config.maxNodeParts;
    rootLeng = config.L;
    orderExp = [&]() -> std::size_t {
        switch (config.prec) {
            case Precision::LOW:    return 8;
            case Precision::MEDIUM: return 17;
            case Precision::HIGH:   return 26;
        } 
        }();
}

void Node::buildTables(const Config& config) {
    tables = Tables(order, config.prec);
    assert(orderExp == tables.quadCoeffs_.size());
}

void Node::buildRotationMats() {
    for (int dir = 0; dir < 14; ++dir) {
        auto X = 
            dir < 8 ? 
            idx2pm(dir) : // vertex directions (for M2M and L2L)
            [dir] { 
            switch (dir-8) {
                case 0: return vec3d(0, 0, 1);
                case 1: return vec3d(0, 0, -1);
                case 2: return vec3d(0, 1, 0);
                case 3: return vec3d(0, -1, 0);
                case 4: return vec3d(1, 0, 0);
                case 5: return vec3d(-1, 0, 0);
            } 
            }(); // cardinal directions (for M2X2L)

        auto R = toSph(X);
        pair2d angles(R[1], R[2]);

        wignerD[dir] = wignerDAlongDir(angles,0);
        wignerDInv[dir] = wignerDAlongDir(angles,1);
        if (dir >=8)
            rotMatR[dir-8] = rotationR(angles);

        // std::cout << dir << '\n' << rotationInvMat[dir][2] << "\n\n";
    }
}

std::vector<matXcd> Node::wignerDAlongDir(const pair2d angles, bool isInv) {
     std::vector<matXcd> mats;
     for (int l = 0; l <= order; ++l)
         mats.push_back(
             isInv ?
             wignerD_l(angles, l) :
             wignerD_l(angles, l).adjoint()
         );

    return mats;
}

// Ylm except constant coefficients and \exp(i m \phi) phase factor
const double Node::legendreLM(const double th, const pair2i lm) {
    auto [l, abs_m] = lm;
    assert(abs_m <= l);

    const auto cos_th = cos(th);
    const auto sin_th = sin(th);
    double legendreSum = 0.0;

    // term is zero for l-k odd
    for (int k = l; k >= abs_m; k -= 2)
        legendreSum += tables.fallingFact_[k][abs_m] * tables.legendreSum_[l][k] * pow(cos_th, k-abs_m);

    return tables.coeffYlm_[l][abs_m] * pow(sin_th, abs_m) * legendreSum;
}

Node::Node(
    const ParticleVec& particles,
    const int branchIdx,
    Node* const base)
    : particles(particles), branchIdx(branchIdx), base(base),
    nodeLeng(base == nullptr ? rootLeng : base->getLeng() / 2.0),
    center(base == nullptr ? zeroVec :
        base->getCenter() + nodeLeng / 2.0 * idx2pm(branchIdx)),
    nodeStat(0), useRot(0)
{
    for (int dir = 0; dir < 6; ++dir) {
        std::vector<vecXcd> expCoeffs_dir;
        for (int k = 0; k < orderExp; ++k)
            expCoeffs_dir.emplace_back( vecXcd::Zero(tables.quadLengs_[k]) );
        expCoeffs[dir] = expCoeffs_dir;
    }
}

std::shared_ptr<Node> const Node::getNeighborGeqSize(const Dir dir) {
    if (isRoot()) return nullptr;

    std::shared_ptr<Node> nbor;

    switch (dir) {
        case U:
            if (branchIdx < 4)
                return base->branches[branchIdx+4];
            nbor = base->getNeighborGeqSize(dir);
            if (nbor == nullptr) return nbor;
            if (nbor->isNodeType<Leaf>()) return nbor;
            return nbor->branches[branchIdx-4];
            break;

        case N:
            if (branchIdx == 0 || branchIdx == 1 || branchIdx == 4 || branchIdx == 5)
                return base->branches[branchIdx+2];
            nbor = base->getNeighborGeqSize(dir);
            if (nbor == nullptr) return nbor;
            if (nbor->isNodeType<Leaf>()) return nbor;
            return nbor->branches[branchIdx-2];
            break;

        case E:
            if (branchIdx == 0 || branchIdx == 2 || branchIdx == 4 || branchIdx == 6)
                return base->branches[branchIdx+1];
            nbor = base->getNeighborGeqSize(dir);
            if (nbor == nullptr) return nbor;
            if (nbor->isNodeType<Leaf>()) return nbor;
            return nbor->branches[branchIdx-1];
            break;

        case W:
            if (branchIdx == 1 || branchIdx == 3 || branchIdx == 5 || branchIdx == 7)
                return base->branches[branchIdx-1];
            nbor = base->getNeighborGeqSize(dir);
            if (nbor == nullptr) return nbor;
            if (nbor->isNodeType<Leaf>()) return nbor;
            return nbor->branches[branchIdx+1];
            break;

        case S:
            if (branchIdx == 2 || branchIdx == 3 || branchIdx == 6 || branchIdx == 7)
                return base->branches[branchIdx-2];
            nbor = base->getNeighborGeqSize(dir);
            if (nbor == nullptr) return nbor;
            if (nbor->isNodeType<Leaf>()) return nbor;
            return nbor->branches[branchIdx+2];
            break;

        case D:
            if (branchIdx >= 4)
                return base->branches[branchIdx-4];
            nbor = base->getNeighborGeqSize(dir);
            if (nbor == nullptr) return nbor;
            if (nbor->isNodeType<Leaf>()) return nbor;
            return nbor->branches[branchIdx+4];
            break;

        case UN:
            if (branchIdx == 0 || branchIdx == 1)
                return base->branches[branchIdx+6];
            if (branchIdx == 6 || branchIdx == 7) {
                nbor = base->getNeighborGeqSize(dir);
                if (nbor == nullptr) return nbor;
                if (nbor->isNodeType<Leaf>()) return nbor;
                return nbor->branches[branchIdx-6];
            } else {
                nbor = (branchIdx == 4 || branchIdx == 5) ?
                    base->getNeighborGeqSize(U) :
                    base->getNeighborGeqSize(N);
                if (nbor == nullptr) return nbor;
                // do not double count neighbor that will be found along a cardinal direction
                if (nbor->isNodeType<Leaf>()) return nullptr;
                return (branchIdx == 2 || branchIdx == 4) ?
                    nbor->branches[6-branchIdx] :
                    nbor->branches[8-branchIdx];
            }
            break;

        case DS:
            if (branchIdx == 6 || branchIdx == 7)
                return base->branches[branchIdx-6];
            if (branchIdx == 0 || branchIdx == 1) {
                nbor = base->getNeighborGeqSize(dir);
                if (nbor == nullptr) return nbor;
                if (nbor->isNodeType<Leaf>()) return nbor;
                return nbor->branches[branchIdx+6];
            } else {
                nbor = (branchIdx == 2 || branchIdx == 3) ?
                    base->getNeighborGeqSize(D) :
                    base->getNeighborGeqSize(S);
                if (nbor == nullptr) return nbor;
                if (nbor->isNodeType<Leaf>()) return nullptr;
                return (branchIdx == 2 || branchIdx == 4) ?
                    nbor->branches[6-branchIdx] :
                    nbor->branches[8-branchIdx];
            }
            break;

        case US:
            if (branchIdx == 2 || branchIdx == 3)
                return base->branches[branchIdx+2];
            if (branchIdx == 4 || branchIdx == 5) {
                nbor = base->getNeighborGeqSize(dir);
                if (nbor == nullptr) return nbor;
                if (nbor->isNodeType<Leaf>()) return nbor;
                return nbor->branches[branchIdx-2];
            } else {
                nbor = (branchIdx == 6 || branchIdx == 7) ?
                    base->getNeighborGeqSize(U) :
                    base->getNeighborGeqSize(S);
                if (nbor == nullptr) return nbor;
                if (nbor->isNodeType<Leaf>()) return nullptr;
                return (branchIdx == 0 || branchIdx == 6) ?
                    nbor->branches[6-branchIdx] :
                    nbor->branches[8-branchIdx];
            }
            break;

        case DN:
            if (branchIdx == 4 || branchIdx == 5)
                return base->branches[branchIdx-2];
            if (branchIdx == 2 || branchIdx == 3) {
                nbor = base->getNeighborGeqSize(dir);
                if (nbor == nullptr) return nbor;
                if (nbor->isNodeType<Leaf>()) return nbor;
                return nbor->branches[branchIdx+2];
            } else {
                nbor = (branchIdx == 0 || branchIdx == 1) ?
                    base->getNeighborGeqSize(D) :
                    base->getNeighborGeqSize(N);
                if (nbor == nullptr) return nbor;
                if (nbor->isNodeType<Leaf>()) return nullptr;
                return (branchIdx == 0 || branchIdx == 6) ?
                    nbor->branches[6-branchIdx] :
                    nbor->branches[8-branchIdx];
            }
            break;

        case UE:
            if (branchIdx == 0 || branchIdx == 2)
                return base->branches[branchIdx+5];
            if (branchIdx == 5 || branchIdx == 7) {
                nbor = base->getNeighborGeqSize(dir);
                if (nbor == nullptr) return nbor;
                if (nbor->isNodeType<Leaf>()) return nbor;
                return nbor->branches[branchIdx-5];
            } else {
                nbor = (branchIdx == 4 || branchIdx == 6) ?
                    base->getNeighborGeqSize(U) :
                    base->getNeighborGeqSize(E);
                if (nbor == nullptr) return nbor;
                if (nbor->isNodeType<Leaf>()) return nullptr;
                return (branchIdx == 1 || branchIdx == 4) ?
                    nbor->branches[5-branchIdx] :
                    nbor->branches[9-branchIdx];
            }
            break;

        case DW:
            if (branchIdx == 5 || branchIdx == 7)
                return base->branches[branchIdx-5];
            if (branchIdx == 0 || branchIdx == 2) {
                nbor = base->getNeighborGeqSize(dir);
                if (nbor == nullptr) return nbor;
                if (nbor->isNodeType<Leaf>()) return nbor;
                return nbor->branches[branchIdx+5];
            } else {
                nbor = (branchIdx == 1 || branchIdx == 3) ?
                    base->getNeighborGeqSize(D) :
                    base->getNeighborGeqSize(W);
                if (nbor == nullptr) return nbor;
                if (nbor->isNodeType<Leaf>()) return nullptr;
                return (branchIdx == 1 || branchIdx == 4) ?
                    nbor->branches[5-branchIdx] :
                    nbor->branches[9-branchIdx];
            }
            break;

        case UW:
            if (branchIdx == 1 || branchIdx == 3)
                return base->branches[branchIdx+3];
            if (branchIdx == 4 || branchIdx == 6) {
                nbor = base->getNeighborGeqSize(dir);
                if (nbor == nullptr) return nbor;
                if (nbor->isNodeType<Leaf>()) return nbor;
                return nbor->branches[branchIdx-3];
            } else {
                nbor = (branchIdx == 5 || branchIdx == 7) ?
                    base->getNeighborGeqSize(U) :
                    base->getNeighborGeqSize(W);
                if (nbor == nullptr) return nbor;
                if (nbor->isNodeType<Leaf>()) return nullptr;
                return (branchIdx == 0 || branchIdx == 5) ?
                    nbor->branches[5-branchIdx] :
                    nbor->branches[9-branchIdx];
            }
            break;

        case DE:
            if (branchIdx == 4 || branchIdx == 6)
                return base->branches[branchIdx-3];
            if (branchIdx == 1 || branchIdx == 3) {
                nbor = base->getNeighborGeqSize(dir);
                if (nbor == nullptr) return nbor;
                if (nbor->isNodeType<Leaf>()) return nbor;
                return nbor->branches[branchIdx+3];
            } else {
                nbor = (branchIdx == 0 || branchIdx == 2) ?
                    base->getNeighborGeqSize(D) :
                    base->getNeighborGeqSize(E);
                if (nbor == nullptr) return nbor;
                if (nbor->isNodeType<Leaf>()) return nullptr;
                return (branchIdx == 0 || branchIdx == 5) ?
                    nbor->branches[5-branchIdx] :
                    nbor->branches[9-branchIdx];
            }
            break;

        case NE:
            if (branchIdx == 0 || branchIdx == 4)
                return base->branches[branchIdx+3];
            if (branchIdx == 3 || branchIdx == 7) {
                nbor = base->getNeighborGeqSize(dir);
                if (nbor == nullptr) return nbor;
                if (nbor->isNodeType<Leaf>()) return nbor;
                return nbor->branches[branchIdx-3];
            } else {
                nbor = (branchIdx == 2 || branchIdx == 6) ?
                    base->getNeighborGeqSize(N) :
                    base->getNeighborGeqSize(E);
                if (nbor == nullptr) return nbor;
                if (nbor->isNodeType<Leaf>()) return nullptr;
                return (branchIdx == 1 || branchIdx == 2) ?
                    nbor->branches[3-branchIdx] :
                    nbor->branches[11-branchIdx];
            }
            break;

        case SW:
            if (branchIdx == 3 || branchIdx == 7)
                return base->branches[branchIdx-3];
            if (branchIdx == 0 || branchIdx == 4) {
                nbor = base->getNeighborGeqSize(dir);
                if (nbor == nullptr) return nbor;
                if (nbor->isNodeType<Leaf>()) return nbor;
                return nbor->branches[branchIdx+3];
            } else {
                nbor = (branchIdx == 1 || branchIdx == 5) ?
                    base->getNeighborGeqSize(S) :
                    base->getNeighborGeqSize(W);
                if (nbor == nullptr) return nbor;
                if (nbor->isNodeType<Leaf>()) return nullptr;
                return (branchIdx == 1 || branchIdx == 2) ?
                    nbor->branches[3-branchIdx] :
                    nbor->branches[11-branchIdx];
            }
            break;

        case NW:
            if (branchIdx == 1 || branchIdx == 5)
                return base->branches[branchIdx+1];
            if (branchIdx == 2 || branchIdx == 6) {
                nbor = base->getNeighborGeqSize(dir);
                if (nbor == nullptr) return nbor;
                if (nbor->isNodeType<Leaf>()) return nbor;
                return nbor->branches[branchIdx-1];
            } else {
                nbor = (branchIdx == 3 || branchIdx == 7) ?
                    base->getNeighborGeqSize(N) :
                    base->getNeighborGeqSize(W);
                if (nbor == nullptr) return nbor;
                if (nbor->isNodeType<Leaf>()) return nullptr;
                return (branchIdx == 0 || branchIdx == 3) ?
                    nbor->branches[3-branchIdx] :
                    nbor->branches[11-branchIdx];
            }
            break;

        case SE:
            if (branchIdx == 2 || branchIdx == 6)
                return base->branches[branchIdx-1];
            if (branchIdx == 1 || branchIdx == 5) {
                nbor = base->getNeighborGeqSize(dir);
                if (nbor == nullptr) return nbor;
                if (nbor->isNodeType<Leaf>()) return nbor;
                return nbor->branches[branchIdx+1];
            } else {
                nbor = (branchIdx == 0 || branchIdx == 4) ?
                    base->getNeighborGeqSize(S) :
                    base->getNeighborGeqSize(E);
                if (nbor == nullptr) return nbor;
                if (nbor->isNodeType<Leaf>()) return nullptr;
                return (branchIdx == 0 || branchIdx == 3) ?
                    nbor->branches[3-branchIdx] :
                    nbor->branches[11-branchIdx];
            }
            break;

        case UNE:
            if (branchIdx == 0)
                return base->branches[7];
            else {
                nbor = [this] {
                    switch (branchIdx) {
                        case 1: return base->getNeighborGeqSize(E);
                        case 2: return base->getNeighborGeqSize(N);
                        case 3: return base->getNeighborGeqSize(NE);
                        case 4: return base->getNeighborGeqSize(U);
                        case 5: return base->getNeighborGeqSize(UE);
                        case 6: return base->getNeighborGeqSize(UN);
                        case 7: return base->getNeighborGeqSize(UNE);
                    } } ();
                    if (nbor == nullptr) return nbor;
                    if (nbor->isNodeType<Leaf>()) return nullptr;
                    return nbor->branches[7-branchIdx];
            }
            break;

        case DSW:
            if (branchIdx == 7)
                return base->branches[0];
            else {
                nbor = [this] {
                    switch (branchIdx) {
                        case 0: return base->getNeighborGeqSize(DSW);
                        case 1: return base->getNeighborGeqSize(DS);
                        case 2: return base->getNeighborGeqSize(DW);
                        case 3: return base->getNeighborGeqSize(D);
                        case 4: return base->getNeighborGeqSize(SW);
                        case 5: return base->getNeighborGeqSize(S);
                        case 6: return base->getNeighborGeqSize(W);
                    } } ();
                    if (nbor == nullptr) return nbor;
                    if (nbor->isNodeType<Leaf>()) return nullptr;
                    return nbor->branches[7-branchIdx];
            }
            break;

        case UNW:
            if (branchIdx == 1)
                return base->branches[6];
            else {
                nbor = [this] {
                    switch (branchIdx) {
                        case 0: return base->getNeighborGeqSize(W);
                        case 2: return base->getNeighborGeqSize(NW);
                        case 3: return base->getNeighborGeqSize(N);
                        case 4: return base->getNeighborGeqSize(UW);
                        case 5: return base->getNeighborGeqSize(U);
                        case 6: return base->getNeighborGeqSize(UNW);
                        case 7: return base->getNeighborGeqSize(UN);
                    } } ();
                    if (nbor == nullptr) return nbor;
                    if (nbor->isNodeType<Leaf>()) return nullptr;
                    return nbor->branches[7-branchIdx];
            }
            break;

        case DSE:
            if (branchIdx == 6)
                return base->branches[1];
            else {
                nbor = [this] {
                    switch (branchIdx) {
                        case 0: return base->getNeighborGeqSize(DS);
                        case 1: return base->getNeighborGeqSize(DSE);
                        case 2: return base->getNeighborGeqSize(D);
                        case 3: return base->getNeighborGeqSize(DE);
                        case 4: return base->getNeighborGeqSize(S);
                        case 5: return base->getNeighborGeqSize(SE);
                        case 7: return base->getNeighborGeqSize(E);
                    } } ();
                    if (nbor == nullptr) return nbor;
                    if (nbor->isNodeType<Leaf>()) return nullptr;
                    return nbor->branches[7-branchIdx];
            }
            break;

        case USE:
            if (branchIdx == 2)
                return base->branches[5];
            else {
                nbor = [this] {
                    switch (branchIdx) {
                        case 0: return base->getNeighborGeqSize(S);
                        case 1: return base->getNeighborGeqSize(SE);
                        case 3: return base->getNeighborGeqSize(E);
                        case 4: return base->getNeighborGeqSize(US);
                        case 5: return base->getNeighborGeqSize(USE);
                        case 6: return base->getNeighborGeqSize(U);
                        case 7: return base->getNeighborGeqSize(UE);
                    } } ();
                    if (nbor == nullptr) return nbor;
                    if (nbor->isNodeType<Leaf>()) return nullptr;
                    return nbor->branches[7-branchIdx];
            }
            break;

        case DNW:
            if (branchIdx == 5)
                return base->branches[2];
            else {
                nbor = [this] {
                    switch (branchIdx) {
                        case 0: return base->getNeighborGeqSize(DW);
                        case 1: return base->getNeighborGeqSize(D);
                        case 2: return base->getNeighborGeqSize(DNW);
                        case 3: return base->getNeighborGeqSize(DN);
                        case 4: return base->getNeighborGeqSize(W);
                        case 6: return base->getNeighborGeqSize(NW);
                        case 7: return base->getNeighborGeqSize(N);
                    } } ();
                    if (nbor == nullptr) return nbor;
                    if (nbor->isNodeType<Leaf>()) return nullptr;
                    return nbor->branches[7-branchIdx];
            }
            break;

        case USW:
            if (branchIdx == 3)
                return base->branches[4];
            else {
                nbor = [this] {
                    switch (branchIdx) {
                        case 0: return base->getNeighborGeqSize(SW);
                        case 1: return base->getNeighborGeqSize(S);
                        case 2: return base->getNeighborGeqSize(W);
                        case 4: return base->getNeighborGeqSize(USW);
                        case 5: return base->getNeighborGeqSize(US);
                        case 6: return base->getNeighborGeqSize(UW);
                        case 7: return base->getNeighborGeqSize(U);
                    } } ();
                    if (nbor == nullptr) return nbor;
                    if (nbor->isNodeType<Leaf>()) return nullptr;
                    return nbor->branches[7-branchIdx];
            }
            break;

        case DNE:
            if (branchIdx == 4)
                return base->branches[3];
            else {
                nbor = [this] {
                    switch (branchIdx) {
                        case 0: return base->getNeighborGeqSize(D);
                        case 1: return base->getNeighborGeqSize(DE);
                        case 2: return base->getNeighborGeqSize(DN);
                        case 3: return base->getNeighborGeqSize(DNE);
                        case 5: return base->getNeighborGeqSize(E);
                        case 6: return base->getNeighborGeqSize(N);
                        case 7: return base->getNeighborGeqSize(NE);
                    } } ();
                    if (nbor == nullptr) return nbor;
                    if (nbor->isNodeType<Leaf>()) return nullptr;
                    return nbor->branches[7-branchIdx];
            }
            break;
        default:
            throw std::runtime_error("Invalid direction");
            break;
    }
}

void Node::buildNearNeighbors() {
    assert(!isRoot());
    for (int i = 0; i < numDir; ++i) {
        Dir dir = static_cast<Dir>(i);
        auto nbor = getNeighborGeqSize(dir);
        if (nbor != nullptr) {
            nbors.push_back(nbor);
            // dirNbors[i] = nbor;
        }
    }
    assert(nbors.size() <= numDir);
}

void Node::buildInteractionList() {
    assert(!isRoot());

    auto baseNbors = base->getNearNeighbors();

    for (const auto& baseNbor : baseNbors)
        if (baseNbor->isNodeType<Leaf>() && !contains<std::shared_ptr<Node>>(nbors, baseNbor))
            iList.push_back(baseNbor);
        else 
            for (const auto& branch : baseNbor->branches) 
                if (!contains<std::shared_ptr<Node>>(nbors, branch))
                    iList.push_back(branch);

    const int maxSize = constPow(6, DIM) - constPow(3, DIM);
    assert(iList.size() <= maxSize);
}

void Node::buildDirectedIList() {
    assert(!isRoot());
    auto contains = [](NodeVec& vec, std::shared_ptr<Node> val) {
        return std::find(vec.begin(), vec.end(), val) != vec.end();
    };

    //auto arrContains = []<int N>(std::array<NodeVec,N>& vecArr, std::shared_ptr<Node> val) {
    //    return std::find(vec.begin(), vec.end(), val) != vec.end();
    //    };

    auto baseNbors = base->getNearNeighbors();

    for (const auto& baseNbor : baseNbors) {
        auto center0 = baseNbor->getCenter();

        if (baseNbor->isNodeType<Leaf>() && !contains(nbors, baseNbor))
            leafIlist.push_back(baseNbor); // list 4
        else {
            for (const auto& branch : baseNbor->branches) {
                auto center1 = branch->getCenter();

                if (center0[2] > center[2]) { 
                    if (!contains(nbors, branch)
                        && (center1[2] - center[2]) > nodeLeng)
                        dirList[0].push_back(branch); // uplist
                } else if (center0[2] < center[2]) { 
                    if (!contains(nbors, branch)
                        && (center1[2] - center[2]) < -nodeLeng)
                        dirList[1].push_back(branch); // downlist
                }

                if (center0[1] > center[1]) { 
                    if (!contains(nbors, branch)
                        && !contains(dirList[0], branch)
                        && !contains(dirList[1], branch)
                        && (center1[1] - center[1]) > nodeLeng)
                        dirList[2].push_back(branch); // northlist
                } else if (center0[1] < center[1]) { 
                    if (!contains(nbors, branch)
                        && !contains(dirList[0], branch)
                        && !contains(dirList[1], branch)
                        && (center1[1] - center[1]) < -nodeLeng)
                        dirList[3].push_back(branch); // southlist
                }

                if (center0[0] > center[0]) { 
                    if (!contains(nbors, branch)
                        && !contains(dirList[0], branch)
                        && !contains(dirList[1], branch)
                        && !contains(dirList[2], branch)
                        && !contains(dirList[3], branch)
                        && (center1[0] - center[0]) > nodeLeng)
                        dirList[4].push_back(branch); // eastlist
                } else if (center0[0] < center[0]) { 
                    if (!contains(nbors, branch)
                        && !contains(dirList[0], branch)
                        && !contains(dirList[1], branch)
                        && !contains(dirList[2], branch)
                        && !contains(dirList[3], branch)
                        && (center1[0] - center[0]) < -nodeLeng)
                        dirList[5].push_back(branch); // westlist
                } else
                    throw std::runtime_error("Invalid interaction node");
            }
        }
    }
    
    std::size_t dirListSize = std::accumulate(
        dirList.begin(), dirList.end(), std::size_t{0},
        [](std::size_t sum, const auto& vec) {
            // std::cout << vec.size() << ' ';
            return sum + vec.size();
        }
    );
    // std::cout << dirListSize << ' ' << leafIlist.size() << ' ' << iList.size() << '\n';
    assert(dirListSize + leafIlist.size() == iList.size());
}

// no rotation matrices
void Node::buildLocalCoeffsFromLeafIlist() {
    for (const auto& iNode : leafIlist) {
        auto mpoleCoeffs( iNode->getMpoleCoeffs() );
        auto dR = toSph(iNode->getCenter() - center);
        double r = dR[0], th = dR[1], ph = dR[2];

        for (int j = 0; j <= order; ++j) {
            for (int k = -j; k <= j; ++k) {
                int k_ = k + j;

                for (int n = 0; n <= order; ++n) {
                    for (int m = -n; m <= n;  ++m) {
                        int m_ = m + n;

                        localCoeffs[j][k_] +=
                            mpoleCoeffs[n][m_] * pow(iu, abs(k-m)-abs(k)-abs(m))
                            * tables.A_[n][m_] * tables.A_[j][k_] / tables.A_[j+n][m-k+j+n]
                            * legendreLM(th, pair2i(j+n, abs(m-k))) * expI(static_cast<double>(m-k)*ph)
                            / ( pm(n) * pow(r, j+n+1) );
                    }
                }
            }
        }
    }
}

// no rotation matrices
/*const std::vector<vecXcd> Node::getShiftedLocalCoeffs(const int branchIdx) {

    std::vector<vecXcd> shiftedCoeffs;
    for (int l = 0; l <= order; ++l)
        shiftedCoeffs.emplace_back(vecXcd::Zero(2*l+1));

    auto dR = toSph(branches[branchIdx]->getCenter() - center);
    double r = dR[0], th = dR[1], ph = dR[2];

    for (int j = 0; j <= order; ++j) {
        for (int k = -j; k <= j; ++k) {
            int k_ = k + j;

            for (int n = j; n <= order; ++n) {
                for (int m = std::max(j-n+k,-n); m <= std::min(n-j+k,n);  ++m) {
                    int m_ = m + n;

                    shiftedCoeffs[j][k_] +=
                        localCoeffs[n][m_] * pow(iu, abs(m)-abs(k)-abs(m-k))
                        * tables.A_[n-j][m-k+n-j] * tables.A_[j][k_] / tables.A_[n][m_]
                        * legendreLM(th, pair2i(n-j, abs(m-k))) * expI(static_cast<double>(m-k)*ph)
                        * pow(r, n-j) / pm(n+j);
                        
                }
            }
        }
    }

    return shiftedCoeffs;
}*/

// rotation matrices
const std::vector<vecXcd> Node::getShiftedLocalCoeffs(const int branchIdx) {
    std::vector<vecXcd> shiftedLocalCoeffs;
    double r = (branches[branchIdx]->getCenter() - center).norm();

    for (int j = 0; j <= order; ++j) {
        shiftedLocalCoeffs.emplace_back(vecXcd::Zero(2*j+1));

        // apply rotation (rotation axis is opposite from M2M)
        localCoeffs[j] = wignerD[7-branchIdx][j] * localCoeffs[j];

        vecXcd shiftedLocalCoeffs_j = vecXcd::Zero(2*j+1);
        for (int k = -j; k <= j; ++k) {
            int k_j = k + j;

            for (int n = j; n <= order; ++n) {
                shiftedLocalCoeffs_j[k_j] +=
                    localCoeffs[n][k+n] *
                    tables.A_[n-j][n-j] * tables.A_[j][k_j] / tables.A_[n][k+n] *
                    pow(r, n-j) / pm(n+j);
                //for (int m = -n; m <= n; ++m) {
                //    if (m != 0) continue;
                //    int m_ = m + n;

                //    shiftedLocalCoeffs_j[k_] +=
                //        localCoeffs[n][m_] * pow(iu, abs(m)-abs(k)-abs(k-m)) *
                //        tables.A[n-j][m_-k_] * tables.A[j][k_] / tables.A[n][m_] *
                //        legendreLM(0, n-j, abs(m-k)) *
                //        pow(r, n-j) / pm(n+j);
                //}
            }
        }

        // apply inverse rotation
        shiftedLocalCoeffs[j] += wignerDInv[7-branchIdx][j] * shiftedLocalCoeffs_j;
    }
    return shiftedLocalCoeffs;
}

const double Node::getDirectPhi(const vec3d& X) {
    double phi = 0;
    for (const auto& src : particles) {
        auto r = (X - src->getPos()).norm();
        if (r > 1.0E-9)
            phi += src->getCharge() / r;
    }
    return phi;
}

const realVec Node::getDirectPhis() {
    realVec phis;

    for (const auto& obs : particles) {
        phis.push_back(getDirectPhi(obs->getPos()));
    }
    return phis;
}

const std::vector<vec3d> Node::getDirectFlds() {
    std::vector<vec3d> flds;

    return flds;
}