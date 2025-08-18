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
    auto wignerDAlongDir = [](const pair2d angles, bool isInv) {
        std::vector<matXcd> mats;
        for (int l = 0; l <= order; ++l)
            mats.push_back(
                isInv ?
                wignerD_l(angles, l) :
                wignerD_l(angles, l).adjoint()
            );
        return mats;
    };

    for (int dir = 0; dir < 14; ++dir) {
        auto X = 
            dir < 8 ? 
            idx2pm(dir) :   // diagonal directions (for M2M and L2L)
            [dir] { 
            switch (dir-8) {
                case 0: return vec3d(0, 0, 1);
                case 1: return vec3d(0, 0, -1);
                case 2: return vec3d(0, 1, 0);
                case 3: return vec3d(0, -1, 0);
                case 4: return vec3d(1, 0, 0);
                case 5: return vec3d(-1, 0, 0);
            } 
            }();            // cardinal directions (for M2X and X2L)

        auto R = toSph(X);
        pair2d angles(R[1], R[2]);

        wignerD[dir] = wignerDAlongDir(angles,0);
        wignerDInv[dir] = wignerDAlongDir(angles,1);
        if (dir >=8)
            rotMatR[dir-8] = rotationR(angles);
    }
}

// theta part of Ylm
const double Node::legendreCos(const double th, const int l, const int abs_m) {
    assert(0 <= abs_m && abs_m <= l);
    const double cos_th = cos(th), sin_th = sin(th);

    double legendreSum = 0.0;
    for (int k = l; k >= abs_m; k -= 2) // expression is zero for l-k odd
        legendreSum +=  
            tables.fallingFact_[k][abs_m] * tables.legendreSum_[l][k] 
            * pow(cos_th, k-abs_m);

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
    // preallocate exp coeffs
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
        if (nbor != nullptr)
            nbors.push_back(nbor);
    }

    assert(nbors.size() <= numDir);
}

void Node::buildInteractionList() {
    assert(!isRoot());
    assert(nbors.size());

    auto notContains = [](NodeVec& vec, std::shared_ptr<Node> ele) {
        return std::find(vec.begin(), vec.end(), ele) == vec.end();
        };

    NodeVec iList;
    auto baseNbors = base->getNearNeighbors();

    for (const auto& baseNbor : baseNbors)
        if (baseNbor->isNodeType<Leaf>() && notContains(nbors, baseNbor))
            leafIlist.push_back(baseNbor); // list 4
        else {
            const auto center0 = base->getCenter();
            const auto maxDist = base->getLeng();
            for (const auto& branch : baseNbor->branches)
                if (notContains(nbors, branch) && 
                    (branch->getCenter()-center0).lpNorm<Eigen::Infinity>() < maxDist )
                    iList.push_back(branch);
        }

    assert(iList.size() <= pow(4, DIM) - pow(3, DIM));

    // assign each interaction node to a dirlist
    // pick minDist \in (nodeLeng, 2.0*nodeLeng) to avoid rounding errors
    const double minDist = 1.5 * nodeLeng;

    for (const auto& iNode : iList) {
        auto center0 = iNode->getCenter();

        if (center0[2] - center[2] > minDist) // uplist
            dirList[0].push_back(iNode);
        else if (center[2] - center0[2] > minDist) // downlist
            dirList[1].push_back(iNode);
        else if (center0[1] - center[1] > minDist) // northlist
            dirList[2].push_back(iNode);
        else if (center[1] - center0[1] > minDist) // southlist
            dirList[3].push_back(iNode);
        else if (center0[0] - center[0] > minDist) // eastlist
            dirList[4].push_back(iNode);
        else if (center[0] - center0[0] > minDist) // westlist
            dirList[5].push_back(iNode);
        else
            throw std::runtime_error("Invalid interaction node");
    }
}

void Node::buildOuterInteractionList() {
    assert(!isRoot());

    const double minDist = nodeLeng;
    for (const auto& nbor : nbors) {
        if (nbor->isNodeType<Leaf>()) continue;
            
        auto nodes = nbor->getBranches();
        for (const auto& node : nodes) {
            auto center0 = node->getCenter();
            if (center0[2] - center[2] > minDist) // uplist
                outerDirList[0].push_back(node);
            else if (center[2] - center0[2] > minDist) // downlist
                outerDirList[1].push_back(node);
            else if (center0[1] - center[1] > minDist) // northlist
                outerDirList[2].push_back(node);
            else if (center[1] - center0[1] > minDist) // southlist
                outerDirList[3].push_back(node);
            else if (center0[0] - center[0] > minDist) // eastlist
                outerDirList[4].push_back(node);
            else if (center[0] - center0[0] > minDist) // westlist
                outerDirList[5].push_back(node);
        }
    }

    auto dirListSize =
        std::accumulate(outerDirList.begin(), outerDirList.end(), size_t{ 0 },
            [](size_t sum, const auto& vec) { return sum + vec.size(); });
    assert(dirListSize <= pow(6, DIM) - pow(4, DIM));
}

void Node::buildLocalCoeffsFromLeafIlist() {
    // if # particles is small, evaluate sol at particles directly
    // if (particles.size() <= order*order) {
        for (const auto& obs : particles)
            for (const auto& iNode : leafIlist) {
                obs->addToPhi(iNode->getDirectPhi(obs->getPos()));
                obs->addToFld(iNode->getDirectFld(obs->getPos()));
            }
        return;
    // }

    // std::cout << "Computing local coeffs from list 4 of size " << leafIlist.size() << '\n';
    // how to precompute rotation matrices for list 4?
    /*for (const auto& iNode : leafIlist) {
        auto mpoleCoeffs = iNode->getMpoleCoeffs();
        auto dR = toSph(iNode->getCenter() - center);
        double r = dR[0], th = dR[1], ph = dR[2];

        for (int j = 0; j <= order; ++j) {
            // double r2jpn = pow(r,j);
            for (int k = -j; k <= j; ++k) {
                int k_j = k + j;

                for (int n = 0; n <= order; ++n) {
                    for (int m = -n; m <= n;  ++m) {
                        int m_n = m + n;

                        localCoeffs[j][k_j] +=
                            mpoleCoeffs[n][m_n] * pow(iu, abs(k-m)-abs(k)-abs(m))
                            * tables.A_[n][m_n] * tables.A_[j][k_j] / tables.A_[j+n][m-k+j+n]
                            * legendreCos(th, j+n, abs(m-k)) * expI(static_cast<double>(m-k)*ph)
                            / ( pm(n) * pow(r, j+n+1) );
                    }
                }
            }
        }
    }*/
}

const std::vector<vecXcd> Node::getShiftedLocalCoeffs(const int branchIdx) const {
    std::vector<vecXcd> shiftedLocalCoeffs, rotatedLocalCoeffs;
    const double r = (branches[branchIdx]->getCenter()-center).norm();

    // apply rotation (rotation axis is opposite from M2M)
    for (int j = 0; j <= order; ++j)
        rotatedLocalCoeffs.push_back(wignerD[7-branchIdx][j] * localCoeffs[j]);

    for (int j = 0; j <= order; ++j) {
        vecXcd shiftedLocalCoeffs_j = vecXcd::Zero(2*j+1);
        for (int k = -j; k <= j; ++k) {
            int k_j = k + j;
            double r2nmj = 1.0;
            for (int n = j; n <= order; ++n) {
                shiftedLocalCoeffs_j[k_j] +=
                    rotatedLocalCoeffs[n][k+n] *
                    tables.A_[n-j][n-j] * tables.A_[j][k_j] / tables.A_[n][k+n]
                    * r2nmj / pm(n+j);
                r2nmj *= r;
            }
        }

        // apply inverse rotation
        shiftedLocalCoeffs.push_back(wignerDInv[7-branchIdx][j] * shiftedLocalCoeffs_j);
    }
    return shiftedLocalCoeffs;
}

// return phi at X due to all particles in this node
const double Node::getDirectPhi(const vec3d& X) {
    double phi = 0;
    for (const auto& src : particles) {
        auto r = (X - src->getPos()).norm();
        if (r > 1.0E-9)
            phi += src->getCharge() / r;
    }
    return phi;
}

// return phi at all particles in this node due to all other particles in this node
const realVec Node::getDirectPhis() {
    realVec phis;

    for (const auto& obs : particles) {
        phis.push_back(getDirectPhi(obs->getPos()));
    }
    return phis;
}

// return fld at X due to all particles in this node
const vec3d Node::getDirectFld(const vec3d& X) {
    vec3d fld = vec3d::Zero();
    for (const auto& src : particles) {
        auto dX = X - src->getPos();
        auto r = dX.norm();
        if (r > 1.0E-9)
            fld += src->getCharge() * dX / pow(r, 3);
    }
    return fld;
}

// return phi at all particles in this node due to all other particles in this node
const std::vector<vec3d> Node::getDirectFlds() {
    std::vector<vec3d> flds;

    for (const auto& obs : particles) {
        flds.push_back(getDirectFld(obs->getPos()));
    }
    return flds;
}