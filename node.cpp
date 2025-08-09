#include "node.h"

using enum Dir;

const int numDir = std::pow(3, DIM) - 1;

int Node::order;
int Node::orderExp;
int Node::maxNodeParts;
double Node::rootLeng;
Tables Node::tables;
std::vector<matXcdVec> Node::rotationMat;
std::vector<matXcdVec> Node::rotationInvMat;

void Node::setNodeParams(const Config& config) {
    order = ceil(-std::log(config.EPS) / std::log(2));
    maxNodeParts = config.maxNodeParts;
    rootLeng = config.L;
}

void Node::buildTables(const Config& config) {
    tables = Tables(order, Precision::LOW); // config.prec;
    orderExp = tables.quadCoeffs_.size();
}

void Node::buildRotationMats() {
    // do cardinal directions

    for (int dir = 18; dir < numDir; ++dir) {
        rotationMat.push_back(rotationMatrixAlongDir(dir, 0));
        rotationInvMat.push_back(rotationMatrixAlongDir(dir, 1));
    }
}

 matXcdVec Node::rotationMatrixAlongDir(int dir, const bool isInv) {
    vec3d R;
    if (dir >= 18) R = toSph(idx2pm(dir-=18));
    else if (dir < 8) {
        // do cardinal direction
    } else
        throw std::runtime_error("Invalid rotation direction");

    pair2d angles(R[1], R[2]);
        
    matXcdVec mats;
    for (int l = 0; l <= order; ++l)
        mats.push_back(isInv ?
            wignerD_l(angles, l) :
            wignerD_l(angles, l).adjoint()
        );

    return mats;
}

// Ylm except constant coefficients and \exp(i m \phi) phase factor
const double Node::legendreLM(const double th, const pair2i& lm) {
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

//const cmplx Node::sphHarmonic(const double th, const double ph, int l, int m) {
//    return legendreLM(th,l,std::abs(m)) * std::exp(iu*static_cast<double>(m)*ph);
//}

Node::Node(
    const ParticleVec& particles, 
    const int branchIdx, 
    Node* const base)
    : particles(particles), branchIdx(branchIdx), base(base),
    nodeLeng(base == nullptr ? rootLeng : base->getLeng() / 2.0),
    center(base == nullptr ? zeroVec :
        base->getCenter() + nodeLeng / 2.0 * idx2pm(branchIdx) ),
    nodeStat(0), useRot(0)
{
};

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
    for (int i = 0; i < numDir; ++i) {
        Dir dir = static_cast<Dir>(i);
        auto nbor = getNeighborGeqSize(dir);
        if (nbor != nullptr)
            nbors.push_back(nbor);
    }
    assert(nbors.size() <= numDir);
}

void Node::buildInteractionList() {
    assert( !isRoot() );

    buildNearNeighbors(); // comment out for iList test
    // base->buildNearNeighbors(); // uncomment for iList test

    auto baseNbors = base->getNearNeighbors();

    for (const auto& baseNbor : baseNbors)
        if (baseNbor->isNodeType<Leaf>() && !contains<std::shared_ptr<Node>>(nbors,baseNbor))
            iList.push_back(baseNbor);
        else 
            for (const auto& branch : baseNbor->branches) 
                if (!contains<std::shared_ptr<Node>>(nbors, branch))
                    iList.push_back(branch);

    const int maxSize = constPow(6, DIM) - constPow(3, DIM);
    // std::cout << iList.size() << '\n';
    assert(iList.size() <= maxSize);
}

// no rotation matrices, no exponential expansions
void Node::buildMpoleToLocalCoeffs() {
    if ( isRoot() ) return;

    buildInteractionList();
    
    for (int l = 0; l <= order; ++l)
        localCoeffs.emplace_back(vecXcd::Zero(2*l+1));

    for (const auto& iNode : iList) {
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

    if (!base->isRoot()) 
        for (int l = 0; l <= order; ++l)
            localCoeffs[l] += (base->getShiftedLocalCoeffs(branchIdx))[l];
    // iList.clear();
}

// no rotation matrices
/*const vecXcdVec Node::getShiftedLocalCoeffs(const int branchIdx) {

    // std::cout << "Shifting local coeffs" << '\n';

    vecXcdVec shiftedCoeffs;
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
                        * tables.A[n-j][m-k+n-j] * tables.A[j][k_] / tables.A[n][m_]
                        * legendreLM(th, n-j, abs(m-k)) * expI(static_cast<double>(m-k)*ph)
                        * pow(r, n-j) / pm(n+j);
                        
                }
            }
        }
    }

    return shiftedCoeffs;
}*/

// rotation matrices
const vecXcdVec Node::getShiftedLocalCoeffs(const int branchIdx) {
    vecXcdVec shiftedCoeffs;
    double r = (branches[branchIdx]->getCenter() - center).norm();

    for (int j = 0; j <= order; ++j) {
        shiftedCoeffs.emplace_back(vecXcd::Zero(2*j+1));

        // apply rotation (rotation axis is opposite from mpole2mpole)
        localCoeffs[j] = rotationMat[7-branchIdx][j] * localCoeffs[j];

        vecXcd shiftedLocalCoeffs_j = vecXcd::Zero(2*j+1);
        for (int k = -j; k <= j; ++k) {
            int k_ = k + j;

            for (int n = j; n <= order; ++n) {
                shiftedLocalCoeffs_j[k_] +=
                    localCoeffs[n][k+n] *
                    tables.A_[n-j][n-j] * tables.A_[j][k_] / tables.A_[n][k+n] *
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
        shiftedCoeffs[j] += rotationInvMat[7-branchIdx][j] * shiftedLocalCoeffs_j;
    }

    return shiftedCoeffs;
}

const double Node::getDirectPhi(const vec3d& X) {
    double phi = 0;
    for (const auto& p : particles) {
        auto r = (X - p->getPos()).norm();
        if (r > 1.0E-9)
            phi += p->getCharge() / r;
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

const vec3dVec Node::getDirectFlds() {
    vec3dVec flds;

    return flds;
}