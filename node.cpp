#include "node.h"

using enum Dir;

const int numDir = std::pow(3, DIM) - 1;

int Node::order;
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

void Node::buildTables() {
    auto binom = [](double x, int k) {
        return fallingFactorial(x, k) / factorial(k);
        };

    for (int l = 0; l <= order; ++l) {
        realVec coeffYlm_l, fallingFact_l, legendreSum_l, A_l;

        for (int m = 0; m <= l; ++m) {
            coeffYlm_l.push_back( coeffYlm(l, m) );
            fallingFact_l.push_back( fallingFactorial(l, m) );
            legendreSum_l.push_back( binom(l, m) * binom((l+m-1)/2.0, l) );
        }

        auto pm_l = pm(l);
        for (int m = -l; m <= l; ++m)
            A_l.push_back(pm_l /
                std::sqrt(static_cast<double>(factorial(l-m)*factorial(l+m))) );

        tables.coeffYlm.push_back(coeffYlm_l);
        tables.fallingFact.push_back(fallingFact_l);
        tables.legendreSum.push_back(legendreSum_l);
        tables.A.push_back(A_l);

        //std::cout << "l = " << l << '\n';
        //for (int m = 0; m <= l; ++m)
        //    std::cout << tables.legendreSum[l][m] << ' ';
        //std::cout << '\n';
    }
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
const double Node::legendreLM(const double th, const int l, const int abs_m) {
    assert(abs_m <= l);

    const auto cos_th = cos(th);
    const auto sin_th = sin(th);
    double legendreSum = 0.0;

    // term is zero for l-k odd
    for (int k = l; k >= abs_m; k -= 2)
        legendreSum += tables.fallingFact[k][abs_m] * tables.legendreSum[l][k] * pow(cos_th, k-abs_m);

    return tables.coeffYlm[l][abs_m] * pow(sin_th, abs_m) * legendreSum;
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

    // buildNearNeighbors(); // comment out for iList test
    base->buildNearNeighbors(); // uncomment for iList test
    auto baseNbors = base->getNearNeighbors();

    for (const auto& baseNbor : baseNbors)
        if (baseNbor->isNodeType<Leaf>() && !contains<std::shared_ptr<Node>>(nbors,baseNbor))
            iList.push_back(baseNbor);
        else 
            for (const auto& branch : baseNbor->branches) 
                if (!contains<std::shared_ptr<Node>>(nbors, branch))
                    iList.push_back(branch);

    const int maxSize = constPow(6, DIM) - constPow(3, DIM);
    assert(iList.size() <= maxSize);
}

void Node::buildMpoleToLocalCoeffs() {
    if ( isRoot() ) return;

    /*buildInteractionList();
    localCoeffs.resize(order+1);

    for (const auto& iNode : iList) {
        auto mpoleCoeffs( iNode->getMpoleCoeffs() );
        auto dz( iNode->getCenter() - center );
        auto mdz2k( vec3d(1,0) );

        vec3dVec innerCoeffs; // innerCoeffs[k] = mpoleCoeffs[k] / (-dz)^k
        for (size_t k = 0; k <= order; ++k) {
            innerCoeffs.push_back(mpoleCoeffs[k] / mdz2k);
            mdz2k *= -dz;
        }

        localCoeffs[0] += innerCoeffs[0] * std::log(-dz);
        for (size_t k = 1; k <= order; ++k)
            localCoeffs[0] += innerCoeffs[k];

        auto dz2l = dz;
        for (size_t l = 1; l <= order; ++l) {
            localCoeffs[l] -= innerCoeffs[0] / (static_cast<double>(l) * dz2l);
            for (size_t k = 1; k <= order; ++k)
                localCoeffs[l] += innerCoeffs[k] / dz2l
                                    * static_cast<double>(binomTable[k+l-1][k-1]);
            dz2l *= dz;
        }
    }

    if (!base->isRoot()) localCoeffs += base->getShiftedLocalCoeffs(center);
    iList.clear();*/
}

const vecXcdVec Node::getShiftedLocalCoeffs(const vec3d z0) {
    auto shiftedCoeffs( localCoeffs );
    //for (size_t j = 0; j <= order - 1; ++j)
    //    for (size_t k = order - j - 1; k <= order - 1; ++k)
    //        shiftedCoeffs[k] += shiftedCoeffs[k+1] * (z0 - center);

    return shiftedCoeffs;
}

const vec3dVec Node::getDirectPhis() {
    vec3dVec phis;

    //for (const auto& obs : particles) {
    //    vec3d phi;
    //    for (const auto& src : particles)
    //        if (src != obs) 
    //            phi -= src->getCharge() * std::log(obs->getPos() - src->getPos());
    //    phis.push_back(phi);
    //}
    return phis;
}

const vec3dVec Node::getDirectFlds() {
    vec3dVec flds;

    //for (const auto& obs : particles) {
    //    vec3d fld;
    //    for (const auto& src : particles) {
    //        if (src != obs) {
    //            auto dz = obs->getPos() - src->getPos();
    //            fld += src->getCharge() * dz / std::norm(dz);
    //        }
    //    }
    //    flds.push_back(fld);
    //}
    return flds;
}