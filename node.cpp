#include "node.h"

using enum Dir;

int Node::order = ceil(-log(Param::EPS) / log(2)); // # terms in multipole expansion
std::vector<std::vector<uint64_t>> Node::binomTable;

void Node::buildBinomTable() {
    for (int n = 0; n <= 2 * order - 1; ++n) {
        std::vector<uint64_t> binomN;
        for (int k = 0; k <= std::min(n,order-1); ++k)
            binomN.emplace_back(binom(n, k));
        binomTable.push_back( binomN );
    }
}

std::shared_ptr<Node> const Node::getNeighborGeqSize(const Dir dir) {
    if ( isRoot() ) return nullptr;

    std::shared_ptr<Node> nbor;

    switch (dir) {
        case N :
            if (branchIdx == 0 || branchIdx == 1)
                return base->branches[branchIdx+2];
            nbor = base->getNeighborGeqSize(dir);
            if (nbor == nullptr)
                return nbor;
            if (nbor->isNodeType<Leaf>())
                return nbor;
            return nbor->branches[branchIdx-2];
            break;

        case E :
            if (branchIdx == 0 || branchIdx == 2)
                return base->branches[branchIdx+1];
            nbor = base->getNeighborGeqSize(dir);
            if (nbor == nullptr)
                return nbor;
            if (nbor->isNodeType<Leaf>())
                return nbor;
            return nbor->branches[branchIdx-1];
            break;

        case W :
            if (branchIdx == 1 || branchIdx == 3 )
                return base->branches[branchIdx-1];
            nbor = base->getNeighborGeqSize(dir);
            if (nbor == nullptr)
                return nbor;
            if (nbor->isNodeType<Leaf>())
                return nbor;
            return nbor->branches[branchIdx+1];
            break;

        case S :
            if (branchIdx == 2 || branchIdx == 3 )
                return base->branches[branchIdx-2];
            nbor = base->getNeighborGeqSize(dir);
            if (nbor == nullptr)
                return nbor;
            if (nbor->isNodeType<Leaf>())
                return nbor;
            return nbor->branches[branchIdx+2];
            break;

        case NE :
            if (branchIdx == 0)
                return base->branches[3];
            if (branchIdx == 3) {
                nbor = base->getNeighborGeqSize(dir);
                if (nbor == nullptr)
                    return nbor;
                if (nbor->isNodeType<Leaf>())
                    return nbor;
                return nbor->branches[0];
            } else {
                nbor = branchIdx == 2 ?
                base->getNeighborGeqSize(N):
                base->getNeighborGeqSize(E);
                if (nbor == nullptr)
                    return nbor;
                if (nbor->isNodeType<Leaf>())
                    // do not double count neighbor that will be found along a cardinal direction
                    return nullptr; 
                return nbor->branches[3-branchIdx];
            }
            break;

        case NW :
            if (branchIdx == 1)
                return base->branches[2];
            if (branchIdx == 2) {
                nbor = base->getNeighborGeqSize(dir);
                if (nbor == nullptr)
                    return nbor;
                if (nbor->isNodeType<Leaf>())
                    return nbor;
                return nbor->branches[1];
            } else {
                nbor = branchIdx == 3 ?
                base->getNeighborGeqSize(N):
                base->getNeighborGeqSize(W);
                if (nbor == nullptr)
                    return nbor;
                if (nbor->isNodeType<Leaf>())
                    return nullptr;
                return nbor->branches[3-branchIdx];
            }
            break;

        case SE :
            if (branchIdx == 2)
                return base->branches[1];
            if (branchIdx == 1) {
                nbor = base->getNeighborGeqSize(dir);
                if (nbor == nullptr)
                    return nbor;
                if (nbor->isNodeType<Leaf>())
                    return nbor;
                return nbor->branches[2];
            } else {
                nbor = branchIdx == 0 ?
                base->getNeighborGeqSize(S):
                base->getNeighborGeqSize(E);
                if (nbor == nullptr)
                    return nbor;
                if (nbor->isNodeType<Leaf>())
                    return nullptr;
                return nbor->branches[3-branchIdx];
            }
            break;

        case SW :
            if (branchIdx == 3)
                return base->branches[0];
            if (branchIdx == 0) {
                nbor = base->getNeighborGeqSize(dir);
                if (nbor == nullptr)
                    return nbor;
                if (nbor->isNodeType<Leaf>())
                    return nbor;
                return nbor->branches[3];
            } else {
                nbor = branchIdx == 1 ?
                base->getNeighborGeqSize(S):
                base->getNeighborGeqSize(W);
                if (nbor == nullptr)
                    return nbor;
                if (nbor->isNodeType<Leaf>())
                    return nullptr;
                return nbor->branches[3-branchIdx];
            }
            break;

        default :
            throw std::runtime_error("Invalid direction");
        break;
    }
}

void Node::buildNearNeighbors() {
    const int nDir = std::pow(3, Param::DIM) - 1;
    for (int i = 0; i < nDir; ++i) {
        Dir dir = static_cast<Dir>(i);
        auto nbor = getNeighborGeqSize(dir);
        if (nbor != nullptr)
            nbors.push_back(nbor);
    }
    assert(nbors.size() <= nDir);
}

void Node::buildInteractionList() {
    assert( !isRoot() );

    auto baseNbors = base->getNearNeighbors();

    for (const auto& baseNbor : baseNbors)
        if (baseNbor->isNodeType<Leaf>() && !contains<std::shared_ptr<Node>>(nbors,baseNbor))
            iList.push_back(baseNbor);
        else 
            for (const auto& branch : baseNbor->branches) 
                if (!contains<std::shared_ptr<Node>>(nbors, branch))
                    iList.push_back(branch);

    assert(iList.size() <= pow(6,Param::DIM) - pow(3,Param::DIM));
}

void Node::buildMpoleToLocalCoeffs() {
    buildNearNeighbors();

    if (!isRoot()) {
        buildInteractionList();
        localCoeffs.resize(order+1);

        for (const auto& iNode : iList) {
            auto mpoleCoeffs( iNode->getMpoleCoeffs() );
            auto dz( iNode->getCenter() - center );
            auto mdz2k( cmplx(1,0) );

            cmplxVec innerCoeffs; // innerCoeffs[k] = mpoleCoeffs[k] / (-dz)^k
            for (size_t k = 0; k <= order; ++k) {
                innerCoeffs.push_back(mpoleCoeffs[k] / mdz2k);
                mdz2k *= -dz;
            }

            localCoeffs[0] += innerCoeffs[0] * std::log(-dz);
            for (size_t k = 1; k <= order; ++k)
                localCoeffs[0] += innerCoeffs[k];

            for (size_t l = 1; l <= order; ++l) {
                cmplx dz2l = std::pow(dz, l);
                localCoeffs[l] -= innerCoeffs[0] / (static_cast<double>(l) * dz2l);
                for (size_t k = 1; k <= order; ++k)
                    localCoeffs[l] += innerCoeffs[k] / dz2l
                                        * static_cast<double>(binomTable[k+l-1][k-1]);
            }
        }

        if (!base->isRoot()) localCoeffs += base->getShiftedLocalCoeffs(center);
        iList.clear();
    }
}

const cmplxVec Node::getShiftedLocalCoeffs(const cmplx z0) {
    assert( !isRoot() );
    
    auto shiftedCoeffs( localCoeffs );
    for (size_t j = 0; j <= order - 1; ++j)
        for (size_t k = order - j - 1; k <= order - 1; ++k)
            shiftedCoeffs[k] += shiftedCoeffs[k+1] * (z0 - center);

    return shiftedCoeffs;
}

