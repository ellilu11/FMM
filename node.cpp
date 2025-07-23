#include "math.h"
#include "node.h"
#include "leaf.h"
#include <cassert>
#include <chrono>
#include <iostream>
#include <random>

using enum Dir;

int Node::P_ = ceil(-log(Param::EPS) / log(2)); // # terms in multipole expansion
std::vector<std::vector<uint64_t>> Node::binomTable;

void Node::buildBinomTable() {
    for (int n = 0; n <= 2 * P_ - 1; ++n) {
        std::vector<uint64_t> binomN;
        for (int k = 0; k <= std::min(n,P_-1); ++k)
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
            // throw exception
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

    // base->buildNearNeighbors(); // uncomment for iListTest() only
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

const cmplxVec Node::getShiftedLocalCoeffs(const cmplx z0) {
    assert( !isRoot() );
    
    auto shiftedCoeffs( localCoeffs );
    for (size_t j = 0; j <= P_ - 1; ++j)
        for (size_t k = P_ - j - 1; k <= P_ - 1; ++k)
            shiftedCoeffs[k] += (z0 - zk) * shiftedCoeffs[k+1];

    return shiftedCoeffs;
}


