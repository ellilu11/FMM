#include "node.h"

using enum Dir;

std::shared_ptr<Node> const Node::getNeighborGeqSize(const Dir dir) {
    if (isRoot()) return nullptr;

    std::shared_ptr<Node> nbor;

    switch (dir) {
        case N:
            if (branchIdx == 0 || branchIdx == 1)
                return base->branches[branchIdx+2];
            nbor = base->getNeighborGeqSize(dir);
            if (nbor == nullptr)
                return nbor;
            if (nbor->isNodeType<Leaf>())
                return nbor;
            return nbor->branches[branchIdx-2];
            break;

        case E:
            if (branchIdx == 0 || branchIdx == 2)
                return base->branches[branchIdx+1];
            nbor = base->getNeighborGeqSize(dir);
            if (nbor == nullptr)
                return nbor;
            if (nbor->isNodeType<Leaf>())
                return nbor;
            return nbor->branches[branchIdx-1];
            break;

        case W:
            if (branchIdx == 1 || branchIdx == 3)
                return base->branches[branchIdx-1];
            nbor = base->getNeighborGeqSize(dir);
            if (nbor == nullptr)
                return nbor;
            if (nbor->isNodeType<Leaf>())
                return nbor;
            return nbor->branches[branchIdx+1];
            break;

        case S:
            if (branchIdx == 2 || branchIdx == 3)
                return base->branches[branchIdx-2];
            nbor = base->getNeighborGeqSize(dir);
            if (nbor == nullptr)
                return nbor;
            if (nbor->isNodeType<Leaf>())
                return nbor;
            return nbor->branches[branchIdx+2];
            break;

        case NE:
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
                    base->getNeighborGeqSize(N) :
                    base->getNeighborGeqSize(E);
                if (nbor == nullptr)
                    return nbor;
                if (nbor->isNodeType<Leaf>())
                    // do not double count neighbor that will be found along a cardinal direction
                    return nullptr;
                return nbor->branches[3-branchIdx];
            }
            break;

        case NW:
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
                    base->getNeighborGeqSize(N) :
                    base->getNeighborGeqSize(W);
                if (nbor == nullptr)
                    return nbor;
                if (nbor->isNodeType<Leaf>())
                    return nullptr;
                return nbor->branches[3-branchIdx];
            }
            break;

        case SE:
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
                    base->getNeighborGeqSize(S) :
                    base->getNeighborGeqSize(E);
                if (nbor == nullptr)
                    return nbor;
                if (nbor->isNodeType<Leaf>())
                    return nullptr;
                return nbor->branches[3-branchIdx];
            }
            break;

        case SW:
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
                    base->getNeighborGeqSize(S) :
                    base->getNeighborGeqSize(W);
                if (nbor == nullptr)
                    return nbor;
                if (nbor->isNodeType<Leaf>())
                    return nullptr;
                return nbor->branches[3-branchIdx];
            }
            break;
    }
}

NodeVec Node::getNeighborsLeqSize(
    const std::shared_ptr<Node>& nborGeqSize, const Dir dir) const
{
    NodeVec nbors{};
    NodeVec queue{ nborGeqSize };
    if (isRoot()) return nbors;

    while (!queue.empty()) {
        auto nbor = queue[0];

        if (nbor->isNodeType<Leaf>())
            nbors.push_back(nbor);
        else {
            switch (dir) {
                case N:
                    queue.push_back(nbor->branches[0]);
                    queue.push_back(nbor->branches[1]);
                    break;

                case E:
                    queue.push_back(nbor->branches[0]);
                    queue.push_back(nbor->branches[2]);
                    break;

                case W:
                    queue.push_back(nbor->branches[1]);
                    queue.push_back(nbor->branches[3]);
                    break;

                case S:
                    queue.push_back(nbor->branches[2]);
                    queue.push_back(nbor->branches[3]);
                    break;

                case NE:
                    queue.push_back(nbor->branches[0]);
                    break;

                case NW:
                    queue.push_back(nbor->branches[1]);
                    break;

                case SE:
                    queue.push_back(nbor->branches[2]);
                    break;

                case SW:
                    queue.push_back(nbor->branches[3]);
                    break;
            }
        }
        queue.erase(queue.begin());
    }
    return nbors;
}