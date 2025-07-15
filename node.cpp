#include "node.h"
#include "leaf.h"
#include <cassert>
#include <iostream>

using enum Dir;

std::shared_ptr<Node> const Node::getNeighborGeqSize(const Dir dir) {
	// if node is root node, it has no neighbors in any direction
	if (base == nullptr)
		return nullptr;

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

		case W:
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
					base->getNeighborGeqSize(N):
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
			}
			else {
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
			}
			else {
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

		default:
			std::cout << "Invalid dir" << std::endl;
			break;
	}
}

std::vector<std::shared_ptr<Node>> const Node::getNearNeighbors() {
	std::vector<std::shared_ptr<Node>> nbors;
	// for (int i = N; i != Last; ++i) {
	for (int i = 0; i < 8; ++i) {
		Dir dir = static_cast<Dir>(i);
		auto nbor = getNeighborGeqSize(dir);
		if (nbor != nullptr)
			nbors.push_back(getNeighborGeqSize(dir));
	}
	assert(nbors.size() <= 8);
	return nbors;
}

void Node::setInteractionList() {
	auto nbors = getNearNeighbors();
	for (const auto& nbor : nbors) // flag near neighbors of this node
		nbor->setNodeStat(1);

	auto baseNbors = base->getNearNeighbors();

	for (const auto& nbor : baseNbors)
		if (nbor->branches.empty() && !nbor->nodeStat)
			iList.push_back(nbor);
		else 
			for (const auto& branch : nbor->branches) 
				if (!branch->nodeStat)
					iList.push_back(branch);

}