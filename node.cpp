#include "node.h"
#include "leaf.h"
#include <iostream>

std::vector<std::shared_ptr<Node>> Node::getNearNeighbors() {
	std::vector<std::shared_ptr<Node>> nbors;
	// for (int i = N; i != Last; ++i) {
	for (int i = 0; i < 4; ++i) {
		Dir dir = static_cast<Dir>(i);
		auto nbor = getNeighborGeqSize(dir);
		if ( nbor != nullptr )
			nbors.push_back( getNeighborGeqSize(dir) );
		// getNeighborLessSize(dir);
	}
	return nbors;
}

std::shared_ptr<Node> Node::getNeighborGeqSize(const Dir dir) {
	// if node is root node, it has no neighbors in any direction
	if (root == nullptr)
		return nullptr;

	std::shared_ptr<Node> nbor;

	switch (dir) {
		case Dir::N :
			if (branchIdx < 2)
				return root->branches[branchIdx + 2];
			break;

		case Dir::E :
			if (~branchIdx % 2)
				return root->branches[branchIdx + 1];
			break;

		case Dir::S :
			if (branchIdx >= 2 )
				return root->branches[branchIdx - 2];
			break;

		case Dir::W :
			if (branchIdx % 2)
				return root->branches[branchIdx - 1];
			break;

		//case Dir::NE : 
		//	if (branchIdx == 0)
		//		return root->branches[3];
		//	break;

		default:
			std::cout << "Invalid dir" << std::endl;
			break;
	}
	
	nbor = root->getNeighborGeqSize(dir);
	if (nbor == nullptr)
		return nbor;
	if (typeid(*nbor) == typeid(Leaf))
		return nbor;

	switch (dir) {
		case Dir::N :
			return nbor->branches[branchIdx - 2];
			break;

		case Dir::E :
			return nbor->branches[branchIdx - 1];
			break;

		case Dir::S :
			return nbor->branches[branchIdx + 2];
			break;

		case Dir::W :
			return nbor->branches[branchIdx + 1];
			break;

		default:
			std::cout << "Invalid dir" << std::endl;
			break;
	}

}