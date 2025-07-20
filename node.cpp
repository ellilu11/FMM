#include "math.h"
#include "node.h"
#include "leaf.h"
#include <cassert>
#include <iostream>
#include <random>

using enum Dir;

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

		default :
			std::cout << "Invalid dir" << std::endl;
			break;
	}
}

void Node::buildNearNeighbors() {
	for (int i = 0; i < 8; ++i) {
		Dir dir = static_cast<Dir>(i);
		auto nbor = getNeighborGeqSize(dir);
		if (nbor != nullptr)
			nbors.push_back(nbor);
	}
	assert(nbors.size() <= 8);
}

void Node::buildInteractionList() {
    assert( !isRoot() );
	//for (const auto& nbor : nbors) // flag near neighbors of this node
	//	nbor->setNodeStat(1);

	// base->buildNearNeighbors(); // uncomment for testing only
	auto baseNbors = base->getNearNeighbors();

	for (const auto& baseNbor : baseNbors)
		if (baseNbor->isNodeType<Leaf>() && !vecContains<std::shared_ptr<Node>>(nbors,baseNbor))
			iList.push_back(baseNbor);
		else 
			for (const auto& branch : baseNbor->branches) 
				if (!vecContains<std::shared_ptr<Node>>(nbors, branch))
					iList.push_back(branch);
}

const cmplxVec Node::shiftBaseLocalCoeffs(const int P) {
    assert( !isRoot() && !base->isRoot() );
    auto shftCoeffs( base->getLocalCoeffs() );

    for (size_t j = 0; j < P - 1; ++j)
        for (size_t k = P - j - 1; k < P - 1; ++k)
            shftCoeffs[k] -= (zk + base->getCenter()) * shftCoeffs[k + 1];

    return shftCoeffs;
}

const cmplx Node::evaluateFfield(const cmplx z, const int P) {
	auto phi( -coeffs[0] * std::log(z) );

	for (size_t k = 1; k < P; ++k)
		phi -= coeffs[k] / std::pow(z, k);

	return phi;
}

const cmplx Node::evaluateFfieldAnl(const cmplx z) {
    cmplx phi;
    for (size_t n = 0; n < psn.size(); ++n)
        phi -= qs[n] * std::log(z - psn[n]);
    return phi;
}

void Node::ffieldTest(const int P, const int Nobs){
    const double R (10.0 * L);

    std::ofstream obsFile, ffFile, ffAnlFile;
    // obsFile.open("obss.txt");
    ffFile.open("out/ff.txt"); //, std::ios::app);
    ffAnlFile.open("out/ffAnl.txt"); // , std::ios::app);

    for (int n = 0; n < Nobs; ++n) {
        const double th = 2 * M_PI * static_cast<double>(n) / static_cast<double>(Nobs);
        const cmplx z(R * cos(th), R * sin(th));

        const auto phi = evaluateFfield(z, P);

        const auto phiAnl = evaluateFfieldAnl(z);

        // obsFile << z << std::endl;
        ffFile << phi.real() << " ";
        // if (P == 1) 
            ffAnlFile << phiAnl.real() << " ";
    }
    ffFile << std::endl;
    ffAnlFile << std::endl;

}

