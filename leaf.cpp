#include "leaf.h"
#include <iostream>

Leaf::Leaf(std::vector<cmplx>& pos,
	std::vector<double>& qs,
	const cmplx zk,
	const double L,
	const int lvl,
	const int branchIdx,
	Node* const base)
	: Node(pos, qs, zk, L, lvl, branchIdx, base)
{
}

void Leaf::buildCoeffs(const int P) {
	//auto sumOverPos =
	//	[this](const double q, const cmplx pos, const int idx) {
	//		return q *
	//			idx == 0 ? 1 : -pow(pos-z0, idx) / static_cast<double>(idx); 
	//	};
	cmplx a_k;
	for (int k = 0; k < P; ++k) {
		// auto a_k = std::accumulate()
		for (size_t i = 0; i < pos.size(); ++i)
			a_k += qs[i] * 
				k == 0 ? 1 : -pow(pos[i]-zk, k) / static_cast<double>(k);
		coeffs.push_back(a_k);
	}
}

std::vector<cmplx> Leaf::shiftLocalCoeffs(const std::vector<cmplx>& coeffs, const int P) {
	return localCoeffs;
}

void Leaf::buildLocalCoeffs(const int P) {

}

void Leaf::printNode(std::ofstream& f) {
	f << zk.real() << " " << zk.imag() << " " << L << " " << nodeStat << std::endl;
}

void Leaf::iListTest() {
	setNodeStat(3);
	setInteractionList();
	auto nbors = getInteractionList();

	for (const auto& nbor : nbors)
		nbor->setNodeStat(2);

	std::ofstream posFile, nodeFile;
	posFile.open("positions.txt");
	nodeFile.open("nodes.txt");

	printPos(posFile);
	printNode(nodeFile);
}