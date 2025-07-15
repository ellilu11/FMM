#ifndef NODE_H
#define NODE_H

#include <complex>
#include <optional>
#include <vector>
#include "math.h"

enum class Dir {
	SW,
	SE,
	NW,
	NE,
	S,
	W,
	E,
	N
};

class Node {
public:
	Node(std::vector<cmplx>& pos,
		 std::vector<double>& qs,
		 const cmplx z0,
		 const double L,
		 const int lvl,
		 const int branchIdx,
		 Node* const base)
		 : pos(pos), qs(qs), z0(z0), L(L), lvl(lvl), branchIdx(branchIdx), base(base), nodeStat(0)
	{
	};

	const cmplx getCenter() const { return z0; }
	const int getLvl() const { return lvl; }
	const std::shared_ptr<Node> getBranches(const size_t idx) const { return branches[idx]; }

	template <typename T>
	bool isNodeType() const { return typeid(*this) == typeid(T); }

	void setNodeStat(int flag) { nodeStat = flag; }

	std::vector<cmplx> getCoeffs() const { return coeffs; }
	std::shared_ptr<Node> const getNeighborGeqSize(const Dir);
	std::vector<std::shared_ptr<Node>> const getNearNeighbors();

	void setInteractionList();
	std::vector<std::shared_ptr<Node>> getInteractionList() const { return iList; }

	void printPos(std::ofstream& f) {
		for (const auto& z : pos)
			f << z.real() << " " << z.imag() << std::endl;
	}

	virtual void buildCoeffs(const int) = 0;
	virtual void buildLocalCoeffs(const int) = 0;
	virtual void printNode(std::ofstream&) = 0;
	// virtual void printCoeffs(std::ofstream&) = 0;

protected:
	const std::vector<double> qs;
	const cmplx z0;
	const double L;
	const int lvl;
	const int branchIdx;
	Node* const base;

	std::vector<std::shared_ptr<Node>> branches;
	std::vector<std::shared_ptr<Node>> iList;

	std::vector<cmplx> pos;
	std::vector<cmplx> coeffs;
	std::vector<cmplx> localCoeffs;

	int nodeStat;
};

#endif