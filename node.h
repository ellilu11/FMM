#ifndef NODE_H
#define NODE_H

#include <vector>
#include <complex>
#include "math.h"

enum class Dir {
	N,
	E,
	S,
	W,
	NE,
	SE,
	SW,
	NW,
	Last
};

class Node {
public:
	Node(std::vector<cmplx>& pos,
		std::vector<double>& qs,
		cmplx z0,
		const double L,
		const int lvl,
		const int branchIdx,
		Node* const root)
		: pos(pos), qs(qs), z0(z0), L(L), lvl(lvl), branchIdx(branchIdx), root(root), nborFlag(0)
	{
	};

	const std::vector<cmplx> getCoeffs() const { return coeffs; }

	const cmplx getCenter() const { return z0; }

	void printPos(std::ofstream& f) {
		for (const auto& z : pos)
			f << z.real() << " " << z.imag() << std::endl;
	};

	std::vector<std::shared_ptr<Node>> getNearNeighbors();
	std::shared_ptr<Node> getNeighborGeqSize(const Dir);

	const int getLvl() const { return lvl; }
	const std::shared_ptr<Node> getBranches(const size_t idx) const { return branches[idx]; }

	void setNborFlag(int flag) { nborFlag = flag; }

	virtual void buildCoeffs(const int) = 0;
	virtual void buildLocalCoeffs(const int) = 0;
	virtual void printNode(std::ofstream&) = 0;
	// virtual void printCoeffs(std::ofstream&) = 0;

protected:
	const std::vector<double> qs;
	const cmplx z0;
	const double L;
	const int lvl;

	std::vector<std::shared_ptr<Node>> branches;
	const int branchIdx;

	std::vector<cmplx> pos;
	std::vector<cmplx> coeffs;
	std::vector<cmplx> localCoeffs;

	Node* const root;

	int nborFlag;
};

#endif