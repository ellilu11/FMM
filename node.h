#pragma once

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
	Node(cmplxVec& psn,
		 std::vector<double>& qs,
		 const cmplx zk,
		 const double L,
		 const int lvl,
		 const int branchIdx,
		 Node* const base)
		 : psn(psn), qs(qs), zk(zk), L(L), lvl(lvl), branchIdx(branchIdx), base(base), nodeStat(0)
	{
	};

    const cmplxVec getPsn() const { return psn;  }
    const std::vector<double> getQs() const { return qs; }
	const cmplx getCenter() const { return zk; }
	const int getLvl() const { return lvl; }
	const std::vector<std::shared_ptr<Node>> getBranches() const { return branches; }
    const bool isRoot() const { return base == nullptr; }

	template <typename T>
	bool isNodeType() const { return typeid(*this) == typeid(T); }

	void setNodeStat(int flag) { nodeStat = flag; }

	cmplxVec getMpoleCoeffs() const { return coeffs; }
	cmplxVec getLocalCoeffs() const { return localCoeffs; }

	void printpsn(std::ofstream& f) {
		for (const auto& z : psn)
			f << z << std::endl;
	}

	std::shared_ptr<Node> const getNeighborGeqSize(const Dir);

	void buildNearNeighbors();
	std::vector<std::shared_ptr<Node>> const getNearNeighbors() { return nbors; }

	void buildInteractionList();
	std::vector<std::shared_ptr<Node>> const getInteractionList()  { return iList; }

    const cmplxVec shiftBaseLocalCoeffs(const int);

	const cmplx evaluateFfield(const cmplx, const int);

    const cmplx evaluateFfieldAnl(const cmplx);

	virtual void buildMpoleCoeffs(const int) = 0;
	virtual void buildLocalCoeffs(const int) = 0;

    virtual void printPhi(std::ofstream&) = 0;
    virtual void printNode(std::ofstream&) = 0;
    virtual void printLocalCoeffs(std::ofstream&) = 0;

	virtual void iListTest() = 0;

    void ffieldTest(const int, const int);



protected:
	const std::vector<double> qs;
	const cmplx zk;
	const double L;
	const int lvl;
	const int branchIdx;
	Node* const base;

    std::vector<std::shared_ptr<Node>> branches;
	std::vector<std::shared_ptr<Node>> nbors;
	std::vector<std::shared_ptr<Node>> iList;

	cmplxVec psn;
	cmplxVec coeffs;
	cmplxVec localCoeffs;

	int nodeStat;
};
