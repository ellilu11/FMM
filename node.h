#pragma once

#include <complex>
#include <iostream>
#include <vector>
#include "math.h"

namespace Param {
    extern const int DIM;
    extern const double EPS;
}

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
        : psn(psn), qs(qs), zk(zk), L_(L), lvl(lvl), branchIdx(branchIdx), base(base), nodeStat(0)
    {
    };

    static void buildBinomTable();
    static const int getP() { return P_; }
    static void setP(const int p) { P_ = p; }

    const cmplxVec getPsn() const { return psn; }
    const std::vector<double> getQs() const { return qs; }
    const cmplx getCenter() const { return zk; }
    const int getLvl() const { return lvl; }
    const std::vector<std::shared_ptr<Node>> getBranches() const { return branches; }
    Node* getBase() const { return base; }
    const bool isRoot() const { return base == nullptr; }

    template <typename T>
    bool isNodeType() const { return typeid(*this) == typeid(T); }

    cmplxVec getMpoleCoeffs() const { return coeffs; }
    cmplxVec getLocalCoeffs() const { return localCoeffs; }

    void printPsn(std::ofstream& f) {
        for (const auto& z : psn)
            f << z << std::endl;
    }

    void setNodeStat(int flag) { nodeStat = flag; }

    std::shared_ptr<Node> const getNeighborGeqSize(const Dir);

    void buildNearNeighbors();
    std::vector<std::shared_ptr<Node>> const getNearNeighbors() { return nbors; }

    void buildInteractionList();
    std::vector<std::shared_ptr<Node>> const getInteractionList()  { return iList; }

    const cmplxVec shiftBaseLocalCoeffs();
    const cmplx evaluateFfield(const cmplx);
    virtual const cmplx evaluateFfieldFromLeaf(const cmplx) = 0;
    const cmplx evalAnalyticFfield(const cmplx);
    void evalAnalyticNfield(std::ofstream&);

    virtual void buildMpoleCoeffs() = 0;
    virtual void resetNode() = 0;

    virtual void buildLocalCoeffs() = 0;
    virtual void printPhi(std::ofstream&) = 0;
    virtual void printNode(std::ofstream&) = 0;
    virtual void printMpoleCoeffs(std::ofstream&) = 0;
    virtual void printLocalCoeffs(std::ofstream&) = 0;

    // tests (move later)
    void ffieldTest(const int);
    virtual void mpoleToLocalTest() = 0;
    void nfieldTest();
    // virtual void iListTest() = 0;

protected:
    static int P_;
    static std::vector<std::vector<uint64_t>> binomTable;

    const std::vector<double> qs;
    const cmplx zk;
    const double L_;
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
