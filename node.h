#pragma once

#include <cassert>
#include <chrono>
#include <iostream>
#include "math.h"
#include "particle.h"

namespace Param {
    extern const int DIM;
    extern const double L;
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
    Node(
        ParticleVec& particles,
        const cmplx center,
        const int lvl,
        const int branchIdx,
        Node* const base)
        : particles(particles), center(center), lvl(lvl), branchIdx(branchIdx), base(base),
        nodeLeng( rootLeng / static_cast<double>(std::pow(2,lvl)) ),
        nodeStat(0)
    {
    };

    static const int getExpansionOrder() { return order; }
    static void setExpansionOrder(const int p) { order = p; }
    static void setMaxLvl(const int lvl) { maxLvl = lvl; }

    ParticleVec getParticles() const { return particles; }
    const cmplx getCenter() const { return center; }
    const int getLvl() const { return lvl; }
    const std::vector<std::shared_ptr<Node>> getBranches() const { return branches; }
    Node* getBase() const { return base; }
    const bool isRoot() const { return base == nullptr; }

    template <typename T>
    bool isNodeType() const { return typeid(*this) == typeid(T); }

    cmplxVec getMpoleCoeffs() const { return coeffs; }
    cmplxVec getLocalCoeffs() const { return localCoeffs; }

    void setNodeStat(int flag) { nodeStat = flag; }

    static void buildBinomTable();
    std::shared_ptr<Node> const getNeighborGeqSize(const Dir);

    void buildNearNeighbors();
    std::vector<std::shared_ptr<Node>> const getNearNeighbors() { return nbors; }

    void buildInteractionList();
    std::vector<std::shared_ptr<Node>> const getInteractionList()  { return iList; }

    void buildMpoleToLocalCoeffs();

    const cmplxVec getShiftedLocalCoeffs(const cmplx);
    const cmplx getDirectPhiFar(const cmplx);
    const cmplx getDirectPhi(const cmplx);
    const cmplxVec getDirectPhis();

    virtual void buildMpoleCoeffs() = 0;
    virtual void resetNode() = 0;

    virtual void buildLocalCoeffs() = 0;
    virtual void printPhis(std::ofstream&) = 0;
    virtual void printNode(std::ofstream&) = 0;

protected:
    static int order;
    static int maxLvl;
    static double rootLeng;
    static std::vector<std::vector<uint64_t>> binomTable;

    ParticleVec particles;
    const cmplx center;
    const int lvl;
    const double nodeLeng;
    const int branchIdx;
    Node* const base;

    std::vector<std::shared_ptr<Node>> branches;
    std::vector<std::shared_ptr<Node>> nbors;
    std::vector<std::shared_ptr<Node>> iList;

    cmplxVec coeffs;
    cmplxVec localCoeffs;

    int nodeStat;
};
