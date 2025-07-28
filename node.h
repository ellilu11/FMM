#pragma once

#include <cassert>
#include <chrono>
#include <iostream>
#include "config.h"
#include "math.h"
#include "particle.h"

extern const int DIM;

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

class Node;

using NodeVec = std::vector<std::shared_ptr<Node>>;

class Node {
public:
    Node(const ParticleVec&, const int, Node* const);

    static const int getExpansionOrder() { return order; }
    static void setExpansionOrder(const int p) { order = p; }

    ParticleVec getParticles() const { return particles; }
    const double getLeng() const { return nodeLeng; }
    const cmplx getCenter() const { return center; }
    Node* getBase() const { return base; }
    const NodeVec getBranches() const { return branches; }
    NodeVec getNearNeighbors() const { return nbors; }
    NodeVec getInteractionList() const { return iList; }
    cmplxVec getMpoleCoeffs() const { return coeffs; }
    cmplxVec getLocalCoeffs() const { return localCoeffs; }

    const bool isRoot() const { return base == nullptr; }

    template <typename T>
    const bool isNodeType() const { return typeid(*this) == typeid(T); }

    static void setNodeParams(const Config&);
    static void buildBinomTable();
    std::shared_ptr<Node> const getNeighborGeqSize(const Dir);
    void buildNearNeighbors();
    void buildInteractionList();
    void buildMpoleToLocalCoeffs();
    const cmplxVec getShiftedLocalCoeffs(const cmplx);
    //const cmplx getDirectPhiFar(const cmplx);
    //const cmplx getDirectPhi(const cmplx);
    const cmplxVec getDirectPhis();
    const cmplxVec getDirectFlds();

    virtual void buildMpoleCoeffs() = 0;
    virtual void buildLocalCoeffs() = 0;
    virtual void resetNode() = 0;
    virtual void printPhis(std::ofstream&) = 0;
    virtual void printNode(std::ofstream&) = 0;

protected:
    static int order;
    static int maxNodeParts;
    static double rootLeng;
    static std::vector<std::vector<uint64_t>> binomTable;

    ParticleVec particles;
    const int branchIdx;
    Node* const base;
    const double nodeLeng;
    const cmplx center;

    NodeVec branches;
    NodeVec nbors;
    NodeVec iList;

    cmplxVec coeffs;
    cmplxVec localCoeffs;
};
