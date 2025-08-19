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
    static const int getExpansionOrder() { return order; }

    static void setExpansionOrder(const int p) { order = p; }

    static const int getNumNodes() { return numNodes; }

    static void setNodeParams(const Config&);

    static void buildBinomTable();

    /* member access and utility functions */
    ParticleVec getParticles() const { return particles; }

    const double getLeng() const { return nodeLeng; }
    
    const cmplx getCenter() const { return center; }
    
    Node* getBase() const { return base; }
    
    const NodeVec getBranches() const { return branches; }
    
    NodeVec getNeighbors() const { return nbors; }
    
    NodeVec getInteractionList() const { return iList; }
    
    cmplxVec getMpoleCoeffs() const { return coeffs; }
    
    cmplxVec getLocalCoeffs() const { return localCoeffs; }

    const bool isRoot() const { return base == nullptr; }

    template <typename T>
    const bool isNodeType() const { return typeid(*this) == typeid(T); }

    /* defined in node.cpp */
    Node(const ParticleVec&, const int, Node* const);

    std::shared_ptr<Node> const getNeighborGeqSize(const Dir);

    void buildNeighbors();

    void buildInteractionList();

    void buildMpoleToLocalCoeffs();

    const cmplxVec getShiftedLocalCoeffs(const cmplx);

    const cmplxVec getDirectPhis();

    const cmplxVec getDirectFlds();

    /* pure abstract */
    virtual void buildLists() = 0;

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
    static int numNodes;

    ParticleVec particles;
    const int branchIdx;
    Node* const base;
    const double nodeLeng;
    const cmplx center;

    NodeVec branches;
    NodeVec nbors;
    NodeVec iList;
    NodeVec leafIlist;

    cmplxVec coeffs;
    cmplxVec localCoeffs;
};
