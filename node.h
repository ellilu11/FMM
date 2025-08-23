#pragma once

#include <cassert>
#include <iostream>
#include <numeric>
#include "config.h"
#include "clock.h"
#include "enum.h"
#include "particle.h"
#include "tables.h"
#include "vec3d.h"

extern ClockTimes t;

constexpr int DIM = 3;
constexpr int numDir = 26;

class Node;

using NodeVec = std::vector<std::shared_ptr<Node>>;

class Node {
public:
    /* static member functions */
    static const int getExpansionOrder() { return order; }

    static void setExpansionOrder(const int p) { order = p; }
    
    static const int getExponentialOrder() { return orderExp; }
    
    static const int getNumNodes() { return numNodes; }

    static void findNearNborPairs();
    
    static void setNodeParams(const Config&);
    
    static void buildTables(const Config&);
    
    static void buildRotationMats();
    
    static const double legendreCos(const double, const int, const int);

    static const double dLegendreCos(const double, const int, const int);

    /* member access and utility functions */
    ParticleVec getParticles() const { return particles; }
    
    const int getBranchIdx() const { return branchIdx; }
    
    const double getLeng() const { return nodeLeng; }
    
    const vec3d getCenter() const { return center; }
    
    Node* getBase() const { return base; }
    
    const NodeVec getBranches() const { return branches; }
    
    NodeVec getNbors() const { return nbors; }

    NodeVec getDirList(const int dir) const { return dirList[dir]; }

    NodeVec getLeafIlist() const { return leafIlist; }
    
    std::vector<vecXcd> getMpoleCoeffs() const { return coeffs; }

    std::vector<vecXcd> getExpCoeffs(const int dirIdx) const { return expCoeffsOut[dirIdx]; }

    std::vector<vecXcd> getLocalCoeffs() const { return localCoeffs; }

    const bool isRoot() const { return base == nullptr; }
    
    template <typename T>
    const bool isNodeType() const { return typeid(*this) == typeid(T); }

    /* definition in node.cpp */
    Node(const ParticleVec&, const int, Node* const);
    
    std::shared_ptr<Node> getNeighborGeqSize(const Dir) const;

    NodeVec getNeighborsLeqSize(const std::shared_ptr<Node>&, const Dir) const;
    
    void assignToDirList(std::array<NodeVec, 6>&, const std::shared_ptr<Node>&, const double);
    
    void buildInteractionList();

    void buildOuterInteractionList();
    
    void pushSelfToNearNonNbors();
   
    void evalLocalCoeffsFromDirList();
    
    const std::vector<vecXcd> getShiftedLocalCoeffs(const int) const;

    void evalLeafIlistSols();

    std::vector<vecXcd> getMpoleToExpCoeffs(const int);

    const std::vector<vecXcd> getMergedExpCoeffs(const int) const;

    void addShiftedExpCoeffsFromBranch(const std::vector<vecXcd>&, const vec3d&, const int);

    void addShiftedExpCoeffs(const std::vector<vecXcd>&, const vec3d&, const int);

    void evalPairSols(const std::shared_ptr<Node>&);

    void evalSelfSols();

    void resetSols();
   
    /* pure virtual */
    virtual std::shared_ptr<Node> getSelf() = 0;
    
    virtual void buildNeighbors() = 0;

    virtual void buildLists() = 0;
    
    virtual void buildMpoleCoeffs() = 0;
    
    virtual void propagateExpCoeffs() = 0;
    
    virtual void buildLocalCoeffs() = 0;
    
    virtual void printNode(std::ofstream&) = 0;

    // ========== Test methods ==========
    int getLvl() { return std::round(std::log(rootLeng/nodeLeng)/std::log(2)); }

    void labelNode(int label_) { label += label_; }

    // definition under test/nodetest.cpp
    void labelNodes();
    const pairSol getDirectSol(const vec3d&, const double = 1.0E-12);
    const cmplx getPhiFromMpole(const vec3d&);
    void ffieldTest(const int, const int, const int);
    void nfieldTest();

    // definition under test/[stemtest.cpp, leaftest.cpp]
    virtual std::shared_ptr<Node> getRandNode(int) = 0;
    virtual const cmplx getPhiFromBranchMpole(const vec3d&, const int) = 0;
    virtual void printMpoleCoeffs(std::ofstream&) = 0;
    virtual void resetNode() = 0;

    // definition under test/exptest.cpp
    const cmplx getPhiFromExp(const vec3d&, const std::vector<vecXcd>&, const int);
    const cmplx getPhiFromLocal(const vec3d&);
    void mpoleToExpToLocalTest();

protected:
    static int order;
    static int orderExp;
    static int maxNodeParts;
    static double rootLeng;
    static int numNodes;
    static Tables tables;
    static std::array<std::vector<matXcd>,14> wignerD;
    static std::array<std::vector<matXcd>,14> wignerDInv;
    static std::array<mat3d,6> rotMatR;

    ParticleVec particles;
    const int branchIdx;
    Node* const base;
    const double nodeLeng;
    const vec3d center;

    NodeVec branches;
    NodeVec nbors; // list 1
    std::array<NodeVec, 6> dirList; // list 2, indexed by direction
    std::array<NodeVec, 6> outerDirList; // intersection of list 2 of all branches
    NodeVec leafIlist; // list 4

    std::vector<vecXcd> coeffs;
    std::array<std::vector<vecXcd>,6> expCoeffsOut;
    std::array<std::vector<vecXcd>,6> expCoeffs;
    std::vector<vecXcd> localCoeffs;

    int label;
};
