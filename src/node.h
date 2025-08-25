#pragma once

#include <cassert>
#include <iostream>
#include <numeric>
#include <queue>
#include "config.h"
#include "clock.h"
#include "particle.h"
#include "tables.h"

extern ClockTimes t;

constexpr int DIM = 3;
constexpr int numDir = 26;

enum class Dir {
    W, E, S, N, D, U,
    SW, SE, NW, NE, DW, DE, UW, UE, DS, DN, US, UN,
    DSW, DSE, DNW, DNE, USW, USE, UNW, UNE
};

class Node;

using NodeVec = std::vector<std::shared_ptr<Node>>;

class Node {

public:
    static const int getExpansionOrder() { return order; }

    static void setExpansionOrder(const int p) { order = p; }
    
    static const int getExponentialOrder() { return orderExp; }
    
    static const int getNumNodes() { return numNodes; }
    
    static void setNodeParams(const Config&);
    
    static void buildTables(const Config&);
    
    static void buildRotationMats();
    
    static double legendreCos(const double, const int, const int);

    static double dLegendreCos(const double, const int, const int);

public:
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
    
    std::vector<vecXcd> getLocalCoeffs() const { return localCoeffs; }

    const bool isRoot() const { return base == nullptr; }
    
    template <typename T>
    const bool isNodeType() const { return typeid(*this) == typeid(T); }

    void resetSols() { for (const auto& p : particles) p->resetSol(); }

public:
    Node(const ParticleVec&, const int, Node* const);
    
    std::shared_ptr<Node> getNeighborGeqSize(const Dir) const;

    NodeVec getNeighborsLeqSize(const std::shared_ptr<Node>&, const Dir) const;
    
    void assignToDirList(std::array<NodeVec, 6>&, const std::shared_ptr<Node>&, const double);
    
    void buildInteractionList();
    
    void pushSelfToNearNonNbors();
   
    void evalExpToLocalCoeffs();
    
    std::vector<vecXcd> getShiftedLocalCoeffs(const int) const;

    void evalLeafIlistSols();

    const std::vector<vecXcd> getMpoleToExpCoeffs(const int) const;
    
    void addShiftedExpCoeffs(const std::vector<vecXcd>&, const vec3d&, const int);

    void evalPairSols(const std::shared_ptr<Node>&);

    void evalSelfSols();
   
    virtual std::shared_ptr<Node> getSelf() = 0;
    
    virtual void buildNeighbors() = 0;

    virtual void buildLists() = 0;
    
    virtual void buildMpoleCoeffs() = 0;
    
    virtual void propagateExpCoeffs() = 0;
    
    virtual void buildLocalCoeffs() = 0;
    
    virtual void printNode(std::ofstream&) = 0;

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
    std::array<NodeVec,6> dirList; // list 2, indexed by direction
    NodeVec leafIlist; // list 4

    std::vector<vecXcd> coeffs;
    std::array<std::vector<vecXcd>,6> expCoeffs;
    std::vector<vecXcd> localCoeffs;

    int label;
};