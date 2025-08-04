#pragma once

#include <cassert>
#include <chrono>
#include <iostream>
#include "config.h"
#include "math.h"
#include "particle.h"
#include "vec3d.h"

extern const int DIM;

enum class Dir {
    W, E, S, N, D, U,
    SW, SE, NW, NE, DW, DE, UW, UE, DS, DN, US, UN,
    DSW, DSE, DNW, DNE, USW, USE, UNW, UNE
};

struct Tables {
    std::vector<realVec> coeffYlm;
    std::vector<realVec> fallingFact;
    std::vector<realVec> legendreSum;
    std::vector<realVec> A;
    // std::vector<realVec> quadWeights;
};

class Node;

using NodeVec = std::vector<std::shared_ptr<Node>>;

class Node {
public:
    Node(const ParticleVec&, const int, Node* const);

    static const int getExpansionOrder() { return order; }
    static void setExpansionOrder(const int p) { order = p; }
    
    static matXcdVec getRotationMatrixAlongDir(int dir) { return rotationMat[dir]; }

    ParticleVec getParticles() const { return particles; }
    const int getBranchIdx() const { return branchIdx; }
    const double getLeng() const { return nodeLeng; }
    const vec3d getCenter() const { return center; }
    Node* getBase() const { return base; }
    const NodeVec getBranches() const { return branches; }
    NodeVec const getNearNeighbors() { return nbors; }
    NodeVec const getInteractionList() { return iList; }
    vecXcdVec getMpoleCoeffs() const { return coeffs; }
    vecXcdVec getLocalCoeffs() const { return localCoeffs; }

    const bool isRoot() const { return base == nullptr; }

    template <typename T>
    const bool isNodeType() const { return typeid(*this) == typeid(T); }

    static void setNodeParams(const Config&);
    static void buildTables();
    static void buildRotationMats();
    static matXcdVec rotationMatrixAlongDir(int, const bool);

    static const double legendreLM(const double, const int, const int);

    std::shared_ptr<Node> const getNeighborGeqSize(const Dir);
    void buildNearNeighbors();
    void buildInteractionList();
    void buildMpoleToLocalCoeffs();
    const vecXcdVec getShiftedLocalCoeffs(const vec3d);

    const vec3dVec getDirectPhis();
    const vec3dVec getDirectFlds();

    virtual void buildMpoleCoeffs() = 0;
    virtual void buildLocalCoeffs() = 0;
    virtual void resetNode() = 0;
    virtual void printPhis(std::ofstream&) = 0;
    virtual void printNode(std::ofstream&) = 0;

    // === Test methods ===
    int getLvl() { return std::round(std::log(rootLeng/nodeLeng)/std::log(2)); }
    void setNodeStat(int stat) { nodeStat = stat; }

    // definition under test/nodetest.cpp
    void setRandNodeStats();
    const double getDirectPhi(const vec3d&);
    const cmplx getPhiFromMpole(const vec3d&);
    void ffieldTest(const int, const int, const int);
    // void nfieldTest();

    // definition under test/[stemtest.cpp, leaftest.cpp]
    virtual std::shared_ptr<Node> getRandNode(int) = 0;
    virtual const cmplx getPhiFromBranchMpole(const vec3d&, const int) = 0;
    virtual void printMpoleCoeffs(std::ofstream&) = 0;

protected:
    static int order;
    static int maxNodeParts;
    static double rootLeng;
    static Tables tables;
    static std::vector<matXcdVec> rotationMat;
    static std::vector<matXcdVec> rotationInvMat;

    ParticleVec particles;
    const int branchIdx;
    Node* const base;
    const double nodeLeng;
    const vec3d center;

    NodeVec branches;
    NodeVec nbors;
    NodeVec iList;

    vecXcdVec coeffs;
    vecXcdVec localCoeffs;

    // test members
    int nodeStat;
};
