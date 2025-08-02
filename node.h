#pragma once

#include <cassert>
#include <chrono>
#include <iostream>
#include "config.h"
#include "vec3d.h"
#include "particle.h"

extern const int DIM;

enum class Dir {
    W, E, S, N, D, U,
    SW, SE, NW, NE, DW, DE, UW, UE, DS, DN, US, UN,
    DSW, DSE, DNW, DNE, USW, USE, UNW, UNE
};

//enum class Dir {
//    D, U, S, N, W, E,
//    DS, US, DN, UN, DW, DE, UW, UE, SW, NW, SE, NE,
//    DSW, USW, DNW, UNW, DSE, USE, DNE, UNE
//};

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
    static matXcdVec rotationMatrixAlongDir(const int);

    static const double legendreLM(const double, const int, const int);

    std::shared_ptr<Node> const getNeighborGeqSize(const Dir);
    void buildNearNeighbors();
    void buildInteractionList();
    void buildMpoleToLocalCoeffs();
    const vecXcdVec getShiftedLocalCoeffs(const vec3d);
    //const vec3d getDirectPhiFar(const vec3d);
    //const vec3d getDirectPhi(const vec3d);
    const vec3dVec getDirectPhis();
    const vec3dVec getDirectFlds();

    virtual void buildMpoleCoeffs() = 0;
    virtual void buildLocalCoeffs() = 0;
    virtual void resetNode() = 0;
    virtual void printPhis(std::ofstream&) = 0;
    virtual void printNode(std::ofstream&) = 0;

    // test methods
    int getLvl() { return std::round(std::log(rootLeng/nodeLeng)/std::log(2)); }
    void setNodeStat(int stat) { nodeStat = stat; }
    virtual std::shared_ptr<Node> getRandNode(int) = 0;
    void setRandNodeStats();

protected:
    static int order;
    static int maxNodeParts;
    static double rootLeng;
    static std::vector<realVec> sphHarmonicTable;
    static std::vector<realVec> fallingFactTable;
    static std::vector<realVec> legendreSumTable;
    static std::vector<realVec> A;
    static std::vector<matXcdVec> rotationMat;

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
