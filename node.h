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

class Node;

using NodeVec = std::vector<std::shared_ptr<Node>>;

class Node {
public:
    Node(ParticleVec&, const int, Node* const);

    static const int getExpansionOrder() { return order; }
    static void setExpansionOrder(const int p) { order = p; }
    static void setMaxLvl(const int lvl) { maxLvl = lvl; }

    ParticleVec getParticles() const { return particles; }
    const int getLvl() const { return lvl; }
    const double getLeng() const { return nodeLeng; }
    const vec3d getCenter() const { return center; }
    Node* getBase() const { return base; }
    const NodeVec getBranches() const { return branches; }
    NodeVec const getNearNeighbors() { return nbors; }
    NodeVec const getInteractionList() { return iList; }
    vec3dVec getMpoleCoeffs() const { return coeffs; }
    vec3dVec getLocalCoeffs() const { return localCoeffs; }

    const bool isRoot() const { return base == nullptr; }

    template <typename T>
    bool isNodeType() const { return typeid(*this) == typeid(T); }

    static void buildBinomTable();
    std::shared_ptr<Node> const getNeighborGeqSize(const Dir);
    void buildNearNeighbors();
    void buildInteractionList();
    void buildMpoleToLocalCoeffs();
    const vec3dVec getShiftedLocalCoeffs(const vec3d);
    //const vec3d getDirectPhiFar(const vec3d);
    //const vec3d getDirectPhi(const vec3d);
    const vec3dVec getDirectPhis();
    const vec3dVec getDirectFlds();

    virtual void buildMpoleCoeffs() = 0;
    virtual void buildLocalCoeffs() = 0;
    virtual void resetNode() = 0;
    virtual void printPhis(std::ofstream&) = 0;
    virtual void printNode(std::ofstream&) = 0;

protected:
    static int order;
    static int maxLvl;
    static double rootLeng;
    static std::vector<std::vector<uint64_t>> binomTable;

    ParticleVec particles;
    const int branchIdx;
    Node* const base;
    const double nodeLeng;
    const vec3d center;
    const int lvl;

    NodeVec branches;
    NodeVec nbors;
    NodeVec iList;

    vec3dVec coeffs;
    vec3dVec localCoeffs;
};
