#pragma once

#include <cassert>
#include <chrono>
#include <iostream>
#include "config.h"
#include "vec3d.h"
#include "particle.h"

extern const int DIM;

enum class Dir {
    DSW, DSE, DNW, DNE, DS, DW, DE, DN, D,
    SW, SE, NW, NE, S, W, E, N,
    USW, USE, UNW, UNE, US, UW, UE, UN, U
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
    const vec3d getCenter() const { return center; }
    Node* getBase() const { return base; }
    const NodeVec getBranches() const { return branches; }
    NodeVec const getNearNeighbors() { return nbors; }
    NodeVec const getInteractionList() { return iList; }
    cmplxVec getMpoleCoeffs() const { return coeffs; }
    vec3dVec getLocalCoeffs() const { return localCoeffs; }

    const bool isRoot() const { return base == nullptr; }

    template <typename T>
    const bool isNodeType() const { return typeid(*this) == typeid(T); }

    static void setNodeParams(const Config&);
    static void buildTables();
    static const double legendreLM(const double, int, int);
    static const cmplx sphHarmonic(const double, const double, int, int);

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
    static int maxNodeParts;
    static double rootLeng;
    static std::vector<realVec> sphHarmonicTable;
    static std::vector<realVec> fallingFactTable;
    static std::vector<realVec> legendreSumTable;

    ParticleVec particles;
    const int branchIdx;
    Node* const base;
    const double nodeLeng;
    const vec3d center;

    NodeVec branches;
    NodeVec nbors;
    NodeVec iList;

    cmplxVec coeffs;
    vec3dVec localCoeffs;
};
