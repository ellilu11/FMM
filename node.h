#pragma once

#include <cassert>
#include <chrono>
#include <iostream>
#include <numeric>
#include "config.h"
#include "enum.h"
#include "particle.h"
#include "tables.h"
#include "vec3d.h"

extern const int DIM;
extern std::chrono::duration<double, std::milli> t_M2X;
extern std::chrono::duration<double, std::milli> t_X2X;
extern std::chrono::duration<double, std::milli> t_X2L;
extern std::chrono::duration<double, std::milli> t_X2L_l4;
extern std::chrono::duration<double, std::milli> t_L2L;
extern std::chrono::duration<double, std::milli> t_L2P;
extern std::chrono::duration<double, std::milli> t_dir;

constexpr int numDir = 26; // std::pow(3, DIM) - 1;

class Node;

using NodeVec = std::vector<std::shared_ptr<Node>>;

class Node {
public:
    static const int getExpansionOrder() { return order; }
    static void setExpansionOrder(const int p) { order = p; }
    static const int getExponentialOrder() { return orderExp; }
    static const int getNumNodes() { return numNodes; }

    ParticleVec getParticles() const { return particles; }
    const int getBranchIdx() const { return branchIdx; }
    const double getLeng() const { return nodeLeng; }
    const vec3d getCenter() const { return center; }
    Node* getBase() const { return base; }
    const NodeVec getBranches() const { return branches; }
    NodeVec const getNearNeighbors() { return nbors; }
    std::array<NodeVec, 6> const getDirList() { return dirList; }
    std::array<NodeVec, 6> const getOuterDirList() { return outerDirList; }
    NodeVec const getLeafIlist() { return leafIlist; }
    std::vector<vecXcd> getMpoleCoeffs() const { return coeffs; }
    std::vector<vecXcd> getExpCoeffs(const int dirIdx) const { 
        return expCoeffsOut[dirIdx]; }
    std::vector<vecXcd> getLocalCoeffs() const { return localCoeffs; }
    const bool isRoot() const { return base == nullptr; }
    template <typename T>
    const bool isNodeType() const { return typeid(*this) == typeid(T); }

    static void setNodeParams(const Config&);
    static void buildTables(const Config&);
    static void buildRotationMats();
    // static std::vector<matXcd> wignerDAlongDir(const pair2d, const bool);
    static const double legendreCos(const double, const int, const int);

    Node(const ParticleVec&, const int, Node* const);
    std::shared_ptr<Node> const getNeighborGeqSize(const Dir);
    void buildNearNeighbors();
    void buildInteractionList();
    void buildOuterInteractionList();

    void buildLocalCoeffsFromDirList();
    void buildLocalCoeffsFromLeafIlist();
    const std::vector<vecXcd> getShiftedLocalCoeffs(const int) const;

    std::vector<vecXcd> getMpoleToExpCoeffs(const int);
    const std::vector<vecXcd> getMergedExpCoeffs(const int) const;
    void addShiftedExpCoeffs(const std::vector<vecXcd>&, const vec3d&, const int);

    const double getDirectPhi(const vec3d&);
    const realVec getDirectPhis();
    const vec3d getDirectFld(const vec3d&);
    const std::vector<vec3d> getDirectFlds();

    virtual void buildLists() = 0;
    virtual void buildMpoleCoeffs() = 0;
    virtual void propagateExpCoeffs() = 0;
    virtual void buildLocalCoeffs() = 0;
    virtual void resetNode() = 0;
    virtual void printPhis(std::ofstream&) = 0;
    virtual void printNode(std::ofstream&) = 0;

    // ========== Test methods ==========
    int getLvl() { return std::round(std::log(rootLeng/nodeLeng)/std::log(2)); }
    void setNodeStat(int stat) { nodeStat = stat; }
    static void setExponentialOrder(const Precision prec) {
        orderExp = [&]() -> std::size_t {
            switch (prec) {
                case Precision::LOW:    return 8;
                case Precision::MEDIUM: return 17;
                case Precision::HIGH:   return 26;
            }
            }();
    }

    // definition under test/nodetest.cpp
    void setRandNodeStats();
    const cmplx getPhiFromMpole(const vec3d&);
    void ffieldTest(const int, const int, const int);
    void nfieldTest();

    // definition under test/[stemtest.cpp, leaftest.cpp]
    virtual std::shared_ptr<Node> getRandNode(int) = 0;
    virtual const cmplx getPhiFromBranchMpole(const vec3d&, const int) = 0;
    virtual void printMpoleCoeffs(std::ofstream&) = 0;

    // definition under test/exptest.cpp
    const cmplx getPhiFromExp(const vec3d&, const std::vector<vecXcd>&, const int);
    const cmplx getPhiFromLocal(const vec3d&);
    void mpoleToExpToLocalTest();

protected:
    static int order;
    static int orderExp;
    static int maxNodeParts;
    static double rootLeng;
    static Tables tables;
    static std::array<std::vector<matXcd>,14> wignerD;
    static std::array<std::vector<matXcd>,14> wignerDInv;
    static std::array<mat3d, 6> rotMatR;
    static int numNodes;

    ParticleVec particles;
    const int branchIdx;
    Node* const base;
    const double nodeLeng;
    const vec3d center;

    NodeVec branches;
    NodeVec nbors; // list 1
    std::array<NodeVec,6> dirList; // list 2, indexed by direction
    std::array<NodeVec, 6> outerDirList; // intersection of list 2 of all branches
    NodeVec leafIlist; // list 4

    std::vector<vecXcd> coeffs;
    std::array<std::vector<vecXcd>,6> expCoeffsOut;
    std::array<std::vector<vecXcd>,6> expCoeffs;
    std::vector<vecXcd> localCoeffs;

    // === Test members ===
    int nodeStat;
    bool useRot;

};
