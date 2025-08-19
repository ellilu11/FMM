#pragma once

#include <iostream>
#include <type_traits>
#include "node.h"
#include "stem.h"

class Leaf final : public Node, public std::enable_shared_from_this<Leaf> {
public:
    Leaf(
        const ParticleVec&,
        const int,
        Stem* const);

    // void getNearNeighborsLessSize();

    realVec getFarPhis();

    std::vector<vec3d> getFarFlds();

    realVec getMidPhis();
    
    template <typename T, typename Func> std::vector<T> getNearSols(Func);
    
    void evaluateSolAtParticles();

    void pushToFarNeighbors(const std::shared_ptr<Node>& node) {
        farNbors.push_back(node);
    }

    const NodeVec getFarNbors() const { return farNbors; }

    const NodeVec getNearNbors(const int dirIdx) const { return nearNbors[dirIdx]; }

    void printNode(std::ofstream& f) override {
        f << center << " " << nodeLeng << " " << label << '\n';
    }

    std::shared_ptr<Node> getSelf() override {
        return shared_from_this();
    }

    void buildNeighbors() override;

    void buildLists() override;

    void buildMpoleCoeffs() override;

    void propagateExpCoeffs() override;

    void buildLocalCoeffs() override;

    // test methods
    std::shared_ptr<Node> getRandNode(int);

    const cmplx getPhiFromBranchMpole(const vec3d&, const int);

    void printMpoleCoeffs(std::ofstream&);

    void resetNode();

private:
    NodeVec farNbors; // list 3
    std::array<NodeVec,numDir> nearNbors; // list 1, indexed by direction
};
