#pragma once

#include <iostream>
#include "node.h"
#include "stem.h"

class Leaf final : public Node, public std::enable_shared_from_this<Leaf> {
public:
    Leaf(
        const ParticleVec&,
        const int,
        Stem* const);

    cmplxVec getFarPhis();

    cmplxVec getFarFlds();

    cmplxVec getNearNonNborPhis();

    cmplxVec getNearNonNborFlds();

    template <typename Func> cmplxVec getNearNborSols(Func);

    void evaluateSolAtParticles();

    void pushToNearNonNbors(const std::shared_ptr<Node>& node) {
        nearNonNbors.push_back(node);
    }

    const NodeVec getNearNonNbors() const { return nearNonNbors; }

    const NodeVec getNearNbors() const { return nearNbors; }

    void printNode(std::ofstream& f) {
        f << center << " " << nodeLeng << " " << label << '\n';
    }

    std::shared_ptr<Node> getSelf() override {
        return shared_from_this();
    }

    void buildNbors() override;

    void buildLists() override;

    void buildMpoleCoeffs() override;

    void buildLocalCoeffs() override;

    /* Test methods */
    std::shared_ptr<Node> getRandNode(int);
    void resetNode();

private:
    NodeVec nearNbors; // list 1
    NodeVec nearNonNbors; // list 3

};
