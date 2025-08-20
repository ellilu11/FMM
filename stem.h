#pragma once

#include "node.h"

class Stem final : public Node, public std::enable_shared_from_this<Stem> {
public:
    Stem(
        const ParticleVec&,
        const int,
        Stem* const);

    void printNode(std::ofstream& f) {
        f << center << " " << nodeLeng << " " << label << '\n';
        for (const auto& branch : branches)
            branch->printNode(f);
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


//private:
//    NodeVec branches;

};
