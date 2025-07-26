#pragma once

#include "node.h"

class Stem final : public Node {
public:
    Stem(ParticleVec&,
        const cmplx,
        const int,
        Stem* const);

    void buildMpoleCoeffs();

    void buildLocalCoeffs();

    void printPhis(std::ofstream& f) {
        for (const auto& branch : branches) 
            branch->printPhis(f);
    }

    void printNode(std::ofstream& f) {
        f << center << " " << nodeLeng << " " << '\n';
        for (const auto& branch : branches)
            branch->printNode(f);
    }

    void resetNode() {
        coeffs.clear();
        localCoeffs.clear();
        nbors.clear();
        iList.clear();
        for (const auto& branch : branches)
            branch->resetNode();
    }

//private:
//    NodeVec branches;

};
