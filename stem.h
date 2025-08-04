#pragma once

#include "node.h"

class Stem final : public Node {
public:
    Stem(
        const ParticleVec&,
        const int,
        Stem* const);

    void buildMpoleCoeffs();
    void buildLocalCoeffs();

    void printPhis(std::ofstream& f) {
        for (const auto& branch : branches) 
            branch->printPhis(f);
    }

    void printNode(std::ofstream& f) {
        f << center << " " << nodeLeng << " " << nodeStat << '\n';
        for (const auto& branch : branches)
            branch->printNode(f);
    }

    // test methods
    std::shared_ptr<Node> getRandNode(int);
    const cmplx getPhiFromBranchMpole(const vec3d&, const int);
    void printMpoleCoeffs(std::ofstream&);
    void resetNode();

};
