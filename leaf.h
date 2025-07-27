#pragma once

#include <iostream>
#include "node.h"
#include "stem.h"

class Leaf final : public Node {
public:
    Leaf(
        const ParticleVec&,
        const int,
        Stem* const);

    void buildMpoleCoeffs();
    void buildLocalCoeffs();

    vec3dVec getFarPhis();
    vec3dVec getFarFlds();
    template <typename Func> vec3dVec getNearSols(Func);
    void evaluateSolAtParticles();

    void printPhis(std::ofstream& f) {
        for (const auto& p : particles)
            f << p->getPhi().real() << ' ';
    }

    void printNode(std::ofstream& f) {
        f << center << " " << nodeLeng << " " << '\n';
    }

    void resetNode() {
        coeffs.clear();
        localCoeffs.clear();
        nbors.clear();
        iList.clear();
    }
};
