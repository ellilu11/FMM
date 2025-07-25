#pragma once

#include <iostream>
#include "node.h"
#include "stem.h"

class Leaf final : public Node {
public:
    Leaf(
        ParticleVec&,
        const cmplx,
        const int,
        const int,
        Stem* const);

    void buildMpoleCoeffs();
    void buildLocalCoeffs();
    cmplxVec getPhisFar();
    cmplxVec getPhisNear();
    cmplxVec getFldsFar();
    cmplxVec getFldsNear();
    void evaluateSolAtParticles();

    const cmplx getFfieldFromLeaf(const cmplx z) {
        return getDirectPhiFar(z);
    }

    void resetNode() {
        coeffs.clear();
        localCoeffs.clear();
        nbors.clear();
        iList.clear();
    }

    //void printPhi(std::ofstream& f) {
    //    for (const auto& phi : phis)
    //        f << phi.real() << ' ';
    //}

    void printNode(std::ofstream& f) {
        f << center << " " << nodeLeng << " " << nodeStat << '\n';
    }

    void printMpoleCoeffs(std::ofstream& f) {
        for (const auto& coeff : coeffs)
            f << coeff << " ";
        f << "\n";
    }

    void printLocalCoeffs(std::ofstream& f) {
        for (const auto& coeff : localCoeffs)
            f << coeff << " ";
        f << "\n";
    }

    void mpoleToLocalTest();
    //void iListTest();

//private:
//    cmplxVec phis;
};
