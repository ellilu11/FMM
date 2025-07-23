#pragma once

#include <vector>
#include <complex>
#include "node.h"
#include "stem.h"

class Leaf final : public Node {
public:
    Leaf(cmplxVec&,
        realVec&,
        const cmplx,
        const double,
        const int,
        const int,
        Stem* const);

    void buildMpoleCoeffs();

    void resetNode() { 
        coeffs.clear();
        localCoeffs.clear();
        nbors.clear();
        iList.clear();
        phis.clear();
    }

    const cmplx getFfieldFromLeaf(const cmplx z) {
        return getFfield(z);
    }

    void buildLocalCoeffs();
    cmplxVec getPhiFarSrc();
    cmplxVec getPhiNearSrc();
    void evaluatePhi();

    cmplxVec getPhi() { return phis; }

    void printPhi(std::ofstream& f) {
        for (const auto& phi : phis)
            f << phi.real() << ' ';
    }

    void printNode(std::ofstream& f) {
        f << zk << " " << L_ << " " << nodeStat << '\n';
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

private:
    cmplxVec phis;
};
