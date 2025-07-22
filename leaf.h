#pragma once

#include <vector>
#include <complex>
#include "node.h"
#include "stem.h"

class Leaf final : public Node {
public:
    Leaf(cmplxVec&,
        std::vector<double>&,
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
    }

    const cmplx evaluateFfieldFromLeaf(const cmplx z) {
        return evaluateFfield(z);
    }

    void buildLocalCoeffs();
    void evaluatePhiLocalExp();
    void evaluatePhiDirect();
    void evaluatePhi();

    void printPhi(std::ofstream& f) {
        for (const auto& phi : phis)
            f << phi << std::endl;
    }

    void printNode(std::ofstream& f) {
        f << zk << " " << L_ << " " << nodeStat << std::endl;
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
    cmplxVec flds;
};
