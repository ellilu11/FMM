#pragma once

#include <vector>
#include <complex>
#include "node.h"

class Stem final : public Node {
public:
    Stem(cmplxVec&,
        std::vector<double>&,
        const cmplx,
        const double,
        const int,
        const int,
        Stem* const);

    void buildMpoleCoeffs();
    const cmplx evaluateFfieldFromLeaf(const cmplx z) {
        cmplx phi;
        for (const auto& branch : branches)
            phi += branch->evaluateFfieldFromLeaf(z);
        return phi;
    };

    void clearMpoleCoeffs() { 
        coeffs.clear(); 
        for (const auto& branch : branches)
            branch->clearMpoleCoeffs();
    }

    void buildLocalCoeffs();

    void printPhi(std::ofstream& f) {
        for (const auto& branch : branches) 
            branch->printPhi(f);
    }

    void printNode(std::ofstream& f) {
        f << zk << " " << L_ << " " << nodeStat << std::endl;
        for (const auto& branch : branches)
            branch->printNode(f);
    }

    void printMpoleCoeffs(std::ofstream& f) {

        for (const auto& coeff : coeffs)
            f << coeff << " ";
        f << "\n";
        for (const auto& branch : branches)
            branch->printMpoleCoeffs(f);
    }

    void printLocalCoeffs(std::ofstream& f) {
        for (const auto& coeff : localCoeffs)
            f << coeff << " ";
        f << "\n";
        for (const auto& branch : branches)
            branch->printLocalCoeffs(f);
    }

    void iListTest();

//private:
//    std::vector<std::shared_ptr<Node>> branches;

};
