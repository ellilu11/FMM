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

	void buildMpoleCoeffs(const int);
	void buildLocalCoeffs(const int);

    void printPhi(std::ofstream& f) {
        for (const auto& branch : branches) 
            branch->printPhi(f);
    }

    void printNode(std::ofstream& f) {
        f << zk << " " << L << " " << nodeStat << std::endl;
        for (const auto& branch : branches)
            branch->printNode(f);
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
