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

	void buildMpoleCoeffs(const int);
	void buildLocalCoeffs(const int);
    void evaluatePhiLocalExp(const int);
    void evaluatePhiDirect(const int);
    void evaluatePhi(const int);

    void printPhi(std::ofstream& f) {
        for (size_t n = 0; n < psn.size(); ++n)
            f << psn[n] << " " << phis[n] << std::endl;
    }

    void printNode(std::ofstream& f) {
        f << zk << " " << L << " " << nodeStat << std::endl;
    }

    void printLocalCoeffs(std::ofstream& f) {
        for (const auto& coeff : localCoeffs)
            f << coeff << " ";
        f << "\n";
    }

	void iListTest();

private:
    cmplxVec phis;
    cmplxVec flds;
};
