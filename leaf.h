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

	void buildCoeffs(const int);
	void buildLocalCoeffs(const int);
    void evaluatePhiLocalExp(const int);
    void evaluatePhiDirect(const int);
    void evaluatePhi(const int);

	void printNode(std::ofstream&);
	void iListTest();

private:
    cmplxVec phis;
};
