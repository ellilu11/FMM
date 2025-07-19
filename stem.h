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

	void buildCoeffs(const int);
	void buildLocalCoeffs(const int);
	void printNode(std::ofstream&);
	void iListTest();

//private:
//    std::vector<std::shared_ptr<Node>> branches;

};
