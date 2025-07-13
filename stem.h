#ifndef STEM_H
#define STEM_H

#include <vector>
#include <complex>
#include "node.h"

class Stem final : public Node {
public:
	Stem(std::vector<cmplx>&,
		std::vector<double>&,
		cmplx,
		const double,
		const int,
		const int,
		Node* const);

	void buildCoeffs(const int);
	void buildLocalCoeffs(const int);
	void printNode(std::ofstream&);

};

#endif
