#ifndef STEM_H
#define STEM_H

#include <vector>
#include <complex>
#include "node.h"

class Stem final : public Node {
public:
	Stem(std::vector<cmplx>&,
		std::vector<double>&,
		const cmplx,
		const double,
		const int,
		const int,
		Node* const);

	void buildCoeffs(const int);
	void buildLocalCoeffs(const int);
	std::vector<cmplx> shiftLocalCoeffs(const std::vector<cmplx>&, const int);
	void printNode(std::ofstream&);
	void iListTest();

};

#endif
