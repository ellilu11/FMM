#ifndef LEAF_H
#define LEAF_H

#include <vector>
#include <complex>
#include "node.h"

class Leaf final : public Node {
public:
	Leaf(std::vector<cmplx>&,
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
