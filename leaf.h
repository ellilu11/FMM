#ifndef LEAF_H
#define LEAF_H

#include <vector>
#include <complex>
#include "node.h"

class Leaf final : public Node {
public:
	Leaf(std::vector<cmplx>&,
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
