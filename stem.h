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
		const int);

	void buildCoeffs();
	void printNode(std::ofstream&);

private:
	std::vector<std::shared_ptr<Node>> branches;

};

#endif
