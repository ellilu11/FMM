#ifndef NODE_H
#define NODE_H

#include <vector>
#include <complex>
#include "math.h"

class Node {
public:
	Node(std::vector<cmplx>& pos,
		std::vector<double>& qs,
		cmplx z0,
		const double L,
		const int lvl,
		const int P)
		: pos(pos), qs(qs), z0(z0), L(L), lvl(lvl), P(P)
	{};

	const std::vector<cmplx> getCoeffs() const { return coeffs; }

	const cmplx getCenter() const { return z0; }

	void printPos(std::ofstream& f) {
		for (const auto& z : pos)
			f << z.real() << " " << z.imag() << std::endl;
	};

	virtual void buildCoeffs() = 0;
	virtual void printNode(std::ofstream&) = 0;

protected:
	const std::vector<double> qs;
	const cmplx z0;
	const double L;
	const int lvl;
	const int P;

	std::vector<cmplx> pos;
	std::vector<cmplx> coeffs;
};

#endif