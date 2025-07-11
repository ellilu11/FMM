#ifndef NODE_H
#define NODE_H

#include <vector>
#include <complex>

class Node {
public:
	Node(std::vector<std::complex<double>>&,
		std::vector<double>&,
		std::complex<double>,
		const double,
		const int,
		const int);

	void buildCoeffs();
	const std::vector<std::complex<double>> getCoeffs() const { return coeffs; }

	void printPos(std::ofstream&);
	void printNode(std::ofstream&);

private:
	const std::vector<double> qs;
	const std::complex<double> z0;
	const double L;
	const int lvl;
	const int P;
	const bool isLeaf; // remove this in favor of derived classes for inner and leaf nodes

	std::vector<std::complex<double>> pos;
	std::vector<std::shared_ptr<Node>> branches;
	std::vector<std::complex<double>> coeffs;

};

#endif