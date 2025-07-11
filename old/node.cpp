#include "node.h"
#include <cassert>
#include <iostream>

Node::Node(std::vector<cmplx> &pos,
		   std::vector<double> &qs,
	       cmplx z0,	
		   const double L,
		   const int lvl,
		   const int P)
	: pos(pos), qs(qs), z0(z0), L(L), lvl(lvl), P(P),
      isLeaf(pos.size()<2 || !lvl) // if node has at least two particles and is not lvl 0, then subdivide node
{
	if (!isLeaf) {
		std::vector<std::vector<cmplx>> branchPos(4);
		std::vector<std::vector<double>> branchQs(4);
		// for (const auto& ele : pos) {
		for (int n = 0; n < pos.size(); ++n){
			int k = (pos[n].real() > z0.real()) + 2 * (pos[n].imag() > z0.imag());
			assert(k < 4);
			branchPos[k].push_back(pos[n]);
			branchQs[k].push_back(qs[n]);
		}
		for (int k = 0; k < 4; ++k) {
			cmplx dz0 = (pow(-1,k%2+1) + iu*pow(-1,k/2+1)) * L/4.0;
			auto branch = std::make_shared<Node>(branchPos[k], branchQs[k], z0 + dz0, L / 2, lvl - 1, P); //
			branches.push_back(branch);
		}
	}
}

void Node::buildCoeffs() {
	if (isLeaf) {
		for (int i = 0; i < pos.size(); ++i) {
			coeffs[0] += qs[i];
			for (int k = 1; k < P; ++k)
				coeffs[k] -= qs[i] * pow(pos[i] - z0, k) / static_cast<double>(k); // pos[i]-z0 or pos[i] ?
		}
	} else {
		// recursively combine coeffs from child nodes into coeffs for parent node
		// only combine coeffs of nodes on the same lvl!
		for (const auto& branch : branches){
			branch->buildCoeffs(); 
			auto branchCoeffs = branch->getCoeffs();
			coeffs[0] += branchCoeffs[0];
			for (int k = 1; k < P; ++k)
				for (int l = 0; l < k; ++l)
					coeffs[k] += branchCoeffs[l] * pow(branch->z0, k - l) * binom(k-1, l-1)
						- branchCoeffs[0] * pow(branch->z0, k) / static_cast<double>(k);
		}
	}
}

void Node::printNode(std::ofstream& f) {
	f << z0.real() << " " << z0.imag() << " " << L << std::endl;
	for (const auto& branch : branches) 
		branch->printNode(f);
}
