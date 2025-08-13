#include "../leaf.h"
#include "../node.h"

using namespace std;

shared_ptr<Node> Leaf::getRandNode(int minLvl) {
    return make_shared<Leaf>(*this);
}

const cmplx Leaf::getPhiFromBranchMpole(const vec3d& R, const int Nlvl) {
    // std::cout << "Leaf" << getLvl() << ' ';
    return getPhiFromMpole(R);
}

void Leaf::printMpoleCoeffs(std::ofstream& f) {
    f << "Leaf " << branchIdx << '\n';
    for (int l = 0; l < coeffs.size(); ++l) {
        for (int m = -l; m <= l; ++m)
        // for (int m = l; m >= -l; --m)
            f << l << ' ' << m << ' ' << coeffs[l][m+l] << ' ';
        f << '\n';
    }
    f << '\n';
}

void Leaf::resetNode() {
    coeffs.clear();
    expCoeffs = {};
    localCoeffs.clear();
    /*nbors.clear();
    iList.clear();
    dirList = {};*/
}