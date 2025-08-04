#include "../stem.h"

using namespace std;

shared_ptr<Node> Stem::getRandNode(int minLvl) {
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution branchIdx(0, 7);

    shared_ptr<Node> node = make_shared<Stem>(*this);

    // while (minLvl ? node->getLvl() > 1 : node->isNodeType<Stem>() )
    while (node->isNodeType<Stem>())
        node = (node->getBranches())[branchIdx(gen)];

    return node;
}

const cmplx Stem::getPhiFromBranchMpole(const vec3d& R, int maxLvl) {
    cmplx phi;
    std::cout << "Branch" << getLvl() << ' ';
    if (getLvl() >= maxLvl) return getPhiFromMpole(R);

    for (const auto& branch : branches)
        phi += branch->getPhiFromBranchMpole(R, maxLvl);
    return phi;
};

void Stem::printMpoleCoeffs(std::ofstream& f) {
    f << "Branch " << branchIdx << '\n';
    for (int l = 0; l < coeffs.size(); ++l) {
        for (int m = -l; m <= l; ++m)
            f << l << ' ' << m << ' ' << coeffs[l][m+l] << ' ';
        f << '\n';
    }
    f << '\n';
        
    for (const auto& branch : branches)
        branch->printMpoleCoeffs(f);
}

void Stem::resetNode() {
    coeffs.clear();
    localCoeffs.clear();
    nbors.clear();
    iList.clear();
    for (const auto& branch : branches)
        branch->resetNode();
}