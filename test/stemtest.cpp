#include "../stem.h"

using namespace std;

shared_ptr<Node> Stem::getRandNode(int maxLvl) {
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution branchIdx(0, 7);

    shared_ptr<Node> node = make_shared<Stem>(*this);

    // while (node->getLvl() <= maxLvl)
    while (node->isNodeType<Stem>())
        node = (node->getBranches())[branchIdx(gen)];

    return node;
}

const cmplx Stem::getPhiFromBranchMpole(const vec3d& R, int maxLvl) {
    cmplx phi;
    // std::cout << "Branch" << getLvl() << ' ';
    if (getLvl() >= maxLvl) return getPhiFromMpole(R);

    for (const auto& branch : branches)
        phi += branch->getPhiFromBranchMpole(R, maxLvl);
    return phi;
};

void Stem::printMpoleCoeffs(std::ofstream& f) {
    f << "Branch " << branchIdx << '\n';
    for (int l = 0; l < coeffs.size(); ++l) {
        for (int m = -l; m <= l; ++m)
        // for (int m = l; m >= -l; --m)
            f << coeffs[l][m+l] << ' ';
        f << '\n';
    }
    f << '\n';
        
    for (const auto& branch : branches)
        branch->printMpoleCoeffs(f);
}

void Stem::resetNode() {
    //coeffs.clear();
    //localCoeffs.clear();
    //nbors.clear();
    //iList.clear();
    //dirList = {};

    for (int dir = 0; dir < 6; ++dir) {
        std::vector<vecXcd> expCoeffs_dir;
        for (int k = 0; k < orderExp; ++k)
            expCoeffs_dir.emplace_back(vecXcd::Zero(tables.quadLengs_[k]));
        expCoeffs[dir] = expCoeffs_dir;
    }

    for (const auto& branch : branches)
        branch->resetNode();
}

