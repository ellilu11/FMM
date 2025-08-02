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

const cmplx Stem::getFfieldFromLeaf(const vec3d& R) {
    cmplx phi;
    for (const auto& branch : branches)
        phi += branch->getFfieldFromLeaf(R);
    return phi;
};

void Stem::resetNode() {
    coeffs.clear();
    localCoeffs.clear();
    nbors.clear();
    iList.clear();
    for (const auto& branch : branches)
        branch->resetNode();
}