#include "../leaf.h"
#include "../node.h"

using namespace std;

shared_ptr<Node> Leaf::getRandNode(int minLvl) {
    return make_shared<Leaf>(*this);
}

const cmplx Leaf::getFfieldFromLeaf(const vec3d& R) {
    return getDirectPhiFromMpole(R);
}

void Leaf::resetNode() {
    coeffs.clear();
    localCoeffs.clear();
    nbors.clear();
    iList.clear();
}