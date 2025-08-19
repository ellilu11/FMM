#include "node.h"

int Node::order; // # terms in multipole expansion
int Node::maxNodeParts;
double Node::rootLeng;
std::vector<std::vector<uint64_t>> Node::binomTable;
int Node::numNodes = 0;

void Node::setNodeParams(const Config& config) {
    order = ceil(-std::log(config.EPS) / std::log(2));
    maxNodeParts = config.maxNodeParts;
    rootLeng = config.L;
}

void Node::buildBinomTable() {
    for (int n = 0; n <= 2 * order - 1; ++n) {
        std::vector<uint64_t> binomN;
        for (int k = 0; k <= std::min(n, order-1); ++k)
            binomN.push_back(binom(n, k));
        binomTable.push_back(binomN);
    }
}

Node::Node(
    const ParticleVec& particles,
    const int branchIdx,
    Node* const base)
    : particles(particles), branchIdx(branchIdx), base(base),
    nodeLeng(base == nullptr ? rootLeng : base->getLeng() / 2.0),
    center(base == nullptr ? 0.0 :
        base->getCenter() +
        cmplx(pow(-1, branchIdx%2+1), pow(-1, branchIdx/2+1)) * nodeLeng / 2.0),
    label(0)
{
    numNodes++;
};

void Node::buildInteractionList() {
    assert( !isRoot() );

    auto notContains = [](NodeVec& vec, std::shared_ptr<Node> val) {
        return std::find(vec.begin(), vec.end(), val) == vec.end();
        };

    for (const auto& baseNbor : base->nbors) {
        if (baseNbor->isNodeType<Leaf>() && notContains(nbors, baseNbor)) {
            // iList.push_back(baseNbor);
            leafIlist.push_back(baseNbor);
            continue;
        }
        for (const auto& branch : baseNbor->branches)
            if (notContains(nbors, branch))
                iList.push_back(branch);
    }

    assert(iList.size() <= pow(6,DIM) - pow(3,DIM));
}

// If leaf is in list 4 of self, self is in list 3 of leaf
void Node::pushSelfToNearNonNbors() {
    if (leafIlist.empty()) return;

    auto self = getSelf(); // call shared_from_this()
    for (const auto& node : leafIlist) {
        auto leaf = dynamic_pointer_cast<Leaf>(node);
        leaf->pushToNearNonNbors(self);
    }
}

void Node::buildMpoleToLocalCoeffs() {
    if ( isRoot() ) return;

    localCoeffs.resize(order+1);

    for (const auto& iNode : iList) {
        auto mpoleCoeffs = iNode->getMpoleCoeffs();
        auto dz = iNode->getCenter() - center;
        auto mdz2k = cmplx(1,0);

        cmplxVec innerCoeffs; // innerCoeffs[k] = mpoleCoeffs[k] / (-dz)^k
        for (size_t k = 0; k <= order; ++k) {
            innerCoeffs.push_back(mpoleCoeffs[k] / mdz2k);
            mdz2k *= -dz;
        }

        localCoeffs[0] += innerCoeffs[0] * std::log(-dz);
        for (size_t k = 1; k <= order; ++k)
            localCoeffs[0] += innerCoeffs[k];

        auto dz2l = dz;
        for (size_t l = 1; l <= order; ++l) {
            localCoeffs[l] -= innerCoeffs[0] / (static_cast<double>(l) * dz2l);
            for (size_t k = 1; k <= order; ++k)
                localCoeffs[l] += innerCoeffs[k] / dz2l
                                    * static_cast<double>(binomTable[k+l-1][k-1]);
            dz2l *= dz;
        }
    }

    if (!base->isRoot()) localCoeffs += base->getShiftedLocalCoeffs(center);
    iList.clear();
}

const cmplxVec Node::getShiftedLocalCoeffs(const cmplx z0) {
    auto shiftedCoeffs( localCoeffs );
    for (size_t j = 0; j <= order - 1; ++j)
        for (size_t k = order - j - 1; k <= order - 1; ++k)
            shiftedCoeffs[k] += shiftedCoeffs[k+1] * (z0 - center);

    return shiftedCoeffs;
}

//const cmplx Node::getDirectPhiFar(const cmplx z) {
//    cmplx phi = -coeffs[0] * std::log(z-center);
//
//    for (size_t k = 1; k < order; ++k)
//        phi -= coeffs[k] / std::pow(z-center, k);
//
//    return phi;
//}
//
//const cmplx Node::getDirectPhi(const cmplx z) {
//    cmplx phi;
//    for (const auto& p : particles)
//        phi -= p->getCharge() * std::log(z - p->getPos());
//    return phi;
//}

const cmplxVec Node::getDirectPhis() {
    cmplxVec phis;

    for (const auto& obs : particles) {
        cmplx phi;
        for (const auto& src : particles)
            if (src != obs) 
                phi -= src->getCharge() * std::log(obs->getPos() - src->getPos());
        phis.push_back(phi);
    }
    return phis;
}

const cmplxVec Node::getDirectFlds() {
    cmplxVec flds;

    for (const auto& obs : particles) {
        cmplx fld;
        for (const auto& src : particles) {
            if (src != obs) {
                auto dz = obs->getPos() - src->getPos();
                fld += src->getCharge() * dz / std::norm(dz);
            }
        }
        flds.push_back(fld);
    }
    return flds;
}