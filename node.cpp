#include "node.h"

int Node::order;
int Node::orderExp;
int Node::maxNodeParts;
double Node::rootLeng;
int Node::numNodes;
Tables Node::tables;
std::array<std::vector<matXcd>,14> Node::wignerD;
std::array<std::vector<matXcd>,14> Node::wignerDInv;
std::array<mat3d,6> Node::rotMatR;


void Node::setNodeParams(const Config& config) {
    order = ceil(-std::log(config.EPS) / std::log(2));
    orderExp = [&]() -> std::size_t {
        switch (config.prec) {
            case Precision::LOW:    return 8;
            case Precision::MEDIUM: return 17;
            case Precision::HIGH:   return 26;
        }
        }();
    maxNodeParts = config.maxNodeParts;
    rootLeng = config.L;
    numNodes = 0;
}

void Node::buildTables(const Config& config) {
    tables = Tables(order, config.prec);
    assert(orderExp == tables.quadCoeffs_.size());
}

void Node::buildRotationMats() {
    auto wignerDAlongDir = [](const pair2d angles, bool isInv) {
        std::vector<matXcd> mats;
        for (int l = 0; l <= order; ++l)
            mats.push_back(
                isInv ?
                wignerD_l(angles, l) :
                wignerD_l(angles, l).adjoint()
            );
        return mats;
    };

    for (int dir = 0; dir < 14; ++dir) {
        auto X = 
            dir < 8 ? 
            idx2pm(dir) : // vertex directions (for M2M and L2L)
            [dir] { 
            switch (dir-8) {
                case 0: return vec3d(0, 0, 1);
                case 1: return vec3d(0, 0, -1);
                case 2: return vec3d(0, 1, 0);
                case 3: return vec3d(0, -1, 0);
                case 4: return vec3d(1, 0, 0);
                case 5: return vec3d(-1, 0, 0);
            } 
            }(); // cardinal directions (for M2X and X2L)

        auto R = toSph(X);
        pair2d angles(R[1], R[2]);

        wignerD[dir] = wignerDAlongDir(angles,0);
        wignerDInv[dir] = wignerDAlongDir(angles,1);
        if (dir >=8)
            rotMatR[dir-8] = rotationR(angles);
    }
}

// theta part of Ylm
const double Node::legendreCos(const double th, const int l, const int abs_m) {
    assert(abs_m <= l);

    const auto cos_th = cos(th);
    const auto sin_th = sin(th);
    double legendreSum = 0.0;

    // term is zero for l-k odd
    for (int k = l; k >= abs_m; k -= 2)
        legendreSum += 
            tables.fallingFact_[k][abs_m] * tables.legendreSum_[l][k] 
            * pow(cos_th, k-abs_m);

    return tables.coeffYlm_[l][abs_m] * pow(sin_th, abs_m) * legendreSum;
}

Node::Node(
    const ParticleVec& particles,
    const int branchIdx,
    Node* const base)
    : particles(particles), branchIdx(branchIdx), base(base),
    nodeLeng(base == nullptr ? rootLeng : base->getLeng() / 2.0),
    center(base == nullptr ? zeroVec :
        base->getCenter() + nodeLeng / 2.0 * idx2pm(branchIdx)),
    label(0), useRot(0)
{
    for (int l = 0; l <= order; ++l)
        localCoeffs.emplace_back(vecXcd::Zero(2*l+1));

    for (int dir = 0; dir < 6; ++dir) {
        std::vector<vecXcd> expCoeffs_dir;
        for (int k = 0; k < orderExp; ++k)
            expCoeffs_dir.emplace_back( vecXcd::Zero(tables.quadLengs_[k]) );
        expCoeffs[dir] = expCoeffs_dir;
    }

    numNodes++;
}

void Node::assignToDirList(
    std::array<NodeVec,6>& list, const std::shared_ptr<Node>& node,
    const double minDist) 
{
    auto center0 = node->center;

    if (center0[2] - center[2] > minDist)      // uplist
        list[0].push_back(node);
    else if (center[2] - center0[2] > minDist) // downlist
        list[1].push_back(node);
    else if (center0[1] - center[1] > minDist) // northlist
        list[2].push_back(node);
    else if (center[1] - center0[1] > minDist) // southlist
        list[3].push_back(node);
    else if (center0[0] - center[0] > minDist) // eastlist
        list[4].push_back(node);
    else if (center[0] - center0[0] > minDist) // westlist
        list[5].push_back(node);
}

void Node::buildInteractionList() {
    assert(!isRoot());
    assert(!nbors.empty());

    auto notContains = [](NodeVec& vec, std::shared_ptr<Node> val) {
        return std::find(vec.begin(), vec.end(), val) == vec.end();
    };

    NodeVec iList;
    for (const auto& baseNbor : base->nbors) {
        if (baseNbor->isNodeType<Leaf>() && notContains(nbors, baseNbor)) {
            leafIlist.push_back(baseNbor); // list 4
            continue;
        }
        for (const auto& branch : baseNbor->branches)
            if (notContains(nbors, branch))
                iList.push_back(branch);
    }

    assert(iList.size() <= pow(6, DIM) - pow(3, DIM));

    // pick minDist \in (nodeLeng, 2.0*nodeLeng) to avoid rounding errors
    const double minDist = 1.5 * nodeLeng;

    for (const auto& node : iList)
        assignToDirList(dirList, node, minDist);

    //auto dirListSize =
    //    std::accumulate(dirList.begin(), dirList.end(), size_t{ 0 },
    //        [](size_t sum, const auto& vec) { return sum + vec.size(); });
    //assert(dirListSize == iList.size());
}

/* If leaf is in list 4 of self, self is in list 3 of leaf */
void Node::pushSelfToNearNonNbors() {
    if (leafIlist.empty()) return;

    // if (isNodeType<Stem>()) std::cout << leafIlist.size() << "\n";

    auto self = getSelf(); // call shared_from_this()
    for (const auto& node : leafIlist) {
        auto leaf = dynamic_pointer_cast<Leaf>(node);
        leaf->pushToNearNonNbors(self);
    }
}

const std::vector<vecXcd> Node::getShiftedLocalCoeffs(const int branchIdx) const {
    std::vector<vecXcd> shiftedLocalCoeffs, rotatedLocalCoeffs;
    const double r = (branches[branchIdx]->center - center).norm();

    // apply rotation (rotation axis is opposite from M2M)
    for (int j = 0; j <= order; ++j)
        rotatedLocalCoeffs.push_back(wignerD[7-branchIdx][j] * localCoeffs[j]);

    for (int j = 0; j <= order; ++j) {
        vecXcd shiftedLocalCoeffs_j = vecXcd::Zero(2*j+1);
        for (int k = -j; k <= j; ++k) {
            int k_j = k + j;
            double r2nmj = 1.0;
            for (int n = j; n <= order; ++n) {
                shiftedLocalCoeffs_j[k_j] +=
                    rotatedLocalCoeffs[n][k+n] *
                    tables.A_[n-j][n-j] * tables.A_[j][k_j] / tables.A_[n][k+n]
                    * r2nmj / pm(n+j);
                r2nmj *= r;
            }
        }

        // apply inverse rotation
        shiftedLocalCoeffs.push_back(wignerDInv[7-branchIdx][j] * shiftedLocalCoeffs_j);
    }
    return shiftedLocalCoeffs;
}

void Node::addToLocalCoeffsFromLeafIlist() {
    /* if # observers is small, evaluate sol there directly */
    if (particles.size() <= order*order) {
        for (const auto& obs : particles)
            for (const auto& iNode : leafIlist) {
                obs->addToPhi(iNode->getDirectPhi(obs->getPos()));
                obs->addToFld(iNode->getDirectFld(obs->getPos()));
            }
        return;
    }

    /* otherwise, add to local expansion due to srcs in list 4 */
    for (const auto& node : leafIlist) {
        for (const auto& src : node->particles){
            auto dR = toSph(src->getPos() - center);
            double r = dR[0], th = dR[1], ph = dR[2];

            double r2lpp = r;
            for (int l = 0; l <= order; ++l) {
                realVec legendreCoeffs;

                for (int m = 0; m <= l; ++m)
                    legendreCoeffs.push_back(legendreCos(th, l, m));

                for (int m = -l; m <= l; ++m){
                    localCoeffs[l][m+l] +=
                        src->getCharge() / r2lpp *
                        legendreCoeffs[abs(-m)] * expI(-m*ph);
                }
                r2lpp *= r;
            }
        }
    }
}

// return phi at X due to all particles in this node
const double Node::getDirectPhi(const vec3d& X) {
    double phi = 0;
    for (const auto& src : particles) {
        auto r = (X - src->getPos()).norm();
        if (r > 1.0E-9)
            phi += src->getCharge() / r;
    }
    return phi;
}

// return phi at all particles in this node due to all other particles in this node
const realVec Node::getDirectPhis() {
    realVec phis;

    for (const auto& obs : particles) {
        phis.push_back(getDirectPhi(obs->getPos()));
    }
    return phis;
}

// return fld at X due to all particles in this node
const vec3d Node::getDirectFld(const vec3d& X) {
    vec3d fld = vec3d::Zero();
    for (const auto& src : particles) {
        auto dX = X - src->getPos();
        auto r = dX.norm();
        if (r > 1.0E-9)
            fld += src->getCharge() * dX / pow(r, 3);
    }
    return fld;
}

// return phi at all particles in this node due to all other particles in this node
const std::vector<vec3d> Node::getDirectFlds() {
    std::vector<vec3d> flds;

    for (const auto& obs : particles) {
        flds.push_back(getDirectFld(obs->getPos()));
    }
    return flds;
}