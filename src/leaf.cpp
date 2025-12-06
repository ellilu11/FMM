#include "leaf.h"

LeafVec Leaf::leaves;

Leaf::Leaf(
    const ParticleVec& particles,
    const int branchIdx,
    Stem* const base)
    : Node(particles, branchIdx, base)
{
}

/* buildNeighbors()
 * Find all neighbor nodes of equal or greater size
 * Also find all neighbor leaves of equal or lesser size (list 1)
 */
void Leaf::buildNeighbors() {
    assert(!isRoot());

    for (int i = 0; i < numDir; ++i) {
        Dir dir = static_cast<Dir>(i);
        auto nbor = getNeighborGeqSize(dir);

        if (nbor != nullptr) {
            nbors.push_back(nbor);
            auto nbors = getNeighborsLeqSize(nbor, dir);
            nearNbors.insert(nearNbors.end(), nbors.begin(), nbors.end());
        }
    }

    assert(nbors.size() <= numDir);
}

/* buildLists()
 * Add self to list of leaves 
 * Find neighbor and interaction lists
 * Add self as near non-neighbor (list 3 node) of any list 4 nodes
 */
void Leaf::buildLists() {
    leaves.push_back(shared_from_this()); 

    if (isRoot()) return;
    
    buildNeighbors();

    buildInteractionList();

    pushSelfToNearNonNbors();
}

/* buildMpoleCoeffs()
 * (P2M) Build mpole expansions from particles in this node  
 */
void Leaf::buildMpoleCoeffs() {
    if (isRoot()) return;

    for (int l = 0; l <= order; ++l)
        coeffs.push_back(vecXcd::Zero(2*l+1));

    for (const auto& src : particles) {
        const auto dR = Math::toSph(src->getPos() - center);
        const double r = dR[0], th = dR[1], ph = dR[2];
        double r2l = 1.0;

        for (int l = 0; l <= order; ++l) {
            realVec legendre_l;
            for (int m = 0; m <= l; ++m)
                legendre_l.push_back(legendreCos(th, l, m));

            for (int m = -l; m <= l; ++m) {
                coeffs[l][m+l] +=
                    src->getCharge() * r2l *
                    legendre_l[abs(-m)] * Math::expI(-m*ph);
            }

            r2l *= r;
        }
    }
}

/* propagateExpCoeffs()
 * (M2X) Convert mpole coeffs into exp coeffs
 * (X2X) Translate exp coeffs to nodes in all dirlists
 */
void Leaf::propagateExpCoeffs() {
    if (isRoot()) return;

    for (int dir = 0; dir < 6; ++dir){
        
        auto expCoeffs = getMpoleToExpCoeffs(dir);

        auto iList = dirList[dir];

        for (const auto& iNode : iList)
            iNode->addShiftedExpCoeffs(expCoeffs, center, dir);
    }
}

/* buildLocalCoeffs()
 * (X2L) Receive incoming exp coeffs and add to local coeffs
 * (P2L) Add contribution from list 4 nodes t to local coeffs
 * (L2L) Shift base local coeffs to center and add to local coeffs
 */
void Leaf::buildLocalCoeffs() {
    if (isRoot()) return;

    evalLeafIlistSols();

    evalExpToLocalCoeffs();
        
    if (!base->isRoot()) {
        auto shiftedLocalCoeffs = base->getShiftedLocalCoeffs(branchIdx);

        for (int l = 0; l <= order; ++l)
            localCoeffs[l] += shiftedLocalCoeffs[l];
    }
}

/* evalFarSols()
 * (L2P) Evaluate sols from local expansion due to far nodes
 */
void Leaf::evalFarSols() {

    for (const auto& obs : particles) {
        const auto dR = Math::toSph(obs->getPos() - center);
        const double r = dR[0], th = dR[1], ph = dR[2];

        cmplx phi(0,0);
        vec3cd fld = vec3cd::Zero();
        double r2l = 1.0;

        for (int l = 0; l <= order; ++l) {
            realVec legendre_l, dLegendre_l;

            for (int m = 0; m <= l; ++m) {
                legendre_l.push_back(legendreCos(th, l, m));
                dLegendre_l.push_back(dLegendreCos(th, l, m));
            }

            for (int m = -l; m <= l; ++m) {
                const size_t abs_m = abs(m);
                const cmplx coeff = localCoeffs[l][m+l] * r2l * Math::expI(m*ph);

                phi += coeff * legendre_l[abs_m];
                fld -= coeff / r *
                    vec3cd(
                        l * legendre_l[abs_m],           // E_r
                        dLegendre_l[abs_m],              // E_th
                        cmplx(0,m) * legendre_l[abs_m]); // E_ph * sin(th)
            }

            r2l *= r;
        }

        // Convert to cartesian components
        fld = Math::matFromSph(th,ph) * fld;

        obs->addToSol(phi.real(), fld.real());
    }
}

/* evalNearNonNborSols()
 * (M2P) Evaluate sols from mpole expansion due to list 3 nodes
 */
void Leaf::evalNearNonNborSols() {

    for (const auto& node : nearNonNbors) {

        const auto srcs = node->getParticles();

        if (srcs.size() <= order*order)
            // Do nothing! Contribution from list 3 node was 
            // already evaluated by Node::evalLeafIlistSols()
            continue;

        // # srcs large, use mpole expansion of list 3 node
        for (const auto& obs : particles) {
            const auto obsPos = obs->getPos();

            cmplx phi(0,0);
            vec3cd fld = vec3cd::Zero();
            vec3cd fld_R = fld;

            const auto srcCoeffs = node->getMpoleCoeffs();

            const auto dR = Math::toSph(obsPos - node->getCenter());
            const double r = dR[0], th = dR[1], ph = dR[2];
            double r2lpp = r;

            for (int l = 0; l <= order; ++l) {
                realVec legendre_l, dLegendre_l;

                for (int m = 0; m <= l; ++m) {
                    legendre_l.push_back(legendreCos(th, l, m));
                    dLegendre_l.push_back(dLegendreCos(th, l, m));
                }

                for (int m = -l; m <= l; ++m) {
                    const size_t abs_m = abs(m);
                    const cmplx coeff = srcCoeffs[l][m+l] / r2lpp * Math::expI(m*ph);

                    phi += coeff * legendre_l[abs_m];
                    fld_R -= coeff / r * 
                        vec3cd(
                            -(l+1) * legendre_l[abs_m],      // E_r
                            dLegendre_l[abs_m],              // E_th
                            cmplx(0,m) * legendre_l[abs_m]); // E_ph * sin(th)
                }

                r2lpp *= r;
            }

            // Convert to cartesian components
            fld += Math::matFromSph(th, ph) * fld_R;

            obs->addToSol(phi.real(), fld.real());
        }
    }
}

/* findNearNborPairs()
 * From list of leaves, find all near neighbor leaf pairs
 */
std::vector<LeafPair> Leaf::findNearNborPairs(){
    std::vector<LeafPair> leafPairs;

    for (const auto& leaf : leaves) {
        for (const auto& nbor : leaf->nearNbors) {
            auto nborLeaf = dynamic_pointer_cast<Leaf>(nbor);
            if (leaf < nborLeaf)
                leafPairs.emplace_back(leaf, nborLeaf);
        }
    }

    return leafPairs;
}

/* evaluateSols()
 * Sum solutions at all particles in all leaves 
 */ 
void Leaf::evaluateSols() {

    for (const auto& leaf : leaves) {
        leaf->evalFarSols();

        leaf->evalNearNonNborSols();

        leaf->evalSelfSols();
    }

    for (const auto& pair : findNearNborPairs()) {
        auto [obsLeaf, srcLeaf] = pair;
        obsLeaf->evalPairSols(srcLeaf);
    }

}