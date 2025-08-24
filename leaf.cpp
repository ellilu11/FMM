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
        const auto dR = toSph(src->getPos() - center);
        const double r = dR[0], th = dR[1], ph = dR[2];
        double r2l = 1.0;

        for (int l = 0; l <= order; ++l) {
            realVec legendreCoeffs;
            for (int m = 0; m <= l; ++m)
                legendreCoeffs.push_back(legendreCos(th, l, m));

            for (int m = -l; m <= l; ++m) {
                coeffs[l][m+l] +=
                    src->getCharge() * r2l *
                    legendreCoeffs[abs(-m)] * expI(-m*ph);
            }
            r2l *= r;
        }
    }
}

/*void Leaf::propagateExpCoeffs() {
    if (isRoot()) return;

    auto start = std::chrono::high_resolution_clock::now();
    for (int dir = 0; dir < 6; ++dir) {
        // for lvl > 1, propagate own exp coeffs to inner dirlist
        if (!base->isRoot())
            for (const auto& node : dirList[dir])
                node->addShiftedExpCoeffs(expCoeffsOut[dir], center, dir);
    }
    t.X2X += std::chrono::high_resolution_clock::now() - start;

    expCoeffsOut = {};
}*/

/* propagateExpCoeffs()
 * (M2X) Convert mpole coeffs into exp coeffs
 * (X2X) Translate exp coeffs to nodes in all dirlists
 */
void Leaf::propagateExpCoeffs() {
    if (isRoot()) return;

    Clock::time_point start;

    for (int dir = 0; dir < 1; ++dir) {
        start = Clock::now();

        if (!base->isRoot())
            for (const auto& node : dirList[dir])
                node->addShiftedExpCoeffs(expCoeffsOut[dir], center, dir);

        t.X2X += Clock::now() - start;
    }

    expCoeffsOut = {};

    // Puzzle: How to combine coordinate rotation to north/south/east/westlist
    // with translation to base?
    // TODO: Remove this section once puzzle solved
    for (int dir = 1; dir < 6; ++dir) {
        start = Clock::now();

        auto expCoeffs = getMpoleToExpCoeffs(dir);

        t.M2X += Clock::now() - start;

        start = Clock::now();

        if (!base->isRoot())
            for (const auto& node : dirList[dir])
                node->addShiftedExpCoeffs(expCoeffs, center, dir);

        t.X2X += Clock::now() - start;

    }
}

/*void Leaf::propagateExpCoeffs() {
    if (isRoot()) return;

    for (int dir = 0; dir < 6; ++dir) {

        auto start = std::chrono::high_resolution_clock::now();

        auto expCoeffs = getMpoleToExpCoeffs(dir);

        t_M2X += std::chrono::high_resolution_clock::now() - start;

        start = std::chrono::high_resolution_clock::now();
        // for lvl > 1, propagate own exp coeffs to inner dirlist
        if (!base->isRoot())
            for (const auto& node : dirList[dir])
                node->addShiftedExpCoeffs(expCoeffs, center, dir);

        t.X2X += std::chrono::high_resolution_clock::now() - start;
    }
 

    expCoeffsOut = {};
}*/

void Leaf::buildLocalCoeffs() {
    if (isRoot()) return;

    auto start = Clock::now();

    evalLeafIlistSols();

    t.P2L += Clock::now() - start;

    start = Clock::now();

    evalExpToLocalCoeffs();

    t.X2L += Clock::now() - start;

    start = Clock::now();

    if (!base->isRoot()) {
        auto shiftedLocalCoeffs = base->getShiftedLocalCoeffs(branchIdx);

        for (int l = 0; l <= order; ++l)
            localCoeffs[l] += shiftedLocalCoeffs[l];
    }

    t.L2L += Clock::now() - start;
}

/* evalFarSols()
 * (L2P) Evaluate sols from local expansion due to far nodes
 */
void Leaf::evalFarSols() {

    for (const auto& obs : particles) {
        const auto dR = toSph(obs->getPos() - center);
        const double r = dR[0], th = dR[1], ph = dR[2];

        cmplx phi(0,0);
        vec3cd fld = vec3cd::Zero();
        double r2l = 1.0;

        for (int l = 0; l <= order; ++l) {
            realVec legendre_l, dLegendre_l;

            //for (int m = 0; m <= l; ++m) {
            //    double legendre_lm = legendreCos(th, l, m);
            //    legendre_l.push_back(legendre_lm);
            //    double dLegendre_prev_m =
            //        (m < l) ? legendre_prev[m] : 0.0;
            //    dLegendre_l.push_back(
            //        dLegendreCos(th, l, m,
            //            legendre_lm, dLegendre_prev_m));
            //}

            for (int m = 0; m <= l; ++m) {
                legendre_l.push_back(legendreCos(th, l, m));
                dLegendre_l.push_back(dLegendreCos(th, l, m));
            }

            for (int m = -l; m <= l; ++m) {
                const size_t abs_m = abs(m);
                const cmplx coeff = localCoeffs[l][m+l] * r2l * expI(m*ph);

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
        fld = matFromSph(th,ph) * fld;

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

            const auto dR = toSph(obsPos - node->getCenter());
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
                    const cmplx coeff = srcCoeffs[l][m+l] / r2lpp * expI(m*ph);

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
            fld += matFromSph(th, ph) * fld_R;

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

    auto start = Clock::now();

    for (const auto& leaf : leaves)
        leaf->evalFarSols();

    t.L2P += Clock::now() - start;

    start = Clock::now();

    for (const auto& leaf : leaves)
        leaf->evalNearNonNborSols();

    t.M2P += Clock::now() - start;

    start = Clock::now();

    for (const auto& pair : findNearNborPairs()) {
        auto [obsLeaf, srcLeaf] = pair;
        obsLeaf->evalPairSols(srcLeaf);
    }

    for (const auto& leaf : leaves)
        leaf->evalSelfSols();

    t.P2P += Clock::now() - start;
}