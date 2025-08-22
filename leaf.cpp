#include "leaf.h"

LeafVec Leaf::leaves;

Leaf::Leaf(
    const ParticleVec& particles,
    const int branchIdx,
    Stem* const base)
    : Node(particles, branchIdx, base)
{
}

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

void Leaf::buildLists() {
    if (isRoot()) return;
    
    buildNeighbors();

    buildInteractionList();

    pushSelfToNearNonNbors();

    leaves.push_back(shared_from_this()); // record self in list of leaves
}

/* buildMpoleCoeffs : 
   Build mpole expansions from particles in this node (P2M) */
void Leaf::buildMpoleCoeffs() {
    for (int l = 0; l <= order; ++l)
        coeffs.emplace_back(vecXcd::Zero(2*l+1));

    for (const auto& src : particles) {
        auto dR = toSph(src->getPos() - center);
        double r = dR[0], th = dR[1], ph = dR[2];
        double r2l = 1.0;

        for (int l = 0; l <= order; ++l) {
            realVec legendreCoeffs;
            for (int m = 0; m <= l; ++m)
                legendreCoeffs.push_back(legendreCos(th,l,m));

            for (int m = -l; m <= l; ++m) {
                coeffs[l][m+l] +=
                    src->getCharge() * r2l *
                    legendreCoeffs[abs(-m)] * expI(-m*ph); 
            }
            r2l *= r;
        }
    }
}

void Leaf::propagateExpCoeffs() {
    if (isRoot()) return;

    auto start = std::chrono::high_resolution_clock::now();
    for (int dir = 0; dir < 6; ++dir) {
        // for lvl > 1, propagate own exp coeffs to inner dirlist
        if (!base->isRoot())
            for (const auto& node : dirList[dir])
                node->addShiftedExpCoeffs(expCoeffsOut[dir], center, dir);
    }
    t_X2X += std::chrono::high_resolution_clock::now() - start;

    expCoeffsOut = {};
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

        t_X2X += std::chrono::high_resolution_clock::now() - start;
    }
 

    expCoeffsOut = {};
}*/

void Leaf::buildLocalCoeffs() {
    if (!isRoot()) {
        auto start = std::chrono::high_resolution_clock::now();

        evalLocalCoeffsFromLeafIlist();

        t_P2L += std::chrono::high_resolution_clock::now() - start;

        start = std::chrono::high_resolution_clock::now();

        evalLocalCoeffsFromDirList();

        t_X2L += std::chrono::high_resolution_clock::now() - start;

        start = std::chrono::high_resolution_clock::now();

        if (!base->isRoot()) {
            auto shiftedLocalCoeffs = base->getShiftedLocalCoeffs(branchIdx);
            for (int l = 0; l <= order; ++l)
                localCoeffs[l] += shiftedLocalCoeffs[l];
        }

        t_L2L += std::chrono::high_resolution_clock::now() - start;
    }
}

/* evalFarSols : 
   Evaluate sols from local expansion due to far nodes (L2P) */
void Leaf::evalFarSols() {

    for (const auto& obs : particles) {
        const auto dR = toSph(obs->getPos() - center);
        const double r = dR[0], th = dR[1], ph = dR[2];
        double r2l = 1.0, r2lmm = 1.0 / r;

        cmplx phi(0,0);
        vec3cd fld = vec3cd::Zero();

        for (int l = 0; l <= order; ++l) {
            realVec legendreCoeffs, dLegendreCoeffs;

            for (int m = 0; m <= l; ++m) {
                legendreCoeffs.push_back(legendreCos(th,l,m));
                dLegendreCoeffs.push_back(dLegendreCos(th,l,m));
            }

            for (int m = -l; m <= l; ++m) {
                const int abs_m = abs(m);

                phi += localCoeffs[l][m+l] * r2l *
                    legendreCoeffs[abs(m)] * expI(m*ph);

                fld -= localCoeffs[l][m+l] * r2lmm * expI(m*ph) *
                    vec3cd(
                        l * legendreCoeffs[abs_m],           // E_r
                        dLegendreCoeffs[abs_m],              // E_th
                        cmplx(0,m) * legendreCoeffs[abs_m]); // E_ph * sin(th)
;
            }
            r2l *= r;
            r2lmm *= r;
        }

        // Convert to cartesian components
        fld = mat3d{
            {  sin(th)*cos(ph),  cos(th)*cos(ph), -sin(ph)/sin(th) },
            {  sin(th)*sin(ph),  cos(th)*sin(ph),  cos(ph)/sin(th) },
            {  cos(th),         -sin(th),          0.0             }
        } * fld;

        obs->addToSol(phi.real(), fld.real());
    }
}

/* evalNearNonNborSols : 
   Get sols directly or from mpole expansion due to near non-neighbor nodes (list 3) */
void Leaf::evalNearNonNborSols() {

    for (const auto& obs : particles) {
        const auto obsPos = obs->getPos();

        cmplx phi(0,0);
        vec3cd fld = vec3cd::Zero();
        vec3cd fld_R = fld;

        for (const auto& node : nearNonNbors) {

            // # srcs small, do direct
            const auto srcs = node->getParticles();
            if (srcs.size() <= order*order) {
                for (const auto& src : srcs) {
                    // assert(obs != src);

                    const auto dX = obsPos - src->getPos();
                    const auto dr = dX.norm();
                    const auto srcPhi = src->getCharge() / dr;

                    phi += srcPhi;
                    fld += srcPhi * dX / (dr*dr);
                }
                continue;
            }

            // # srcs large, use mpole expansion of src node
            const auto srcCoeffs = node->getMpoleCoeffs();

            const auto dR = toSph(obsPos - node->getCenter());
            const double r = dR[0], th = dR[1], ph = dR[2];
            double r2lpp = r, r2lp2 = r*r;

            for (int l = 0; l <= order; ++l) {
                realVec legendreCoeffs, dLegendreCoeffs;

                for (int m = 0; m <= l; ++m) {
                    legendreCoeffs.push_back(legendreCos(th, l, m));
                    dLegendreCoeffs.push_back(dLegendreCos(th, l, m));
                }

                for (int m = -l; m <= l; ++m) {
                    const int abs_m = abs(m);

                    phi += srcCoeffs[l][m+l] / r2lpp *
                        legendreCoeffs[abs(m)] * expI(m*ph);

                    fld_R -=
                        srcCoeffs[l][m+l] / r2lp2 * expI(m*ph) *
                        vec3cd(
                            -(l+1) * legendreCoeffs[abs_m],      // E_r
                            dLegendreCoeffs[abs_m],              // E_th
                            cmplx(0,m) * legendreCoeffs[abs_m]); // E_ph * sin(th)
                }
                r2lpp *= r;
                r2lp2 *= r;
            }

            // Convert to cartesian components
            fld += mat3d{
                {  sin(th)*cos(ph),  cos(th)*cos(ph), -sin(ph)/sin(th) },
                {  sin(th)*sin(ph),  cos(th)*sin(ph),  cos(ph)/sin(th) },
                {  cos(th),         -sin(th),          0.0             }
            } * fld_R;
        }

        obs->addToSol(phi.real(), fld.real());
    }
}

/* findNearNborPairs :
   From list of leaves, find all near neighbor pairs */
std::vector<LeafPair> Leaf::findNearNborPairs(){
    std::vector<LeafPair> leafPairs;

    for (const auto& leaf : leaves) {
        for (const auto& nbor : leaf->getNearNbors()) {
            auto nborLeaf = dynamic_pointer_cast<Leaf>(nbor);
            if (leaf < nbor)
                leafPairs.emplace_back(std::make_pair(leaf, nborLeaf));
        }
    }
    //std::cout << "# leaf pairs: " << leafPairs.size() << '\n';

    return leafPairs;
}

/* evaluateSols: 
   Sum solutions at all particles in all leaf nodes */ 
void Leaf::evaluateSols() {

    auto start = std::chrono::high_resolution_clock::now();

    for (const auto& leaf : leaves) {
        leaf->evalFarSols();
        leaf->evalNearNonNborSols();
    }

    t_L2P += std::chrono::high_resolution_clock::now() - start;

    start = std::chrono::high_resolution_clock::now();

    for (const auto& pair : findNearNborPairs()) {
        auto [obsNode, srcNode] = pair;
        const bool evalAtSrcs = 1;
        obsNode->evalDirectSols(srcNode, evalAtSrcs);
    }

    for (const auto& leaf : leaves)
        leaf->evalSelfSols();

    t_dir += std::chrono::high_resolution_clock::now() - start;
}

