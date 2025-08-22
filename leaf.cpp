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
    // add self to list of leaves
    leaves.push_back(shared_from_this()); 

    if (isRoot()) return;
    
    buildNeighbors();

    buildInteractionList();

    pushSelfToNearNonNbors();
}

/* buildMpoleCoeffs : 
   Build mpole expansions from particles in this node (P2M) */
void Leaf::buildMpoleCoeffs() {
    if (isRoot()) return;

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

    for (int dir = 0; dir < 6; ++dir){
        auto start = Clock::now();
        
        auto expCoeffs = getMpoleToExpCoeffs(dir);
        
        t.M2X += Clock::now() - start;

        auto iList = dirList[dir];

        start = Clock::now();
        
        for (const auto& iNode : iList)
            iNode->addShiftedExpCoeffs(expCoeffs, center, dir);
        
        t.X2X += Clock::now() - start;
    }
}

void Leaf::buildLocalCoeffs() {
    if (isRoot()) return;
    
    auto start = Clock::now();

    evalLocalCoeffsFromLeafIlist();

    t.P2L += Clock::now() - start;

    start = Clock::now();

    evalLocalCoeffsFromDirList();

    t.X2L += Clock::now() - start;

    start = Clock::now();
        
    if (!base->isRoot()) {
        auto shiftedLocalCoeffs = base->getShiftedLocalCoeffs(branchIdx);
        for (int l = 0; l <= order; ++l)
            localCoeffs[l] += shiftedLocalCoeffs[l];
    }
        
    t.L2L += Clock::now() - start;
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
   Get sols from mpole expansion due to near non-neighbor nodes (list 3) */
void Leaf::evalNearNonNborSols() {

    for (const auto& obs : particles) {
        const auto obsPos = obs->getPos();

        cmplx phi(0,0);
        vec3cd fld = vec3cd::Zero();
        vec3cd fld_R = fld;

        for (const auto& node : nearNonNbors) {

            const auto srcs = node->getParticles();
            if (srcs.size() <= order*order)
                // do nothing since contribution from list 3 node was 
                // calculated by Node::evalLocalCoeffsFromLeafIlist()
                continue;

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
   From list of leaves, find all near neighbor leaf pairs */
std::vector<LeafPair> Leaf::findNearNborPairs(){
    std::vector<LeafPair> leafPairs;

    for (const auto& leaf : leaves) {
        for (const auto& nbor : leaf->nearNbors) {
            auto nborLeaf = dynamic_pointer_cast<Leaf>(nbor);
            if (leaf < nborLeaf)
                leafPairs.emplace_back(std::make_pair(leaf, nborLeaf));
        }
    }

    return leafPairs;
}

/* evaluateSols: 
   Sum solutions at all particles in all leaves */ 
void Leaf::evaluateSols() {

    auto start = Clock::now();

    for (const auto& leaf : leaves) {
        leaf->evalFarSols();
        leaf->evalNearNonNborSols();
    }

    t.L2P += Clock::now() - start;

    start = Clock::now();

    for (const auto& pair : findNearNborPairs()) {
        auto [obsLeaf, srcLeaf] = pair;
        const bool evalAtSrcs = 1;
        obsLeaf->evalPairSols(srcLeaf);
    }

    for (const auto& leaf : leaves)
        leaf->evalSelfSols();

    t.P2P += Clock::now() - start;
}

