#include "leaf.h"

Leaf::Leaf(
    const ParticleVec& particles,
    const int branchIdx,
    Stem* const base)
    : Node(particles, branchIdx, base)
{
}

void Leaf::buildNbors() {
    assert(!isRoot());

    for (int i = 0; i < numDir; ++i) {
        Dir dir = static_cast<Dir>(i);
        auto nbor = getNeighborGeqSize(dir);
        if (nbor != nullptr) {
            nbors.push_back(nbor);
            auto nbors = getNeighborsLeqSize(nbor, dir);
            nearNbors.reserve(nearNbors.size() + nbors.size());
            nearNbors.insert(nearNbors.end(), nbors.begin(), nbors.end());
        }
    }
    assert(nbors.size() <= numDir);
}

void Leaf::buildLists() {
    if (isRoot()) return;
    
    buildNbors();
    buildInteractionList();
    pushSelfToNearNonNbors();
}

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
            for (const auto& iNode : dirList[dir])
                iNode->addShiftedExpCoeffs(expCoeffsOut[dir], center, dir);
    }
    t_X2X += std::chrono::high_resolution_clock::now() - start;

    expCoeffsOut = {};
}

void Leaf::buildLocalCoeffs() {

    if (!isRoot()) {
        auto start = std::chrono::high_resolution_clock::now();

        addToLocalCoeffsFromLeafIlist();

        t_P2L += std::chrono::high_resolution_clock::now() - start;

        start = std::chrono::high_resolution_clock::now();

        addToLocalCoeffsFromDirList();

        t_X2L += std::chrono::high_resolution_clock::now() - start;

        start = std::chrono::high_resolution_clock::now();

        if (!base->isRoot()) {
            auto shiftedLocalCoeffs = base->getShiftedLocalCoeffs(branchIdx);
            for (int l = 0; l <= order; ++l)
                localCoeffs[l] += shiftedLocalCoeffs[l];
        }

        t_L2L += std::chrono::high_resolution_clock::now() - start;
    }

    evaluateSolAtParticles();
}

/* from local expansion */
solVec Leaf::getFarSols() const {
    solVec sols;

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

        sols.push_back(pairSol(phi.real(), fld.real()));
    }

    return sols;
}

/* from particles in near non-neighbors (list 3) */
solVec Leaf::getNearNonNborSols() const {
    solVec sols;

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
                    assert(obs != src);

                    const auto dX = obsPos - src->getPos();
                    const auto dist = dX.norm();
                    const auto charge = src->getCharge();

                    phi += charge / dist;
                    fld += charge * dX / pow(dist, 3);
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

        sols.push_back(pairSol(phi.real(), fld.real()));
    }

    return sols;
}

/* from particles in near neighbors (list 1) */
solVec Leaf::getNearNborSols() const {
    solVec sols;

    for (const auto& obs : particles) {
        double phi = 0.0;
        vec3d fld = vec3d::Zero();

        auto obsPos = obs->getPos();

        for (const auto& node : nearNbors) {
            auto srcs = node->getParticles();

            for (const auto& src : node->getParticles()) {
                assert(obs != src);

                vec3d dX = obs->getPos() - src->getPos();
                auto dist = dX.norm();
                auto charge = src->getCharge();

                phi += charge / dist;
                fld += charge * dX / pow(dist, 3);
            }
        }
        sols.push_back(pairSol(phi,fld));
    }
    return sols;
}

/* from other particles in this node */
solVec Leaf::getSelfSols() const {
    solVec sols;

    for (const auto& obs : particles) {
        double phi = 0.0;
        vec3d fld = vec3d::Zero();

        auto obsPos = obs->getPos();

        // due to other particles in this node (implement reciprocity later)
        for (const auto& src : particles) {
            if (src == obs) continue;

            vec3d dX = obs->getPos() - src->getPos();
            auto dist = dX.norm();
            auto charge = src->getCharge();

            phi += charge / dist;
            fld += charge * dX / pow(dist, 3);

        }
        sols.push_back(pairSol(phi,fld));
    }
    return sols;
}

// pass calculated phi and fld to particles
void Leaf::evaluateSolAtParticles() {
    if (isRoot()) return; // fix later

    auto start = std::chrono::high_resolution_clock::now();

    auto sols = 
        getFarSols()
        + getNearNonNborSols()
        ;

    t_L2P += std::chrono::high_resolution_clock::now() - start;

    start = std::chrono::high_resolution_clock::now();

    sols = sols
        + getNearNborSols()
        + getSelfSols();

    t_dir += std::chrono::high_resolution_clock::now() - start;

    for (size_t n = 0; n < particles.size(); ++n) {
        auto p = particles[n];
        p->addFromSol(sols[n]);
    }
}

/*
void Leaf::evaluateSolAtParticles() {
    if (isRoot()) return; // fix later

    auto start = std::chrono::high_resolution_clock::now();

    auto phis = 
        getFarPhis() 
        + getNearNonNborPhis()
        ;

    auto flds =
        getFarFlds()
        + getNearNonNborFlds()
        ;

    t_L2P += std::chrono::high_resolution_clock::now() - start;

    start = std::chrono::high_resolution_clock::now();

    phis = phis 
        + getNearNborSols<double>(
        [](vec3d X) { return 1.0 / X.norm(); });

    flds = flds 
        + getNearNborSols<vec3d>(
        [](vec3d X) { return X / pow(X.norm(), 3); });

    t_dir += std::chrono::high_resolution_clock::now() - start;

    // for (auto [p, phi, fld] : std::views::zip(particles, phis, flds)) {
    for (size_t n = 0; n < particles.size(); ++n) {
        auto p = particles[n];
        p->addToPhi(phis[n]);
        p->addToFld(flds[n]);
    }
}*/