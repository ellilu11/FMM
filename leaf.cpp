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

        buildLocalCoeffsFromLeafIlist();

        t_X2L_l4 += std::chrono::high_resolution_clock::now() - start;

        start = std::chrono::high_resolution_clock::now();

        buildLocalCoeffsFromDirList();

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

realVec Leaf::getFarPhis() const {
    realVec phis;
    for (const auto& obs : particles) {
        auto dR = toSph(obs->getPos() - center);
        double r = dR[0], th = dR[1], ph = dR[2];

        // using Horner's scheme
        /*cmplxVec polyCoeffs;
        for (int l = 0; l <= order; ++l) {
            realVec legendreCoeffs;
            cmplx polyCoeffs_l;
            for (int m = 0; m <= l; ++m)
                legendreCoeffs.push_back(legendreCos(th, pair2i(l,m)));

            for (int m = -l; m <= l; ++m)
                polyCoeffs_l += localCoeffs[l][m+l] *
                    legendreCoeffs[std::abs(m)] * expI(static_cast<double>(m)*ph);
            polyCoeffs.push_back(polyCoeffs_l);
        }
        phis.push_back(evaluatePoly<cmplx>(polyCoeffs,r).real());*/

        cmplx phi(0, 0);
        double r2l = 1.0;
        for (int l = 0; l <= order; ++l) {
            realVec legendreCoeffs;
            for (int m = 0; m <= l; ++m)
                legendreCoeffs.push_back(legendreCos(th,l,m));

            for (int m = -l; m <= l; ++m)
                phi += localCoeffs[l][m+l] * r2l *
                        legendreCoeffs[abs(m)] * expI(m*ph);
            r2l *= r;
        }
        phis.push_back(phi.real());
    }

    return phis;
}

std::vector<vec3d> Leaf::getFarFlds() const {
    std::vector<vec3d> flds;
    for (const auto& obs : particles) {
        const auto dR = toSph(obs->getPos() - center);
        const double r = dR[0], th = dR[1], ph = dR[2];

        vec3cd fld_R = vec3cd::Zero();
        double r2lmm = 1.0 / r;
        for (int l = 0; l <= order; ++l) {
            realVec legendreCoeffs, dLegendreCoeffs;
            for (int m = 0; m <= l; ++m) {
                legendreCoeffs.push_back(legendreCos(th, l, m));
                dLegendreCoeffs.push_back(dLegendreCos(th, l, m));
            }

            for (int m = -l; m <= l; ++m) {
                const int abs_m = abs(m);
                fld_R -=
                    localCoeffs[l][m+l] * r2lmm * expI(m*ph) *
                    vec3cd(
                        l * legendreCoeffs[abs_m],             // E_r
                        dLegendreCoeffs[abs_m],                // E_th
                        cmplx(0,m) * legendreCoeffs[abs_m]);   // E_ph * sin(th)
            }
            r2lmm *= r;
        }

        // Convert to cartesian components
        vec3cd fld_X = mat3d{
            {  sin(th)*cos(ph),  cos(th)*cos(ph), -sin(ph)/sin(th) },
            {  sin(th)*sin(ph),  cos(th)*sin(ph),  cos(ph)/sin(th) },
            {  cos(th),         -sin(th),          0.0             }
        } * fld_R;

        flds.push_back(fld_X.real());
    }
    return flds;
}

// due to particles in near non-neighbors (list 3)
realVec Leaf::getNearNonNborPhis() const {
    realVec phis;

    for (const auto& obs : particles) {
        cmplx phi(0, 0);
        auto obsPos = obs->getPos();

        for (const auto& node : nearNonNbors) {
            auto srcs = node->getParticles();

            // # srcs small, do direct
            if (srcs.size() <= order*order) {
                for (const auto& src : srcs)
                    phi += src->getCharge() / (obsPos - src->getPos()).norm();
                continue;
            }

            // # srcs large, use mpole expansion of src node
            auto dR = toSph(obsPos - node->getCenter());
            auto srcCoeffs = node->getMpoleCoeffs();

            double r = dR[0], th = dR[1], ph = dR[2];

            double r2lpp = r;
            for (int l = 0; l <= order; ++l) {
                realVec legendreCoeffs;
                for (int m = 0; m <= l; ++m)
                    legendreCoeffs.push_back(legendreCos(th, l, m));

                for (int m = -l; m <= l; ++m) {
                    int m_ = m + l;
                    phi += srcCoeffs[l][m_] / r2lpp *
                        legendreCoeffs[abs(m)] * expI(m*ph);

                }
                r2lpp *= r;
            }
        }
        phis.push_back(phi.real());
    }
    return phis;
}

std::vector<vec3d> Leaf::getNearNonNborFlds() const {
    std::vector<vec3d> flds;

    for (const auto& obs : particles) {
        vec3cd fld_X = vec3cd::Zero();
        vec3cd fld_R = fld_X;
        auto obsPos = obs->getPos();

        for (const auto& node : nearNonNbors) {
            auto srcs = node->getParticles();

            // # srcs small, do direct
            if (srcs.size() <= order*order) {
                for (const auto& src : srcs) {
                    auto dX = obsPos - src->getPos();
                    fld_X += src->getCharge() * dX / pow(dX.norm(), 3);
                }
                continue;
            }

            // # srcs large, use mpole expansion of src node
            auto srcCoeffs = node->getMpoleCoeffs();

            auto dR = toSph(obsPos - node->getCenter());
            double r = dR[0], th = dR[1], ph = dR[2];

            double r2lp2 = r*r;
            for (int l = 0; l <= order; ++l) {
                realVec legendreCoeffs, dLegendreCoeffs;
                for (int m = 0; m <= l; ++m) {
                    legendreCoeffs.push_back(legendreCos(th, l, m));
                    dLegendreCoeffs.push_back(dLegendreCos(th, l, m));
                }

                for (int m = -l; m <= l; ++m) {
                    const int abs_m = abs(m);
                    fld_R -=
                        srcCoeffs[l][m+l] / r2lp2 * expI(m*ph) *
                        vec3cd(
                            -(l+1) * legendreCoeffs[abs_m],        // E_r
                            dLegendreCoeffs[abs_m],                // E_th
                            cmplx(0,m) * legendreCoeffs[abs_m]);   // E_ph * sin(th)
                }
                r2lp2 *= r;
            }

            // Convert to cartesian components
            fld_X = mat3d{
                {  sin(th)*cos(ph),  cos(th)*cos(ph), -sin(ph)/sin(th) },
                {  sin(th)*sin(ph),  cos(th)*sin(ph),  cos(ph)/sin(th) },
                {  cos(th),         -sin(th),          0.0             }
            } * fld_R;
        }

        flds.push_back(fld_X.real());
    }
    return flds;
}

template <typename T, typename Func>
std::vector<T> Leaf::getNearNborSols(Func kernel) {
    std::vector<T> sols;

    for (const auto& obs : particles) {
        T sol;
        if constexpr (std::is_same_v<T, double>)
            sol = 0.0;
        else
            sol = vec3d::Zero();

        auto obsPos = obs->getPos();

        /* due to other particles in this node */
        for (const auto& src : particles)
            if (src != obs)
                sol += src->getCharge() * kernel(obs->getPos() - src->getPos());

        /* due to particles in all neighbors */
        //for (const auto& node : nbors) {
        //    auto srcs = node->getParticles();
        //    for (const auto& src : srcs)
        //        sol += src->getCharge() * kernel(obsPos - src->getPos());
        //}

        /* due to particles in near neighbors (list 1) */
        for (const auto& node : nearNbors) {
            auto srcs = node->getParticles();
            for (const auto& src : srcs)
                sol += src->getCharge() * kernel(obsPos - src->getPos());
        }

        sols.push_back(sol);
    }
    return sols;
}

// pass calculated phi and fld to particles
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
}