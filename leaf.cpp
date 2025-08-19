#include "leaf.h"

Leaf::Leaf(
    const ParticleVec& particles,
    const int branchIdx,
    Stem* const base)
    : Node(particles, branchIdx, base)
{
}

/*
void Leaf::getFarNeighbors() {
    assert(!isRoot());

    const double minDist = nodeLeng;
    for (const auto& nbor : nbors) {
        if (nbor->isNodeType<Leaf>()) continue;

        auto nodes = nbor->getBranches();
        for (const auto& node : nodes)
            
    }
}*/

void Leaf::buildLists() {
    if (!isRoot()) {
        buildNearNeighbors();
        buildInteractionList();
    }
}

void Leaf::buildMpoleCoeffs() {
    for (int l = 0; l <= order; ++l)
        coeffs.emplace_back(vecXcd::Zero(2*l+1));

    for (const auto& src : particles) {
        auto dR = toSph(src->getPos() - center);
        double r = dR[0], th = dR[1], ph = dR[2];
        double r2l = 1.0;

        for (int l = 0; l <= order; ++l) {
            realVec legendreCosCoeffs;
            for (int m = 0; m <= l; ++m)
                legendreCosCoeffs.push_back(legendreCos(th,l,m));

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
        auto start = std::chrono::high_resolution_clock::now();
        auto expCoeffs = getMpoleToExpCoeffs(dir);
        t_M2X += std::chrono::high_resolution_clock::now() - start;

        auto iList = dirList[dir];

        start = std::chrono::high_resolution_clock::now();
        for (const auto& iNode : iList)
            iNode->addShiftedExpCoeffs(expCoeffs, center, dir);
        t_X2X += std::chrono::high_resolution_clock::now() - start;
    }
}

void Leaf::buildLocalCoeffs() {
    auto start = std::chrono::high_resolution_clock::now();
    if (!isRoot()) {
        buildLocalCoeffsFromLeafIlist();
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

    start = std::chrono::high_resolution_clock::now();
    evaluateSolAtParticles();
    t_direct += std::chrono::high_resolution_clock::now() - start;
}

realVec Leaf::getFarPhis() {
    realVec phis;
    for (const auto& obs : particles) {
        auto dR = toSph(obs->getPos() - center);
        double r = dR[0], th = dR[1], ph = dR[2];

        // using Horner's scheme
        /*cmplxVec polyCoeffs;
        for (int l = 0; l <= order; ++l) {
            realVec legendreCosCoeffs;
            cmplx polyCoeffs_l;
            for (int m = 0; m <= l; ++m)
                legendreCosCoeffs.push_back(legendreCos(th, pair2i(l,m)));

            for (int m = -l; m <= l; ++m)
                polyCoeffs_l += localCoeffs[l][m+l] *
                    legendreCosCoeffs[std::abs(m)] * expI(static_cast<double>(m)*ph);
            polyCoeffs.push_back(polyCoeffs_l);
        }
        phis.push_back(evaluatePoly<cmplx>(polyCoeffs,r).real());*/

        cmplx phi(0, 0);
        double r2l = 1.0;
        for (int l = 0; l <= order; ++l) {
            realVec legendreCosCoeffs;
            for (int m = 0; m <= l; ++m)
                legendreCosCoeffs.push_back(legendreCos(th,l,m));

            for (int m = -l; m <= l; ++m)
                phi += localCoeffs[l][m+l] * r2l *
                        legendreCoeffs[std::abs(m)] * expI(m*ph);
            r2l *= r;
        }
        phis.push_back(phi.real());
    }

    return phis;
}

std::vector<vec3d> Leaf::getFarFlds() {
    auto dLegendreCos = [](const double th, const int l, const int abs_m) {
        if (!l) return 0.0;
        assert(0 <= abs_m && abs_m <= l);
        return
            l / tan(th) * legendreCos(th, l, abs_m) -
            (abs_m <= (l-1) ?
                (l + abs_m) / sin(th) * legendreCos(th, l-1, abs_m)
                * tables.fracCoeffYlm_[l][abs_m] : // sqrt((l-abs_m)/static_cast<double>(l+abs_m))
                0.0);
        };

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
                    localCoeffs[l][m+l] * r2lmm * expI(static_cast<double>(m)*ph) *
                    vec3cd(
                        l * legendreCoeffs[abs_m],             // E_r
                        dLegendreCoeffs[abs_m],                // E_th
                        cmplx(0, m) * legendreCoeffs[abs_m]); // E_ph * sin(th)
            }
            r2lmm *= r;
        }

        // Convert to cartesian components
        vec3cd fld_X = mat3d{
            {  sin(th)*cos(ph),  cos(th)*cos(ph), -sin(ph)/sin(th) },
            {  sin(th)*sin(ph),  cos(th)*sin(ph),  cos(ph)/sin(th) },
            {  cos(th),         -sin(th),          0.0             }
        } *fld_R;

        flds.push_back(fld_X.real());
    }
    return flds;
}

template <typename T, typename Func>
std::vector<T> Leaf::getNearSols(Func kernel) {
    std::vector<T> sols;

    for (const auto& obs : particles) {
        T sol;
        if constexpr (std::is_same_v<T, double>)
            sol = 0.0;
        else
            sol = vec3d::Zero();

        auto obsPos = obs->getPos();

        // due to other particles in this node (implement reciprocity later)
        for (const auto& src : particles)
            if (src != obs)
                sol += src->getCharge() * kernel(obsPos - src->getPos());

        // due to particles in neighboring nodes (implement reciprocity much later)
        for (const auto& nbor : nbors) {
            auto srcsNbor = nbor->getParticles();
            for (const auto& src : srcsNbor)
                sol += src->getCharge() * kernel(obsPos - src->getPos());
        }
        sols.push_back(sol);
    }
    return sols;
}

// pass calculated phi and fld to particles
void Leaf::evaluateSolAtParticles() {
    if (isRoot()) return;

    auto phis = getFarPhis() + getNearSols<double>( 
        [](vec3d X) { return 1.0 / X.norm(); } 
    );

    auto flds = getFarFlds() + getNearSols<vec3d>(
        [](vec3d X) { return X / pow(X.norm(), 3); }
    );

    // for (auto [p, phi, fld] : std::views::zip(particles, phis, flds)) {
    for (size_t n = 0; n < particles.size(); ++n) {
        auto p = particles[n];
        p->addToPhi(phis[n]);
        p->addToFld(flds[n]);
    }
}