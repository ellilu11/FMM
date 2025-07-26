#include "leaf.h"

Leaf::Leaf(
    ParticleVec& particles,
    const cmplx center,
    const int branchIdx,
    Stem* const base)
    : Node(particles, center, branchIdx, base)
{
}

void Leaf::buildMpoleCoeffs() {
    for (int k = 0; k <= order; ++k) {
        cmplx a_k;
        for (const auto& src : particles)
            a_k += src->getCharge() * // care with operator precedence here
                    ( k == 0 ? 1.0 : -pow(src->getPos()-center, k) / static_cast<double>(k) );
        coeffs.push_back(a_k);
    }
}

void Leaf::buildLocalCoeffs() {
    buildMpoleToLocalCoeffs();
    evaluateSolAtParticles();
}

cmplxVec Leaf::getFarPhis() {
    //cmplxVec phis(particles.size());
    //auto evaluateLocalExp = [this](std::shared_ptr<Particle> p) {
    //    return -evaluatePoly<cmplx>(localCoeffs, p->getPos()-center);
    //    };
    //std::transform(particles.begin(), particles.end(), phis.begin(), evaluateLocalExp);

    cmplxVec phis;
    for (const auto& obs : particles)
        phis.push_back( -evaluatePoly<cmplx>(localCoeffs, obs->getPos()-center) );

    return phis;
}

cmplxVec Leaf::getFarFlds() {
    cmplxVec flds, dcoeffs;

    for (size_t l = 1; l <= order; ++l)
        dcoeffs.push_back(static_cast<double>(l) * localCoeffs[l]);

    for (const auto& obs : particles) {
        auto dphi = -evaluatePoly<cmplx>(dcoeffs, obs->getPos()-center);
        flds.push_back( cmplx(-dphi.real(), dphi.imag()) ); // care with the minus signs
    }

    return flds;
}

template <typename Func>
cmplxVec Leaf::getNearSols(Func kernel) {
    cmplxVec sols;

    for (const auto& obs : particles) {
        cmplx sol;
        auto obsPos = obs->getPos();

        // due to other particles in this node (apply reciprocity later)
        for (const auto& src : particles) {
            if (src != obs)
                sol += src->getCharge() * kernel(obsPos - src->getPos());
        }

        // due to particles in neighboring nodes (apply reciprocity much later)
        for (const auto& nbor : nbors) {
            auto srcsNbor = nbor->getParticles();
            for (const auto& src : srcsNbor)
                sol += src->getCharge() * kernel(obsPos - src->getPos());
        }
        sols.push_back(sol);
    }
    return sols;
}

void Leaf::evaluateSolAtParticles() {
    if (isRoot()) return;

    auto phis = getFarPhis() + getNearSols(
        [](cmplx z) { return -std::log(z); }
    );
    auto flds = getFarFlds() + getNearSols(
        [](cmplx z) { return z / std::norm(z); }
    );

    // for (auto [p, phi, fld] : std::views::zip(particles, phis, flds)) {
    for (size_t n = 0; n < particles.size(); ++n) {
        // auto p = particles[n];
        particles[n]->setPhi(phis[n]);
        particles[n]->setFld(flds[n]);
    }
}

