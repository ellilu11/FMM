#include "leaf.h"

Leaf::Leaf(
    ParticleVec& particles,
    const cmplx center,
    const int lvl,
    const int branchIdx,
    Stem* const base)
    : Node(particles, center, lvl, branchIdx, base)
{
}

void Leaf::buildMpoleCoeffs() {
    for (int k = 0; k <= order; ++k) {
        cmplx a_k;
        for (const auto& src : particles)
            a_k += src->getCharge() *
            k == 0 ? 1 : -pow(src->getPos()-center, k) / static_cast<double>(k);
        coeffs.push_back(a_k);
    }
}

void Leaf::buildLocalCoeffs() {
    buildMpoleToLocalCoeffs();
    evaluateSolAtParticles();
}

cmplxVec Leaf::getPhisFar() {
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

cmplxVec Leaf::getFldsFar() {
    cmplxVec flds, dcoeffs;

    for (size_t l = 1; l <= order; ++l)
        dcoeffs.push_back(static_cast<double>(l) * localCoeffs[l]);

    for (const auto& obs : particles) {
        auto dphi = -evaluatePoly<cmplx>(dcoeffs, obs->getPos()-center);
        flds.push_back( cmplx(-dphi.real(), dphi.imag()) ); // care with the minus signs
    }

    return flds;
}

cmplxVec Leaf::getPhisNear() {
    cmplxVec phis;

    for (const auto& obs : particles) {
        cmplx phi;
        auto obsPos = obs->getPos();

        // phi due to other particles in this node (apply reciprocity later)
        for (const auto& src : particles)
            if (src != obs) phi -= src->getCharge() * std::log(obsPos - src->getPos());

        // phi due to particles in neighboring nodes
        for (const auto& nbor : nbors) {
            auto srcsNbor = nbor->getParticles();
            for (const auto& src : srcsNbor)
                phi -= src->getCharge() * std::log(obsPos - src->getPos());
        }
        phis.push_back(phi);
    }
    return phis;
}

void Leaf::evaluateSolAtParticles() {
    auto phis = getPhisFar() + getPhisNear();
    auto flds = getFldsFar(); // + getFldsNear()

    // any way to avoid the indexing?
    for (size_t n = 0; n < particles.size(); ++n) {
        particles[n]->setPhi(phis[n]);
        particles[n]->setFld(flds[n]);
    }
}

