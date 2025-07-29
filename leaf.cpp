#include "leaf.h"

Leaf::Leaf(
    const ParticleVec& particles,
    const int branchIdx,
    Stem* const base)
    : Node(particles, branchIdx, base)
{
}

void Leaf::buildMpoleCoeffs() {
    for (int l = 0; l <= order; ++l) {
        cmplxVec coeffsL;
        for (int m = -l; m <= l; ++m) {
            cmplx a_lm;
            for (const auto& src : particles) {
                auto [r, th, ph] = src->getPos();
                a_lm += src->getCharge() * pow(r, l) * sphHarmonic(th, ph, l, -m);
            }
            coeffsL.push_back(a_lm);
        }
        coeffs.push_back(coeffsL);
    }
}

void Leaf::buildLocalCoeffs() {
    buildMpoleToLocalCoeffs();
    evaluateSolAtParticles();
}

vec3dVec Leaf::getFarPhis() {
    //vec3dVec phis(particles.size());
    //auto evaluateLocalExp = [this](std::shared_ptr<Particle> p) {
    //    return -evaluatePoly<vec3d>(localCoeffs, p->getPos()-center);
    //    };
    //std::transform(particles.begin(), particles.end(), phis.begin(), evaluateLocalExp);

    vec3dVec phis;
    //for (const auto& obs : particles)
    //    phis.push_back( -evaluatePoly<vec3d>(localCoeffs, obs->getPos()-center) );

    return phis;
}

vec3dVec Leaf::getFarFlds() {
    vec3dVec flds, dcoeffs;

    //for (size_t l = 1; l <= order; ++l)
    //    dcoeffs.push_back(static_cast<double>(l) * localCoeffs[l]);

    //for (const auto& obs : particles) {
    //    auto dphi = -evaluatePoly<vec3d>(dcoeffs, obs->getPos()-center);
    //    flds.push_back( vec3d(-dphi.real(), dphi.imag()) ); // care with minus signs
    //}

    return flds;
}

template <typename Func>
vec3dVec Leaf::getNearSols(Func kernel) {
    vec3dVec sols;

    for (const auto& obs : particles) {
        vec3d sol;
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

    /*
    auto phis = getFarPhis() + getNearSols(
        [](vec3d X) { return 1 / abs(X); }
    );
    auto flds = getFarFlds() + getNearSols(
        [](vec3d X) { return X / pow(abs(X),3); }
    );

    // for (auto [p, phi, fld] : std::views::zip(particles, phis, flds)) {
    for (size_t n = 0; n < particles.size(); ++n) {
        auto p = particles[n];
        p->setPhi(phis[n]);
        p->setFld(flds[n]);
    }*/
}

