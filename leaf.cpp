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
    if (!isRoot()) {
        buildNbors();
        buildInteractionList();
        pushSelfToNearNonNbors();
    }
}

void Leaf::buildMpoleCoeffs() {
    for (int k = 0; k <= order; ++k) {
        cmplx a_k;
        for (const auto& src : particles)
            a_k += src->getCharge() * // care with operator precedence
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
        flds.push_back( cmplx(-dphi.real(), dphi.imag()) ); // care with minus signs
    }

    return flds;
}

cmplxVec Leaf::getNearNonNborPhis() {
    cmplxVec phis;

    for (const auto& obs : particles) {
        cmplx phi(0,0);
        auto obsPos = obs->getPos();

        for (const auto& node : nearNonNbors) {
            auto dz = obsPos - node->getCenter();
            auto srcCoeffs = node->getMpoleCoeffs();

            phi -= srcCoeffs[0] * log(dz);

            auto dz2k = dz;
            for (int k = 1; k <= order; ++k) {
                phi -= srcCoeffs[k] / dz2k;
                dz2k *= dz;
            }
        }
        phis.push_back(phi);
    }

    return phis;
}

cmplxVec Leaf::getNearNonNborFlds() {
    cmplxVec flds;

    for (const auto& obs : particles) {
        cmplx fld(0, 0);
        auto obsPos = obs->getPos();

        for (const auto& node : nearNonNbors) {
            auto dz = obsPos - node->getCenter();
            auto srcCoeffs = node->getMpoleCoeffs();

            auto b_0 = srcCoeffs[0] / dz;
            fld += cmplx( b_0.real(), -b_0.imag() ) ;

            auto dz2kpp = dz*dz;
            for (int k = 1; k <= order; ++k) {
                auto b_k = static_cast<double>(k) * srcCoeffs[k] / dz2kpp;
                fld += cmplx(-b_k.real(), b_k.imag());
                dz2kpp *= dz;
            }
        }
        flds.push_back(fld);
    }

    return flds;
}

template <typename Func>
cmplxVec Leaf::getNearNborSols(Func kernel) {
    cmplxVec sols;

    for (const auto& obs : particles) {
        cmplx sol(0,0);
        auto obsPos = obs->getPos();

        // due to other particles in this node (implement reciprocity later)
        for (const auto& src : particles)
            if (src != obs)
                sol += src->getCharge() * kernel(obsPos - src->getPos());

        // due to particles in neighboring nodes (implement reciprocity much later)
        // for (const auto& nbor : nbors) {
        for (const auto& nbor : nearNbors) {
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

    auto phis = 
        getFarPhis() + 
        getNearNonNborPhis() + 
        getNearNborSols(
        [](cmplx z) { return -std::log(z); }
        );

    auto flds = 
        getFarFlds() + 
        getNearNonNborFlds() +
        getNearNborSols(
        [](cmplx z) { return z / std::norm(z); }
        );

    // for (auto [p, phi, fld] : std::views::zip(particles, phis, flds)) {
    for (size_t n = 0; n < particles.size(); ++n) {
        auto p = particles[n];
        p->setPhi(phis[n]);
        p->setFld(flds[n]);
    }
}

