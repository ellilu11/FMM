#include "leaf.h"

Leaf::Leaf(
    const ParticleVec& particles,
    const int branchIdx,
    Stem* const base)
    : Node(particles, branchIdx, base)
{
    // std::cout << "Leaf has " << particles.size() << " particles\n";
}

void Leaf::buildLists() {
    if (!isRoot()) {
        buildNearNeighbors();
        buildInteractionList();
        buildDirectedIList();
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
            realVec legendreLMCoeffs;
            for (int m = 0; m <= l; ++m)
                legendreLMCoeffs.push_back(legendreLM(th, pair2i(l,m)));

            for (int m = -l; m <= l; ++m) {
                coeffs[l][m+l] +=
                    src->getCharge() * r2l *
                    legendreLMCoeffs[std::abs(-m)] * expI(static_cast<double>(-m)*ph)
                    // * (m < 0 ? pm(m) : 1.0);
                    ; 
            }
            r2l *= r;
        }
    }
}

void Leaf::propagateExpCoeffs() {
    if (isRoot()) return;

    for (int dir = 0; dir < 6; ++dir){
        auto expCoeffs = getMpoleToExpCoeffs(dir);
        auto iList = dirList[dir];
        for (const auto& iNode : iList)
            iNode->buildShiftedExpCoeffs(expCoeffs, center, dir);
    }
}

void Leaf::buildLocalCoeffs() {
    if (!isRoot()) {
        buildLocalCoeffsFromLeafIlist();
        buildLocalCoeffsFromDirList();
    }

    evaluateSolAtParticles();
}

realVec Leaf::getFarPhis() {
    realVec phis;
    for (const auto& obs : particles) {
        auto dR = toSph(obs->getPos() - center);

        double r = dR[0], th = dR[1], ph = dR[2];
        cmplx phi(0, 0);

        for (int l = 0; l <= order; ++l) {
            realVec legendreLMCoeffs;
            for (int m = 0; m <= l; ++m)
                legendreLMCoeffs.push_back(legendreLM(th, pair2i(l,m)));

            for (int m = -l; m <= l; ++m)
                phi += localCoeffs[l][m+l] * pow(r, l) *
                    legendreLMCoeffs[std::abs(m)] * expI(static_cast<double>(m)*ph);
        }
        phis.push_back(phi.real());
    }

    return phis;
}

std::vector<vec3d> Leaf::getFarFlds() {
    std::vector<vec3d> flds, dcoeffs;

    //for (size_t l = 1; l <= order; ++l)
    //    dcoeffs.push_back(static_cast<double>(l) * localCoeffs[l]);

    //for (const auto& obs : particles) {
    //    auto dphi = -evaluatePoly<vec3d>(dcoeffs, obs->getPos()-center);
    //    flds.push_back( vec3d(-dphi.real(), dphi.imag()) ); // care with minus signs
    //}

    return flds;
}

template <typename Func>
realVec Leaf::getNearPhis(Func kernel) {
    realVec sols;

    for (const auto& obs : particles) {
        double sol = 0;
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

//template <typename T, typename Func>
//std::vector<T> Leaf::getNearSols(Func kernel) {
//    std::vector<T> sols;
//
//    for (const auto& obs : particles) {
//        T sol{};
//        auto obsPos = obs->getPos();
//
//        // due to other particles in this node (implement reciprocity later)
//        for (const auto& src : particles)
//            if (src != obs)
//                sol += src->getCharge() * kernel(obsPos - src->getPos());
//
//        // due to particles in neighboring nodes (implement reciprocity much later)
//        for (const auto& nbor : nbors) {
//            auto srcsNbor = nbor->getParticles();
//            for (const auto& src : srcsNbor)
//                sol += src->getCharge() * kernel(obsPos - src->getPos());
//        }
//        sols.push_back(sol);
//    }
//    return sols;
//}

// pass calculated phi and fld to particles
void Leaf::evaluateSolAtParticles() {
    if (isRoot()) return;

    auto phis = getFarPhis() 
        // + getNearPhis( [](vec3d X) { return 1.0 / X.norm(); } )
        ;
    //auto flds = getFarFlds() + getNearSols<vec3d>(
    //    [](vec3d X) { return X / pow(abs(X),3); }
    //);

    // for (auto [p, phi, fld] : std::views::zip(particles, phis, flds)) {
    for (size_t n = 0; n < particles.size(); ++n) {
        auto p = particles[n];
        p->setPhi(phis[n]);
        // p->setFld(flds[n]);
    }
}