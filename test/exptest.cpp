#include "../node.h"

using namespace std;

const cmplx Node::getPhiFromExp(const vec3d& X, const std::vector<vecXcd>& expCoeffs, const int dir) {
    const vec3d dX = rotMatR[dir] * (X - center);
    const double dx = dX[0], dy = dX[1], dz = dX[2];

    const auto idX = round(Eigen::Array3d(dX)/nodeLeng);
    const int idx = idX[0], idy = idX[1], idz = idX[2];
    // const size_t l = idx + 7*idy + 49*idz - 74; // = (idx+3) + (idy+3)*7 + (idz-2)*49;
    // std::cout << dir << ' ' << idx << ' ' << idy << ' ' << idz << '\n';

    //assert(0 <= l && l < 98);
    //assert(1 <= idz && idz <= 4);
    //assert(sqrt(idx*idx + idy*idy) <= 4.0*sqrt(2.0));
    // std::cout << dx / nodeLeng << ' ' << dy / nodeLeng << ' ' << dz / nodeLeng << '\n';

    cmplx phi(0,0);

    //for (size_t k = 0; k < orderExp; ++k)
    //    for (size_t j = 0; j < tables.quadLengs_[k]; ++j)
    //        phi += expCoeffs[k][j] * tables.exps_[k][j][l];

    for (int k = 0; k < orderExp; ++k) {
        auto M_k = tables.quadLengs_[k];
        auto [lmd_k, w_k] = tables.quadCoeffs_[k];
        double l_k = lmd_k / nodeLeng;
        for (int j = 0; j < M_k; ++j) {
            double a_kj = tables.alphas_[k][j];
            assert(a_kj == 2.0*PI*(j+1)/static_cast<double>(M_k));
            phi += expCoeffs[k][j]
                * exp(-l_k * dz + iu*l_k * (dx*cos(a_kj) + dy*sin(a_kj)));
        }
    }
    return phi;
}

const cmplx Node::getPhiFromLocal(const vec3d& X) {
    auto dR = toSph(X - center);

    double r = dR[0], th = dR[1], ph = dR[2];
    cmplx phi(0, 0);

    for (int l = 0; l <= order; ++l) {
        realVec legendreCoeffs;
        for (int m = 0; m <= l; ++m)
            legendreCoeffs.push_back(legendreCos(th, l, m));

        for (int m = -l; m <= l; ++m)
            phi += localCoeffs[l][m+l] * pow(r, l) *
            legendreCoeffs[std::abs(m)] * expI(static_cast<double>(m)*ph);
    }
    return phi;
}

void Node::mpoleToExpToLocalTest() {
    ofstream outFile, outAnlFile, coeffsFile, obsFile;
    outFile.open("out/local.txt");
    outAnlFile.open("out/localAnl.txt");
    coeffsFile.open("out/mpolecoeffs.txt");
    obsFile.open("config/obss.txt");

    outFile << setprecision(15) << scientific;
    outAnlFile << setprecision(15) << scientific;

    auto revDir = [](int dir) {
        return dir - 2*(dir%2) + 1;
        };

    cout << " Computing upward pass...\n";
    buildMpoleCoeffs();

    // Select a node at random for each trial
    constexpr int ntrials = 1;
    for (int trial = 0; trial < ntrials; ++trial) {
        const int maxLvl = 1;
        auto node = getRandNode(maxLvl);
        while (!node->particles.size())
            node = getRandNode(maxLvl);
        node->printMpoleCoeffs(coeffsFile);
        cout << " Trial " << trial << ", # Particles in this node: " << node->getParticles().size() << '\n';

        // Node is source node, interaction list are target nodes
        for (int dir = 0; dir < 6; ++dir){
            auto expCoeffs = node->getMpoleToExpCoeffs(dir);
            auto iList = (node->getDirList())[dir];
            cout << " # INodes: " << iList.size() << '\n';

            // from exp coeffs
            /*for (const auto& iNode : iList) {
                for (const auto& obs : iNode->getParticles()) {
                    // outFile << node->getPhiFromMpole(obs->getPos()).real() << ' ';
                    outFile << node->getPhiFromExp(obs->getPos(), node->getExpCoeffs(dir), dir).real() << ' ';

                    // analytic
                    double phiAnl = 0;
                    for (const auto& src : node->particles)
                        phiAnl += src->getCharge() / (obs->getPos() - src->getPos()).norm();
                    outAnlFile << phiAnl << ' ';
                }
            }*/

            // from local coeffs
            for (const auto& iNode : iList) {
                iNode->addShiftedExpCoeffs(expCoeffs, node->getCenter(), dir);
                iNode->buildLocalCoeffsFromDirList();
                for (const auto& obs : iNode->getParticles()) {
                    outFile << iNode->getPhiFromLocal(obs->getPos()).real() << ' ';

                    // analytic
                    double phiAnl = 0;
                    for (const auto& src : node->getParticles())
                        phiAnl += src->getCharge() / (obs->getPos() - src->getPos()).norm();
                    outAnlFile << phiAnl << ' ';
                }
            }
        }

        // Node is target node, interaction list are source nodes
        /*for (int dir = 0; dir < 1; ++dir) {
            auto iList = (node->getDirList())[revDir(dir)];
            // cout << " # INodes: " << iList.size() << '\n';

            for (const auto& iNode : iList) {
                auto expCoeffs = iNode->getMpoleToExpCoeffs(dir);
                node->addShiftedExpCoeffs(expCoeffs, iNode->getCenter(), dir);
            }
        }
        node->buildLocalCoeffsFromDirList();

        for (const auto& obs : node->getParticles()) {
            outFile << node->getPhiFromLocal(obs->getPos()).real() << ' ';

            // analytic
            double phiAnl = 0;
            for (int dir = 0; dir < 1; ++dir) {
                auto iList = (node->getDirList())[revDir(dir)];
                for (const auto& iNode : iList)
                    for (const auto& src : iNode->getParticles())
                        phiAnl += src->getCharge() / (obs->getPos() - src->getPos()).norm();
            }
            outAnlFile << phiAnl << ' ';
        }
        // reset exp and local coeffs of this node
        node->resetNode();*/
    }
}