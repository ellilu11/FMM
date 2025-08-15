#include "../node.h"

using namespace std;

const cmplx Node::getPhiFromExp(const vec3d& X, const std::vector<vecXcd>& expCoeffs, const int dirIdx) {
    vec3d dX = rotMatR[dirIdx] * (X - center);
    double dx = dX[0], dy = dX[1], dz = dX[2];
    assert(nodeLeng <= dz && dz <= 4.0*nodeLeng);
    assert(sqrt(dx*dx + dy*dy) <= 4.0*sqrt(2.0)*nodeLeng);
    // std::cout << dx / nodeLeng << ' ' << dy / nodeLeng << ' ' << dz / nodeLeng << '\n';

    cmplx phi(0,0);
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
        realVec legendreLMCoeffs;
        for (int m = 0; m <= l; ++m)
            legendreLMCoeffs.push_back(legendreLM(th, pair2i(l, m)));

        for (int m = -l; m <= l; ++m)
            phi += localCoeffs[l][m+l] * pow(r, l) *
            legendreLMCoeffs[std::abs(m)] * expI(static_cast<double>(m)*ph);
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

    auto revDir = [](int dirIdx) {
        return dirIdx - 2*(dirIdx%2) + 1;
        };

    cout << " Computing upward pass...\n";
    buildMpoleCoeffs();

    // Select a node at random for each trial
    int ntrials = 50;
    for (int trial = 0; trial < ntrials; ++trial) {

        const int maxLvl = 1;
        auto node = getRandNode(maxLvl);
        while (!node->getParticles().size())
            node = getRandNode(maxLvl);
        node->printMpoleCoeffs(coeffsFile);
        cout << trial << ", # Particles in this node: " << node->getParticles().size() << '\n';

        // Node is source node, interaction list are target nodes
        /*for (int dirIdx = 5; dirIdx < 6; ++dirIdx){
            auto expCoeffs = node->getMpoleToExpCoeffs(dirIdx);
            auto iList = (node->getDirList())[dirIdx];
            cout << " # INodes: " << iList.size() << '\n';

            for (const auto& iNode : iList) {
                // from exp coeffs
                //for (const auto& obs : iNode->getParticles()) {
                //    // outFile << node->getPhiFromMpole(obs->getPos()).real() << ' ';
                //    outFile << node->getPhiFromExp(obs->getPos(), expCoeffs, dirIdx).real() << ' ';
                //    obsFile << obs->getPos() << '\n';
                iNode->addShiftedExpCoeffs(expCoeffs, node->getCenter(), dirIdx);
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
        }*/

        // Node is target node, interaction list are source nodes


        for (int dirIdx = 0; dirIdx < 6; ++dirIdx) {
            auto iList = (node->getDirList())[revDir(dirIdx)];
            // cout << " # INodes: " << iList.size() << '\n';

            for (const auto& iNode : iList) {
                auto expCoeffs = iNode->getMpoleToExpCoeffs(dirIdx);
                node->addShiftedExpCoeffs(expCoeffs, iNode->getCenter(), dirIdx);
            }
        }
        node->buildLocalCoeffsFromDirList();

        for (const auto& obs : node->getParticles()) {
            outFile << node->getPhiFromLocal(obs->getPos()).real() << ' ';

            // analytic
            double phiAnl = 0;
            for (int dirIdx = 0; dirIdx < 6; ++dirIdx) {
                auto iList = (node->getDirList())[revDir(dirIdx)];
                for (const auto& iNode : iList)
                    for (const auto& src : iNode->getParticles())
                        phiAnl += src->getCharge() / (obs->getPos() - src->getPos()).norm();
            }
            outAnlFile << phiAnl << ' ';
        }
        // reset exp and local coeffs of this node
        node->resetNode();
    }
}