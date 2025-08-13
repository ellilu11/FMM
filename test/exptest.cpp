#include "../node.h"

using namespace std;

const cmplx Node::getPhiFromExp(const vec3d& X, const std::vector<vecXcd>& expCoeffs, const int dirIdx) {
    vec3d dX = X - center;
    // vec3d dX = rotMatR[dirIdx] * (X - center);
    double dx = dX[0], dy = dX[1], dz = dX[2];
    assert(nodeLeng <= dz && dz <= 4.0*nodeLeng);
    assert(sqrt(dx*dx + dy*dy) <= 4.0*sqrt(2.0)*nodeLeng);
    // std::cout << dx / nodeLeng << ' ' << dy / nodeLeng << ' ' << dz / nodeLeng << '\n';

    cmplx phi(0,0);
    for (int k = 0; k < orderExp; ++k) {
        auto M_k = tables.quadLengs_[k];
        auto [l_k, w_k] = tables.quadCoeffs_[k];
        for (int j = 0; j < M_k; ++j) {
            double a_kj = tables.alphas_[k][j];
            phi += expCoeffs[k][j]
                * exp(-l_k / nodeLeng * dz)
                * expI(l_k / nodeLeng * (dx*cos(a_kj) + dy*sin(a_kj)));
        }
    }

    return phi;
}

const cmplx Node::getPhiFromLocal(const vec3d& X) {
    auto dR = toSph(X);

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
    ofstream outFile, outAnlFile;
    outFile.open("out/local.txt");
    outAnlFile.open("out/localAnl.txt");

    outFile << setprecision(15) << scientific;
    outAnlFile << setprecision(15) << scientific;

    cout << " Computing upward pass...\n";
    buildMpoleCoeffs();

    // Select a node at random
    const int maxLvl = 1;
    auto node = getRandNode(maxLvl);
    while (!node->getParticles().size())
        node = getRandNode(maxLvl);
    cout << " This node has " << node->getParticles().size() << " sources\n";

    const int dirIdx = 0;
    auto iList = (node->getDirList())[dirIdx]; // get all nodes in dirlist of source node
    cout << " This ilist has " << iList.size() << " nodes\n";

    for (int prec = 0; prec < 3; ++prec) {
        Node::setExponentialOrder(prec);
        cout << " Exponential order: " << orderExp << '\n';

        auto expCoeffs = node->getMpoleToExpCoeffs(dirIdx);

        // using exponential coeffs
        for (const auto& iNode : iList) {
            for (const auto& obs : iNode->getParticles())
                // outFile << node->getPhiFromMpole(obs->getPos()).real() << ' ';
                outFile << node->getPhiFromExp(obs->getPos(), expCoeffs, dirIdx).real() << ' ';
        }

        // using local coeffs
        //for (const auto& iNode : iList) {
        //    iNode->buildShiftedExpCoeffs(expCoeffs, node->getCenter(), dirIdx);
        //    iNode->buildLocalCoeffsFromDirList();
        //    for (const auto& obs : iNode->getParticles())
        //        outFile << iNode->getPhiFromLocal(obs->getPos()-node->getCenter()).real() << ' ';
        //}

        outFile << '\n';
    }

    // analytic
    for (const auto& iNode : iList) {
        for (const auto& obs : iNode->getParticles()) {
            double phiAnl = 0;
            for (const auto& src : node->getParticles())
                // if (src != obs)
                phiAnl += src->getCharge() / (obs->getPos() - src->getPos()).norm();
            outAnlFile << phiAnl << ' ';
        }
    }
    outAnlFile << '\n';
}