#include "../node.h"

using namespace std;

void Node::setRandNodeStats() {
    auto node = getRandNode(0);
    node->setNodeStat(3);

    node->buildNearNeighbors();
    for (const auto& nbor: node->getNearNeighbors())
        nbor->setNodeStat(1);

    node->buildInteractionList();
    for (const auto& iNode : node->getInteractionList())
        iNode->setNodeStat(2);
}

const double Node::getDirectPhi(const vec3d& X) {
    double phi = 0;
    for (const auto& p : particles)
        phi += p->getCharge() / (X - p->getPos()).norm();
    return phi;
}

const cmplx Node::getPhiFromMpole(const vec3d& X) {
    auto dR = toSph(X - center);
    double r = dR[0], th = dR[1], ph = dR[2];
    cmplx phi(0, 0);

    for (int l = 0; l <= order; ++l) {
        realVec legendreLMCoeffs;
        for (int m = 0; m <= l; ++m)
            legendreLMCoeffs.push_back(legendreLM(th, l, m));

        for (int m = -l; m <= l; ++m) {
            int m_ = m + l;
            phi += coeffs[l][m_] / pow(r, l+1) *
                legendreLMCoeffs[std::abs(m)] * expI(static_cast<double>(m)*ph);
        }
    }

    return phi;
}

void Node::ffieldTest(const int Nr, const int Nth, const int Nph) {
    std::ofstream obsFile, outFile, outAnlFile, coeffsFile;
    obsFile.open("config/obss.txt");
    outFile.open("out/ff.txt");
    outAnlFile.open("out/ffAnl.txt");
    coeffsFile.open("out/mpolecoeffs.txt");

    //outFile << setprecision(9);
    //outAnlFile << setprecision(9);

    vec3dVec obss;
    for (int ir = 0; ir < Nr; ++ir){
        double r = 5.0*(ir+1.0)*rootLeng;
        for (int ith = 0; ith < Nth; ++ith) {
            double th = PI * ith / static_cast<double>(Nth);
            for (int iph = 0; iph < Nph; ++iph) {
                double ph = TAU * iph / static_cast<double>(Nph);
                auto obs = toCart(vec3d(r, th, ph));
                obss.push_back(obs);
                obsFile << obs << '\n';
            }
        }
    }

    const int order = Node::getExpansionOrder();
    for (int p = 1; p <= order; ++p) {
        Node::setExpansionOrder(p);

        cout << " Computing upward pass...   (" << "Expansion order: " << p << ")\n";
        auto start = chrono::high_resolution_clock::now();
        
        buildMpoleCoeffs();

        auto end = chrono::high_resolution_clock::now();
        chrono::duration<double, milli> duration_ms = end - start;
        cout << "   Elapsed time: " << duration_ms.count() << " ms\n";

        for (const auto& obs : obss) {
            // auto phi = getPhiFromMpole(obs);
            auto phi = getPhiFromBranchMpole(obs,0);
            outFile << phi.real() << " ";
        }
        outFile << '\n';
        if (p < order) resetNode();
    }

    printMpoleCoeffs(coeffsFile);

    for (const auto& obs : obss) 
        outAnlFile << getDirectPhi(obs) << " ";
    outAnlFile << '\n';

}