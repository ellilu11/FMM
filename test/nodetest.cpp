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

const double Node::getDirectPhi(const vec3d& R) {
    double phi = 0;
    auto X = toCart(R);
    for (const auto& p : particles)
        phi += p->getCharge() / (X - p->getPos()).norm();
    return phi;
}

const cmplx Node::getDirectPhiFromMpole(const vec3d& R) {
    auto r = R[0], th = R[1], ph = R[2];
    cmplx phi(0, 0);

    for (int l = 0; l <= order; ++l) {
        realVec legendreLMCoeffs;
        for (int m = 0; m <= l; ++m)
            legendreLMCoeffs.push_back(legendreLM(th, l, m));

        for (int m = -l; m <= l; ++m) {
            int m_ = m+l;
            phi += coeffs[l][m_] / pow(r, l+1) *
                legendreLMCoeffs[std::abs(m)] * expI(static_cast<double>(-m)*ph);
        }
    }

    return phi;
}

void Node::ffieldTest(const pair<int,int>& Nangles) {
    auto [Nth, Nph] = Nangles;
    const double r(2.0*rootLeng);

    std::ofstream obsFile, outFile, outAnlFile;
    obsFile.open("config/obss.txt");
    outFile.open("out/ff.txt");
    outAnlFile.open("out/ffAnl.txt");

    vec3dVec obss; // observer positions in spherical coordinates
    for (int ith = 0; ith < Nth; ++ith) {
        for (int iph = 0; iph < Nph; ++iph) {
            const double th = PI * static_cast<double>(ith) / static_cast<double>(Nth);
            const double ph = TAU * static_cast<double>(iph) / static_cast<double>(Nph);
            vec3d obs(r, th, ph);
            // vec3d obs(r*sin(th)*cos(ph), r*sin(th)*sin(ph), r*cos(th));
            obss.push_back(obs);
            obsFile << obs << '\n';
        }
    }

    const int order = Node::getExpansionOrder();
    for (int p = order; p <= order; ++p) {
        Node::setExpansionOrder(p);

        cout << " Computing upward pass...   (" << "Expansion order: " << p << ")\n";
        auto start = chrono::high_resolution_clock::now();
        
        buildMpoleCoeffs();

        auto end = chrono::high_resolution_clock::now();
        chrono::duration<double, milli> duration_ms = end - start;
        cout << "   Elapsed time: " << duration_ms.count() << " ms\n";

        for (const auto& obs : obss) {
            // auto phi = evaluateFfield(obs);
            auto phi = getFfieldFromLeaf(obs);
            outFile << phi.real() << "\n";
        }
        outFile << '\n';
        if (p < order) resetNode();
    }

    for (const auto& obs : obss) 
        outAnlFile << getDirectPhi(obs) << "\n";
    outAnlFile << '\n';

}