#include <cassert>
#include <chrono>
#include <iostream>
#include <random>
#include "../math.h"
#include "../node.h"
#include "../leaf.h"

const cmplx Node::getFfield(const cmplx z) {
    cmplx phi = -coeffs[0] * std::log(z-zk);

    for (size_t k = 1; k < P_; ++k)
        phi -= coeffs[k] / std::pow(z-zk, k);

    return phi;
}

const cmplx Node::getAnalyticField(const cmplx z) {
    cmplx phi;
    for (size_t n = 0; n < psn.size(); ++n)
        phi -= qs[n] * std::log(z - psn[n]);
    return phi;
}

const cmplxVec Node::getAnalyticNfields() {
    cmplxVec phis;

    for (size_t obs = 0; obs < psn.size(); ++obs) {
        cmplx phi;
        for (size_t src = 0; src < psn.size(); ++src)
            if (src != obs) phi -= qs[src] * std::log(psn[obs] - psn[src]);
        phis.push_back(phi);
    }
    return phis;
}

/*
void Node::ffieldTest(const int Nobs) {
    const double R(2.0*L_);

    std::ofstream obsFile, outFile, outAnlFile;
    obsFile.open("out/obss.txt");
    outFile.open("out/ff.txt");
    outAnlFile.open("out/ffAnl.txt");

    cmplxVec obss;
    for (int n = 0; n < Nobs; ++n) {
        const double th = 2 * 3.1415927 * static_cast<double>(n) / static_cast<double>(Nobs);
        obss.emplace_back(cmplx(R * cos(th), R * sin(th)));
        obsFile << obss[n] << '\n';
    }

    const int P = Node::getP();
    for (int p = 1; p <= P; ++p) {
        Node::setP(p);
        buildMpoleCoeffs();

        for (const auto& obs : obss) {
            auto phi = evaluateFfield(obs);
            // auto phi = evaluateFfieldFromLeaf(obs);
            outFile << phi.real() << " ";
        }
        outFile << '\n';
        if (p < P) resetNode();
    }

    for (const auto& obs : obss) {
        auto phi = evalAnalyticField(obs);
        outAnlFile << phi.real() << " ";
    }
    outAnlFile << '\n';
}*/

void Node::nfieldTest() {
    using namespace std;

    ofstream outFile, outAnlFile;
    outFile.open("out/nf.txt");
    outAnlFile.open("out/nfAnl.txt");

    outFile << setprecision(15) << scientific;
    outAnlFile << setprecision(15) << scientific;

    const int P = Node::getP();
    for (int p = 1; p <= P; ++p) {
        Node::setP(p);
        Node::buildBinomTable();

        cout << " Computing upward pass...   (" << " P = " << p << " )\n";
        auto start = chrono::high_resolution_clock::now();

        buildMpoleCoeffs();

        auto end = chrono::high_resolution_clock::now();
        chrono::duration<double, milli> duration_ms = end - start;
        cout << "   Elapsed time: " << duration_ms.count() << " ms\n";

        cout << " Computing downward pass... (" << " P = " << P_ << " )\n";
        start = chrono::high_resolution_clock::now();

        buildLocalCoeffs();

        end = chrono::high_resolution_clock::now();
        duration_ms = end - start;
        cout << "   Elapsed time: " << duration_ms.count() << " ms\n";

        printPhi(outFile);
        outFile << '\n';

        if (p < P) {
            resetNode();
            binomTable.clear();
        }
    }

    cout << " Computing pairwise..." << endl;
    auto start = chrono::high_resolution_clock::now();

    auto phis = getAnalyticNfields();
    for (const auto& phi : phis)
        outAnlFile << phi.real() << '\n';

    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double, milli> duration_ms = end - start;
    cout << "   Elapsed time: " << duration_ms.count() << " ms\n";
}

