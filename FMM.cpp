#include <iostream>
#include <chrono>
#include <cmath>
#include <fstream>
#include <random>
#include "fmm.h"

using namespace std;

namespace Param {
    extern constexpr int    DIM     = 2;
    extern constexpr double L       = 10.0;
    extern constexpr double EPS     = 1.0E-1;
}

enum class Mode {
    READ_FROM_FILE,
    GEN_RNG
};

int main(int argc, char *argv[])
{
    // ==================== Make particles ==================== //
    int Nsrcs;
    Mode mode = Mode::READ_FROM_FILE;
    ParticleVec srcs;
    ofstream srcFile;

    switch (mode) {
        case Mode::READ_FROM_FILE :
             srcs = import_particles("config/srcs.txt");
             Nsrcs = srcs.size();
             break;

        case Mode::GEN_RNG : 
            Nsrcs = 5000;
            srcs = makeRNGParticles<uniform_real_distribution<double>>
                    (Nsrcs, -Param::L/2, Param::L/2);
            srcFile.open("config/srcs.txt");
            for (const auto& src : srcs) srcFile << *src;
            break;
        //case Mode::GEN_GAUSSIAN :
        //    srcs =
        //        makeRNGParticles<normal_distribution<double>, uniform_real_distribution<double>>
        //            (Nsrcs, 0.0, 0.1*Param::L);
        //    break;
        default : 
            throw std::runtime_error("Invalid mode");
    }

    // ==================== Partition domain ==================== //
    cout << " Partitioning domain...    (" << " N = " << Nsrcs << " )\n";
    auto start = chrono::high_resolution_clock::now();

    const int Nlvl = ceil(log(Nsrcs) / log(4.0));
    shared_ptr<Node> root;
    if (Nsrcs > 1)
        root = make_shared<Stem>(srcs, 0, Param::L, Nlvl, 0, nullptr);
    else
        root = make_shared<Leaf>(srcs, 0, Param::L, Nlvl, 0, nullptr);
    
    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double, milli> duration_ms = end - start;
    cout << "   Elapsed time: " << duration_ms.count() << " ms\n";

    std::ofstream nodeFile;
    nodeFile.open("out/nodes.txt");
    root->printNode(nodeFile);

    // ==================== Upward pass ==================== //
    const int P = Node::getP();
    cout << " Computing upward pass...   (" << " Order = " << P << " )\n";
    start = chrono::high_resolution_clock::now();

    Node::buildBinomTable();
    root->buildMpoleCoeffs();

    end = chrono::high_resolution_clock::now();
    duration_ms = end - start;
    cout << "   Elapsed time: " << duration_ms.count() << " ms\n";

    // constexpr int NOBS = 1000;
    // root->ffieldTest(NOBS);
    // std::ofstream mpoleCoeffFile;
    // mpoleCoeffFile.open("out/mpolecoeffs.txt");
    // root->printMpoleCoeffs(mpoleCoeffFile);

    // ==================== Downward pass ==================== //
    cout << " Computing downward pass... (" << " Order = " << P << " )\n";
    start = chrono::high_resolution_clock::now();

    root->buildLocalCoeffs();

    end = chrono::high_resolution_clock::now();
    duration_ms = end - start;
    cout << "   Elapsed time: " << duration_ms.count() << " ms\n";
    
    ofstream outFile;
    outFile.open("out/phi.txt");
    root->printPhi(outFile);

    // root->mpoleToLocalTest();
    // root->nfieldTest();
    // std::ofstream localCoeffFile;
    // localCoeffFile.open("out/localcoeffs.txt");
    // root->printLocalCoeffs(localCoeffFile);

    // ==================== Compute pairwise ==================== //
    cout << " Computing pairwise..." << endl;
    start = chrono::high_resolution_clock::now();

    auto phisAnl = root->getAnalyticNfields();

    end = chrono::high_resolution_clock::now();
    duration_ms = end - start;
    cout << "   Elapsed time: " << duration_ms.count() << " ms\n";

    ofstream outAnlFile;
    outAnlFile.open("out/phiAnl.txt");
    for (const auto& phi : phisAnl)
        outAnlFile << phi.real() << ' ';
    outAnlFile << '\n';

    return 0;
}