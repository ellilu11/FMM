#include <chrono>
#include <fstream>
#include <iostream>
#include "fmm.h"

using namespace std;

namespace Param {
    extern constexpr int    DIM     = 2;
    extern constexpr double L       = 10.0;
    extern constexpr double EPS     = 1.0E-1;
    extern constexpr int maxPartsPerNode = 5;
}

int main(int argc, char *argv[])
{
    // ==================== Make particles ==================== //
    enum class Mode {
        READ,
        GEN
    };

    constexpr Mode mode = Mode::READ;
    ParticleVec srcs;
    int Nsrcs;
    ofstream srcFile;

    switch (mode) {
        case Mode::READ :
             srcs = importParticles("config/uniform_minusQ.txt");
             Nsrcs = srcs.size();
             break;

        case Mode::GEN : 
            Nsrcs = 5000;
            srcs = makeRNGParticles<uniform_real_distribution<double>>
                    (Nsrcs, -Param::L/2, Param::L/2);

            srcFile.open("config/uniform.txt");
            for (const auto& src : srcs) srcFile << *src;
            break;

        default : 
            throw std::runtime_error("Invalid mode");
    }

    // ==================== Partition domain ==================== //
    //const int Nlvl = ceil(log(Nsrcs) / log(4.0));
    //Node::setMaxLvl(Nlvl);
    cout << " Partitioning domain...     (" << " Nsrcs = " << Nsrcs << " )\n";
    auto start = chrono::high_resolution_clock::now();

    shared_ptr<Node> root;
    if (Nsrcs > 1)
        root = make_shared<Stem>(srcs, 0, nullptr);
    else
        root = make_shared<Leaf>(srcs, 0, nullptr);
    
    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double, milli> duration_ms = end - start;
    cout << "   Elapsed time: " << duration_ms.count() << " ms\n";
    std::ofstream nodeFile("out/nodes.txt");
    root->printNode(nodeFile);

    // ==================== Upward pass ==================== //
    const int order = Node::getExpansionOrder();
    cout << " Computing upward pass...   (" << " Order = " << order << " )\n";
    start = chrono::high_resolution_clock::now();

    Node::buildBinomTable();
    root->buildMpoleCoeffs();

    end = chrono::high_resolution_clock::now();
    duration_ms = end - start;
    cout << "   Elapsed time: " << duration_ms.count() << " ms\n";

    // ==================== Downward pass ==================== //
    cout << " Computing downward pass... (" << " Order = " << order << " )\n";
    start = chrono::high_resolution_clock::now();

    root->buildLocalCoeffs();

    end = chrono::high_resolution_clock::now();
    duration_ms = end - start;
    cout << "   Elapsed time: " << duration_ms.count() << " ms\n";

    printSols(srcs, "out/phi.txt", "out/fld.txt");

    // ==================== Compute pairwise ==================== //
    cout << " Computing pairwise..." << endl;
    start = chrono::high_resolution_clock::now();

    auto phisAnl = root->getDirectPhis();

    end = chrono::high_resolution_clock::now();
    duration_ms = end - start;
    cout << "   Elapsed time: " << duration_ms.count() << " ms\n";

    ofstream phiAnlFile("out/phiAnl.txt");
    for (const auto& phi : phisAnl)
        phiAnlFile << phi.real() << '\n';

    auto fldsAnl = root->getDirectFlds();

    ofstream fldAnlFile("out/fldAnl.txt");
    for (const auto& fld : fldsAnl)
        fldAnlFile << fld << '\n';

    return 0;

}