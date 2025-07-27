#include <chrono>
#include <fstream>
#include <iostream>
#include "config.h"
#include "fmm.h"

using namespace std;

extern constexpr int DIM = 2;

int main(int argc, char *argv[])
{
    Config config("config/config.txt");

    // ==================== Make particles ==================== //
    const auto fname = makeFname(config);
    ParticleVec srcs;
    int Nsrcs;
    
    switch (config.mode) {
        case Mode::READ :
            srcs = importParticles(fname);
            Nsrcs = srcs.size();
            break;

        case Mode::GEN: {
            Nsrcs = config.nsrcs;
            ofstream srcFile(fname);
            if (config.dist == Dist::UNIFORM)
                srcs = makeParticles<uniform_real_distribution<double>>
                        (Nsrcs, -config.L/2, config.L/2, config.dist, config.cdist);
            else
                srcs = makeParticles<uniform_real_distribution<double>>
                        (Nsrcs, 0, 1, config.dist, config.cdist);

            for (const auto& src : srcs) srcFile << *src;
            break;
        }
        default : 
            throw std::runtime_error("Invalid mode");
    }

    cout << " Mode:         " << (config.mode == Mode::READ ? "READ" : "GEN") << '\n';
    cout << " Src file:     " << fname << '\n';
    cout << " Nsrcs:        " << Nsrcs << '\n';
    cout << " Root length:  " << config.L << '\n';
    cout << " EPS:          " << config.EPS << '\n';
    cout << " maxNodeParts: " << config.maxNodeParts << '\n' << '\n';

    // ==================== Partition domain ==================== //
    cout << " Setting up domain...\n";
    auto start = chrono::high_resolution_clock::now();

    Node::setNodeParams(config);
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
    cout << " Computing upward pass...   (" << " Expansion order: " << order << " )\n";
    start = chrono::high_resolution_clock::now();

    Node::buildBinomTable();
    root->buildMpoleCoeffs();

    end = chrono::high_resolution_clock::now();
    duration_ms = end - start;
    cout << "   Elapsed time: " << duration_ms.count() << " ms\n";

    // ==================== Downward pass ==================== //
    cout << " Computing downward pass...\n";
    start = chrono::high_resolution_clock::now();

    root->buildLocalCoeffs();

    end = chrono::high_resolution_clock::now();
    duration_ms = end - start;
    cout << "   Elapsed time: " << duration_ms.count() << " ms\n";

    printSols(srcs, "out/phi.txt", "out/fld.txt");

    // ==================== Compute direct ==================== //
    if (!config.evalDirect) return 0;
    cout << " Computing direct..." << endl;
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