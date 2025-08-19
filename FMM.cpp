#include <chrono>
#include <fstream>
#include <iostream>
#include "config.h"
#include "fmm.h"

using namespace std;

extern constexpr int DIM = 2;

int main(int argc, char *argv[])
{
    Config config("config/config2D.txt");

    // ==================== Make particles ==================== //
    const auto fname = makeFname(config);
    // const auto fname = "config/srcs.txt";
    ParticleVec srcs;
    int Nsrcs;
    
    switch (config.mode) {
        case Mode::READ :
            srcs = importParticles(fname);
            Nsrcs = srcs.size();
            break;

        case Mode::WRITE : {
            srcs = makeParticles<uniform_real_distribution<double>>(config);
            Nsrcs = config.nsrcs;

            ofstream srcFile(fname);
            for (const auto& src : srcs) srcFile << *src;
            break;
        }
        default : 
            throw std::runtime_error("Invalid mode");
    }

    Node::setNodeParams(config);
    const int order = Node::getExpansionOrder();

    cout << " Mode:            " << (config.mode == Mode::READ ? "READ" : "WRITE") << '\n';
    cout << " Source file:     " << fname << '\n';
    cout << " # sources:       " << Nsrcs << '\n';
    cout << " Root length:     " << config.L << '\n';
    cout << " Error tol.:      " << config.EPS << '\n';
    cout << " Expansion order: " << order << '\n';
    cout << " Max node parts:  " << config.maxNodeParts << '\n' << '\n';

    // ==================== Set up domain ==================== //
    cout << " Setting up domain...\n";
    auto start = chrono::high_resolution_clock::now();
    auto start_ = start;

    shared_ptr<Node> root;
    if (Nsrcs > config.maxNodeParts)
        root = make_shared<Stem>(srcs, 0, nullptr);
    else
        root = make_shared<Leaf>(srcs, 0, nullptr);
    
    root->buildLists();

    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double, milli> duration_ms = end - start;
    cout << "   # Nodes: " << Node::getNumNodes() << '\n';
    cout << "   Elapsed time: " << duration_ms.count() << " ms\n";

    // ==================== Tests ==================== //
    root->labelNodes();

    std::ofstream nodeFile("out/nodes.txt");
    root->printNode(nodeFile);

    return 0;

    // ==================== Upward pass ==================== //
    cout << " Computing upward pass...\n";
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
    chrono::duration<double, milli> fmm_duration_ms = end - start_;
    cout << "   Elapsed time: " << duration_ms.count() << " ms\n";

    printSols(srcs, "out/phi.txt", "out/fld.txt");

    cout << " FMM total elapsed time: " << fmm_duration_ms.count() << " ms\n\n";

    // ==================== Compute direct ==================== //
    if (!config.evalDirect) return 0;
    cout << " Computing direct..." << endl;
    start = chrono::high_resolution_clock::now();

    auto phisAnl = root->getDirectPhis();

    end = chrono::high_resolution_clock::now();
    duration_ms = end - start;
    cout << "   Elapsed time: " << duration_ms.count() << " ms\n";

    ofstream phiAnlFile("out/phiAnl.txt");
    phiAnlFile << setprecision(15) << scientific;
    for (const auto& phi : phisAnl)
        phiAnlFile << phi.real() << '\n';

    cout << " Computing direct fld..." << endl;
    start = chrono::high_resolution_clock::now();

    auto fldsAnl = root->getDirectFlds();

    end = chrono::high_resolution_clock::now();
    duration_ms = end - start;
    cout << "   Elapsed time: " << duration_ms.count() << " ms\n";

    ofstream fldAnlFile("out/fldAnl.txt");
    fldAnlFile << setprecision(15) << scientific;
    for (const auto& fld : fldsAnl)
        fldAnlFile << fld << '\n';

    return 0;
}