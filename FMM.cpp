#include <chrono>
#include <fstream>
#include <iostream>
#include "clock.h"
#include "config.h"
#include "fmm.h"

using namespace std;

extern auto t = ClockTimes();

int main() {
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

    cout << " Mode:              " << (config.mode == Mode::READ ? "READ" : "WRITE") << '\n';
    cout << " Source file:       " << fname << '\n';
    cout << " # Sources:         " << Nsrcs << '\n';
    cout << " Root length:       " << config.L << '\n';
    cout << " Error tol.:        " << config.EPS << '\n';
    cout << " Expansion order:   " << order << '\n';
    cout << " Exponential order: " << Node::getExponentialOrder() << '\n';
    cout << " Max node parts:    " << config.maxNodeParts << '\n' << '\n';

    // ==================== Build tables ===================== //
    auto fmm_start = Clock::now();
    
    Node::buildTables(config);

    Node::buildRotationMats();

    // ==================== Set up domain ==================== //
    cout << " Setting up domain...\n";
    auto start = Clock::now();

    shared_ptr<Node> root;
    if (Nsrcs > config.maxNodeParts)
        root = make_shared<Stem>(srcs, 0, nullptr);
    else
        root = make_shared<Leaf>(srcs, 0, nullptr);

    root->buildLists();

    auto end = Clock::now();
    Time duration_ms = end - start;
    cout << "   # Nodes: " << Node::getNumNodes() << '\n';
    cout << "   # Leaves: " << Leaf::getNumLeaves() << '\n';
    cout << "   Elapsed time: " << duration_ms.count() << " ms\n";

    // ==================== Tests ==================== //
    std::ofstream nodeFile("out/nodes.txt");
    root->printNode(nodeFile);

    // ==============================================
    // root->ffieldTest(1,10,10);
    // ==============================================   
    // root->mpoleToExpToLocalTest();
    // ==============================================   
    // root->nfieldTest();

    // return 0;
    // ==================== Upward pass ==================== //
    cout << " Computing upward pass...\n";
    start = Clock::now();

    root->buildMpoleCoeffs();

    end = Clock::now();
    duration_ms = end - start;
    cout << "   Elapsed time: " << duration_ms.count() << " ms\n";

    // ==================== Downward pass ==================== //
    cout << " Computing downward pass...\n";
    start = Clock::now();

    root->propagateExpCoeffs();

    root->buildLocalCoeffs();

    end = Clock::now();
    duration_ms = end - start;

    cout << "   Elapsed time: " << duration_ms.count() << " ms\n";
    cout << "   Elapsed time (M2X): " << t.M2X.count() << " ms\n";
    cout << "   Elapsed time (X2X): " << t.X2X.count() << " ms\n";
    cout << "   Elapsed time (X2L): " << t.X2L.count() << " ms\n";
    cout << "   Elapsed time (P2L): " << t.P2L.count() << " ms\n";
    cout << "   Elapsed time (L2L): " << t.L2L.count() << " ms\n";

    // ================== Evaluate solutions ================= //
    cout << " Evaluating solutions...\n";
    start = Clock::now();

    Leaf::evaluateSols();

    end = Clock::now();
    duration_ms = end - start;
    Time fmm_duration_ms = end - fmm_start;

    cout << "   Elapsed time: " << duration_ms.count() << " ms\n";
    cout << "   Elapsed time (L2P): " << t.L2P.count() << " ms\n";
    cout << "   Elapsed time (M2P): " << t.M2P.count() << " ms\n";
    cout << "   Elapsed time (P2P): " << t.P2P.count() << " ms\n";
    cout << " FMM total elapsed time: " << fmm_duration_ms.count() << " ms\n";

    printSols(srcs, "out/phi.txt", "out/fld.txt");

    // ================== Compute direct ===================== //
    if (!config.evalDirect) return 0;

    root->resetSols();

    cout << "\n Computing direct...\n";
    start = Clock::now();

    root->evalSelfSols();

    end = Clock::now();
    duration_ms = end - start;
    cout << "   Elapsed time: " << duration_ms.count() << " ms\n";

    printSols(srcs, "out/phiAnl.txt", "out/fldAnl.txt");

    return 0;
}