#include <chrono>
#include <fstream>
#include <iostream>
#include "config.h"
#include "fmm.h"

using namespace std;

extern constexpr int DIM = 3;

extern std::chrono::duration<double, std::milli> t_M2X{ 0 };
extern std::chrono::duration<double, std::milli> t_X2X{ 0 };
extern std::chrono::duration<double, std::milli> t_X2L{ 0 };
extern std::chrono::duration<double, std::milli> t_P2L{ 0 };
extern std::chrono::duration<double, std::milli> t_L2L{ 0 };
extern std::chrono::duration<double, std::milli> t_L2P{ 0 };
extern std::chrono::duration<double, std::milli> t_dir{ 0 };

int main() {
    Config config("config/config.txt");

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

    cout << " Mode:              " << (config.mode == Mode::READ ? "READ" : "WRITE") << '\n';
    cout << " Source file:       " << fname << '\n';
    cout << " # sources:         " << Nsrcs << '\n';
    cout << " Root length:       " << config.L << '\n';
    cout << " Error tol.:        " << config.EPS << '\n';
    cout << " Expansion order:   " << order << '\n';
    cout << " Exponential order: " << Node::getExponentialOrder() << '\n';
    cout << " Max node parts:    " << config.maxNodeParts << '\n' << '\n';

    // ==================== Build tables ===================== //
    auto fmm_start = chrono::high_resolution_clock::now();
    Node::buildTables(config);
    Node::buildRotationMats();

    // ==================== Set up domain ==================== //
    cout << " Setting up domain...\n";
    auto start = chrono::high_resolution_clock::now();

    // std::shared_ptr<Node> root = std::make_shared<Stem>(srcs, 0, nullptr);
    shared_ptr<Node> root;
    if (Nsrcs > config.maxNodeParts)
        root = make_shared<Stem>(srcs, 0, nullptr);
    else
        root = make_shared<Leaf>(srcs, 0, nullptr);

    root->buildLists();

    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double, milli> duration_ms = end - start;
    cout << "   # Nodes: " << Node::getNumNodes() << '\n';
    cout << "   # Leaves: " << Leaf::getNumLeaves() << '\n';
    cout << "   Elapsed time: " << duration_ms.count() << " ms\n";

    // ==================== Tests ==================== //
    //int ntrials = 1000;

    //for (int trial = 0; trial < ntrials; ++trial) {
    //    cout << "Trial # " << trial << '\n';
    //    root->resetNode();
    //    root->labelNodes();
    //}
    // ==============================================

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
    start = chrono::high_resolution_clock::now();

    root->buildMpoleCoeffs();

    end = chrono::high_resolution_clock::now();
    duration_ms = end - start;
    cout << "   Elapsed time: " << duration_ms.count() << " ms\n";

    // ==================== Downward pass ==================== //
    cout << " Propagating exponential coeffs...\n";
    start = chrono::high_resolution_clock::now();

    root->propagateExpCoeffs();

    end = chrono::high_resolution_clock::now();
    duration_ms = end - start;
    cout << "   Elapsed time: " << duration_ms.count() << " ms\n";
    cout << "   Elapsed time (M2X): " << t_M2X.count() << " ms\n";
    cout << "   Elapsed time (X2X): " << t_X2X.count() << " ms\n";

    cout << " Computing downward pass...\n";
    start = chrono::high_resolution_clock::now();

    root->buildLocalCoeffs();

    end = chrono::high_resolution_clock::now();
    duration_ms = end - start;

    cout << "   Elapsed time: " << duration_ms.count() << " ms\n";
    cout << "   Elapsed time (X2L): " << t_X2L.count() << " ms\n";
    cout << "   Elapsed time (P2L): " << t_P2L.count() << " ms\n";
    cout << "   Elapsed time (L2L): " << t_L2L.count() << " ms\n";

    // ================== Evaluate solutions ================= //
    cout << " Evaluating solutions...\n";
    start = chrono::high_resolution_clock::now();

    Leaf::evaluateSols();

    end = chrono::high_resolution_clock::now();
    duration_ms = end - start;
    chrono::duration<double, milli> fmm_duration_ms = end - fmm_start;

    cout << "   Elapsed time: " << duration_ms.count() << " ms\n";
    cout << "   Elapsed time (L2P): " << t_L2P.count() << " ms\n";
    cout << "   Elapsed time (Direct): " << t_dir.count() << " ms\n";
    cout << " FMM total elapsed time: " << fmm_duration_ms.count() << " ms\n\n";

    printSols(srcs, "out/phi.txt", "out/fld.txt");

    // ================== Compute direct ===================== //
    if (!config.evalDirect) return 0;

    root->resetSols();

    cout << " Computing direct..." << endl;
    start = chrono::high_resolution_clock::now();

    root->evalSelfSols();

    end = chrono::high_resolution_clock::now();
    duration_ms = end - start;
    cout << "   Elapsed time: " << duration_ms.count() << " ms\n";

    printSols(srcs, "out/phiAnl.txt", "out/fldAnl.txt");

    return 0;
}