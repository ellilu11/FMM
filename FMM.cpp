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
extern std::chrono::duration<double, std::milli> t_L2L{ 0 };
extern std::chrono::duration<double, std::milli> t_direct{ 0 };

int main()
{
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

        case Mode::GEN : {
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

    cout << " Mode:              " << (config.mode == Mode::READ ? "READ" : "GEN") << '\n';
    cout << " Source file:       " << fname << '\n';
    cout << " # sources:         " << Nsrcs << '\n';
    cout << " Root length:       " << config.L << '\n';
    cout << " Error tol.:        " << config.EPS << '\n';
    cout << " Expansion order:   " << order << '\n';
    cout << " Exponential order: " << Node::getExponentialOrder() << '\n';
    cout << " Max node parts:    " << config.maxNodeParts << '\n' << '\n';

    // ==================== Build tables ==================== //
    cout << " Building tables..\n";
    auto start = chrono::high_resolution_clock::now();

    Node::buildTables(config);
    Node::buildRotationMats();

    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double, milli> duration_ms = end - start;
    cout << "   Elapsed time: " << duration_ms.count() << " ms\n";

    // ==================== Set up domain ==================== //
    cout << " Setting up domain...\n";
    start = chrono::high_resolution_clock::now();

    shared_ptr<Node> root;
    // root = make_shared<Stem>(srcs, 0, nullptr);
    if (Nsrcs > config.maxNodeParts)
        root = make_shared<Stem>(srcs, 0, nullptr);
    else
        root = make_shared<Leaf>(srcs, 0, nullptr);

    root->buildLists();

    end = chrono::high_resolution_clock::now();
    duration_ms = end - start;
    cout << "   Elapsed time: " << duration_ms.count() << " ms\n";

    // ==================== Tests ==================== //
    // root->setRandNodeStats();
    std::ofstream nodeFile("out/nodes.txt");
    root->printNode(nodeFile);

    //const pair2d angles(PI, 0);
    //for (int l = 0; l < 4; ++l)
    //    cout << wignerD_l(angles, l) << "\n\n" << wignerD_l(angles, l).inverse() << "\n\n";

    // ==============================================
    // root->ffieldTest(1,10,10);
    // ==============================================   
    // root->mpoleToExpToLocalTest();
    // ==============================================   
    // root->nfieldTest();
    //
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
    cout << "   Elapsed time (L2L): " << t_L2L.count() << " ms\n";
    cout << "   Elapsed time (direct): " << t_direct.count() << " ms\n";

    printSols(srcs, "out/phi.txt", "out/fld.txt");

    // ==================== Compute direct ==================== //
    if (!config.evalDirect) return 0;
    cout << " Computing direct phi..." << endl;
    start = chrono::high_resolution_clock::now();

    auto phisAnl = root->getDirectPhis();

    end = chrono::high_resolution_clock::now();
    duration_ms = end - start;
    cout << "   Elapsed time: " << duration_ms.count() << " ms\n";

    ofstream phiAnlFile("out/phiAnl.txt");
    for (const auto& phi : phisAnl)
        phiAnlFile << phi << '\n';

    cout << " Computing direct fld..." << endl;
    start = chrono::high_resolution_clock::now();

    auto fldsAnl = root->getDirectFlds();

    end = chrono::high_resolution_clock::now();
    duration_ms = end - start;
    cout << "   Elapsed time: " << duration_ms.count() << " ms\n";

    ofstream fldAnlFile("out/fldAnl.txt");
    for (const auto& fld : fldsAnl)
        fldAnlFile << fld << '\n';

    return 0;
}