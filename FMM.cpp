#include <chrono>
#include <fstream>
#include <iostream>
#include "config.h"
#include "fmm.h"

using namespace std;

extern constexpr int DIM = 3;

int main(int argc, char *argv[])
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

    cout << " Mode:           " << (config.mode == Mode::READ ? "READ" : "GEN") << '\n';
    cout << " Source file:    " << fname << '\n';
    cout << " # sources:      " << Nsrcs << '\n';
    cout << " Root length:    " << config.L << '\n';
    cout << " Error tol.:     " << config.EPS << '\n';
    cout << " Max node parts: " << config.maxNodeParts << '\n' << '\n';

    // ==================== Set up domain ==================== //
    cout << " Setting up domain...\n";
    auto start = chrono::high_resolution_clock::now();

    Node::setNodeParams(config);
    shared_ptr<Node> root;
    if (Nsrcs > config.maxNodeParts)
        root = make_shared<Stem>(srcs, 0, nullptr);
    else
        root = make_shared<Leaf>(srcs, 0, nullptr);
    
    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double, milli> duration_ms = end - start;
    cout << "   Elapsed time: " << duration_ms.count() << " ms\n";

    // ==================== Tests ==================== //
    root->setRandNodeStats();
    std::ofstream nodeFile("out/nodes.txt");
    root->printNode(nodeFile);

    Node::setNodeParams(config);
    cout << " Expansion order: " << Node::getExpansionOrder() << "\n";

    cout << " Building rotation matrices...\n";
    start = chrono::high_resolution_clock::now();

    Node::buildRotationMats();
    int l = 1;
    for (int dir = 0; dir < 1; ++dir) {
        auto mat = Node::getRotationMatrixAlongDir(dir)[l];
        // cout << mat << '\n' << mat.adjoint() << '\n' << '\n';
        // cout << mat * mat.adjoint() << '\n' << '\n';
    }
    // cout << rotationMatrix(pair2d(0, PI), 1) << '\n';

    end = chrono::high_resolution_clock::now();
    duration_ms = end - start;
    cout << "   Elapsed time: " << duration_ms.count() << " ms\n";

    //const int order = Node::getExpansionOrder();
    //Node::buildTables();

    //double th = PI/4;
    //double ph = 0;

    //for (int l = 0; l < order; ++l) {
    //    for (int m = -l; m <= l; ++m) {
    //        auto Ylm = Node::legendreLM(th, l, abs(m));
    //        cout << (abs(Ylm) > 1.0E-9 ? Ylm : 0.0 ) << ' ';
    //    }
    //    cout << '\n';
    //}

    // ==================== Upward pass ==================== //
    const int order = Node::getExpansionOrder();
    cout << " Computing upward pass...   (" << " Expansion order: " << order << ")\n";
    start = chrono::high_resolution_clock::now();

    Node::buildTables();
    root->buildMpoleCoeffs();

    end = chrono::high_resolution_clock::now();
    duration_ms = end - start;
    cout << "   Elapsed time: " << duration_ms.count() << " ms\n";

    return 0;

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
        phiAnlFile << phi << '\n';

    auto fldsAnl = root->getDirectFlds();

    ofstream fldAnlFile("out/fldAnl.txt");
    for (const auto& fld : fldsAnl)
        fldAnlFile << fld << '\n';

    return 0;
}