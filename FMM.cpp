#include "FMM.h"

using namespace std;

namespace Param {
    extern constexpr int    DIM     = 2;
    extern constexpr double L       = 10.0;
    extern constexpr double EPS     = 1.0E-1;
}

enum class Mode {
    READ_FROM_FILE,
    GEN_UNIFORM,
    GEN_GAUSSIAN
};

template <class T, class U = T>
ParticleVec makeRNGParticles(const int N, const double param0, const double param1) {

    ParticleVec particles;
    constexpr double Q = 1.0;
    constexpr double M = 1.0;

    random_device rd;
    mt19937 gen(rd());

    T rand0(param0, param1);
    U rand1(param0, param1);
    // T rand0(param0, param1);
    // U rand1(0, 2.0 * 3.1415927);

    for (int n = 0; n < N; ++n) {
        //if constexpr (is_integral<T>::uniform_real_distribution<double>)
        //    z(rand0(gen), imag(gen));
        //else if constexpr (is_integral<T>::normal_distribution<double>) {
        //    const auto R = abs(rand0(gen));
        //    z(R * cos(rand1(gen)), R * sin(rand1(gen)));
        //}
        //else
        //    throw std::runtime_error("Invalid probability distribution");

        cmplx z(rand0(gen), rand1(gen));
        // const auto R = abs(rand0(gen));
        // cmplx z(R * cos(rand1(gen)), R * sin(rand1(gen)));
        particles.emplace_back( make_shared<Particle>(z, Q, M) );
    }
    return particles;
}

int main(int argc, char *argv[])
{
    // ==================== Make particles ==================== //
    int Nsrcs = 10000;
    Mode mode = Mode::GEN_UNIFORM;
    ParticleVec srcs;

    switch (mode) {
        // case Mode::READ_FROM_FILE :
            // srcs = import_particles("out/srcs.txt");
            // break;

        case Mode::GEN_UNIFORM : 
            srcs = makeRNGParticles<uniform_real_distribution<double>>
                    (Nsrcs, -Param::L/2, Param::L/2);
            break;
        case Mode::GEN_GAUSSIAN :
            srcs =
                makeRNGParticles<normal_distribution<double>, uniform_real_distribution<double>>
                    (Nsrcs, 0.0, 0.1*Param::L);
            break;
        default : 
            throw std::runtime_error("Invalid mode");
    }
    ofstream srcFile;
    srcFile.open("out/srcs.txt");
    for (const auto& src : srcs) srcFile << *src;

    // ==================== Partition domain ==================== //
    cout << " Partitioning domain..." << endl;
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
    cout << " Computing upward pass...   (" << " P = " << P << " )\n";
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
    cout << " Computing downward pass... (" << " P = " << P << " )\n";
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