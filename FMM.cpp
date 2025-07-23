#include "FMM.h"

using namespace std;

namespace Param {
    extern constexpr int    DIM     = 2;
    extern constexpr double L       = 10.0;
    extern constexpr double EPS     = 1.0E-3;
}

int main(int argc, char *argv[])
{

    // ==================== Populate domain ==================== //
    random_device rd;
    mt19937 gen(rd());

    uniform_real_distribution<double> real(-Param::L / 2, Param::L / 2);
    uniform_real_distribution<double> imag(-Param::L / 2, Param::L / 2);

    // uniform_real_distribution<double> real(0, Param::L / 2);
    // uniform_real_distribution<double> imag(0, Param::L / 2);

    // normal_distribution<double> u(0, 0.2 * Param::L);
    // uniform_real_distribution<double> th(0, 2.0 * M_PI);

    constexpr int N = 1000;
    constexpr double Q = 1.0;
    const int Nlvl = ceil(log(N) / log(4.0)) - 1;

    vector<cmplx> psn;
    vector<double> qs;

    for (int n = 0; n < N; ++n) {
        cmplx z(real(gen), imag(gen));
        // const auto R = std::abs(u(gen));
        // cmplx z(R * std::cos(th(gen)), R * std::sin(th(gen)));
        psn.push_back(z);
        qs.push_back(Q);
    }

    // ==================== Partition domain ==================== //
    cout << " Partitioning domain..." << endl;
    auto start = chrono::high_resolution_clock::now();

    shared_ptr<Node> root;
    if (N > 1)
        root = make_shared<Stem>(psn, qs, 0, Param::L, Nlvl, 0, nullptr);
    else
        root = make_shared<Leaf>(psn, qs, 0, Param::L, Nlvl, 0, nullptr);
    
    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double, milli> duration_ms = end - start;
    cout << "   Elapsed time: " << duration_ms.count() << " ms\n";

    // ==================== Upward pass ==================== //
    cout << " Computing upward pass..." << endl;
    start = chrono::high_resolution_clock::now();

    Node::buildBinomTable();
    root->buildMpoleCoeffs();

    // constexpr int NOBS = 1000;
    // root->ffieldTest(NOBS);

    end = chrono::high_resolution_clock::now();
    duration_ms = end - start;
    cout << "   Elapsed time: " << duration_ms.count() << " ms\n";

    std::ofstream mpoleCoeffFile;
    mpoleCoeffFile.open("out/mpolecoeffs.txt");
    root->printMpoleCoeffs(mpoleCoeffFile);

    // ==================== Downward pass ==================== //
    //cout << " Computing downward pass..." << endl;
    //start = chrono::high_resolution_clock::now();

    // root->buildLocalCoeffs();

    // root->mpoleToLocalTest();
    root->nfieldTest();

    //end = chrono::high_resolution_clock::now();
    //duration_ms = end - start;
    //cout << "   Elapsed time: " << duration_ms.count() << " ms\n";
    
    std::ofstream localCoeffFile;
    localCoeffFile.open("out/localcoeffs.txt");
    root->printLocalCoeffs(localCoeffFile);

    return 0;
}