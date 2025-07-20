#define _USE_MATH_DEFINES
#include <chrono>
#include <cmath>
#include <fstream>
#include <random>
#include "node.h"
#include "node.cpp"
#include "stem.cpp"
#include "leaf.cpp"

using namespace std;

namespace Param {
    extern const int DIM = 2;
    constexpr double L = 10.0;
    constexpr double EPS = 1.0E-10;
    const int P = ceil(-log(EPS) / log(2)); // # terms in multipole expansion
}

int main(int argc, char *argv[])
{
    // ==================== Populate domain ==================== //
    random_device rd;
    mt19937 gen(rd());

    uniform_real_distribution<double> real(0, Param::L / 2);
    uniform_real_distribution<double> imag(0, Param::L / 2);

    //normal_distribution<double> u(0, 0.2 * L);
    //uniform_real_distribution<double> th(0, 2.0 * M_PI);

    constexpr int N = 100;
    constexpr double Q = 1.0;
    const int Nlvl = ceil(log(N) / log(4.0));

    vector<cmplx> psn;
    vector<double> qs;

    for (int n = 0; n < N; ++n) {
        cmplx z(real(gen), imag(gen));
        //const auto R = std::abs(u(gen));
        //cmplx z(R * std::cos(th(gen)), R * std::sin(th(gen)));
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

    root->iListTest(); // Test interaction list finding

    // ==================== Upward pass ==================== //
    cout << " Computing upward pass..." << endl;
    start = chrono::high_resolution_clock::now();

    root->buildMpoleCoeffs(Param::P);

    end = chrono::high_resolution_clock::now();
    duration_ms = end - start;
    cout << "   Elapsed time: " << duration_ms.count() << " ms\n";

    std::ofstream mpoleCoeffFile;
    mpoleCoeffFile.open("out/mpolecoeffs.txt");
    root->printMpoleCoeffs(mpoleCoeffFile);

    /*constexpr int NOBS = 1000;
    for (int p = P; p <= P; ++p) {
        root->buildMpoleCoeffs(p);
        root->ffieldTest(p, NOBS);
    }*/

    // ==================== Downward pass ==================== //
    cout << " Computing downward pass..." << endl;
    start = chrono::high_resolution_clock::now();

    root->buildLocalCoeffs(Param::P);

    end = chrono::high_resolution_clock::now();
    duration_ms = end - start;
    cout << "   Elapsed time: " << duration_ms.count() << " ms\n";

    std::ofstream phiFile;
    phiFile.open("out/nf.txt");
    root->printPhi(phiFile);

    std::ofstream localCoeffFile;
    localCoeffFile.open("out/localcoeffs.txt");
    root->printLocalCoeffs(localCoeffFile);

    return 0;
}