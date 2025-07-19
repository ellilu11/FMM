#define _USE_MATH_DEFINES
#include <cmath>
#include <fstream>
#include <random>
#include "node.h"
#include "node.cpp"
#include "stem.cpp"
#include "leaf.cpp"

using namespace std;

/* namespace Param {

}*/

int main()
{
    // Define computational domain
    constexpr double L = 10.0;
    constexpr double EPS = 1.0E-10;
    const int P = ceil(-log(EPS) / log(2)); // # terms in multipole expansion

    random_device rd;
    mt19937 gen(rd());

    uniform_real_distribution<double> real(-L / 2, L / 2);
    uniform_real_distribution<double> imag(-L / 2, L / 2);

    //normal_distribution<double> u(0, 0.2 * L);
    //uniform_real_distribution<double> th(0, 2.0 * M_PI);

    // Populate domain
    constexpr int N = 1000;
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

    // Refine domain
    shared_ptr<Node> root;
    if (N > 1)
        root = make_shared<Stem>(psn, qs, 0, L, Nlvl, 0, nullptr);
    else
        root = make_shared<Leaf>(psn, qs, 0, L, Nlvl, 0, nullptr);

    // Tests
    root->iListTest();

    constexpr int NOBS = 1000;
    for (int p = P; p <= P; ++p) {
        // Upward pass: Aggregate multipole expansions
        root->buildCoeffs(p);
        root->ffieldTest(p, NOBS);
    }

    // Downward pass
    // root->buildLocalCoeffs(P);

	return 0;
}