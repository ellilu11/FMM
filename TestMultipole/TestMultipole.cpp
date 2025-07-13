#include "pch.h"
#include "CppUnitTest.h"

#include <complex>
#include <fstream>
#include <memory>
#include <random>
#include <vector>
#include "../node.h"
#include "../stem.h"
#include "../stem.cpp"
#include "../leaf.h"
#include "../leaf.cpp"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace TestMultipole
{
	TEST_CLASS(TestMultipole)
	{
	public:
		
		TEST_METHOD_INITIALIZE(initObs)
		{
			const double L = 10.0;
			const double c = 5.0;
			const int Nobs = 10;
			std::vector<cmplx> posObs;
			for (int n = 0; n < Nobs; ++n) {
				double r = static_cast<double>(n) / static_cast<double>(Nobs);
				posObs.push_back(c * L * (1.0 + iu * (r - 1 / 2)));
			}

		}

		TEST_METHOD_INITIALIZE(initSrc)
		{
			constexpr int N = 1000;
			constexpr double Q = 1.0;
			constexpr double eps = 1.0E-6;
			const double L = 10.0;

			const int Nlvl = ceil(log(N) / log(4.0));
			const int P = ceil(-log(eps) / log(2));

			std::random_device rd;
			std::mt19937 gen(rd());

			std::uniform_real_distribution<double> real(-L / 2, L / 2);
			std::uniform_real_distribution<double> imag(-L / 2, L / 2);

			std::vector<cmplx> pos;
			std::vector<double> qs;

			for (int n = 0; n < N; ++n) {
				cmplx z(real(gen), imag(gen));
				pos.push_back(z);
				qs.push_back(Q);
			}
			std::shared_ptr<Node> master;
			if (N > 1 && Nlvl)
				master = std::make_shared<Stem>(pos, qs, 0, L, Nlvl, P);
			else
				master = std::make_shared<Leaf>(pos, qs, 0, L, Nlvl, P);

			master->buildCoeffs();

		}

		TEST_METHOD(diff_farfield)
		{

		}
	};
}
