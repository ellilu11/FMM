#include "particle.h"

using namespace std;

ParticleVec importParticles(const std::string& fname) {
    ifstream inFile(fname);
    string line;
    ParticleVec particles;

    while (getline(inFile, line)) {
        istringstream iss(line);
        
        Particle p;
        if (iss >> p) 
            particles.emplace_back(make_shared<Particle>(p));
        else
            throw std::runtime_error("Unable to parse line");
    }
    return particles;
}

template <class T, class U = T>
ParticleVec makeRNGParticles(const int N, const double param0, const double param1) {
    ParticleVec particles;
    constexpr double e = 1.0;
    constexpr double M = 1.0;

    random_device rd;
    mt19937 gen(rd());

    T rand0(param0, param1);
    U rand1(param0, param1);

    uniform_int_distribution pm(0, 0);

    for (int n = 0; n < N; ++n) {
        cmplx z(rand0(gen), rand1(gen));
        //const double R = sqrt(-2 * log(rand0(gen)));
        //const double th = 2.0 * 3.1415927 * rand1(gen);
        //cmplx z(R * cos(th), R * sin(th));
        particles.emplace_back(make_shared<Particle>(z, pow(-1,pm(gen))*e, M));
    }
    return particles;
}

void printSols(ParticleVec& particles, const std::string& pname, const std::string& fname) {
    ofstream phiFile(pname);
    ofstream fldFile(fname);

    for (const auto& p : particles) {
        p->printPhi(phiFile);
        p->printFld(fldFile);
    }
}