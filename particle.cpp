#include "particle.h"

using namespace std;

ParticleVec importParticles(const std::string& fname) {
    ifstream inFile(fname);
    if (!inFile) throw std::runtime_error("Unable to find file");
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

template <class dist0, class dist1 = dist0>
ParticleVec makeParticles(const Config& config) 
{
    ParticleVec particles;
    constexpr double e = 1.0;
    constexpr double M = 1.0;

    random_device rd;
    mt19937 gen(rd());

    T rand0(param0, param1);
    U rand1(param0, param1);
    V rand2(param0, param1);

    switch (config.dist) {
        case Dist::UNIFORM:
            rand0 = dist0(-config.L/2, config.L/2);
            rand1 = dist1(-config.L/2, config.L/2);
            break;
        case Dist::GAUSSIAN:
            rand0 = dist0(0, 1);
            rand1 = dist1(0, 1);
            break;
        default:
            throw std::runtime_error("Invalid distribution");
    }

    for (int n = 0; n < config.nsrcs; ++n) {
        cmplx z = [&] {
            switch (config.dist) {
                case Dist::UNIFORM:
                    return cmplx(rand0(gen), rand1(gen));
                case Dist::GAUSSIAN: {
                    const double R = sqrt(-2 * log(rand0(gen)));
                    const double th = 2.0 * 3.1415927 * rand1(gen);
                    return cmplx(R * cos(th), R * sin(th));
                }
            } } ();

        int pm = [&] () -> int {
            switch (config.cdist) {
                case ChargeDist::PLUS:  return 1;
                case ChargeDist::MINUS: return 0;
                case ChargeDist::DIP:   
                    return z.real() > 0;
                case ChargeDist::QUAD: {
                    auto idx = bools2Idx(z > 0);
                    return (idx == 0 || idx == 3);
                }
                case ChargeDist::RAND: {
                    uniform_int_distribution randi(0, 1);
                    return randi(gen);
                }
            } } ();
        pm = 2*pm-1;
        particles.emplace_back(make_shared<Particle>(z, pm*e, M));
    }
    return particles;
}

void printSols(ParticleVec& particles, 
    const std::string& pname, 
    const std::string& fname) {
    ofstream phiFile(pname);
    ofstream fldFile(fname);

    for (const auto& p : particles) {
        p->printPhi(phiFile);
        p->printFld(fldFile);
    }
}