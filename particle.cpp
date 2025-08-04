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

template <class dist0, class dist1 = dist0, class dist2 = dist0>
ParticleVec makeParticles(const Config& config) 
{
    ParticleVec particles;
    constexpr double e = 1.0;
    constexpr double M = 1.0;

    random_device rd;
    mt19937 gen(rd());

    dist0 rand0;
    dist1 rand1;
    dist2 rand2;

    switch (config.dist) {
        case Dist::UNIFORM:
            rand0 = dist0(-config.L/2, config.L/2);
            rand1 = dist1(-config.L/2, config.L/2);
            rand2 = dist2(-config.L/2, config.L/2);
            //rand0 = dist0(0.0, config.L/2);
            //rand1 = dist1(0.0, config.L/2);
            //rand2 = dist2(0.0, config.L/2);
            break;
        case Dist::GAUSSIAN:
            rand0 = dist0(0, 1);
            rand1 = dist1(0, 1);
            break;
        default:
            throw std::runtime_error("Invalid distribution");
    }

    for (int n = 0; n < config.nsrcs; ++n) {
        vec3d X;
        X = [&] {
            switch (config.dist) {
                case Dist::UNIFORM:
                    return vec3d{ rand0(gen), rand1(gen), rand2(gen) };
                case Dist::GAUSSIAN: {
                    const double r = sqrt(-2 * log(rand0(gen))); // fix
                    const double th = 2.0 * 3.1415927 * rand1(gen); // fix
                    const double ph = rand2(gen);
                    return toCart(vec3d(r,th,ph));
                }
            } } ();
        
        auto x = X[0], y = X[1], z = X[2];
        auto idx = bools2Idx(X > zeroVec );
        int pm = [&] () -> int {
            switch (config.cdist) {
                case ChargeDist::PLUS:  return 1;
                case ChargeDist::MINUS: return 0;
                case ChargeDist::DIP:   return x > 0;
                case ChargeDist::QUAD:
                    return (idx == 0 || idx == 3 || idx == 4 || idx == 7);
                case ChargeDist::OCT:
                    return (idx == 0 || idx == 3 || idx == 5 || idx == 6);
                case ChargeDist::RAND: {
                    uniform_int_distribution randi(0, 1);
                    return randi(gen);
                }
            } } () * 2 - 1;
        particles.emplace_back(make_shared<Particle>(X, pm*e, M));
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