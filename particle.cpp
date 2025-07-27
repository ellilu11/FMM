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

template <class T, class U = T>
ParticleVec makeParticles(const int N, 
    const double param0, const double param1, 
    const Dist dist, const ChargeDist cdist) 
{
    ParticleVec particles;
    constexpr double e = 1.0;
    constexpr double M = 1.0;

    random_device rd;
    mt19937 gen(rd());

    T rand0(param0, param1);
    U rand1(param0, param1);

    cmplx z;
    int pm;

    for (int n = 0; n < N; ++n) {
        switch (dist) {
            case Dist::UNIFORM : 
                z = cmplx(rand0(gen), rand1(gen));
                break;
            case Dist::GAUSSIAN: {
                const double R = sqrt(-2 * log(rand0(gen)));
                const double th = 2.0 * 3.1415927 * rand1(gen);
                z = cmplx(R * cos(th), R * sin(th));
                break;
            }
            default :
                throw std::runtime_error("Invalid distribution");
                break;
        }

        switch (cdist) {
            case ChargeDist::PLUS : pm = 1;            break;
            case ChargeDist::MINUS: pm = 0;            break;
            case ChargeDist::DIP:   pm = z.real() > 0; break;
            case ChargeDist::QUAD: {
                auto idx = bools2Idx(z > 0);
                pm = (idx == 0 || idx == 3);
                break;
            }
            case ChargeDist::RAND: {
                uniform_int_distribution randi(0,1);
                pm = randi(gen);
                break;
            }
            default :
                throw std::runtime_error("Invalid charge distribution");
                break;
        }
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