#include "particle.h"

using namespace std;

ParticleVec importParticles(const std::filesystem::path& fpath) {
    ifstream inFile(fpath);
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

    dist0 rand0(0, 1);
    dist1 rand1(0, 1);
    dist2 rand2(0, 1);

    if (config.dist == Dist::UNIFORM) {
        double lim = config.L/2.0;
        rand0 = dist0(-lim, lim);
        rand1 = dist1(-lim, lim);
        rand2 = dist2(-lim, lim);
    }

    for (int n = 0; n < config.nsrcs; ++n) {
        vec3d X = [&] {
            double r, th, ph, z;

            switch (config.dist) {
                case Dist::UNIFORM:
                    return vec3d( rand0(gen), rand1(gen), rand2(gen) );

                // 2D Gaussian at z = 0; TODO: 3D Gaussian
                case Dist::GAUSSIAN: { 
                    r = sqrt(-2.0 * log(rand0(gen))); 
                    th = PI / 2.0;
                    ph = 2.0 * PI * rand1(gen);
                    return Math::fromSph(vec3d(r,th,ph));
                }

                case Dist::SPHERE: {
                    r = 0.90 * (config.L / 2.0);
                    th = acos(2.0*rand0(gen) - 1.0);
                    ph = 2.0 * PI * rand1(gen);
                    return Math::fromSph(vec3d(r,th,ph));
                }

                case Dist::CYLINDER: {
                    r = 0.45 * (config.L / 2.0);
                    ph = 2.0 * PI * rand1(gen);
                    z = 0.90 * config.L * (rand1(gen) - 0.5);
                    return Math::fromCyl(vec3d(r, ph, z));
                }
            } 
            }();
        
        auto x = X[0], y = X[1], z = X[2];
        auto idx = Math::bools2Idx(X > zeroVec );
        int pm = [&] () -> int {
            switch (config.qdist) {
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
            } 
            }() * 2 - 1;
        particles.emplace_back(make_shared<Particle>(X, pm*e, M));
    }
    return particles;
}

void printSols(ParticleVec& particles, 
    const std::string& pname, const std::string& fname) 
{
    namespace fs = std::filesystem;
    fs::path dir = "out";
    std::error_code ec;

    if (fs::create_directory(dir, ec))
        std::cout << " Created directory " << dir.generic_string() << "/\n";
    else if (ec)
        std::cerr << " Error creating directory " << ec.message() << "\n";

    ofstream phiFile(dir/pname);
    ofstream fldFile(dir/fname);

    phiFile << setprecision(15) << scientific;
    fldFile << setprecision(15) << scientific;

    for (const auto& p : particles) {
        p->printPhi(phiFile);
        p->printFld(fldFile);
    }
}