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
    constexpr double Q = 1.0;
    constexpr double M = 1.0;

    random_device rd;
    mt19937 gen(rd());

    T rand0(param0, param1);
    U rand1(param0, param1);
    // T rand0(param0, param1);
    // U rand1(0, 2.0 * M_PI);

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
        particles.emplace_back(make_shared<Particle>(z, Q, M));
    }
    return particles;
}

//realVec getCharges(ParticleVec& particles) {
//    realVec charges;
//    for (const auto& p : particles)
//        charges.push_back(p->getCharge());
//    return charges;
//}