#pragma once

#include <complex>
#include <vector>
#include "math.h"

class Particle;

using ParticleVec = std::vector<std::shared_ptr<Particle>>;

class Particle {
public :
    Particle(const cmplx z, const double q, const double m)
        : pos(z), charge(q), mass(m)
    {
    };

    cmplx getPos() const { return pos; }
    const double getCharge() const { return charge; }
    const double getMass() const { return mass; }

    friend std::ostream& operator<<(std::ostream& out, Particle& p) {
        out << p.pos << " " << p.charge << " " << p.mass << '\n';
        return out;
    }

    friend std::istream& operator>>(std::istream& in, Particle& p) {
        in >> p.pos >> p.charge >> p.mass;
        return in;
    }

private :
    cmplx pos;
    double charge;
    double mass;
    // cmplx vel;
};