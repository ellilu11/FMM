#pragma once

#include "fmm.h"

class Particle;

using ParticleVec = std::vector<std::shared_ptr<Particle>>;

class Particle {
public :
    Particle() = default;
    Particle(const cmplx z, const double q, const double m)
        : pos(z), charge(q), mass(m)
    {
    };

    cmplx getPos() const { return pos; }
    const double getCharge() const { return charge; }
    const double getMass() const { return mass; }

    friend std::ostream& operator<<(std::ostream& os, Particle& p) {
        os << p.pos << " " << p.charge << " " << p.mass << '\n';
        return os;
    }

    friend std::istream& operator>>(std::istream& is, Particle& p) {
        is >> p.pos >> p.charge >> p.mass;
        return is;
    }

private :
    cmplx pos;
    double charge;
    double mass;
    // cmplx vel;
};

