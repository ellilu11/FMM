#pragma once

#include <random>
#include "config.h"
#include "math.h"

class Particle;

using ParticleVec = std::vector<std::shared_ptr<Particle>>;

class Particle {
public :
    Particle() = default;
    Particle(const vec3d z, const double q, const double m)
        : pos(z), charge(q), mass(m)
    {
    };

    vec3d getPos() const { return pos; }
    double getCharge() const { return charge; }
    double getMass() const { return mass; }

    const vec3d getPhi() const { return phi; }
    void setPhi(vec3d phi_) { phi = phi_; }
    void printPhi(std::ofstream& f) const {
        f << phi.real() << '\n';
    }

    void setFld(vec3d fld_) { fld = fld_; }
    void printFld(std::ofstream& f) const {
        f << fld << '\n';
    }

    friend std::ostream& operator<<(std::ostream& os, Particle& p) {
        os << p.pos << " " << p.charge << " " << p.mass << '\n';
        return os;
    }

    friend std::istream& operator>>(std::istream& is, Particle& p) {
        is >> p.pos >> p.charge >> p.mass;
        return is;
    }

private :
    vec3d pos;
    // vec3d vel;
    vec3d phi;
    vec3d fld;

    double charge;
    double mass;
   
};

