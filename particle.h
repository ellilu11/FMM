#pragma once

#include <random>
#include "config.h"
#include "math.h"

class Particle;

using ParticleVec = std::vector<std::shared_ptr<Particle>>;

class Particle {
public :
    Particle() = default;
    Particle(const vec3d X, const vec3d R, const double q, const double m)
        : pos(X), sph(R), charge(q), mass(m)
    {
    };

    vec3d getPos() const { return pos; }
    vec3d getSph() const { return sph; }
    double getCharge() const { return charge; }
    double getMass() const { return mass; }

    const double getPhi() const { return phi; }
    void setPhi(double phi_) { phi = phi_; }
    void printPhi(std::ofstream& f) const {
        f << phi << '\n';
    }

    void setFld(vec3d fld_) { fld = fld_; }
    void printFld(std::ofstream& f) const {
        f << fld << '\n';
    }

    friend std::ostream& operator<<(std::ostream& os, Particle& p) {
        os << p.pos << " " << " " << p.charge << " " << p.mass << '\n';
        return os;
    }

    friend std::istream& operator>>(std::istream& is, Particle& p) {
        is >> p.pos >> p.charge >> p.mass;
        return is;
    }

private :
    vec3d pos;
    vec3d sph;
    // vec3d vel;
    double phi;
    vec3d fld;

    double charge;
    double mass;
   
};

