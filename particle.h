#pragma once

#include <random>
#include "config.h"
#include "math.h"

class Particle;

using ParticleVec = std::vector<std::shared_ptr<Particle>>;

class Particle {
public :
    Particle()
        : phi(0.0), fld(zeroVec) {};
    Particle(const vec3d X, const double q, const double m)
        : pos(X), charge(q), mass(m), phi(0.0), fld(zeroVec) {};

    vec3d getPos() const { return pos; }

    double getCharge() const { return charge; }

    double getMass() const { return mass; }

    const double getPhi() const { return phi; }

    //void addToPhi(double phi_) { phi += phi_; }

    //void addToFld(vec3d fld_) { fld += fld_; }

    void addFromSol(const pairSol& sol) { 
        phi += sol.first;
        fld += sol.second;
    }

    void printPhi(std::ofstream& f) const {
        // f << phi << ' ';
        f << phi << '\n';
    }

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
    double phi;
    vec3d fld;

    double charge;
    double mass;
   
};

