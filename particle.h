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

    const cmplx getPhi() const { return phi; }
    void setPhi(cmplx phi_) { phi = phi_; }
    void printPhi(std::ofstream& f) const {
        f << phi.real() << ' ';
    }

    void setFld(cmplx fld_) { fld = fld_; }
    void printFld(std::ofstream& f) const {
        f << fld.real() << ' ' << fld.imag() << ' ';
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
    cmplx pos;
    // cmplx vel;
    cmplx phi;
    cmplx fld;

    double charge;
    double mass;
   
};

