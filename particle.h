#pragma once

#include "math.h"

using particleVec = std::vector<Particle>;

class Particle {
public :
    const cmplxVec getPsn() const { return psn; }
    const realVec getQs() const { return qs; }

private :
    const realVec ms;
    const realVec qs;
    cmplxVec psn;
    // cmplxVec vel;

};