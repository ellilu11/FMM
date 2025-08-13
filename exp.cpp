#include "node.h"

const std::vector<vecXcd> Node::getMpoleToExpCoeffs(const int dirIdx) {
    std::vector<vecXcd> rotatedCoeffs, expCoeffs;

    // apply rotation
    for (int l = 0; l <= order; ++l) {
        if (dirIdx == 0)
            rotatedCoeffs.push_back(coeffs[l]);
        else 
            rotatedCoeffs.push_back(wignerD[dirIdx+8][l] * coeffs[l]);
    }

    for (int k = 0; k < orderExp; ++k) {
        auto M_k = tables.quadLengs_[k];
        auto [l_k, w_k] = tables.quadCoeffs_[k];
        vecXcd innerCoeffs = vecXcd::Zero(2*order+1);

        for (int m = -order; m <= order; ++m){
            for (int l = abs(m); l <= order; ++l) {
                int m_l = m+l;
                innerCoeffs[m+order] += 
                    rotatedCoeffs[l][m_l]
                    // / std::sqrt(static_cast<double>(factorial(l-m)*factorial(l+m)))
                    * tables.Aexp_[l][m_l] 
                    * pow(l_k/nodeLeng, l);
            }
            // innerCoeffs[m_p] *= pow(-iu, abs(m));
        }

        expCoeffs.emplace_back(vecXcd::Zero(M_k));
        for (int j = 0; j < M_k; ++j) {
            for (int m = -order; m <= order; ++m)
                expCoeffs[k][j] +=
                pow(-iu, abs(m)) *
                expI(static_cast<double>(m)*tables.alphas_[k][j])
                * innerCoeffs[m+order];
            expCoeffs[k][j] *= w_k / (nodeLeng * static_cast<double>(M_k));
        }

        // expCoeffs[k] *= w_k / (nodeLeng * static_cast<double>(M_k));
    }
    return expCoeffs;
}

void Node::buildShiftedExpCoeffs(
    const std::vector<vecXcd>& srcExpCoeffs, const vec3d& center0, const int dirIdx) {

    vec3d dX = rotMatR[dirIdx] * (center - center0); // rotate dX so this center is in uplist of center0
    double dx = dX[0], dy = dX[1], dz = dX[2];
    //if (dirIdx == 3)
    //    std::cout << nodeLeng0 << ' ' << nodeLeng << ' ' 
    //              << dx / nodeLeng0 << ' ' << dy / nodeLeng0 << ' ' << dz / nodeLeng0 << '\n';

    // assert(nodeLeng == nodeLeng0);
    assert(nodeLeng <= dz && dz <= 4.0*nodeLeng);
    assert(sqrt(dx*dx + dy*dy) <= 4.0*sqrt(2.0)*nodeLeng);
    
    for (int k = 0; k < orderExp; ++k) {
        auto M_k = tables.quadLengs_[k];
        auto [l_k, w_k] = tables.quadCoeffs_[k];

        for (int j = 0; j < M_k; ++j) {
            double a_kj = tables.alphas_[k][j];
            expCoeffs[dirIdx][k][j] +=
                srcExpCoeffs[k][j]
                * exp( -l_k / nodeLeng * dz )
                * expI( l_k / nodeLeng * (dx*cos(a_kj) + dy*sin(a_kj)) );
        }
    }
}

// exponential expansions
void Node::buildLocalCoeffsFromDirList() {
    assert(!isRoot());

    for (int dirIdx = 0; dirIdx < 6; ++dirIdx) {
        std::vector<vecXcd> rotatedLocalCoeffs;

        for (int k = 0; k < orderExp; ++k) {
            auto [l_k, w_k] = tables.quadCoeffs_[k];
            auto M_k = tables.quadLengs_[k];
            vecXcd innerCoeffs = vecXcd::Zero(2*order+1);

            for (int m = -order; m <= order; ++m) {
                // int m_p = m+order;
                for (int j = 0; j < M_k; ++j) {
                    innerCoeffs[m+order] +=
                        expCoeffs[dirIdx][k][j]
                        * expI(m*tables.alphas_[k][j]);
                }
                innerCoeffs[m+order] *= pow(-iu, abs(m));
            }

            for (int l = 0; l <= order; ++l) {
                rotatedLocalCoeffs.emplace_back(vecXcd::Zero(2*l+1));
                for (int m = -l; m <= l; ++m) {
                    rotatedLocalCoeffs[l][m+l] +=
                        tables.Aexp_[l][m+l] 
                        * pow(-l_k / nodeLeng, l)
                        * innerCoeffs[m+l];
                }
            }
        }

        // apply inverse rotation
        for (int l = 0; l <= order; ++l)
            localCoeffs[l] += wignerDInv[dirIdx+8][l] * rotatedLocalCoeffs[l];
    }

    // comment out for mpole2local test
    //if (!base->isRoot())
    //    for (int l = 0; l <= order; ++l)
    //        localCoeffs[l] += (base->getShiftedLocalCoeffs(branchIdx))[l];


    //for (int l = 0; l <= order; ++l)
    //    std::cout << localCoeffs[l].transpose() << '\n';
}