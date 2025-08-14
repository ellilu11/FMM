#include "node.h"

const std::vector<vecXcd> Node::getMpoleToExpCoeffs(const int dirIdx) {
    std::vector<vecXcd> rotatedCoeffs, expCoeffs;

    // apply rotation
    for (int l = 0; l <= order; ++l)
        rotatedCoeffs.push_back(
            dirIdx ?
            wignerD[dirIdx+8][l] * coeffs[l] :
            coeffs[l]
        );

    for (int k = 0; k < orderExp; ++k) {
        auto M_k = tables.quadLengs_[k];
        auto [lmd_k, w_k] = tables.quadCoeffs_[k];
        double l_k = lmd_k / nodeLeng;

        // slow way
        /*expCoeffs.emplace_back(vecXcd::Zero(M_k));
        for (int j = 0; j < M_k; ++j) {
            assert(tables.alphas_[k][j] == 2.0*PI*(j+1)/static_cast<double>(M_k));

            for (int m = -order; m <= order; ++m) {
                for (int l = abs(m); l <= order; ++l) {
                    int m_l = m+l;
                    assert(tables.Aexp_[l][m_l] ==
                        1.0 / std::sqrt(static_cast<double>(factorial(l-m)*factorial(l+m))));

                    expCoeffs[k][j] +=
                        rotatedCoeffs[l][m_l]
                        * tables.Aexp_[l][m_l]
                        * pow(l_k, l)
                        * pow(iu, abs(m))
                        * exp(iu*static_cast<double>(m)*tables.alphas_[k][j]);
                }
            }
            expCoeffs[k][j] *= w_k / nodeLeng / static_cast<double>(M_k);
        }*/

        // fast way
        vecXcd innerCoeffs = vecXcd::Zero(2*order+1);
        for (int m = -order; m <= order; ++m){
            for (int l = abs(m); l <= order; ++l) {
                int m_l = m+l;
                innerCoeffs[m+order] += 
                    rotatedCoeffs[l][m_l] 
                    * tables.Aexp_[l][m_l] 
                    * pow(l_k, l);
            }
            innerCoeffs[m+order] *= pow(iu, abs(m));
        }

        expCoeffs.emplace_back(vecXcd::Zero(M_k));
        for (int j = 0; j < M_k; ++j) {
            for (int m = -order; m <= order; ++m)
                expCoeffs[k][j] +=
                expI(static_cast<double>(m)*tables.alphas_[k][j])
                * innerCoeffs[m+order];
        }
        expCoeffs[k] *= w_k / (nodeLeng * static_cast<double>(M_k));
        
    }
    return expCoeffs;
}

void Node::buildShiftedExpCoeffs(
    const std::vector<vecXcd>& srcExpCoeffs, const vec3d& center0, const int dirIdx) {
    // rotate dX so this center is in uplist of center0
    vec3d dX = rotMatR[dirIdx] * (center - center0); 
    double dx = dX[0], dy = dX[1], dz = dX[2];

    assert(nodeLeng <= dz && dz <= 4.0*nodeLeng);
    assert(sqrt(dx*dx + dy*dy) <= 4.0*sqrt(2.0)*nodeLeng);
    
    for (int k = 0; k < orderExp; ++k) {
        auto M_k = tables.quadLengs_[k];
        auto [lmd_k, w_k] = tables.quadCoeffs_[k];
        double l_k = lmd_k / nodeLeng;

        for (int j = 0; j < M_k; ++j) {
            double a_kj = tables.alphas_[k][j];
            expCoeffs[dirIdx][k][j] += srcExpCoeffs[k][j]
                * exp(-l_k * dz + iu*l_k * (dx*cos(a_kj) + dy*sin(a_kj)));
        }
    }
}

void Node::buildLocalCoeffsFromDirList() {
    assert(!isRoot());

    // move to constructor later
    for (int l = 0; l <= order; ++l)
        localCoeffs.emplace_back(vecXcd::Zero(2*l+1));

    for (int dirIdx = 0; dirIdx < 6; ++dirIdx) {
        std::vector<vecXcd> rotatedLocalCoeffs;

        for (int k = 0; k < orderExp; ++k) {
            auto [lmd_k, w_k] = tables.quadCoeffs_[k];
            auto M_k = tables.quadLengs_[k];
            double l_k = lmd_k / nodeLeng;

            // slow way
            //for (int l = 0; l <= order; ++l) {
            //    rotatedLocalCoeffs.emplace_back(vecXcd::Zero(2*l+1));
            //    for (int m = -l; m <= l; ++m) {
            //        int m_l = m+l;
            //        for (int j = 0; j < M_k; ++j)
            //            rotatedLocalCoeffs[l][m_l] +=
            //                tables.Aexp_[l][m_l]
            //                * pow(-l_k, l)
            //                * expCoeffs[dirIdx][k][j]
            //                * expI(-m*tables.alphas_[k][j])
            //                * pow(iu, abs(m))
            //                ;
            //    }
            //}

            // fast way
            vecXcd innerCoeffs = vecXcd::Zero(2*order+1);
            for (int m = -order; m <= order; ++m) {
                // int m_p = m+order;
                for (int j = 0; j < M_k; ++j)
                    innerCoeffs[m+order] +=
                        expCoeffs[dirIdx][k][j]
                        * expI(-m*tables.alphas_[k][j]);
                innerCoeffs[m+order] *= pow(iu, abs(m));
            }

            for (int l = 0; l <= order; ++l) {
                rotatedLocalCoeffs.emplace_back(vecXcd::Zero(2*l+1));
                for (int m = -l; m <= l; ++m) {
                    int m_l = m+l;
                    rotatedLocalCoeffs[l][m_l] +=
                        tables.Aexp_[l][m_l] 
                        * pow(-l_k, l)
                        * innerCoeffs[m+order];
                }
            }
        }

        // apply inverse rotation
        for (int l = 0; l <= order; ++l)
            localCoeffs[l] +=
                ( dirIdx ?
                  wignerDInv[dirIdx+8][l] * rotatedLocalCoeffs[l] :
                  rotatedLocalCoeffs[l] );
    }

    // comment out for mpole2local test
    //if (!base->isRoot())
    //    for (int l = 0; l <= order; ++l)
    //        localCoeffs[l] += (base->getShiftedLocalCoeffs(branchIdx))[l];


    //for (int l = 0; l <= order; ++l)
    //    std::cout << localCoeffs[l].transpose() << '\n';
}