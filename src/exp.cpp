#include "node.h"

/* getMpoleToExpCoeffs(dirIdx)
 * (M2X) Convert mpole coeffs into outgoing exp coeffs along direction dirIdx
 * dirIdx : direction of outgoing exp coeffs ( 0 = up, 1 = down, ...)
 */
const std::vector<vecXcd> Node::getMpoleToExpCoeffs(const int dirIdx) const {
    std::vector<vecXcd> rotatedCoeffs, expCoeffs;

    // apply rotation
    for (int l = 0; l <= order; ++l)
        rotatedCoeffs.push_back( dirIdx ?
            wignerD[dirIdx+8][l] * coeffs[l] :
            coeffs[l] );

    for (int k = 0; k < orderExp; ++k) {
        auto M_k = tables.quadLengs_[k];
        double l_k = tables.quadCoeffs_[k].first / nodeLeng;
        double coeff_k = tables.quadCoeffs_[k].second / (nodeLeng * M_k);

        vecXcd innerCoeffs = vecXcd::Zero(2*order+1);

        for (int m = -order; m <= order; ++m) {
            int abs_m = abs(m);
            double l_k2l = pow(l_k, abs_m);
            int m_p = m+order;

            for (int l = abs_m; l <= order; ++l) {
                int m_l = m+l;
                innerCoeffs[m_p] +=
                    rotatedCoeffs[l][m_l]
                    * tables.Aexp_[l][m_l] * l_k2l;
                l_k2l *= l_k;
            }

            innerCoeffs[m_p] *= Math::powI(abs_m); // plus sign
        }

        expCoeffs.push_back(vecXcd::Zero(M_k));

        for (int j = 0; j < M_k; ++j) {
            for (int m_p = 0; m_p <= 2*order; ++m_p)
                expCoeffs[k][j] +=
                    innerCoeffs[m_p]
                    * tables.expI_alphas_[k][j][m_p];
        }

        expCoeffs[k] *= coeff_k;
    }

    return expCoeffs;
}

/* addShiftedExpCoeffs(srcExpCoeffs, center0, dirIdx)
 * (X2X) Shift outgoing exp coeffs srcExpCoeffs from direction dirIdx
 * and source node at center0 to center and add to incoming exp coeffs
 * srcExpCoeffs : outgoing exp coeffs
 * center0      : center of source node
 * dirIdx       : direction of outgoing exp coeffs ( 0 = up, 1 = down, ...)
 */
void Node::addShiftedExpCoeffs(
    const std::vector<vecXcd>& srcExpCoeffs, const vec3d& center0, const int dirIdx) 
{
    // rotate dX so this center is in uplist of center0
    const auto& dX = rotMatR[dirIdx] * (center - center0);
    const auto& idX = round(Eigen::Array3d(dX)/nodeLeng);

    const size_t l = idX[0] + 7*idX[1] + 49*idX[2] - 74;
    // assert(0 <= l && l < 98);

    for (size_t k = 0; k < orderExp; ++k)
        for (size_t j = 0; j < tables.quadLengs_[k]; ++j)
            expCoeffs[dirIdx][k][j] +=
                srcExpCoeffs[k][j] * tables.exps_[k][j][l];
}

/* evalExpToLocalCoeffs()
 * (X2L) Convert incoming exp coeffs from all 6 directions into 
 * local coeffs and sum them
 */
void Node::evalExpToLocalCoeffs() {
    assert(!isRoot());

    for (int dirIdx = 0; dirIdx < 6; ++dirIdx) {
        std::vector<vecXcd> rotatedLocalCoeffs;

        for (int k = 0; k < orderExp; ++k) {
            auto M_k = tables.quadLengs_[k];
            double l_k = tables.quadCoeffs_[k].first / nodeLeng;

            vecXcd innerCoeffs = vecXcd::Zero(2*order+1);

            for (int m = -order; m <= order; ++m) {
                int m_p = m+order;

                for (int j = 0; j < M_k; ++j)
                    innerCoeffs[m_p] +=
                        expCoeffs[dirIdx][k][j]
                        * conj(tables.expI_alphas_[k][j][m_p]); // conj

                innerCoeffs[m_p] *= Math::powI(abs(m)); // plus sign
            }

            double ml_k2l = 1.0;
            for (int l = 0; l <= order; ++l) {
                rotatedLocalCoeffs.push_back(vecXcd::Zero(2*l+1));

                for (int m = -l; m <= l; ++m) {
                    int m_l = m+l;
                    rotatedLocalCoeffs[l][m_l] +=
                        innerCoeffs[m+order]
                        * tables.Aexp_[l][m_l]
                        * ml_k2l;
                }

                ml_k2l *= -l_k;
            }
        }

        // apply inverse rotation
        for (int l = 0; l <= order; ++l)
            localCoeffs[l] +=
            (dirIdx ?
                wignerDInv[dirIdx+8][l] * rotatedLocalCoeffs[l] :
                rotatedLocalCoeffs[l]);
    }
}