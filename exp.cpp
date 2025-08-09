#include "node.h"

const vecXcdVec Node::getMpoleToExpCoeffsFromNode() {
    vecXcdVec expCoeffs;

    for (int k = 0; k < orderExp; ++k) {
        auto M_k = tables.quadLengs_[k];
        auto [w_k, l_k] = tables.quadCoeffs_[k];
        double coeff_k = w_k / (nodeLeng * static_cast<double>(M_k));

        expCoeffs.emplace_back(vecXcd::Zero(M_k-1));
        for (int j = 1; j <= M_k; ++j) {
            int j_ = j - 1;
            double a_jk = 2.0*PI*j/M_k;
            
            for (int m = -order; m <= order; ++m) {
                int abs_m = abs(m);
                cmplx coeff_m = coeff_k * pow(-iu, abs_m) * expI(static_cast<double>(m)*a_jk);

                for (int n = abs_m; n <= order; ++n) {
                    int m_ = m+n;
                    expCoeffs[k][j_] += coeff_m *
                        coeffs[n][m_] * tables.A_[n][m_] / pm(n) * pow(l_k/nodeLeng, n);

                }
            }
        }
    }
    return expCoeffs;
}