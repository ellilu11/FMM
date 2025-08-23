#include "node.h"

std::vector<vecXcd> Node::getMpoleToExpCoeffs(const int dirIdx) {
    std::vector<vecXcd> rotatedCoeffs, expCoeffs;

    // apply rotation
    for (int l = 0; l <= order; ++l)
        rotatedCoeffs.push_back(dirIdx ?
            wignerD[dirIdx+8][l] * coeffs[l] :
            coeffs[l]);

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

            innerCoeffs[m_p] *= powI(abs_m); // plus sign
        }

        expCoeffs.emplace_back(vecXcd::Zero(M_k));

        for (int j = 0; j < M_k; ++j) {
            for (int m_p = 0; m_p <= 2*order; ++m_p)
                expCoeffs[k][j] +=
                innerCoeffs[m_p]
                * tables.expI_alphas_[k][j][m_p];
        }

        expCoeffs[k] *= coeff_k;
    }

    expCoeffsOut[dirIdx] = expCoeffs;

    return expCoeffs;
}

/*const std::vector<vecXcd> Node::getMergedExpCoeffs(const int dirIdx) const {
    std::vector<vecXcd> mergedCoeffs;
    for (int k = 0; k < orderExp; ++k)
        mergedCoeffs.emplace_back(vecXcd::Zero(tables.quadLengs_[k]));

    for (const auto& branch : branches) {
        auto expCoeffs = branch->getMpoleToExpCoeffs(dirIdx);
        const auto dX = rotMatR[dirIdx] * (center - branch->getCenter());
        // const auto dX = center - branch->getCenter();
        const double dx = dX[0], dy = dX[1], dz = dX[2];

        for (int k = 0; k < orderExp; ++k) {
            const double l_k = tables.quadCoeffs_[k].first / nodeLeng;
            for (int j = 0; j < tables.quadLengs_[k]; ++j) {
                const double a_kj = tables.alphas_[k][j];
                mergedCoeffs[k][j] +=
                    expCoeffs[k][j]
                    * exp(l_k / 4.0
                        * cmplx(-1.0*dz,
                            dx*cos(a_kj) + dy*sin(a_kj)));
            }
        }
    }
    return mergedCoeffs;
}*/

/*void Node::addShiftedExpCoeffs(
    const std::vector<vecXcd>& srcExpCoeffs, const vec3d& center0, const int dirIdx) {
    // rotate dX so this center is in uplist of center0
    const auto dX = rotMatR[dirIdx] * (center - center0);
    const double dx = dX[0], dy = dX[1], dz = dX[2];
    // std::cout << dx / nodeLeng << ' ' << dy / nodeLeng << ' ' << dz / nodeLeng << '\n';

    for (int k = 0; k < orderExp; ++k) {
        const double l_k = tables.quadCoeffs_[k].first / nodeLeng;
        for (int j = 0; j < tables.quadLengs_[k]; ++j) {
            const double a_kj = tables.alphas_[k][j];
            expCoeffs[dirIdx][k][j] +=
                srcExpCoeffs[k][j]
                * exp(-l_k * dz + iu*l_k * (dx*cos(a_kj) + dy*sin(a_kj)));
        }
    }
}*/

const std::vector<vecXcd> Node::getMergedExpCoeffs(const int dirIdx) const {
    std::vector<vecXcd> mergedCoeffs;
    for (int k = 0; k < orderExp; ++k)
        mergedCoeffs.emplace_back(vecXcd::Zero(tables.quadLengs_[k]));

    for (const auto& branch : branches) {
        // Construct exp coeffs of each branch here
        auto expCoeffs = branch->getMpoleToExpCoeffs(dirIdx);
        auto branchIdx = branch->getBranchIdx();

        // Shift branch exp coeffs to center of this node and merge them
        for (int k = 0; k < orderExp; ++k)
            for (int j = 0; j < tables.quadLengs_[k]; ++j)
                mergedCoeffs[k][j] +=
                    expCoeffs[k][j] * tables.expsMerge_[k][j][branchIdx];
    }
    return mergedCoeffs;
}

void Node::addShiftedExpCoeffsFromBranch(
    const std::vector<vecXcd>& srcExpCoeffs, const vec3d& center0, const int dirIdx) 
{
    // rotate dX so this center is in uplist of center0
    const auto dX = rotMatR[dirIdx] * (center - center0);
    const auto idX = round(Eigen::Array3d(dX)/nodeLeng);
    const size_t l = idX[0] + 7*idX[1] + 49*idX[2] - 74;

    assert(0 <= l && l < 98);
    //assert(idz == 2 || idz == 3);
    //assert(sqrt(idx*idx + idy*idy) <= 4.0*sqrt(2.0));

    // shift to dirlist: [idX, idY] \in {-3, -2, -1, 0, 1, 2, 3}, idZ \in {2, 3} 
    for (size_t k = 0; k < orderExp; ++k)
        for (size_t j = 0; j < tables.quadLengs_[k]; ++j)
            expCoeffs[dirIdx][k][j] +=
                srcExpCoeffs[k][j] * tables.exps_[k][j][l];
}

void Node::addShiftedExpCoeffs(
    const std::vector<vecXcd>& srcExpCoeffs, const vec3d& center0, const int dirIdx) 
{
    // rotate dX so this center is in uplist of center0 (center of source node)
    const auto idX = rotMatR[dirIdx] * (center - center0) / nodeLeng;
    const size_t idz = round(2.0*idX[2]); 
    // assert(idz == 4 || idz == 5 );

    size_t l;

    // shift to inner ilist: [idX, idY] \in {-2, -1, 0, 1, 2}, idZ = 2 
    if (idz == 4) {
        l = round(idX[0]+2.0) + 5*round(idX[1]+2.0);
        assert(0 <= l && l < 25);

        for (size_t k = 0; k < orderExp; ++k)
            for (size_t j = 0; j < tables.quadLengs_[k]; ++j)
                expCoeffs[dirIdx][k][j] +=
                    srcExpCoeffs[k][j] * tables.expsInner_[k][j][l];

    // shift to outer ilist: [idX, idY] \in {-2.5, -1.5, -0.5, 0.5, 1.5, 2.5}, idZ = 2.5 
    } else {
        l = round(idX[0]+2.5) + 6*round(idX[1]+2.5);
        assert(0 <= l && l < 36);

        for (size_t k = 0; k < orderExp; ++k)
            for (size_t j = 0; j < tables.quadLengs_[k]; ++j)
                expCoeffs[dirIdx][k][j] +=
                srcExpCoeffs[k][j] * tables.expsOuter_[k][j][l];
    }
}

void Node::evalLocalCoeffsFromDirList() {
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

                innerCoeffs[m_p] *= powI(abs(m)); // plus sign
            }

            double ml_k2l = 1.0;
            for (int l = 0; l <= order; ++l) {
                rotatedLocalCoeffs.emplace_back(vecXcd::Zero(2*l+1));

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