#pragma once

#include "enum.h"
#include "math.h"

struct Tables {
    Tables() = default;
    Tables(const int order, const Precision prec) {
        buildYlmTables(order);
        buildQuadTables(prec);
        buildExpTables(order);
    }

    void buildYlmTables(const int);
    void buildQuadTables(const Precision);
    void buildExpTables(const int);

    // Ylm tables
    std::vector<realVec> coeffYlm_;
    std::vector<realVec> fallingFact_;
    std::vector<realVec> legendreSum_;
    std::vector<realVec> fracCoeffYlm_;
    std::vector<realVec> A_;
    std::vector<realVec> Aexp_;

    // quad tables
    std::vector<pair2d> quadCoeffs_;
    std::vector<int> quadLengs_;

    // exp tables
    std::vector<realVec> alphas_;
    std::vector<std::vector<cmplxVec>> expI_alphas_;
    std::vector<std::vector<std::array<cmplx,25>>> expsInner_;
    std::vector<std::vector<std::array<cmplx,36>>> expsOuter_;
    std::vector<std::vector<std::array<cmplx,98>>> exps_;
    std::vector<std::vector<std::array<cmplx,8>>> expsMerge_;
};

void Tables::buildYlmTables(const int order) {
    auto binom = [](double x, int k) {
        return fallingFactorial(x, k) / factorial(k);
        };

    for (int l = 0; l <= 2*order; ++l) {
        realVec coeffYlm_l, fallingFact_l, legendreSum_l, fracCoeffYlm_l, A_l, Aexp_l;

        for (int m = 0; m <= l; ++m) {
            coeffYlm_l.push_back(coeffYlm(l, m));
            fallingFact_l.push_back(fallingFactorial(l, m));
            legendreSum_l.push_back(binom(l, m) * binom((l+m-1)/2.0, l));
            fracCoeffYlm_l.push_back(sqrt((l-m)/static_cast<double>(l+m)));
        }

        auto pm_l = pm(l);
        for (int m = -l; m <= l; ++m) {
            A_l.push_back(pm_l /
                std::sqrt(static_cast<double>(factorial(l-m)*factorial(l+m))));
            Aexp_l.push_back(1.0 /
                std::sqrt(static_cast<double>(factorial(l-m)*factorial(l+m))));
        }

        coeffYlm_.push_back(coeffYlm_l);
        fallingFact_.push_back(fallingFact_l);
        legendreSum_.push_back(legendreSum_l);
        fracCoeffYlm_.push_back(fracCoeffYlm_l);
        A_.push_back(A_l);
        Aexp_.push_back(Aexp_l);
    }
}

void Tables::buildQuadTables(const Precision prec) {
    switch (prec) {
        case Precision::LOW:
            quadCoeffs_ = {
                {0.10934746769000, 0.27107502662774},
                {0.51769741015341, 0.52769158843946},
                {1.13306591611192, 0.69151504413879},
                {1.88135015110740, 0.79834400406452},
                {2.71785409601205, 0.87164160121354},
                {3.61650274907449, 0.92643839116924},
                {4.56271053303821, 0.97294622259483},
                {5.54900885348528, 1.02413865844686}
            };
            quadLengs_ = { 4,8,16,16,24,24,8,4 };
            break;

        case Precision::MEDIUM:
            quadCoeffs_ = {
                {0.05599002531749, 0.14239483712194},
                {0.28485138101968, 0.31017671029271},
                {0.66535367065853, 0.44557516683709},
                {1.16667904805296, 0.55303383994159},
                {1.76443027413431, 0.63944903363523},
                {2.44029832236380, 0.70997911214019},
                {3.18032180991515, 0.76828253949732},
                {3.97371715777193, 0.81713201141707},
                {4.81216799410634, 0.85872191623337},
                {5.68932314511487, 0.89480789582390},
                {6.60040479444377, 0.92680189417317},
                {7.54190497469911, 0.95586282708096},
                {8.51136569298099, 0.98299145008230},
                {9.50723242759128, 1.00913395385703},
                {10.52874809650967, 1.03531774600508},
                {11.57587019602884, 1.06318427913963},
                {12.65078163968520, 1.10232109521088}
            };
            quadLengs_ = { 8,8,16,16,24,32,32,32,48,48,48,48,48,48,48,8,4 };
            break;

        case Precision::HIGH:
            quadCoeffs_ = {
                {0.03705701953816, 0.09473396337900},
                {0.19219683859955, 0.21384206006426},
                {0.46045971214897, 0.32031528543989},
                {0.82805130101422, 0.41254929390710},
                {1.28121229944787, 0.49176691815621},
                {1.80792019276297, 0.55998309037174},
                {2.39814728074333, 0.61909314036708},
                {3.04359012306582, 0.67064351982741},
                {3.73732742924096, 0.71586567032066},
                {4.47354768940212, 0.75576118553096},
                {5.24735518169467, 0.79116885492295},
                {6.05462948620944, 0.82280556212477},
                {6.89191648795972, 0.85129012269433},
                {7.75633860708838, 0.87715909928110},
                {8.64551915195994, 0.90087981520398},
                {9.55751929613924, 0.92286282936149},
                {10.49078760616705, 0.94347471535979},
                {11.44412262341269, 0.96305166489156},
                {12.41664955395045, 0.98191478773737},
                {13.40781311788324, 1.00038891281291},
                {14.41739038894472, 1.01882849188686},
                {15.44553016867884, 1.03765781507554},
                {16.49282861241170, 1.05744113465683},
                {17.56045648926099, 1.07903824697122},
                {18.65046484106274, 1.10434337868208},
                {19.76847686619416, 1.14488166506896}
            };
            quadLengs_ = { 8,16,16,16,24,32,32,32,48,48,48,64,64,
                      64,64,72,72,80,80,88,88,88,88,72,32,4 };
            break;
    }

    assert(quadCoeffs_.size() == quadLengs_.size());
}

void Tables::buildExpTables(const int order) {
    for (int k = 0; k < quadCoeffs_.size(); ++k) {
        double M_k = quadLengs_[k];
        realVec alphas_k; //
        std::vector<cmplxVec> expI_alphas_k;
        std::vector<std::array<cmplx,25>> expsInner_k;
        std::vector<std::array<cmplx,36>> expsOuter_k;
        std::vector<std::array<cmplx,98>> exps_k;
        std::vector<std::array<cmplx,8>> expsMerge_k;

        for (int j = 0; j < M_k; ++j) {
            double alpha_kj = 2.0 * PI * (j+1) / static_cast<double>(M_k);
            alphas_k.push_back(alpha_kj); //

            cmplxVec expI_alphas_kj;
            for (int m = -order; m <= order; ++m)
                expI_alphas_kj.push_back(expI(m*alpha_kj));
            expI_alphas_k.push_back(expI_alphas_kj);

            std::array<cmplx, 98> exps_kj;
            size_t l = 0;
            for (int dz = 2; dz <= 3; ++dz)
                for (int dy = -3; dy <= 3; ++dy)
                    for (int dx = -3; dx <= 3; ++dx) {
                        exps_kj[l++] =
                            exp(quadCoeffs_[k].first
                                * cmplx(-1.0*dz,
                                    dx*cos(alpha_kj) + dy*sin(alpha_kj)));
                    }
            assert(l == 98);
            exps_k.push_back(exps_kj);

            std::array<cmplx,25> expsInner_kj;
            size_t m = 0;
            constexpr double dzInner = 2.0;
            for (int dy = -2; dy <= 2; dy++)
                for (int dx = -2; dx <= 2; dx++) {
                    expsInner_kj[m++] =
                        exp(quadCoeffs_[k].first
                            * cmplx(-dzInner,
                                dx*cos(alpha_kj) + dy*sin(alpha_kj)));
                }
            assert(m == 25);
            expsInner_k.push_back(expsInner_kj);

            std::array<cmplx,36> expsOuter_kj;
            size_t n = 0;
            constexpr double dzOuter = 2.5;
            for (double dy = -2.5; dy <= 2.5; dy += 1.0)
                for (double dx = -2.5; dx <= 2.5; dx += 1.0) {
                    expsOuter_kj[n++] =
                        exp(quadCoeffs_[k].first
                            * cmplx(-dzOuter,
                                dx*cos(alpha_kj) + dy*sin(alpha_kj)));
                }
            assert(n == 36);
            expsOuter_k.push_back(expsOuter_kj);

            std::array<cmplx, 8> expsMerge_kj;
            for (int dir = 0; dir < 8; ++dir){
                auto dX = -idx2pm(dir);
                expsMerge_kj[dir] = // 1.0;
                    exp(quadCoeffs_[k].first / 4.0 // 4.0
                        * cmplx(-1.0*dX[2],
                            dX[0]*cos(alpha_kj) + dX[1]*sin(alpha_kj)));
            }
            expsMerge_k.push_back(expsMerge_kj);
        }
        alphas_.push_back(alphas_k); //
        expI_alphas_.push_back(expI_alphas_k);
        expsInner_.push_back(expsInner_k);
        expsOuter_.push_back(expsOuter_k);
        exps_.push_back(exps_k);
        expsMerge_.push_back(expsMerge_k);
    }
}