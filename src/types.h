#pragma once

#include <Eigen/Dense>

using realVec = std::vector<double>;
using cmplx = std::complex<double>;
using cmplxVec = std::vector<cmplx>;

using pair2i = std::pair<int, int>;
using pair2d = std::pair<double, double>;

using vec3d = Eigen::Vector3d;
using vec3cd = Eigen::Vector3cd;
using vecXcd = Eigen::VectorXcd;

using mat3d = Eigen::Matrix3d;
using matXcd = Eigen::MatrixXcd;

constexpr cmplx iu(0, 1);
const double PI = std::acos(-1.0);
const vec3d zeroVec = vec3d::Zero();