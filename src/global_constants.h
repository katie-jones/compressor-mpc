#ifndef GLOBAL_CONSTANTS_H
#define GLOBAL_CONSTANTS_H

#include <Eigen/Eigen>

template <size_t N>
using Vec = Eigen::Matrix<double, N, 1>;

const double pi = 3.14159265358979323846;

const double speed_sound = 340;

#endif

