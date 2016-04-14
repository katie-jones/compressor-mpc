#ifndef COMPRESSOR_DEFS_H
#define COMPRESSOR_DEFS_H

#include <cmath>
#include <Eigen/Eigen>
#include <boost/numeric/odeint/algebra/vector_space_algebra.hpp>

const int n_states = 5;
const int n_inputs = 5;

template <size_t N>
using Vec = Eigen::Matrix<double, N, 1>;

typedef Eigen::Array<double, n_states, 1> comp_state;
typedef Eigen::Array<double, n_inputs, 1> comp_input;

const double pi = 3.14159265358979323846;

// Define norm of Eigen::Array
namespace boost {
namespace numeric {
namespace odeint {
template <>
struct vector_space_norm_inf<comp_state> {
  typedef double result_type;
  double operator()(comp_state x) const {
    double absval = 0;
    for (int i = 0; i < n_states; i++) absval += x[i] * x[i];
    return absval;
  }
};
}
}
}

#endif
