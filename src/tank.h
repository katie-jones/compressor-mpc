#ifndef TANK_H
#define TANK_H

#include <Eigen/Eigen>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/algebra/vector_space_algebra.hpp>

#include "global_constants.h"

/*
 * Class containing compressor tank
 */
class Tank {
 public:
  constexpr static int n_states = 1;
  constexpr static int n_inputs = 1;

  typedef Eigen::Array<double, n_inputs, 1> TankInput;
  typedef Eigen::Array<double, n_states, 1> TankState;

  struct Params {
    double pout;
    double volume;
    Vec<8> D;
    double m_out_c;
    Params();
    Params(const Params &x);
  };

  Tank(Params params) : params(params) {}
  Tank() : params(Params()) {}

  TankState GetDerivative(const TankState x, const TankInput u,
                          const double mass_flow_compressors) const;
  static inline TankInput GetDefaultInput() {
    return ((TankInput() << 0.7).finished());
  }
  static inline TankState GetDefaultState() {
    return ((TankState() << 1.12).finished());
  }

  const Params params;
};

#endif
