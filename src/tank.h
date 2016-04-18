#ifndef TANK_H
#define TANK_H

#include <Eigen/Eigen>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/algebra/vector_space_algebra.hpp>
#include <iostream>

#include "global_constants.h"

/*
 * Class containing compressor tank
 */
class Tank {
 public:
  typedef double TankInput;
  typedef double TankState;

  struct Params {
    double pout;
    double volume;
    Vec<8> D;
    double m_out_c;
    Params();
    Params(const Params &x);
  };

  constexpr static int n_states = 1;
  constexpr static int n_inputs = 1;
  constexpr static TankInput default_input = 0.7;
  constexpr static TankState default_state = 1.12;

  Tank(Params params = Params()) : params(params) {}

  TankState GetDerivative(const TankState x, const TankInput u,
                          const double mass_flow_compressors) const;

  const Params params;

};

#endif
