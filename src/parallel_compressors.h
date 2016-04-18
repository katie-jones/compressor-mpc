#ifndef PARALLEL_COMPRESSORS_H
#define PARALLEL_COMPRESSORS_H

#include <Eigen/Eigen>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/algebra/vector_space_algebra.hpp>
#include <iostream>

#include "global_constants.h"
#include "child_compressor.h"
#include "tank.h"

/*
 * Class containing system of parallel compressors
 * Overloaded () operator to give state derivative
 */
class ParallelCompressors {
 public:
  constexpr static int n_compressors = 2;
  constexpr static int n_states =
      n_compressors * Comp::n_states + Tank::n_states;
  constexpr static int n_inputs =
      n_compressors * Comp::n_inputs + Tank::n_inputs;

  typedef Vec<n_inputs> SysInput;
  typedef Vec<n_states> SysState;

  const static SysInput default_input;
  const static SysState default_state;

  void operator()(const SysState &x_in, SysState &dxdt,
                  const double /* t */) const {
    dxdt = GetDerivative(x_in, u);
  }

  ParallelCompressors(SysState x = default_state, SysInput u = default_input)
      : x(x), u(u) {}

  SysState GetDerivative(const SysState x, const SysInput u) const;

  SysInput u;
  SysState x;
  Comp comps[n_compressors];
  Tank tank;

 private:
  constexpr static int n_comp_states = Comp::n_states;
  constexpr static int n_tank_states = Tank::n_states;
  constexpr static int n_comp_inputs = Comp::n_inputs;
  constexpr static int n_tank_inputs = Tank::n_inputs;

  typedef boost::numeric::odeint::runge_kutta_dopri5<
      SysState, double, SysState, double,
      boost::numeric::odeint::vector_space_algebra> Dopri5Stepper;
  typedef boost::numeric::odeint::controlled_runge_kutta<Dopri5Stepper>
      ControlledStepper;
};

#endif

