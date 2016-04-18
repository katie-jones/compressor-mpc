#ifndef PARALLEL_COMPRESSORS_H
#define PARALLEL_COMPRESSORS_H

#include <Eigen/Eigen>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/algebra/vector_space_algebra.hpp>
#include <iostream>

#include "global_constants.h"
#include "compressor.h"
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

  void operator()(const SysState &x_in, SysState &dxdt,
                  const double /* t */) const {
    dxdt = GetDerivative(x_in, u);
  }

  ParallelCompressors(SysState x = GetDefaultState(),
                      SysInput u = GetDefaultInput())
      : x(x), u(u) {}

  ParallelCompressors(Comp comps_in[], Tank tank = Tank(),
                      SysState x = GetDefaultState(),
                      SysInput u = GetDefaultInput())
      : x(x), u(u), tank(tank) {
    for (int i = 0; i < n_compressors; i++) {
      comps[i] = comps_in[i];
    }
  }

  ParallelCompressors(Comp comp, Tank tank = Tank(),
                      SysState x = GetDefaultState(),
                      SysInput u = GetDefaultInput())
      : x(x), u(u), tank(tank) {
    for (int i = 0; i < n_compressors; i++) {
      comps[i] = Comp(comp);
    }
  }

  SysState GetDerivative(const SysState x_in, const SysInput u_in) const;

  const static inline SysState GetDefaultState() {
    return ((SysState() << Comp::GetDefaultState().replicate(n_compressors, 1),
             Tank::default_state).finished());
  }

  const static inline SysInput GetDefaultInput() {
    return ((SysInput() << Comp::GetDefaultInput().replicate(n_compressors, 1),
             Tank::default_input).finished());
  }

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
