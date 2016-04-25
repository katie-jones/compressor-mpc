#ifndef PARALLEL_COMPRESSORS_H
#define PARALLEL_COMPRESSORS_H

#include <Eigen/Eigen>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/algebra/vector_space_algebra.hpp>
#include <iostream>

#include "global_constants.h"
#include "compressor.h"
#include "tank.h"
#include "global_constants.h"

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

  typedef Eigen::Array<double,n_inputs,1> SysInput;
  typedef Eigen::Array<double,n_states,1> SysState;

  typedef void (*IntegrationCallbackPtr)(const SysState, const double);

  void operator()(const SysState &x_in, SysState &dxdt,
                  const double /* t */) const {
    dxdt = GetDerivative(x_in, u);
    // dxdt = 0*x_in;
  }

  const static inline SysState GetDefaultState() {
    return ((SysState() << Comp::GetDefaultState().replicate(n_compressors, 1),
             Tank::GetDefaultState()).finished());
  }

  const static inline SysInput GetDefaultInput() {
    SysInput uout;
    uout << Comp::GetDefaultInput().replicate(n_compressors, 1),Tank::GetDefaultInput();

    // return ((SysInput() << Comp::GetDefaultInput().replicate(n_compressors, 1),
             // Tank::GetDefaultInput()).finished());
    // return SysInput();
    return uout;
    
  }

  // ParallelCompressors(SysState x = GetDefaultState(),
                      // SysInput u = GetDefaultInput())
      // : x(x), u(u) {}
  ParallelCompressors() : x(GetDefaultState()), u(GetDefaultInput()) {}

  ParallelCompressors(Comp comps_in[], Tank tank = Tank(),
                      SysState x = GetDefaultState(),
                      SysInput u = GetDefaultInput())
      : x(x), u(u), tank(tank) {
    for (int i = 0; i < n_compressors; i++) comps[i] = Comp(comps_in[i]);
  }

  ParallelCompressors(Comp comp, Tank tank = Tank(),
                      SysState x = GetDefaultState(),
                      SysInput u = GetDefaultInput())
      : x(x), u(u), tank(tank) {
    for (int i = 0; i < n_compressors; i++) comps[i] = Comp(comp);
  }

  SysState GetDerivative(const SysState x_in, const SysInput u_in) const;


  SysInput u;
  SysState x;
  Comp comps[n_compressors];
  Tank tank;

  friend void IntegrateSystem(ParallelCompressors compsys, const double t0,
                                  const double tf, const double dt,
                                  IntegrationCallbackPtr callback,
                                  const double rel_error = 1e-6,
                                  const double abs_error = 1e-6) {
    ControlledStepper stepper =
        make_controlled(rel_error, abs_error, Dopri5Stepper());
    integrate_const(stepper, compsys, compsys.x, t0, tf, dt, callback);
  }

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

namespace boost {
namespace numeric {
namespace odeint {
template <>
struct vector_space_norm_inf<ParallelCompressors::SysState> {
  typedef double result_type;
  double operator()(ParallelCompressors::SysState x) const {
    double absval = 0;
    
    for (int i = 0; i < ParallelCompressors::n_states; i++) absval += x[i] * x[i];
    return sqrt(absval);
  }
};
}
}
}

#endif
