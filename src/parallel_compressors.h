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
  constexpr static int n_compressors = 2;  ///< Number of compressors in system.

  /// Total number of system states.
  constexpr static int n_states =
      n_compressors * Comp::n_states + Tank::n_states;

  /// Total number of system inputs.
  constexpr static int n_inputs =
      n_compressors * Comp::n_inputs + Tank::n_inputs;

  /// Type describing a valid input to the system.
  typedef Eigen::Array<double, n_inputs, 1> SysInput;

  /// Type describing a valid state of the system.
  typedef Eigen::Array<double, n_states, 1> SysState;

  /// Type of function used as callback for the IntegrateSystem function.
  typedef void (*IntegrationCallbackPtr)(const SysState, const double);

  /// Gives the default state of the system.
  const static inline SysState GetDefaultState() {
    const Comp::CompressorState x =
        ((Comp::CompressorState() << 0.916, 1.145, 0.152, 440, 0).finished());
    return ((SysState() << x.replicate(n_compressors, 1),
             Tank::GetDefaultState()).finished());
  }

  /// Gives the default inputs to the system.
  const static inline SysInput GetDefaultInput() {
    SysInput uout;
    const Comp::CompressorInput u =
        ((Comp::CompressorInput() << 0.304, 0.43, 1.0, 0).finished());

    return ((SysInput() << u.replicate(n_compressors, 1),
             Tank::GetDefaultInput()).finished());
  }

  /**
   * Assign initial state and inputs to system.
   * Use default configurations for the tank and compressors.
   */
  ParallelCompressors(SysState x = GetDefaultState(),
                      SysInput u = GetDefaultInput())
      : x(x), u(u), tank(Tank()) {
    for (int i = 0; i < n_compressors; i++) comps[i] = Comp();
  }

  /**
   * Initialize compressors with pre-defined array.
   * Each compressor can be unique. Values for the tank, state and inputs can
   * also be provided.
   */
  ParallelCompressors(Comp comps_in[], Tank tank = Tank(),
                      SysState x = GetDefaultState(),
                      SysInput u = GetDefaultInput())
      : x(x), u(u), tank(tank) {
    for (int i = 0; i < n_compressors; i++) comps[i] = Comp(comps_in[i]);
  }

  /**
   * Initialize all compressors to be identical.
   * All compressors have the same coefficients and flow characteristics. Values
   * for the tank, state and inputs can also be provided.
   */
  ParallelCompressors(Comp comp, Tank tank = Tank(),
                      SysState x = GetDefaultState(),
                      SysInput u = GetDefaultInput())
      : x(x), u(u), tank(tank) {
    for (int i = 0; i < n_compressors; i++) comps[i] = Comp(comp);
  }

  /**
   * Calculate the derivative of the system at the provided operating point.
   * Use the provided x_in and u_in values instead of the member variables.
   */
  SysState GetDerivative(const SysState x_in, const SysInput u_in) const;

  /**
   * Calculate the derivative of the system.
   * Override the () operator to calculate the derivative using the provided
   * state x_in and the inputs given in the member variable u. Return the
   * derivative in the second argument.
   */
  void operator()(const SysState &x_in, SysState &dxdt,
                  const double /* t */) const {
    dxdt = GetDerivative(x_in, u);
  }

  SysInput u;                 ///< Current input values to the system.
  SysState x;                 ///< Current state of the system.
  Comp comps[n_compressors];  ///< Array of compressors contained in the system.
  Tank tank;                  ///< Discharge tank at output of system.

  /**
   * Integrate system compsys from t0 to tf with timestep dt.
   * Integrate the system compsys from its current state and input values at t0
   * until time tf, using timestep dt. After each timestep the function provided
   * in callback is called where for example the system inputs and other
   * parameters can be changed. The integration method used is a variable-step
   * Dormand-Prince algorithm which reduces the step size until the relative and
   * absolute error are with the bounds specified in max_rel_error and
   * max_abs_error respectively.
   */
  friend void IntegrateSystem(ParallelCompressors compsys, const double t0,
                              const double tf, const double dt,
                              IntegrationCallbackPtr callback,
                              const double max_rel_error = 1e-6,
                              const double max_abs_error = 1e-6) {
    ControlledStepper stepper =
        make_controlled(max_rel_error, max_abs_error, Dopri5Stepper());
    integrate_const(stepper, compsys, compsys.x, t0, tf, dt, callback);
  }

 private:
  constexpr static int n_comp_states = Comp::n_states;
  constexpr static int n_tank_states = Tank::n_states;
  constexpr static int n_comp_inputs = Comp::n_inputs;
  constexpr static int n_tank_inputs = Tank::n_inputs;

  // Type of stepper used to integrate
  typedef boost::numeric::odeint::runge_kutta_dopri5<
      SysState, double, SysState, double,
      boost::numeric::odeint::vector_space_algebra> Dopri5Stepper;

  // Type of stepper used for variable step size integration
  typedef boost::numeric::odeint::controlled_runge_kutta<Dopri5Stepper>
      ControlledStepper;
};

// Define vector_space_norm_inf for the state used in order for odeint to work
namespace boost {
namespace numeric {
namespace odeint {
template <>
struct vector_space_norm_inf<ParallelCompressors::SysState> {
  typedef double result_type;
  double operator()(ParallelCompressors::SysState x) const {
    double absval = 0;

    for (int i = 0; i < ParallelCompressors::n_states; i++)
      absval += x[i] * x[i];
    return sqrt(absval);
  }
};
}
}
}

#endif
