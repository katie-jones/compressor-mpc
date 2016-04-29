#ifndef COMPRESSOR_SIMULATION_H
#define COMPRESSOR_SIMULATION_H

#include <Eigen/Eigen>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/algebra/vector_space_algebra.hpp>
#include <iostream>

#include "global_constants.h"
#include "compressor.h"

/**
 * Child class of compressor implementing methods required for simulations.
 * Contains member variables x (state) and u (inputs), overloaded () operator to
 * give state derivative. Also defines a friend function to integrate the system
 * over a given time range.
 */
class CompressorSimulation : public Comp {
 public:
   /// Function pointer used for callback when integrating system.
  typedef void (*IntegrationCallbackPtr)(const CompressorState, const double);

  /// Default outlet pressure.
  constexpr static double default_pout = 1.0;

  /// Provide initial member variable values or use defaults.
  CompressorSimulation(CompressorState x = GetDefaultState(),
                       CompressorInput u = GetDefaultInput(),
                       double pout = default_pout,
                       Coefficients coeffs = Coefficients(),
                       FlowConstants flow_constants = FlowConstants())
      : x(x), u(u), pout(pout), Comp(coeffs, flow_constants) {}

  /**
   * Calculate the derivative of the system.
   * Override the () operator to calculate the derivative using the provided
   * state x_in and the inputs given in the member variable u. Return the
   * derivative in the second argument.
   */
  void operator()(const CompressorState &x_in, CompressorState &dxdt,
                  const double /* t */) const {
    dxdt = GetDerivative(x_in, u, pout);
  }

  /// Return output of system using current state (x)
  inline CompressorOutput GetOutput() { return Comp::GetOutput(x); }

  /**
   * Integrate compressor comp from t0 to tf with timestep dt.
   * Integrate the compressor comp from its current state and input values at t0
   * until time tf, using timestep dt. After each timestep the function provided
   * in callback is called where, for example, the system inputs and other
   * parameters can be changed. The integration method used is a variable-step
   * Dormand-Prince algorithm which reduces the step size until the relative and
   * absolute error are with the bounds specified in max_rel_error and
   * max_abs_error respectively.
   */
  friend void IntegrateCompressor(CompressorSimulation comp, const double t0,
                                  const double tf, const double dt,
                                  IntegrationCallbackPtr callback,
                                  const double max_rel_error = 1e-6,
                                  const double max_abs_error = 1e-6) {
    ControlledStepper stepper =
        make_controlled(max_rel_error, max_abs_error, Dopri5Stepper());
    integrate_const(stepper, comp, comp.x, t0, tf, dt, callback);
  }

  CompressorInput u; ///< Current inputs to compressor.
  CompressorState x; ///< Current state of compressor.
  double pout; ///< Pressure at the outlet of compressor.

 private:
  // Type of stepper used to integrate
  typedef boost::numeric::odeint::runge_kutta_dopri5<
      CompressorState, double, CompressorState, double,
      boost::numeric::odeint::vector_space_algebra> Dopri5Stepper;

  // Type of stepper used for variable step size integration
  typedef boost::numeric::odeint::controlled_runge_kutta<Dopri5Stepper>
      ControlledStepper;
};

// Define norm of Eigen::Array for odeint to work
namespace boost {
namespace numeric {
namespace odeint {
template <>
struct vector_space_norm_inf<Comp::CompressorState> {
  typedef double result_type;
  double operator()(Comp::CompressorState x) const {
    double absval = 0;
    for (int i = 0; i < Comp::n_states; i++) absval += x[i] * x[i];
    return sqrt(absval);
  }
};
}
}
}

#endif
