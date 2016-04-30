#ifndef SIMULATION_SYSTEM_H
#define SIMULATION_SYSTEM_H

#include <Eigen/Eigen>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/algebra/vector_space_algebra.hpp>

#include "dynamic_system.h"

/**
 * Abstract class describing a dynamic system used for simulation.
 * Template with parameters: n_states (number of states), n_inputs (number of
 * inputs), n_outputs (number of outputs). Includes the current system state and
 * inputs as member variables. Overrides the () operator to define the
 * derivative. Defines a friend function to integrate the system for a given
 * time range.
 */
template <int n_states, int n_inputs, int n_outputs>
class SimulationSystem : public DynamicSystem<n_states, n_inputs, n_outputs> {
 protected:
  typedef typename DynamicSystem<n_states, n_inputs, n_outputs>::State State;
  typedef typename DynamicSystem<n_states, n_inputs, n_outputs>::Input Input;
  typedef typename DynamicSystem<n_states, n_inputs, n_outputs>::Output Output;

  // Type of stepper used to integrate
  typedef boost::numeric::odeint::runge_kutta_dopri5<
      State, double, State, double,
      boost::numeric::odeint::vector_space_algebra> Dopri5Stepper;

  // Type of stepper used for variable step size integration
  typedef boost::numeric::odeint::controlled_runge_kutta<Dopri5Stepper>
      ControlledStepper;

  Input u;
  State x;

 public:
  /// Function pointer used for callback when integrating system.
  typedef void (*IntegrationCallbackPtr)(const State, const double);

  Input GetCurrentInput() { return u; }
  State GetCurrentState() { return x; }

  /**
   * Update input values, considering constraints.
   * If one or more of the entries in u_in violates the system constraints, the
   * maximum/minimum value that respects the constraints is used.
   */
  void SetInput(Input u_in) {
    for (int i = 0; i < n_inputs; i++) {
      u[i] = std::fmin(
          std::fmin(std::fmax(std::fmax(u_in[i], this->lower_input_constraint[i]),
                              u[i] - this->lower_input_rate_constraint[i]),
                    this->upper_input_constraint[i]),
          u[i] + this->upper_input_rate_constraint[i]);
    }
  }

  void SetState(State x_in) { x = x_in; }

  /// Return output at current system state x
  inline Output GetOutput() { return GetOutput(x); }

  /**
   * Calculate the derivative of the system.
   * Override the () operator to calculate the derivative using the provided
   * state x_in and the inputs given in the member variable u. Return the
   * derivative in the second argument.
   */
  void operator()(const State &x_in, State &dxdt, const double) const {
    dxdt = GetDerivative(x_in, u);
  }

  /**
   * Integrate system from t0 to tf with timestep dt.
   * Integrate the system from its current state and input values at t0
   * until time tf, using timestep dt. After each timestep the function provided
   * in callback is called where, for example, the system inputs and other
   * parameters can be changed. The integration method used is a variable-step
   * Dormand-Prince algorithm which reduces the step size until the relative and
   * absolute error are with the bounds specified in max_rel_error and
   * max_abs_error respectively.
   */
  friend void IntegrateSystem(SimulationSystem *sys, const double t0,
                              const double tf, const double dt,
                              IntegrationCallbackPtr callback,
                              const double max_rel_error = 1e-6,
                              const double max_abs_error = 1e-6) {
    ControlledStepper stepper =
        make_controlled(max_rel_error, max_abs_error, Dopri5Stepper());
    integrate_const(stepper, sys, sys->x, t0, tf, dt, callback);
  }
};

// Define vector_space_norm_inf for the state used in order for odeint to work
// namespace boost {
// namespace numeric {
// namespace odeint {
// template <>
// struct vector_space_norm_inf<DynamicSystem::State> {
  // typedef double result_type;
  // double operator()(DynamicSystem::State x) const {
    // double absval = 0;

    // for (int i = 0; i < DynamicSystem::n_states; i++) absval += x[i] * x[i];
    // return sqrt(absval);
  // }
// };
// }
// }
// }
#endif
