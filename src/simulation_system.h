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
template <int n_states, int n_inputs, int n_outputs, int n_control_inputs>
class SimulationSystem
    : public virtual DynamicSystem<n_states, n_inputs, n_outputs,
                                   n_control_inputs> {
 public:
  /// Vector of indices for a ControlInput
  typedef std::array<int, n_control_inputs> ControlInputIndex;

 protected:
  typedef typename DynamicSystem<n_states, n_inputs, n_outputs,
                                 n_control_inputs>::State State;
  typedef typename DynamicSystem<n_states, n_inputs, n_outputs,
                                 n_control_inputs>::Input Input;
  typedef typename DynamicSystem<n_states, n_inputs, n_outputs,
                                 n_control_inputs>::Output Output;
  typedef typename DynamicSystem<n_states, n_inputs, n_outputs,
                                 n_control_inputs>::ControlInput ControlInput;

  // Type of stepper used to integrate
  typedef boost::numeric::odeint::runge_kutta_dopri5<
      State, double, State, double,
      boost::numeric::odeint::vector_space_algebra> Dopri5Stepper;

  // Type of stepper used for variable step size integration
  typedef boost::numeric::odeint::controlled_runge_kutta<Dopri5Stepper>
      ControlledStepper;

  State x_;
  Input u_offset_;
  Input u_;

  // index such that ControlInput[i] -> Input[control_input_index_[i]]
  const ControlInputIndex control_input_index_;

  SimulationSystem(Input u_offset, ControlInputIndex input_index,State x_in, 
                   ControlInput u_init = ControlInput::Zero())
      : x_(x_in),
        u_offset_(u_offset),
        control_input_index_(input_index),
        u_(GetPlantInput(u_init)) {}

 public:
  /// Function pointer used for callback when integrating system.
  typedef void (*IntegrationCallbackPtr)(const State, const double);

  ControlInput GetCurrentInput() { return u_; }
  State GetCurrentState() { return x_; }

  /// Update input offset.
  void SetOffset(const Input& u_in) { u_offset_ = u_in; }

  /// Set input from controller.
  void SetInput(const ControlInput u) { u_ = GetPlantInput(u); }
  // void SetInput(const Input &u) { u_ = u; }

  /// Update current state
  void SetState(State x_in) { x_ = x_in; }

  /// Return output at current system state x
  Output GetOutput() const { return GetOutput(x_); }

  /// output plant input based on control input and offset
  const Input GetPlantInput(const ControlInput& u_control) const {
    Input u = u_offset_;
    for (int i = 0; i < n_control_inputs; i++) {
      u(control_input_index_[i]) += u_control(i);
    }
    return u;
  }

  /**
   * Calculate the derivative of the system.
   * Override the () operator to calculate the derivative using the provided
   * state x_in and the inputs given in the member variable u. Return the
   * derivative in the second argument.
   */
  void operator()(const State& x_in, State& dxdt, const double) const {
    dxdt = this->GetDerivative(x_in, u_);
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
  void Integrate(const double t0, const double tf, const double dt,
                 IntegrationCallbackPtr callback,
                 const double max_rel_error = 1e-6,
                 const double max_abs_error = 1e-6) {
    ControlledStepper stepper =
        make_controlled(max_rel_error, max_abs_error, Dopri5Stepper());
    boost::numeric::odeint::integrate_const(stepper, std::ref(*this), x_, t0,
                                            tf, dt, callback);
  }
};

#endif
