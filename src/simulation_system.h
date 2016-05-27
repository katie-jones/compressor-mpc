#ifndef SIMULATION_SYSTEM_H
#define SIMULATION_SYSTEM_H

#include <Eigen/Eigen>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/algebra/vector_space_algebra.hpp>

#include "time_delay.h"

/**
 * Class describing a dynamic system used for simulation.
 * Includes the current system state and
 * inputs as member variables. Overrides the () operator to define the
 * derivative. Defines a function to integrate the system for a given
 * time range.
 */
template <class System, int n_delay_states>
class SimulationSystem {
 public:
  /// Vector of indices for a ControlInput
  typedef std::array<int, System::n_control_inputs> ControlInputIndex;

 private:
  typedef typename System::State State;
  typedef typename System::Input Input;
  typedef typename System::Output Output;
  typedef typename System::ControlInput ControlInput;

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
  System* p_sys_;  // system to simulate
  TimeDelay<System::n_control_inputs, n_delay_states> delayed_inputs_;

  // index such that ControlInput[i] -> Input[control_input_index_[i]]
  const ControlInputIndex control_input_index_;
  const ControlInputIndex n_delay_;  // time delays of each input

 public:
  /// Function pointer used for callback when integrating system.
  typedef void (*IntegrationCallbackPtr)(const State, const double);

  SimulationSystem(System* p_sys, Input u_offset, ControlInputIndex input_index,
                   ControlInputIndex n_delay, State x_in,
                   ControlInput u_init = ControlInput::Zero())
      : p_sys_(p_sys),
        x_(x_in),
        u_offset_(u_offset),
        control_input_index_(input_index),
        n_delay_(n_delay),
        delayed_inputs_(
            TimeDelay<System::n_control_inputs, n_delay_states>(n_delay)),
        u_(GetPlantInput(u_init)) {}

  ControlInput GetCurrentInput() { return u_; }
  State GetCurrentState() { return x_; }

  /// Update input offset.
  void SetOffset(const Input& u_in) { u_offset_ = u_in; }

  /// Set input from controller.
  void SetInput(const ControlInput u) {
    u_ = GetPlantInput(delayed_inputs_.GetDelayedInput(u));
  }

  /// Update current state
  void SetState(State x_in) { x_ = x_in; }

  /// Return output at current system state x
  Output GetOutput() const { return GetOutput(x_); }

  /// output plant input based on control input and offset
  const Input GetPlantInput(const ControlInput& u_control) const {
    Input u = u_offset_;
    for (int i = 0; i < System::n_control_inputs; i++) {
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
    dxdt = p_sys_->GetDerivative(x_in, u_);
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
