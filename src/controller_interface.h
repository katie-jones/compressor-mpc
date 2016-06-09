#ifndef CONTROLLER_INTERFACE_H
#define CONTROLLER_INTERFACE_H

#include <Eigen/Eigen>
#include <Eigen/SparseCore>
#include "dynamic_system.h"
#include "aug_lin_sys.h"
#include "observer.h"
#include "mpc_exceptions.h"
#include "qpOASES.hpp"
#include "input_constraints.h"

template <class System, int p>
class ControllerInterface {
 protected:
  static constexpr int n_control_inputs = System::n_control_inputs;
  static constexpr int n_inputs = System::n_inputs;
  static constexpr int n_outputs = System::n_outputs;
  static constexpr int n_states = System::n_states;

 public:
  /// State of the dynamic system
  typedef Eigen::Matrix<double, n_states, 1> State;

  /// Input to dynamic system
  typedef Eigen::Matrix<double, n_inputs, 1> Input;

  /// Controlled inputs of dynamic system
  typedef Eigen::Matrix<double, n_control_inputs, 1> ControlInput;

  /// Output of dynamic system
  typedef Eigen::Matrix<double, n_outputs, 1> Output;

  /// Output of dynamic system for p time steps
  typedef Eigen::Matrix<double, n_outputs * p, 1> OutputPrediction;

  /// Vector of indices for a ControlInput
  typedef std::array<int, n_control_inputs> ControlInputIndex;

  /// Constructor
  ControllerInterface(const Input& u_offset, const ControlInputIndex& n_delay,
                      const ControlInputIndex& control_input_index)
      : u_offset_(u_offset),
        n_delay_(n_delay),
        control_input_index_(control_input_index) {}

  /**
   * Compute the next input to apply to the system.
   * Linearizes the system about current state estimate and finds the optimal
   * input value by solving a QP it generates using the MPC formulation.
   */
  virtual const ControlInput GetNextInput(const Output& y) = 0;

  /// output plant input based on control input and offset
  const Input GetPlantInput(const ControlInput& u_control) const {
    Input u = u_offset_;
    for (int i = 0; i < n_control_inputs; i++) {
      u(control_input_index_[i]) += u_control(i);
    }
    return u;
  }

  /// output control input based on plant input and offset
  const ControlInput GetControlInput(const Input& u) const {
    ControlInput u_control;
    for (int i = 0; i < n_control_inputs; i++) {
      u_control(i) =
          u(control_input_index_[i]) - u_offset_(control_input_index_[i]);
    }
    return u_control;
  }

 protected:
  const Input u_offset_;             // offset applied to control input
  const ControlInputIndex n_delay_;  // delay states per input
  const ControlInputIndex
      control_input_index_;  // index such that ControlInput[i] ->
                             // Input[control_input_index_[i]]
};

#endif
