#ifndef OBSERVER_H
#define OBSERVER_H

#include <Eigen/Eigen>

/**
 * Observer for an augmented linear system of type AugLinSys.
 * Includes functions to observe system a priori and a posteriori based on the
 * linearized system and a given ObserverMatrix. Keeps track of a differential
 * augmented state.
 */
template <class AugLinSys>
class Observer {
 public:
   /// Total number of delayed states
  constexpr static int n_delay_states = AugLinSys::n_delay_states;
/// Number of states in non-augmented system
  constexpr static int n_states = AugLinSys::n_states;
  /// Number of control inputs
  constexpr static int n_control_inputs = AugLinSys::n_control_inputs;
  /// Number of disturbance/integrator states
  constexpr static int n_disturbance_states = AugLinSys::n_disturbance_states;
  /// Number of observable states
  constexpr static int n_obs_states = n_states + n_disturbance_states;
  /// Number of total states (augmented + real)
  constexpr static int n_total_states = n_obs_states + n_delay_states;

  /// Observer matrix of dynamic system
  typedef Eigen::Matrix<double, n_obs_states, AugLinSys::n_outputs,
                        Eigen::RowMajor>
      ObserverMatrix;

 public:
  /// Augmented state of system
  using AugmentedState = typename AugLinSys::AugmentedState;

  /// State of system
  using State = typename AugLinSys::State;
  /// Output of system
  using Output = typename AugLinSys::Output;
  /// Input to system
  using Input = typename AugLinSys::Input;
  /// Control input to system
  using ControlInput = typename AugLinSys::ControlInput;

 protected:
  const ObserverMatrix M_;        // observer matrix used
  const AugLinSys* p_auglinsys_;  // current augmented linearization
  Output y_old_;                  // past output
  AugmentedState dx_aug_;         // differential augmented state

 public:
  /// Constructor
  Observer(const ObserverMatrix& M, const Output& y_init,
           const ControlInput& u_init = AugLinSys::ControlInput::Zero(),
           const AugmentedState dx_init = AugmentedState::Zero())
      : M_(M), y_old_(y_init), dx_aug_(dx_init) {}

  /// Set the pointer to the linearized system to be used
  void InitializeSystem(const AugLinSys* p_auglinsys) {
    p_auglinsys_ = p_auglinsys;
  }

  /// Apply observer a posteriori
  State ObserveAPosteriori(const Output& y_in);

  /// Apply observer a priori (prediction)
  void ObserveAPriori(const ControlInput& du_in, const ControlInput& u_old);

  /// Get augmented state
  AugmentedState GetStateEstimate() const { return dx_aug_; }

  /// Get previous output
  Output GetPreviousOutput() const { return y_old_; }

  /// Set initial output
  void SetIntialOutput(const Output& y_init) { y_old_ = y_init; }

  /// Set initial augmented state
  void SetInitialAugmentedState(const AugmentedState dx) { dx_aug_ = dx; }
};

#endif
