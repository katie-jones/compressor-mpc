#ifndef OBSERVER_H
#define OBSERVER_H

#include <Eigen/Eigen>

template <class AugLinSys>
class Observer {
 public:
  constexpr static int n_delay_states = AugLinSys::n_delay_states;
  constexpr static int n_states = AugLinSys::n_states;
  constexpr static int n_control_inputs = AugLinSys::n_control_inputs;
  constexpr static int n_disturbance_states = AugLinSys::n_disturbance_states;
  constexpr static int n_obs_states = n_states + n_disturbance_states;
  constexpr static int n_total_states = n_obs_states + n_delay_states;

  /// Observer of dynamic system
  typedef Eigen::Matrix<double, n_obs_states, AugLinSys::n_outputs>
      ObserverMatrix;

 protected:
  /// Augmented state of system
  using AugmentedState = typename AugLinSys::AugmentedState;

  /// State of system
  using State = typename AugLinSys::State;
  using Output = typename AugLinSys::Output;
  using Input = typename AugLinSys::Input;
  using ControlInput = typename AugLinSys::ControlInput;

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

  void InitializeSystem(const AugLinSys* p_auglinsys) {
    p_auglinsys_ = p_auglinsys;
  }

  /// apply observer
  State ObserveAPosteriori(const Output& y_in);

  /// Generate state prediction
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
