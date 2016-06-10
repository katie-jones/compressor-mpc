#ifndef OBSERVER_H
#define OBSERVER_H

#include <Eigen/Eigen>
#include "aug_lin_sys.h"

template <class System, typename Delays, int n_disturbance_states>
class Observer {
 public:
  constexpr static int n_delay_states = Delays::GetSum();
  constexpr static int n_states = System::n_states;
  constexpr static int n_control_inputs = System::n_control_inputs;
  constexpr static int n_obs_states = System::n_states + n_disturbance_states;
  constexpr static int n_total_states = n_obs_states + n_delay_states;

  /// Observer of dynamic system
  typedef Eigen::Matrix<double, n_obs_states, System::n_outputs> ObserverMatrix;

 protected:
  /// Augmented state of system
  using AugmentedState =
      typename AugmentedLinearizedSystem<System, Delays,
                                         n_disturbance_states>::AugmentedState;

  /// State of system
  using State = typename AugmentedLinearizedSystem<System, Delays,
                                                   n_disturbance_states>::State;

  const ObserverMatrix M_;  // observer matrix used
  const AugmentedLinearizedSystem<System, Delays, n_disturbance_states>*
      p_auglinsys_;                // current augmented linearization
  typename System::Output y_old_;  // past output
  AugmentedState dx_aug_;          // differential augmented state

 public:
  /// Constructor
  Observer(const ObserverMatrix& M, const typename System::Output& y_init,
           const typename System::ControlInput& u_init =
               System::ControlInput::Zero(),
           const AugmentedState dx_init = AugmentedState::Zero())
      : M_(M), y_old_(y_init), dx_aug_(dx_init) {}

  void InitializeSystem(const AugmentedLinearizedSystem<
      System, Delays, n_disturbance_states>* p_auglinsys) {
    p_auglinsys_ = p_auglinsys;
  }

  /// apply observer
  State ObserveAPosteriori(const typename System::Output& y_in);

  /// Generate state prediction
  void ObserveAPriori(const typename System::ControlInput& du_in,
                      const typename System::ControlInput& u_old);

  /// Get augmented state
  AugmentedState GetStateEstimate() const { return dx_aug_; }

  /// Get previous output
  typename System::Output GetPreviousOutput() const { return y_old_; }

  /// Set initial output
  void SetIntialOutput(const typename System::Output& y_init) {
    y_old_ = y_init;
  }

  /// Set initial augmented state
  void SetInitialAugmentedState(const AugmentedState dx) { dx_aug_ = dx; }
};

#endif
