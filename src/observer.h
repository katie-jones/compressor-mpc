#ifndef OBSERVER_H
#define OBSERVER_H

#include <Eigen/Eigen>
#include "aug_lin_sys.h"

template <class System, int n_delay_states, int n_disturbance_states>
class Observer {
 public:
  constexpr static int n_states = System::n_states;
  constexpr static int n_control_inputs = System::n_control_inputs;
  constexpr static int n_obs_states = System::n_states + n_disturbance_states;
  constexpr static int n_total_states = n_obs_states + n_delay_states;

  /// Observer of dynamic system
  typedef Eigen::Matrix<double, n_obs_states, System::n_outputs> ObserverMatrix;

  /// Augmented state of system
  typedef Eigen::Matrix<double, n_total_states, 1> AugmentedState;

  /// State of system
  typedef Eigen::Matrix<double, n_states, 1> State;

  /// Constructor
  Observer(const ObserverMatrix& M, const typename System::Output& y_init,
           const typename System::ControlInput& u_init =
               System::ControlInput::Zero(),
           const AugmentedState dx_init = AugmentedState::Zero())
      : M_(M), y_old_(y_init), u_old_(u_init), dx_aug_(dx_init) {}

  void InitializeSystem(const AugmentedLinearizedSystem<
      System, n_delay_states, n_disturbance_states>* p_auglinsys) {
    p_auglinsys_ = p_auglinsys;
  }

  /// apply observer
  State ObserveAPosteriori(const typename System::Output& y_in);

  /// Generate state prediction
  void ObserveAPriori(const typename System::ControlInput& u_in);

  /// Get augmented state
  AugmentedState GetStateEstimate() const { return dx_aug_; }

  /// Get previous output
  typename System::Output GetPreviousOutput() const { return y_old_; }

  /// Set initial output
  void SetIntialOutput(const typename System::Output& y_init) {
    y_old_ = y_init;
  }

 private:
  const ObserverMatrix M_;  // observer matrix used
  const AugmentedLinearizedSystem<System, n_delay_states, n_disturbance_states>*
      p_auglinsys_;                      // current augmented linearization
  typename System::Output y_old_;        // past output
  typename System::ControlInput u_old_;  // past input
  AugmentedState dx_aug_;                // differential augmented state
};

#endif
