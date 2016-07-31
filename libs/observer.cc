#include "observer.h"

/*
 * Predict next state of system based on input value to be applied
 */
template <class AugLinSys>
void Observer<AugLinSys>::ObserveAPriori(
    const ControlInput& du_in, const ControlInput& u_old) {
  ControlInput du = ControlInput::Zero();
  AugmentedState dx = dx_aug_;

  dx.template head<n_states>().setZero();
  int index_delay_states = n_obs_states;

  for (int i = 0; i < n_control_inputs; i++) {
    if (p_auglinsys_->n_delay_[i] == 0) {
      du(i) = du_in(i);
    } else {
      du(i) = u_old[i] + du_in(i);
      dx(index_delay_states) -= u_old[i];
      index_delay_states += p_auglinsys_->n_delay_[i];
    }
  }

  dx_aug_ = p_auglinsys_->B * du + p_auglinsys_->A * dx;
  dx_aug_.template head<n_states>() += p_auglinsys_->f;
}

/*
 * Observe system given new output values from plant
 */
template <class AugLinSys>
typename Observer<AugLinSys>::State
Observer<AugLinSys>::ObserveAPosteriori(
    const Output& y_in) {
  // apply observer to non-delay states
  dx_aug_.template head<n_obs_states>() =
      dx_aug_.template head<n_obs_states>() +
      M_ *
          (y_in - y_old_ -
           static_cast<Output>(
               p_auglinsys_->C * (dx_aug_.template head<n_obs_states>())))
              .matrix();

  // store input value
  y_old_ = y_in;
  return dx_aug_.template head<n_states>();
}

#include "observer_list.h"
