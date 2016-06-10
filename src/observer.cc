#include "observer.h"

/*
 * Predict next state of system based on input value to be applied
 */
template <class System, typename Delays, int n_total_states>
void Observer<System, Delays, n_total_states>::ObserveAPriori(
    const typename System::ControlInput& du_in, const typename System::ControlInput& u_old) {
  typename System::ControlInput du = System::ControlInput::Zero();
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
template <class System, typename Delays, int n_total_states>
typename Observer<System, Delays, n_total_states>::State
Observer<System, Delays, n_total_states>::ObserveAPosteriori(
    const typename System::Output& y_in) {
  // apply observer to non-delay states
  dx_aug_.template head<n_obs_states>() =
      dx_aug_.template head<n_obs_states>() +
      M_ *
          (y_in - y_old_ -
           static_cast<typename System::Output>(
               p_auglinsys_->C * (dx_aug_.template head<n_obs_states>())))
              .matrix();

  // store input value
  y_old_ = y_in;
  return dx_aug_.template head<n_states>();
}

#include "observer_list.h"
