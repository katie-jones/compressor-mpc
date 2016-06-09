#include "mpc_controller.h"

#include <boost/timer/timer.hpp>

/*
 * Constructor
 */
template <class System, int n_delay_states, int n_disturbance_states, int p,
          int m>
MpcController<System, n_delay_states, n_disturbance_states, p, m>::
    MpcController(
        const AugmentedLinearizedSystem<System, n_delay_states,
                                        n_disturbance_states>& sys,
        const Observer<System, n_delay_states, n_disturbance_states>& observer,
        const Input& u_offset, const OutputPrediction& y_ref,
        const ControlInputIndex& input_delay,
        const ControlInputIndex& control_input_index,
        const UWeightType& u_weight, const YWeightType& y_weight,
        const InputConstraints<n_control_inputs>& constraints)
    : auglinsys_(sys),
      observer_(observer),
      ControllerInterface<System, p>(u_offset, input_delay,
                                     control_input_index),
      MpcQpSolver<n_total_states, n_outputs, n_control_inputs, p, m>(
          constraints, y_ref, u_weight, y_weight) {
  // Check number of delay states
  int sum_delay = 0;
  for (int i = 0; i < n_control_inputs; i++) sum_delay += n_delay_[i];
  if (sum_delay != n_delay_states) {
    throw delay_states_wrong();
  }
  observer_.InitializeSystem(&auglinsys_);
}

/*
 * Solve QP and return optimal input to apply
 */
template <class System, int n_delay_states, int n_disturbance_states, int p,
          int m>
const typename MpcController<System, n_delay_states, n_disturbance_states, p,
                             m>::ControlInput
MpcController<System, n_delay_states, n_disturbance_states, p, m>::GetNextInput(
    const Output& y, std::ofstream& cpu_time_out) {
  boost::timer::cpu_timer integrate_timer;
  x_ += observer_.ObserveAPosteriori(y);
  auglinsys_.Update(x_, this->GetPlantInput(u_old_));
  const Prediction pred = auglinsys_.GeneratePrediction(p, m);

  AugmentedState delta_x0;
  delta_x0 << auglinsys_.GetDerivative(),
      observer_.GetStateEstimate().template tail<n_aug_states>();

  int index_delay_states = n_obs_states;
  for (int i = 0; i < n_control_inputs; i++) {
    if (n_delay_[i] != 0) {
      for (int j = 0; j < n_delay_[i]; j++) {
        delta_x0(index_delay_states + j) -= u_old_[i];
      }
      index_delay_states += n_delay_[i];
    }
  }

  // Solve QP and get optimal input
  const QP qp = this->GenerateQP(pred, delta_x0, n_aug_states,
                                 observer_.GetPreviousOutput());
  const ControlInputPrediction usol = this->SolveQP(qp, u_old_);

  // Update state estimation based on input
  observer_.ObserveAPriori(usol.template head<n_control_inputs>(), u_old_);

  // Update previous solution
  u_old_ += usol.template head<n_control_inputs>();

  // Measure time
  boost::timer::cpu_times int_elapsed = integrate_timer.elapsed();
  boost::timer::nanosecond_type elapsed_ns(int_elapsed.wall);
  cpu_time_out << elapsed_ns << std::endl;

  return u_old_;
}

/*
 * Initialize initial states and QP
 */
template <class System, int n_delay_states, int n_disturbance_states, int p,
          int m>
void MpcController<System, n_delay_states, n_disturbance_states, p,
                   m>::SetInitialState(const State& x_init,
                                       const Output& y_init,
                                       const ControlInput& u_init) {
  x_ = x_init;
  u_old_ = u_init;

  observer_.SetIntialOutput(y_init);

  AugmentedState delta_x0;
  delta_x0 << auglinsys_.GetDerivative(),
      observer_.GetStateEstimate().template tail<n_aug_states>();

  int index_delay_states = n_obs_states;
  for (int i = 0; i < n_control_inputs; i++) {
    if (n_delay_[i] != 0) {
      for (int j = 0; j < n_delay_[i]; j++) {
        delta_x0(index_delay_states + j) -= u_old_[i];
      }
      index_delay_states += n_delay_[i];
    }
  }

  const Prediction pred = auglinsys_.GeneratePrediction(p, m);
  const QP qp = this->GenerateQP(pred, delta_x0, n_aug_states,
                                 observer_.GetPreviousOutput());
  this->InitializeQPProblem(qp, u_old_);
}

#include "mpc_controller_list.h"
