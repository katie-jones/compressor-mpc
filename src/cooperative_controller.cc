#include "cooperative_controller.h"
#include <boost/timer/timer.hpp>

/*
 * Constructor
 */
template <class System, int n_delay_states, int n_disturbance_states, int p,
          int m, int n_controllers>
CooperativeController<System, n_delay_states, n_disturbance_states, p, m,
                      n_controllers>::
    CooperativeController(
        const AugmentedLinearizedSystem<System, n_delay_states,
                                        n_disturbance_states>& sys,
        const Observer<System, n_delay_states, n_disturbance_states>& observer,
        const Input& u_offset, const OutputPrediction& y_ref,
        const ControlInputIndex& input_delay,
        const ControlInputIndex& control_input_index,
        const int n_solver_iterations,
        const InputConstraints<n_control_inputs>& constraints,
        const UWeightType& u_weight, const YWeightType& y_weight)
    : ControllerInterface<System, p>(y_ref, u_offset, input_delay,
                                     control_input_index),
      auglinsys_(sys),
      n_solver_iterations_(n_solver_iterations),
      observer_(observer) {
  // Check number of delay states
  int sum_delay = 0;
  for (int i = 0; i < n_control_inputs; i++) sum_delay += n_delay_[i];
  if (sum_delay != n_delay_states) {
    throw delay_states_wrong();
  }

  // Initialize observer pointer
  observer_.InitializeSystem(&auglinsys_);

  // Intialize subsolvers with correct constraints and initial states
  sub_solvers_.reserve(n_controllers);
  InputConstraints<n_sub_control_inputs> sub_constraints;
  for (int i = 0; i < n_controllers; i++) {
    sub_constraints.lower_bound =
        constraints.lower_bound.template segment<n_sub_control_inputs>(
            i * n_sub_control_inputs);
    sub_constraints.upper_bound =
        constraints.upper_bound.template segment<n_sub_control_inputs>(
            i * n_sub_control_inputs);
    sub_constraints.lower_rate_bound =
        constraints.lower_rate_bound.template segment<n_sub_control_inputs>(
            i * n_sub_control_inputs);
    sub_constraints.upper_rate_bound =
        constraints.upper_rate_bound.template segment<n_sub_control_inputs>(
            i * n_sub_control_inputs);

    sub_solvers_.emplace_back(
        &y_ref_, i, sub_constraints,
        u_weight.template block<n_sub_control_inputs, n_sub_control_inputs>(
            i * n_sub_control_inputs, i * n_sub_control_inputs),
        y_weight);
  }
}

/*
 * Set initial state, output and inputs of system and initialize QP
 */
template <class System, int n_delay_states, int n_disturbance_states, int p,
          int m, int n_controllers>
void CooperativeController<System, n_delay_states, n_disturbance_states, p, m,
                           n_controllers>::SetInitialState(const State& x_init,
                                                           const Output& y_init,
                                                           const ControlInput&
                                                               u_init) {
  x_ = x_init;
  observer_.SetIntialOutput(y_init);
  u_old_ = u_init;

  observer_.SetIntialOutput(y_init);

  AugmentedState delta_x0;
  delta_x0 << auglinsys_.GetF(),
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

  auglinsys_.Update(x_init, this->GetPlantInput(u_init));
  const Prediction pred = auglinsys_.GeneratePrediction(p, m);

  QP qp;
  Eigen::MatrixXd sub_Su[n_controllers];
  SplitSuMatrix(sub_Su, pred.Su);
  for (int i = 0; i < n_controllers; i++) {
    sub_solvers_[i].GenerateDistributedQP(qp, sub_Su[i], pred.Sx, pred.Sf,
                                          delta_x0, n_aug_states, y_init);
    sub_solvers_[i].InitializeQPProblem(
        qp, u_old_.template segment<n_sub_control_inputs>(
                i * n_sub_control_inputs));
  }
}

/**
 * Compute the next input to apply to the system.
 * Linearizes the system about current state estimate and finds the optimal
 * input value by solving a QP it generates using the MPC formulation.
 */
template <class System, int n_delay_states, int n_disturbance_states, int p,
          int m, int n_controllers>
const typename CooperativeController<System, n_delay_states,
                                     n_disturbance_states, p, m,
                                     n_controllers>::ControlInput
CooperativeController<System, n_delay_states, n_disturbance_states, p, m,
                      n_controllers>::GetNextInput(const Output& y,
                                                   std::ofstream&
                                                       cpu_time_out) {
  boost::timer::cpu_timer integrate_timer;

  x_ += observer_.ObserveAPosteriori(y);
  auglinsys_.Update(x_, this->GetPlantInput(u_old_));

  const Prediction pred = auglinsys_.GeneratePrediction(p, m);

  AugmentedState delta_x0;
  delta_x0 << auglinsys_.GetF(),
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

  QP qp[n_controllers];
  Eigen::MatrixXd sub_Su[n_controllers];

  SplitSuMatrix(sub_Su, pred.Su);  // assign values to sub_Su

  SubControlInputPrediction du_out[n_controllers];
  for (int i = 0; i < n_controllers; i++) {
    sub_solvers_[i].GenerateDistributedQP(qp[i], sub_Su[i], pred.Sx, pred.Sf,
                                          delta_x0, n_aug_states, y);
    if (i == 0) integrate_timer.stop();
  }

  for (int n_iter = 0; n_iter < n_solver_iterations_; n_iter++) {
    // Re-solve QP with new inputs
    for (int i = 0; i < n_controllers; i++) {
      if (i == 0) integrate_timer.resume();
      sub_solvers_[i].UpdateAndSolveQP(
          qp[i], du_out[i], u_old_.template segment<n_sub_control_inputs>(
                                i * n_sub_control_inputs),
          sub_Su, du_prev_);
      if (i == 0) integrate_timer.stop();
    }

    // Update du_prev for next iteration
    for (int i = 0; i < n_controllers; i++) {
      if (i == 0) integrate_timer.resume();
      du_prev_.template segment<m* n_sub_control_inputs>(
          i * m * n_sub_control_inputs) = du_out[i];
      if (i == 0) integrate_timer.stop();
    }
  }

  // Add new du to u_old_
  ControlInput du_applied;
  for (int i = 0; i < n_controllers; i++) {
    if (i == 0) integrate_timer.resume();
    du_applied.template segment<n_sub_control_inputs>(i *
                                                      n_sub_control_inputs) =
        du_prev_.template segment<n_sub_control_inputs>(i * m *
                                                        n_sub_control_inputs);

    // Initialize du_prev_ for next iteration using predictions
    if (m > 1) {
      du_prev_.template segment<(m - 1)* n_sub_control_inputs>(
          i * m * n_sub_control_inputs) =
          du_prev_.template segment<(m - 1) * n_sub_control_inputs>(
              (i * m + 1) * n_sub_control_inputs);

      // Repeat last prediction
      du_prev_.template segment<n_sub_control_inputs>(((i + 1) * m - 1) *
                                                      n_sub_control_inputs) =
          du_prev_.template segment<n_sub_control_inputs>(((i + 1) * m - 2) *
                                                          n_sub_control_inputs);
    } else {
      du_prev_.template segment<n_sub_control_inputs>(i * n_sub_control_inputs)
          .setZero();
    }
    if (i == 0) integrate_timer.stop();
  }

  boost::timer::cpu_times int_elapsed = integrate_timer.elapsed();
  // boost::timer::nanosecond_type elapsed_ns(int_elapsed.system +
  // int_elapsed.user);
  boost::timer::nanosecond_type elapsed_ns(int_elapsed.wall);
  cpu_time_out << elapsed_ns << std::endl;

  u_old_ += du_applied;
  observer_.ObserveAPriori(du_applied);

  return u_old_;
}

#include "cooperative_controller_list.h"
