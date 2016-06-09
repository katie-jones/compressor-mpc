#include "noncooperative_controller.h"
#include <boost/timer/timer.hpp>

/*
 * Constructor
 */
template <class System, int n_delay_states, int n_disturbance_states, int p,
          int m, int n_controllers, int n_sub_outputs>
NonCooperativeController<System, n_delay_states, n_disturbance_states, p, m,
                         n_controllers, n_sub_outputs>::
    NonCooperativeController(
        const AugmentedLinearizedSystem<System, n_delay_states,
                                        n_disturbance_states>& sys,
        const Observer<System, n_delay_states, n_disturbance_states>& observer,
        const Input& u_offset, const OutputPrediction& y_ref,
        const ControlInputIndex& input_delay,
        const ControlInputIndex& control_input_index,
        const int n_solver_iterations,
        const Eigen::Array<int, n_controllers, n_sub_outputs> sub_output_index,
        const InputConstraints<n_control_inputs>& constraints,
        const UWeightType& u_weight, const YWeightType& y_weight)
    : ControllerInterface<System, p>(u_offset, input_delay,
                                     control_input_index),
      auglinsys_(sys),
      n_solver_iterations_(n_solver_iterations),
      sub_output_index_(sub_output_index),
      observer_(observer) {
  // make array of yweights all mapped to the same data and use them to
  // initialize
  std::array<YWeightType, n_controllers> y_weights;
  for (int i = 0; i < n_controllers; i++) {
    y_weights[i] = Eigen::Map<const YWeightType>(y_weight.data());
  }
  InitializeArguments(y_ref, y_weights, constraints, u_weight);
}

/*
 * Constructor with multiple y weights
 */
template <class System, int n_delay_states, int n_disturbance_states, int p,
          int m, int n_controllers, int n_sub_outputs>
NonCooperativeController<System, n_delay_states, n_disturbance_states, p, m,
                         n_controllers, n_sub_outputs>::
    NonCooperativeController(
        const AugmentedLinearizedSystem<System, n_delay_states,
                                        n_disturbance_states>& sys,
        const Observer<System, n_delay_states, n_disturbance_states>& observer,
        const Input& u_offset, const OutputPrediction& y_ref,
        const ControlInputIndex& input_delay,
        const ControlInputIndex& control_input_index,
        const int n_solver_iterations,
        const Eigen::Array<int, n_controllers, n_sub_outputs> sub_output_index,
        const std::array<YWeightType, n_controllers>& y_weights,
        const InputConstraints<n_control_inputs>& constraints,
        const UWeightType& u_weight)
    : ControllerInterface<System, p>(u_offset, input_delay,
                                     control_input_index),
      auglinsys_(sys),
      n_solver_iterations_(n_solver_iterations),
      sub_output_index_(sub_output_index),
      observer_(observer) {
  InitializeArguments(y_ref, y_weights, constraints, u_weight);
}

/*
 * Initialize member variables (called by constructors)
 */
template <class System, int n_delay_states, int n_disturbance_states, int p,
          int m, int n_controllers, int n_sub_outputs>
void NonCooperativeController<System, n_delay_states, n_disturbance_states, p,
                              m, n_controllers, n_sub_outputs>::
    InitializeArguments(const OutputPrediction& y_ref,
                        const std::array<YWeightType, n_controllers>& y_weights,
                        const InputConstraints<n_control_inputs>& constraints,
                        const UWeightType& u_weight) {
  // Check number of delay states
  int sum_delay = 0;
  for (int i = 0; i < n_control_inputs; i++) sum_delay += n_delay_[i];
  if (sum_delay != n_delay_states) {
    throw delay_states_wrong();
  }

  // Initialize observer pointer
  observer_.InitializeSystem(&auglinsys_);

  // Initialize previous QP solutions
  for (int i = 0; i < n_controllers; i++) {
    du_prev_[i] = SubControlInputPrediction::Zero();
  }

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
        i, sub_constraints, SubOutputPrediction::Zero(),
        u_weight.template block<n_sub_control_inputs, n_sub_control_inputs>(
            i * n_sub_control_inputs, i * n_sub_control_inputs),
        y_weights[i]);
  }

  // Split output reference
  SetOutputReference(y_ref);

}

/*
 * Set initial state, output and inputs of system and initialize QP
 */
template <class System, int n_delay_states, int n_disturbance_states, int p,
          int m, int n_controllers, int n_sub_outputs>
void NonCooperativeController<
    System, n_delay_states, n_disturbance_states, p, m, n_controllers,
    n_sub_outputs>::SetInitialState(const State& x_init, const Output& y_init,
                                    const ControlInput& u_init) {
  x_ = x_init;
  u_old_ = u_init;
  observer_.SetIntialOutput(y_init);
  auglinsys_.Update(x_init, this->GetPlantInput(u_init));

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

  auglinsys_.Update(x_init, this->GetPlantInput(u_init));

  QP qp;

  // Split output
  SubOutput y_subs[n_controllers];
  SplitOutput(y_subs, y_init);

  Prediction pred;
  Eigen::MatrixXd sub_Su(p * n_sub_outputs, m * n_sub_control_inputs);
  int* output_index_data = sub_output_index_.data();
  for (int i = 0; i < n_controllers; i++) {
    pred = auglinsys_.template GenerateSubPrediction<n_sub_outputs>(
        p, m, output_index_data);

    for (int j = 0; j < m; j++) {
      sub_Su.template block<p * n_sub_outputs, n_sub_control_inputs>(
          0, j * n_sub_control_inputs) =
          pred.Su.template block<p * n_sub_outputs, n_sub_control_inputs>(
              0, (j * n_controllers + i) * n_sub_control_inputs);
    }

    sub_solvers_[i].GenerateDistributedQP(qp, sub_Su, pred.Sx, pred.Sf,
                                          delta_x0, n_aug_states, y_subs[i]);
    sub_solvers_[i].InitializeQPProblem(
        qp, u_old_.template segment<n_sub_control_inputs>(
                i * n_sub_control_inputs));

    output_index_data += n_sub_outputs;
  }
}

/**
 * Compute the next input to apply to the system.
 * Linearizes the system about current state estimate and finds the optimal
 * input value by solving a QP it generates using the MPC formulation.
 */
template <class System, int n_delay_states, int n_disturbance_states, int p,
          int m, int n_controllers, int n_sub_outputs>
const typename NonCooperativeController<
    System, n_delay_states, n_disturbance_states, p, m, n_controllers,
    n_sub_outputs>::ControlInput
NonCooperativeController<
    System, n_delay_states, n_disturbance_states, p, m, n_controllers,
    n_sub_outputs>::GetNextInput(const Output& y, std::ofstream& cpu_time_out) {
  boost::timer::cpu_timer integrate_timer;

  x_ += observer_.ObserveAPosteriori(y);
  auglinsys_.Update(x_, this->GetPlantInput(u_old_));

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

  QP qp[n_controllers];

  // Split prediction matrices
  Eigen::MatrixXd sub_Su[n_controllers];
  SubControlInputPrediction du_out[n_controllers];

  // Split output
  SubOutput y_subs[n_controllers];
  SplitOutput(y_subs, y);

  Prediction pred;
  int* output_index_data = sub_output_index_.data();
  for (int i = 0; i < n_controllers; i++) {
    pred = auglinsys_.template GenerateSubPrediction<n_sub_outputs>(
        p, m, output_index_data);
    sub_Su[i] = Eigen::MatrixXd::Zero(p * n_sub_outputs, m * n_control_inputs);

    // Rearrange columns so they're ordered first by controller number, then by
    // prediction number (pred.Su is the other way around)
    for (int k = 0; k < n_controllers; k++) {
      for (int j = 0; j < m; j++) {
        sub_Su[i].template block<p * n_sub_outputs, n_sub_control_inputs>(
            0, (k * m + j) * n_sub_control_inputs) =
            pred.Su.template block<p * n_sub_outputs, n_sub_control_inputs>(
                0, (j * n_controllers + k) * n_sub_control_inputs);
      }
    }

    sub_solvers_[i].GenerateDistributedQP(
        qp[i],
        sub_Su[i].template block<p * n_sub_outputs, m * n_sub_control_inputs>(
            0, i * m * n_sub_control_inputs),
        pred.Sx, pred.Sf, delta_x0, n_aug_states, y_subs[i]);

    if (i == 0) integrate_timer.stop();

    output_index_data += n_sub_outputs;
  }

  for (int n_iter = 0; n_iter < n_solver_iterations_; n_iter++) {
    // Re-solve QP with new inputs
    for (int i = 0; i < n_controllers; i++) {
      if (i == 0) integrate_timer.resume();
      sub_solvers_[i].UpdateAndSolveQP(
          qp[i], du_out[i], u_old_.template segment<n_sub_control_inputs>(
                                i * n_sub_control_inputs),
          sub_Su[i], du_prev_);
      if (i == 0) integrate_timer.stop();
    }

    // Update du_prev for next iteration
    for (int i = 0; i < n_controllers; i++) {
      if (i == 0) integrate_timer.resume();
      du_prev_[i] = du_out[i];
      if (i == 0) integrate_timer.stop();
    }
  }

  // Add new du to u_old_
  ControlInput du_applied;
  for (int i = 0; i < n_controllers; i++) {
    if (i == 0) integrate_timer.resume();
    du_applied.template segment<n_sub_control_inputs>(i *
                                                      n_sub_control_inputs) =
        du_prev_[i].template head<n_sub_control_inputs>();

    // Initialize du_prev_ for next iteration using predictions
    if (m > 1) {
      du_prev_[i].template head<(m - 1)* n_sub_control_inputs>() =
          du_prev_[i].template tail<(m - 1) * n_sub_control_inputs>();

      // Repeat last prediction
      du_prev_[i].template tail<n_sub_control_inputs>() =
          du_prev_[i].template segment<n_sub_control_inputs>(
              (m - 2) * n_sub_control_inputs);
    } else {
      du_prev_[i].setZero();
    }
    if (i == 0) integrate_timer.stop();
  }

  observer_.ObserveAPriori(du_applied, u_old_);

  boost::timer::cpu_times int_elapsed = integrate_timer.elapsed();
  boost::timer::nanosecond_type elapsed_ns(int_elapsed.wall);
  cpu_time_out << elapsed_ns << std::endl;

  u_old_ += du_applied;

  return u_old_;
}

#include "noncooperative_controller_list.h"
