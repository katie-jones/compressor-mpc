#include "cooperative_controller.h"

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
        const OutputPrediction& y_ref, const ControlInputIndex& input_delay,
        const ControlInputIndex& control_input_index, const Input& u_offset,
        const InputConstraints<n_control_inputs>& constraints,
        const UWeightType& u_weight, const YWeightType& y_weight)
    : ControllerInterface<System, p>(y_ref, u_offset, input_delay,
                                     control_input_index),
      auglinsys_(sys),
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

    sub_solvers_[i] = SubSolver(
        y_ref_, SubControlInput::Zero(), sub_constraints,
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
  u_old_ = u_init.template replicate<m, 1>();

  for (int i = 0; i < n_controllers; i++) {
    sub_solvers_[i].SetInitialInput(
        u_init.template segment<n_control_inputs>(i * n_control_inputs));
  }

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

  const Prediction pred = auglinsys_.GeneratePrediction(p, m);
  Prediction sub_pred;
  sub_pred.Sx = pred.Sx;
  sub_pred.Sf = pred.Sf;

  // Initialize for system 0
  Eigen::MatrixXd Su_other(p * n_outputs,
                           (n_controllers - 1) * n_sub_control_inputs);
  Su_other = pred.Su.template block<p * n_outputs,
                                    (n_controllers - 1) * n_sub_control_inputs>(
      0, n_sub_control_inputs);

  for (int i = 0; i < n_controllers; i++) {
    sub_pred.Su = pred.Su.template block<p * n_outputs, n_sub_control_inputs>(
        0, i * n_sub_control_inputs);
    for (int j = 0; j < i; j++) {
      Su_other.template block<p * n_outputs, n_sub_control_inputs>(
          0, j * n_sub_control_inputs) =
          pred.Su.template block<p * n_outputs,
                                 (n_controllers - 1) * n_sub_control_inputs>(
              0, j * n_sub_control_inputs);
    }

    const QP qp = sub_solvers_[i].GenerateQP(sub_pred, delta_x0, n_aug_states,
                                             observer_.GetPreviousOutput());
    sub_solvers_[i]->InitializeQPProblem(qp);
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
                      n_controllers>::GetNextInput(const Output& y) {
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

  // Generate and solve QPs for each subsystem
  QP qp;
  ControlInputPrediction usol;
  Prediction sub_pred(Eigen::MatrixXd(), pred.Sx, pred.Sf);
  for (int i = 0; i < n_controllers; i++) {
    sub_pred.Su = static_cast<Eigen::MatrixXd>(
        pred.Su.template block<p * n_outputs, n_sub_control_inputs>(
            0, i * n_sub_control_inputs));

    // generate initial QP
    qp = sub_solvers_[i].GenerateDistributedQP(sub_pred, delta_x0, n_aug_states,
                                               y);
    // Add cross terms from other systems
    for (int j = 0; j < n_controllers; j++) {
      if (i != j) {
        sub_solvers_[i].ApplyOtherInput(
            qp, du_prev_.template segment<m * n_sub_control_inputs>(
                    j * m * n_sub_control_inputs),
            static_cast<Eigen::MatrixXd>(
                pred.Su.template block<p * n_outputs, n_sub_control_inputs>(
                    0, i * n_sub_control_inputs)));
      }
    }
    const SubControlInputPrediction sub_usol = sub_solvers_[i].SolveQP(qp);

    // store solution and update u_old_
    usol.template segment<m* n_sub_control_inputs>(
        i * m * n_sub_control_inputs) = sub_usol;
    u_old_.template segment<n_sub_control_inputs>(i * n_sub_control_inputs) =
        sub_usol.template head<n_sub_control_inputs>();
  }

  observer_.ObserveAPriori(u_old_);
  du_prev_ = usol;

  return u_old_;
}
