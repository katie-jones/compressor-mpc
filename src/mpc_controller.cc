#include "mpc_controller.h"

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
        const ControlInputIndex& input_delay,
        const ControlInputIndex& control_input_index,
        const UWeightType& u_weight, const YWeightType& y_weight,
        const InputConstraints& constraints, const Input& u_offset)
    : auglinsys_(sys),
      observer_(observer),
      u_offset_(u_offset),
      Ain_(GetConstraintMatrix()),
      n_delay_(input_delay),
      control_input_index_(control_input_index),
      u_constraints_(constraints),
      qp_problem_(
          qpOASES::SQProblem(m * n_control_inputs, m * n_control_inputs)) {
  int sum_delay = 0;
  for (int i = 0; i < n_control_inputs; i++) sum_delay += n_delay_[i];
  if (sum_delay != n_delay_states) {
    throw delay_states_wrong();
  }
  observer_.InitializeSystem(&auglinsys_);

  y_weight_ = Eigen::SparseMatrix<double>(p * n_outputs, p * n_outputs);
  const Eigen::Matrix<double, p * n_outputs, 1> reserve_values =
      Eigen::Matrix<double, p * n_outputs, 1>::Constant(n_outputs);
  y_weight_.reserve(reserve_values);
  for (int i = 0; i < p; i++) {
    for (int j = 0; j < n_outputs * n_outputs; j++) {
      y_weight_.insert(i * n_outputs + j % n_outputs,
                           i * n_outputs + j / n_outputs) =
          y_weight(j % n_outputs, j / n_outputs);
    }
  }

  u_weight_.setZero();
  for (int i = 0; i < m; i++) {
    u_weight_.template block<n_control_inputs, n_control_inputs>(
        i * n_control_inputs, i * n_control_inputs) = u_weight;
  }
}

/*
 * Generate QP matrices from MPC formulation
 */
template <class System, int n_delay_states, int n_disturbance_states, int p,
          int m>
const typename MpcController<System, n_delay_states, n_disturbance_states, p,
                             m>::QP
MpcController<System, n_delay_states, n_disturbance_states, p, m>::GenerateQP()
    const {
  QP qp;
  const typename AugmentedLinearizedSystem<
      System, n_delay_states, n_disturbance_states>::Prediction pred =
      auglinsys_.GeneratePrediction(p, m);

  const OutputPrediction dy_ref =
      y_ref_ - observer_.GetPreviousOutput().template replicate<p, 1>();

  Eigen::Matrix<double, n_aug_states, 1> delta_x0 =
      observer_.GetStateEstimate().template tail<n_aug_states>();
  // delta_x0.template head<n_states>().setZero();

  int index_delay_states = n_aug_states - n_delay_states;
  for (int i = 0; i < n_control_inputs; i++) {
    if (n_delay_[i] != 0) {
      for (int j = 0; j < n_delay_[i]; j++) {
        delta_x0(index_delay_states + j) -= u_old_[i];
      }
      index_delay_states += n_delay_[i];
    }
  }

  qp.H = pred.Su.transpose() * y_weight_ * pred.Su + u_weight_;

  qp.f =
      auglinsys_.GetF().transpose() * pred.Sf.transpose() * y_weight_ * pred.Su -
      dy_ref.transpose() * y_weight_ * pred.Su +
      delta_x0.transpose() * pred.Sx.transpose() * y_weight_ * pred.Su;

  return qp;
}

/*
 * Solve QP and obtain next input to apply
 */
template <class System, int n_delay_states, int n_disturbance_states, int p,
          int m>
const typename MpcController<System, n_delay_states, n_disturbance_states, p,
                             m>::ControlInput
MpcController<System, n_delay_states, n_disturbance_states, p, m>::SolveQP(
    const QP& qp) {
  int n_wsr = n_wsr_max;

  // Replicate constraints for number of moves
  const Eigen::Matrix<double, m * n_control_inputs, 1> lb =
      (u_constraints_.lower_bound - u_old_).template replicate<m, 1>();
  const Eigen::Matrix<double, m * n_control_inputs, 1> ub =
      (u_constraints_.upper_bound - u_old_).template replicate<m, 1>();
  const Eigen::Matrix<double, m * n_control_inputs, 1> lbA =
      u_constraints_.lower_rate_bound.template replicate<m, 1>();
  const Eigen::Matrix<double, m * n_control_inputs, 1> ubA =
      u_constraints_.upper_rate_bound.template replicate<m, 1>();

  qpOASES::returnValue status =
      qp_problem_.hotstart(qp.H.data(), qp.f.data(), Ain_.data(), lb.data(),
                           ub.data(), lbA.data(), ubA.data(), n_wsr, NULL);

  if (status != qpOASES::SUCCESSFUL_RETURN) {
    // QP not solved, return zeros
    return ControlInput::Zero();
  }

  ControlInput u_solution;
  qp_problem_.getPrimalSolution(u_solution.data());
  return u_solution;
}

/*
 * Calculate a plant input based on control input and offset
 */
template <class System, int n_delay_states, int n_disturbance_states, int p,
          int m>
const typename MpcController<System, n_delay_states, n_disturbance_states, p,
                             m>::Input
MpcController<System, n_delay_states, n_disturbance_states, p,
              m>::GetPlantInput(const ControlInput& u_control) const {
  Input u = u_offset_;
  for (int i = 0; i < n_control_inputs; i++) {
    u(control_input_index_[i]) += u_control(i);
  }
  return u;
}

/*
 * Calculate a control input based on plant input and offset
 */
template <class System, int n_delay_states, int n_disturbance_states, int p,
          int m>
const typename MpcController<System, n_delay_states, n_disturbance_states, p,
                             m>::ControlInput
MpcController<System, n_delay_states, n_disturbance_states, p,
              m>::GetControlInput(const Input& u) const {
  ControlInput u_control;
  for (int i = 0; i < n_control_inputs; i++) {
    u_control(i) =
        u(control_input_index_[i]) - u_offset_(control_input_index_[i]);
  }
  return u_control;
}

/*
 * Solve QP and return optimal input to apply
 */
template <class System, int n_delay_states, int n_disturbance_states, int p,
          int m>
const typename MpcController<System, n_delay_states, n_disturbance_states, p,
                             m>::ControlInput
MpcController<System, n_delay_states, n_disturbance_states, p, m>::GetNextInput(
    const Output& y) {
  x_ += observer_.ObserveAPosteriori(y);
  y_old_ = y;
  auglinsys_.Update(x_, GetPlantInput(u_old_));
  const QP qp = GenerateQP();
  const ControlInput usol = SolveQP(qp);
  observer_.ObserveAPriori(usol);

  u_old_ += usol;
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
                                       const ControlInput& u_init,
                                       const AugmentedState& dx_init) {
  x_ = x_init;
  u_old_ = u_init;
  // y_old_ = sys_.GetOutput(x_init);
  y_old_ = y_init;

  auglinsys_.Update(x_, GetPlantInput(u_old_));
  const QP qp = GenerateQP();

  // Replicate constraints for number of moves
  const Eigen::Matrix<double, m * n_control_inputs, 1> lb =
      (u_constraints_.lower_bound - u_old_).template replicate<m, 1>();
  const Eigen::Matrix<double, m * n_control_inputs, 1> ub =
      (u_constraints_.upper_bound - u_old_).template replicate<m, 1>();
  const Eigen::Matrix<double, m * n_control_inputs, 1> lbA =
      u_constraints_.lower_rate_bound.template replicate<m, 1>();
  const Eigen::Matrix<double, m * n_control_inputs, 1> ubA =
      u_constraints_.upper_rate_bound.template replicate<m, 1>();
  int n_wsr = n_wsr_max;

  qp_problem_.setPrintLevel(qpOASES::PrintLevel::PL_LOW);  // only print errors
  qp_problem_.init(qp.H.data(), qp.f.data(), Ain_.data(), lb.data(), ub.data(),
                   lbA.data(), ubA.data(), n_wsr, NULL);
}

/*
 * Default values for input constraints
 */
template <class System, int n_delay_states, int n_disturbance_states, int p,
          int m>
MpcController<System, n_delay_states, n_disturbance_states, p,
              m>::InputConstraints::InputConstraints()
    : upper_bound(ControlInput::Constant(std::nan(""))),
      lower_bound(ControlInput::Constant(-std::nan(""))),
      use_rate_constraints(false) {}

#include "mpc_controller_list.h"
