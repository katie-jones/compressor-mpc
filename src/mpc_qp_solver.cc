#include "mpc_qp_solver.h"

template <int n_total_states, int n_outputs, int n_control_inputs, int p, int m>
MpcQpSolver<n_total_states, n_outputs, n_control_inputs, p, m>::MpcQpSolver(
    const InputConstraints<n_control_inputs>& u_constraints,
    const OutputPrediction y_ref, const UWeightType& u_weight,
    const YWeightType& y_weight)
    : y_ref_(y_ref),
      u_constraints_(u_constraints),
      Ain_(GetConstraintMatrix()),
      qp_problem_(
          qpOASES::SQProblem(m * n_control_inputs, m * n_control_inputs)) {
  SetWeights(u_weight, y_weight);
}

/*
 * Generate QP matrices given y_pred_weight
 */
template <int n_total_states, int n_outputs, int n_control_inputs, int p, int m>
const typename MpcQpSolver<n_total_states, n_outputs, n_control_inputs, p,
                           m>::QP
MpcQpSolver<n_total_states, n_outputs, n_control_inputs, p, m>::GenerateQP(
    const Eigen::MatrixXd& Su, const Eigen::MatrixXd& Sx,
    const Eigen::MatrixXd& Sf, const AugmentedState& delta_x0,
    const int n_aug_states, const Output& y_prev,
    const Eigen::MatrixXd& y_pred_weight) const {
  QP qp;

  const OutputPrediction dy_ref = y_ref_ - y_prev.template replicate<p, 1>();

  qp.H = Su.transpose() * y_pred_weight + u_weight_;

  qp.f =
      delta_x0.head(n_total_states - n_aug_states).transpose() *
          Sf.transpose() * y_pred_weight -
      dy_ref.transpose() * y_pred_weight +
      delta_x0.tail(n_aug_states).transpose() * Sx.transpose() * y_pred_weight;

  return qp;
}

/*
 * Solve QP and obtain next input to apply
 */
template <int n_total_states, int n_outputs, int n_control_inputs, int p, int m>
const typename MpcQpSolver<n_total_states, n_outputs, n_control_inputs, p,
                           m>::ControlInputPrediction
MpcQpSolver<n_total_states, n_outputs, n_control_inputs, p, m>::SolveQP(
    const QP& qp, const ControlInput& u_old) {
  int n_wsr = n_wsr_max;

  // Replicate constraints for number of moves
  const Eigen::Matrix<double, m * n_control_inputs, 1> lb =
      (u_constraints_.lower_bound - u_old).template replicate<m, 1>();
  const Eigen::Matrix<double, m * n_control_inputs, 1> ub =
      (u_constraints_.upper_bound - u_old).template replicate<m, 1>();
  const Eigen::Matrix<double, m * n_control_inputs, 1> lbA =
      u_constraints_.lower_rate_bound.template replicate<m, 1>();
  const Eigen::Matrix<double, m * n_control_inputs, 1> ubA =
      u_constraints_.upper_rate_bound.template replicate<m, 1>();

  qpOASES::returnValue status =
      qp_problem_.hotstart(qp.H.data(), qp.f.data(), Ain_.data(), lb.data(),
                           ub.data(), lbA.data(), ubA.data(), n_wsr, NULL);

  if (status != qpOASES::SUCCESSFUL_RETURN) {
    // QP not solved, return zeros
    return ControlInputPrediction::Zero();
  }

  ControlInputPrediction u_solution;
  qp_problem_.getPrimalSolution(u_solution.data());

  return u_solution;
}

/*
 * Initialize QP
 */
template <int n_total_states, int n_outputs, int n_control_inputs, int p, int m>
void MpcQpSolver<n_total_states, n_outputs, n_control_inputs, p,
                 m>::InitializeQPProblem(const QP& qp,
                                         const ControlInput& u_old) {
  int n_wsr = n_wsr_max;

  // Replicate constraints for number of moves
  const Eigen::Matrix<double, m * n_control_inputs, 1> lb =
      (u_constraints_.lower_bound - u_old).template replicate<m, 1>();
  const Eigen::Matrix<double, m * n_control_inputs, 1> ub =
      (u_constraints_.upper_bound - u_old).template replicate<m, 1>();
  const Eigen::Matrix<double, m * n_control_inputs, 1> lbA =
      u_constraints_.lower_rate_bound.template replicate<m, 1>();
  const Eigen::Matrix<double, m * n_control_inputs, 1> ubA =
      u_constraints_.upper_rate_bound.template replicate<m, 1>();

#ifndef DEBUG
  qp_problem_.setPrintLevel(qpOASES::PrintLevel::PL_LOW);  // only print errors
#endif
  qp_problem_.init(qp.H.data(), qp.f.data(), Ain_.data(), lb.data(), ub.data(),
                   lbA.data(), ubA.data(), n_wsr, NULL);
}

#include "mpc_qp_solver_list.h"
