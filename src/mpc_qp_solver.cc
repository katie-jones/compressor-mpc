#include "mpc_qp_solver.h"

template <int n_total_states, int n_outputs, int n_control_inputs, int p, int m>
MpcQpSolver<n_total_states, n_outputs, n_control_inputs, p, m>::MpcQpSolver(
    const OutputPrediction* y_ref, const ControlInput& u_init,
    const InputConstraints<n_control_inputs>& u_constraints,
    const UWeightType& u_weight, const YWeightType& y_weight)
    : p_y_ref_(y_ref),
      u_constraints_(u_constraints),
      Ain_(GetConstraintMatrix()),
      qp_problem_(
          qpOASES::SQProblem(m * n_control_inputs, m * n_control_inputs)) {
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
template <int n_total_states, int n_outputs, int n_control_inputs, int p, int m>
const typename MpcQpSolver<n_total_states, n_outputs, n_control_inputs, p,
                           m>::QP
MpcQpSolver<n_total_states, n_outputs, n_control_inputs, p, m>::GenerateQP(
    const Prediction& pred, const AugmentedState& delta_x0,
    const int n_aug_states, const Output& y_prev) const {
  QP qp;

  const OutputPrediction dy_ref = *p_y_ref_ - y_prev.template replicate<p, 1>();

  qp.H = pred.Su.transpose() * y_weight_ * pred.Su + u_weight_;

  qp.f = delta_x0.head(n_total_states - n_aug_states).transpose() *
             pred.Sf.transpose() * y_weight_ * pred.Su -
         dy_ref.transpose() * y_weight_ * pred.Su +
         delta_x0.tail(n_aug_states).transpose() * pred.Sx.transpose() *
             y_weight_ * pred.Su;

  return qp;
}

/*
 * Generate QP matrices from MPC formulation
 */
template <int n_total_states, int n_outputs, int n_control_inputs, int p, int m>
const typename MpcQpSolver<n_total_states, n_outputs, n_control_inputs, p,
                           m>::QP
MpcQpSolver<n_total_states, n_outputs, n_control_inputs, p, m>::GenerateQP(
    const Prediction& pred, const AugmentedState& delta_x0,
    const int n_aug_states, const Output& y_prev,
    const Eigen::MatrixXd y_pred_weight) const {
  QP qp;

  const OutputPrediction dy_ref = *p_y_ref_ - y_prev.template replicate<p, 1>();

  qp.H = pred.Su.transpose() * y_weight_ * pred.Su + u_weight_;

  qp.f = delta_x0.head(n_total_states - n_aug_states).transpose() *
             pred.Sf.transpose() * y_weight_ * pred.Su -
         dy_ref.transpose() * y_weight_ * pred.Su +
         delta_x0.tail(n_aug_states).transpose() * pred.Sx.transpose() *
             y_weight_ * pred.Su;

  return qp;
}

/*
 * Solve QP and obtain next input to apply
 */
template <int n_total_states, int n_outputs, int n_control_inputs, int p, int m>
const typename MpcQpSolver<n_total_states, n_outputs, n_control_inputs, p,
                           m>::ControlInputPrediction
MpcQpSolver<n_total_states, n_outputs, n_control_inputs, p, m>::SolveQP(
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
    return ControlInputPrediction::Zero();
  }

  ControlInputPrediction u_solution;
  qp_problem_.getPrimalSolution(u_solution.data());

  u_old_ += u_solution.template head<n_control_inputs>();

  return u_solution;
}

/*
 * Initialize QP
 */
template <int n_total_states, int n_outputs, int n_control_inputs, int p, int m>
void MpcQpSolver<n_total_states, n_outputs, n_control_inputs, p,
                 m>::InitializeQPProblem(const QP& qp) {
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

  qp_problem_.setPrintLevel(qpOASES::PrintLevel::PL_LOW);  // only print errors
  qp_problem_.init(qp.H.data(), qp.f.data(), Ain_.data(), lb.data(), ub.data(),
                   lbA.data(), ubA.data(), n_wsr, NULL);
}

#include "mpc_qp_solver_list.h"
