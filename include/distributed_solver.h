#ifndef DISTRIBUTED_SOLVER_H
#define DISTRIBUTED_SOLVER_H

#include <Eigen/Eigen>
#include "mpc_qp_solver.h"

/**
 * Specialization of MpcQpSolver for distributed MPC controllers. Adds the
 * option of updating and re-solving a QP based on the solutions of other
 * subcontrollers.
 */
template <int n_total_states, int n_outputs, int n_control_inputs, int p, int m>
class DistributedSolver
    : public MpcQpSolver<n_total_states, n_outputs, n_control_inputs, p, m> {
  using MpcQpSolver<n_total_states, n_outputs, n_control_inputs, p,
                    m>::u_weight_;
  using MpcQpSolver<n_total_states, n_outputs, n_control_inputs, p,
                    m>::y_weight_;
  using MpcQpSolver<n_total_states, n_outputs, n_control_inputs, p,
                    m>::u_constraints_;
  using MpcQpSolver<n_total_states, n_outputs, n_control_inputs, p, m>::Ain_;
  using MpcQpSolver<n_total_states, n_outputs, n_control_inputs, p,
                    m>::qp_problem_;

 public:
  /// Augmented state of system
  using AugmentedState =
      typename MpcQpSolver<n_total_states, n_outputs, n_control_inputs, p,
                           m>::AugmentedState;
  /// Matrix of input weights
  using UWeightType = typename MpcQpSolver<n_total_states, n_outputs,
                                           n_control_inputs, p, m>::UWeightType;
  /// Matrix of output weights
  using YWeightType = typename MpcQpSolver<n_total_states, n_outputs,
                                           n_control_inputs, p, m>::YWeightType;
  /// QP to be solved
  using QP = typename MpcQpSolver<n_total_states, n_outputs, n_control_inputs,
                                  p, m>::QP;
  /// Control input to solve for
  using ControlInput =
      typename MpcQpSolver<n_total_states, n_outputs, n_control_inputs, p,
                           m>::ControlInput;
  /// Prediction of control input over move horizon
  using ControlInputPrediction =
      typename MpcQpSolver<n_total_states, n_outputs, n_control_inputs, p,
                           m>::ControlInputPrediction;
  /// Output of system
  using Output = typename MpcQpSolver<n_total_states, n_outputs,
                                      n_control_inputs, p, m>::Output;
  /// Prediction of output over prediction horizon
  using OutputPrediction =
      typename MpcQpSolver<n_total_states, n_outputs, n_control_inputs, p,
                           m>::OutputPrediction;

  /// Number of control inputs
  static constexpr int n_control_inputs_ = n_control_inputs;

  /// Constructor
  DistributedSolver(const int index,
                    const InputConstraints<n_control_inputs>& u_constraints,
                    const OutputPrediction y_ref = OutputPrediction::Zero(),
                    const UWeightType& u_weight = UWeightType::Identity(),
                    const YWeightType& y_weight = YWeightType::Identity())
      : MpcQpSolver<n_total_states, n_outputs, n_control_inputs, p, m>(
            u_constraints, y_ref, u_weight, y_weight),
        index_(index) {}

  /// Return solution of QP generated based on prediction matrices given
  void UpdateAndSolveQP(QP* qp, ControlInputPrediction* du_out,
                        const ControlInput u_old, const Eigen::MatrixXd Su_full,
                        const double* du_full);

  /// Return solution of QP generated based on prediction matrices given
  void UpdateAndSolveQP(QP* qp, ControlInputPrediction* du_out,
                        const ControlInput u_old,
                        const Eigen::MatrixXd Su_other,
                        const Eigen::VectorXd du_other) {
    ApplyOtherInput(qp, du_other.data(), Su_other);
    *du_out = this->SolveQP(*qp, u_old);
  }

  /// generate QP matrices based on linearization
  void GenerateDistributedQP(QP* qp, const Eigen::MatrixXd& Su,
                             const Eigen::MatrixXd& Sx,
                             const Eigen::MatrixXd& Sf,
                             const AugmentedState& delta_x0,
                             const int n_aug_states, const Output& y_prev) {
    // Update prediction weight
    y_pred_weight_ = (y_weight_ * Su);

    // Get original QP
    *qp = this->GenerateQP(Su, Sx, Sf, delta_x0, n_aug_states, y_prev,
                           y_pred_weight_);
  }

 private:
  /// Add effect of another systems input on QP
  void ApplyOtherInput(QP* qp, const double* du_other,
                       const Eigen::MatrixXd& Su_other) {
    const Eigen::VectorXd du_other_vec =
        Eigen::Map<const Eigen::VectorXd>(du_other, Su_other.cols(), 1);
    qp->f += (du_other_vec.transpose() * Su_other.transpose()) * y_pred_weight_;
  }

  Eigen::MatrixXd y_pred_weight_;  // = y_weight_ * Su
  const int index_;                // index of current solver
};

template <int n_total_states, int n_outputs, int n_control_inputs, int p, int m>
void DistributedSolver<n_total_states, n_outputs, n_control_inputs, p,
                       m>::UpdateAndSolveQP(QP* qp,
                                            ControlInputPrediction* du_out,
                                            const ControlInput u_old,
                                            const Eigen::MatrixXd Su_full,
                                            const double* du_full) {
  // Add effect of other subcontrollers' inputs
  ApplyOtherInput(qp, du_full, Su_full);

  // Solve QP and return solution
  *du_out = this->SolveQP(*qp, u_old);
}

#endif
