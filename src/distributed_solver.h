#ifndef DISTRIBUTED_SOLVER_H
#define DISTRIBUTED_SOLVER_H

#include <Eigen/Eigen>
#include "mpc_qp_solver.h"

template <int n_total_states, int n_outputs, int n_control_inputs, int p, int m,
          int n_controllers>
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
  using AugmentedState =
      typename MpcQpSolver<n_total_states, n_outputs, n_control_inputs, p,
                           m>::AugmentedState;
  using UWeightType = typename MpcQpSolver<n_total_states, n_outputs,
                                           n_control_inputs, p, m>::UWeightType;
  using YWeightType = typename MpcQpSolver<n_total_states, n_outputs,
                                           n_control_inputs, p, m>::YWeightType;
  using QP = typename MpcQpSolver<n_total_states, n_outputs, n_control_inputs,
                                  p, m>::QP;
  using ControlInput =
      typename MpcQpSolver<n_total_states, n_outputs, n_control_inputs, p,
                           m>::ControlInput;
  using ControlInputPrediction =
      typename MpcQpSolver<n_total_states, n_outputs, n_control_inputs, p,
                           m>::ControlInputPrediction;
  using Output = typename MpcQpSolver<n_total_states, n_outputs,
                                      n_control_inputs, p, m>::Output;
  using OutputPrediction =
      typename MpcQpSolver<n_total_states, n_outputs, n_control_inputs, p,
                           m>::OutputPrediction;

  using TotalControlInput =
      Eigen::Matrix<double, n_controllers * m * n_control_inputs, 1>;
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
  void UpdateAndSolveQP(QP* qp, ControlInputPrediction& du_out,
                        const ControlInput u_old, const Eigen::MatrixXd Su_full,
                        const ControlInputPrediction du_full[n_controllers]);
  
  /// Return solution of QP generated based on prediction matrices given
  void UpdateAndSolveQP(QP* qp, ControlInputPrediction* du_out,
                        const ControlInput u_old, const Eigen::MatrixXd Su_other,
                        const Eigen::VectorXd du_other) {
    ApplyOtherInput(qp, Su_other, du_other);
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
  template <typename Derived>
  void ApplyOtherInput(QP* qp, const Eigen::MatrixBase<Derived>& du_other,
                       const Eigen::MatrixXd& Su_other) {
    qp->f += (du_other.transpose() * Su_other.transpose()) * y_pred_weight_;
  }

  Eigen::MatrixXd y_pred_weight_;  // = y_weight_ * Su
  const int index_;                // index of current solver
};

template <int n_total_states, int n_outputs, int n_control_inputs, int p, int m,
          int n_controllers>
void DistributedSolver<n_total_states, n_outputs, n_control_inputs, p, m,
                       n_controllers>::
    UpdateAndSolveQP(QP* qp, ControlInputPrediction& du_out,
                     const ControlInput u_old, const Eigen::MatrixXd Su_full,
                     const ControlInputPrediction du_full[n_controllers]) {
  // Add effect of other subcontrollers' inputs
  for (int i = 0; i < n_controllers; i++) {
    if (i != index_) {
      ApplyOtherInput(
          qp, du_full[i],
          Su_full.template block<p * n_outputs, m * n_control_inputs>(
              0, i * m * n_control_inputs));
    }
  }

  // Solve QP and return solution
  du_out = this->SolveQP(*qp, u_old);
}

#endif
