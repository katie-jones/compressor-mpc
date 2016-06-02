#ifndef DISTRIBUTED_SOLVER_H
#define DISTRIBUTED_SOLVER_H

#include <Eigen/Eigen>
#include "mpc_qp_solver.h"

template <int n_total_states, int n_outputs, int n_control_inputs, int p, int m,
          int n_controllers>
class DistributedSolver
    : public MpcQpSolver<n_total_states, n_outputs, n_control_inputs, p, m> {
  using MpcQpSolver<n_total_states, n_outputs, n_control_inputs, p, m>::u_old_;
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
  DistributedSolver(const OutputPrediction* p_y_ref, const int index,
                    const ControlInput& u_init = ControlInput::Zero(),
                    const InputConstraints<n_control_inputs>& u_constraints =
                        InputConstraints<n_control_inputs>(),
                    const UWeightType& u_weight = UWeightType::Identity(),
                    const YWeightType& y_weight = YWeightType::Identity())
      : MpcQpSolver<n_total_states, n_outputs, n_control_inputs, p, m>(
            p_y_ref, u_init, u_constraints, u_weight, y_weight),
        index_(index) {}

  /// Return solution of QP generated based on prediction matrices given
  void UpdateAndSolveQP(QP& qp, ControlInputPrediction& du_out,
                        const Eigen::MatrixXd Su[],
                        const TotalControlInput& du_full);

  /// generate QP matrices based on linearization
  void GenerateDistributedQP(QP& qp, const Eigen::MatrixXd& Su,
                             const Eigen::MatrixXd& Sx,
                             const Eigen::MatrixXd& Sf,
                             const AugmentedState& delta_x0,
                             const int n_aug_states, const Output& y_prev) {
    // Update prediction weight
    y_pred_weight_ = (y_weight_ * Su);

    // Get original QP
    qp = this->GenerateQP(Su, Sx, Sf, delta_x0, n_aug_states, y_prev,
                          y_pred_weight_);
  }

 private:
  /// Add effect of another systems input on QP
  void ApplyOtherInput(QP& qp, const ControlInputPrediction& du_other,
                       const Eigen::MatrixXd& Su_other) {
    qp.f += (du_other.transpose() * Su_other.transpose()) * y_pred_weight_;
  }

  Eigen::MatrixXd y_pred_weight_;  // = y_weight_ * Su
  const int index_;                // index of current solver
};

template <int n_total_states, int n_outputs, int n_control_inputs, int p, int m,
          int n_controllers>
void DistributedSolver<
    n_total_states, n_outputs, n_control_inputs, p, m,
    n_controllers>::UpdateAndSolveQP(QP& qp, ControlInputPrediction& du_out,
                                     const Eigen::MatrixXd Su[],
                                     const TotalControlInput& du_full) {
  // Add effect of other subcontrollers' inputs
  for (int i = 0; i < n_controllers; i++) {
    if (i != index_) {
      ApplyOtherInput(qp, du_full.template segment<m * n_control_inputs>(
                              i * m * n_control_inputs),
                      Su[i]);
    }
  }

  // Solve QP and return solution
  du_out = this->SolveQP(qp);
  u_old_ = du_out.template head<n_control_inputs>();
}

#endif
