#ifndef DISTRIBUTED_SOLVER_H
#define DISTRIBUTED_SOLVER_H

#include <Eigen/Eigen>
#include "mpc_qp_solver.h"

template <int n_total_states, int n_outputs, int n_control_inputs, int p, int m,
          int n_other_controllers>
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

  using OtherControlInput =
      Eigen::Matrix<double, n_other_controllers * n_control_inputs, 1>;

  /// Constructor
  DistributedSolver(const OutputPrediction* y_ref,
                    const ControlInput& u_init = ControlInput::Zero(),
                    const InputConstraints<n_control_inputs>& u_constraints =
                        InputConstraints<n_control_inputs>(),
                    const UWeightType& u_weight = UWeightType::Identity(),
                    const YWeightType& y_weight = YWeightType::Identity());

  /// generate QP matrices based on linearization
  QP GenerateDistributedQP(const Prediction& pred,
                           const AugmentedState& delta_x0,
                           const int n_aug_states, const Output& y_prev) const {
    y_pred_weight_ = y_weight_ * pred.Su;
    return GenerateQP(pred, delta_x0, n_aug_states, y_prev, y_pred_weight_);
  }

  /// Add effect of another systems input on QP
  void ApplyOtherInput(QP& qp, const OtherControlInput& du_other,
                       const Eigen::MatrixXd& Su_other) {
    qp.f += (du_other.transpose() * Su_other.transpose()) * y_pred_weight_;
  }

 private:
  Eigen::MatrixXd y_pred_weight_;  // = y_weight_ * Su
};

#endif
