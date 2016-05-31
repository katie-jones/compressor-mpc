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

  // Add effect of other inputs to QP
  void ApplyOtherInputs(QP& qp, const OtherControlInput& du_other,
                        const Eigen::MatrixXd& Su_other,
                        const Eigen::MatrixXd& Su) {
    Eigen::MatrixXd weight_du_other(Su_other.cols(), Su.cols());
    weight_du_other = Su_other.transpose() * y_weight_ * Su;
    qp.f += du_other.transpose() * weight_du_other;
  }
};

#endif
