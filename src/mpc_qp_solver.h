#ifndef MPC_QP_SOLVER_H
#define MPC_QP_SOLVER_H

#include <Eigen/Eigen>
#include "qpOASES.hpp"
#include "prediction.h"
#include "input_constraints.h"
#include "controller_interface.h"

template <int n_total_states, int n_outputs, int n_control_inputs, int p, int m>
class MpcQpSolver {
 private:
  static constexpr int n_wsr_max = 10;  // max working set recalculations

 public:
  using AugmentedState = Eigen::Matrix<double, n_total_states, 1>;
  using Output = Eigen::Matrix<double, n_outputs, 1>;
  using ControlInput = Eigen::Matrix<double, n_control_inputs, 1>;
  using OutputPrediction = Eigen::Matrix<double, p * n_outputs, 1>;
  using ControlInputPrediction = Eigen::Matrix<double, m * n_control_inputs, 1>;

  /// Matrix of input weight terms
  typedef Eigen::Matrix<double, n_control_inputs, n_control_inputs> UWeightType;

  /// Matrix of output weight terms
  typedef Eigen::Matrix<double, n_outputs, n_outputs> YWeightType;

  /// Structure containing QP problem to solve
  struct QP {
    Eigen::Matrix<double, m * n_control_inputs, m * n_control_inputs,
                  Eigen::RowMajor> H;
    Eigen::Matrix<double, m * n_control_inputs, 1> f;
  };

  /// Constructor
  MpcQpSolver(const OutputPrediction* y_ref,
              const ControlInput& u_init = ControlInput::Zero(),
              const InputConstraints<n_control_inputs>& u_constraints =
                  InputConstraints<n_control_inputs>(),
              const UWeightType& u_weight = UWeightType::Identity(),
              const YWeightType& y_weight = YWeightType::Identity());

  /// Initialize QP Problem so hotstart method can be used
  void InitializeQPProblem(const QP& qp);

  /// generate QP matrices based on linearization
  const QP GenerateQP(const Prediction& pred, const AugmentedState& delta_x0,
                      const int n_aug_states, const Output& y_prev) const;

  /// use QPoases to solve QP
  const ControlInputPrediction SolveQP(const QP& qp);
  
  /// Set initial input
  void SetInitialInput(const ControlInput& u_init) { u_old_ = u_init; }

 protected:
  // output rate constraint matrix
  static const inline std::array<double,
                                 m * n_control_inputs * m * n_control_inputs>
  GetConstraintMatrix() {
    std::array<double, m * n_control_inputs * m * n_control_inputs> A;
    for (int i = 0; i < m * m * n_control_inputs * n_control_inputs; i++)
      A[i] = 0;

    for (int i = 0; i < m * n_control_inputs; i++) {
      A[i + m * n_control_inputs * i] = 1;
      if (i % (m * n_control_inputs) >= n_control_inputs) {
        A[i - n_control_inputs + m * n_control_inputs * i] = -1;
      }
    }
    return A;
  }

  const OutputPrediction* p_y_ref_;  // reference trajectory
  ControlInput u_old_;               // past input
  Eigen::Matrix<double, m * n_control_inputs, m * n_control_inputs> u_weight_;
  Eigen::SparseMatrix<double> y_weight_;
  const InputConstraints<n_control_inputs>
      u_constraints_;  // input constraints of system
  const std::array<double, m * n_control_inputs * m * n_control_inputs>
      Ain_;                        // matrix used for rate constraints
  qpOASES::SQProblem qp_problem_;  // qp problem to be solved using qpoases
};

#endif
