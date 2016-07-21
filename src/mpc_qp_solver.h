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
  double objective_function_value_;

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
  MpcQpSolver(const InputConstraints<n_control_inputs>& u_constraints,
              const OutputPrediction y_ref = OutputPrediction::Zero(),
              const UWeightType& u_weight = UWeightType::Identity(),
              const YWeightType& y_weight = YWeightType::Identity());

  /// Initialize QP Problem so hotstart method can be used
  void InitializeQPProblem(const QP& qp, const ControlInput& u_old);

  /// Set QP weights
  void SetWeights(const UWeightType& uwt, const YWeightType& ywt) {
    y_weight_ = Eigen::SparseMatrix<double>(p * n_outputs, p * n_outputs);
    const Eigen::Matrix<double, p * n_outputs, 1> reserve_values =
        Eigen::Matrix<double, p * n_outputs, 1>::Constant(n_outputs);
    y_weight_.reserve(reserve_values);
    for (int i = 0; i < p; i++) {
      for (int j = 0; j < n_outputs * n_outputs; j++) {
        y_weight_.insert(i * n_outputs + j % n_outputs,
                         i * n_outputs + j / n_outputs) =
            ywt(j % n_outputs, j / n_outputs);
      }
    }

    u_weight_.setZero();
    for (int i = 0; i < m; i++) {
      u_weight_.template block<n_control_inputs, n_control_inputs>(
          i * n_control_inputs, i * n_control_inputs) = uwt;
    }
  }

  /// generate QP matrices based on linearization
  QP GenerateQP(const Prediction& pred, const AugmentedState& delta_x0,
                const int n_aug_states, const Output& y_prev) const {
    Eigen::MatrixXd y_pred_weight = y_weight_ * pred.Su;
    return GenerateQP(pred.Su, pred.Sx, pred.Sf, delta_x0, n_aug_states, y_prev,
                      y_pred_weight);
  }

  /// use QPoases to solve QP
  const ControlInputPrediction SolveQP(const QP& qp, const ControlInput& u_old);

  /// Set output reference to use
  void SetOutputReference(const OutputPrediction& y_ref) { y_ref_ = y_ref; }

  /// Get current output reference
  OutputPrediction GetOutputReference() const { return y_ref_; }

  /// Get objective function value
  double GetObjVal() const { return objective_function_value_; }

 protected:
  // generate QP matrices based on linearization, given the y prediction
  // weight
  const QP GenerateQP(const Eigen::MatrixXd& Su, const Eigen::MatrixXd& Sx,
                      const Eigen::MatrixXd& Sf, const AugmentedState& delta_x0,
                      const int n_aug_states, const Output& y_prev,
                      const Eigen::MatrixXd& y_pred_weight) const;

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

  OutputPrediction y_ref_;  // reference trajectory
  Eigen::Matrix<double, m * n_control_inputs, m * n_control_inputs> u_weight_;
  Eigen::SparseMatrix<double> y_weight_;
  const InputConstraints<n_control_inputs>
      u_constraints_;  // input constraints of system
  const std::array<double, m * n_control_inputs * m * n_control_inputs>
      Ain_;                        // matrix used for rate constraints
  qpOASES::SQProblem qp_problem_;  // qp problem to be solved using qpoases
};

#endif
