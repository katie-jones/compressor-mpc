#ifndef MPC_QP_SOLVER_H
#define MPC_QP_SOLVER_H

#include <Eigen/Eigen>
#include "qpOASES.hpp"
#include "prediction.h"
#include "input_constraints.h"
#include "controller_interface.h"

template <class System, int n_delay_states, int n_disturbance_states, int p,
          int m>
class MpcQpSolver {
 private:
  static constexpr int n_control_inputs = System::n_control_inputs;
  static constexpr int n_inputs = System::n_inputs;
  static constexpr int n_outputs = System::n_outputs;
  static constexpr int n_states = System::n_states;
  static constexpr int n_aug_states = n_delay_states + n_disturbance_states;
  static constexpr int n_obs_states = n_states + n_disturbance_states;

  static constexpr int n_wsr_max = 10;  // max working set recalculations

 public:
  /// Number of states in augmented system
  static constexpr int n_total_states = n_aug_states + n_states;

 protected:
  using State = typename ControllerInterface<System, p>::State;
  using Output = typename ControllerInterface<System, p>::Output;
  using Input = typename ControllerInterface<System, p>::Input;
  using ControlInput = typename ControllerInterface<System, p>::ControlInput;
  using ControlInputIndex =
      typename ControllerInterface<System, p>::ControlInputIndex;
  using OutputPrediction =
      typename ControllerInterface<System, p>::OutputPrediction;

  using AugmentedState = Eigen::Matrix<double, n_total_states, 1>;

 public:
  /// Matrix of input weight terms
  typedef Eigen::Matrix<double, n_control_inputs, n_control_inputs> UWeightType;

  /// Matrix of output weight terms
  typedef Eigen::Matrix<double, n_outputs, n_outputs> YWeightType;

  /// Constructor
  MpcQpSolver(const OutputPrediction* y_ref, const ControlInputIndex& n_delay,
              const ControlInput& u_init = ControlInput::Zero(),
              const InputConstraints<n_control_inputs>& u_constraints =
                  InputConstraints<n_control_inputs>(),
              const UWeightType& u_weight = UWeightType::Identity(),
              const YWeightType& y_weight = YWeightType::Identity());

 protected:
  // Structure containing QP problem to solve
  struct QP {
    Eigen::Matrix<double, m * n_control_inputs, m * n_control_inputs,
                  Eigen::RowMajor> H;
    Eigen::Matrix<double, m * n_control_inputs, 1> f;
  };

  // generate QP matrices based on linearization
  const QP GenerateQP(const Prediction& pred,
                      Eigen::Matrix<double, n_aug_states, 1>& delta_x0,
                      const State& xdot, const Output& y_prev) const;

  // use QPoases to solve QP
  const ControlInput SolveQP(const QP& qp);

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
  const ControlInputIndex n_delay_;  // delay states per input
  const InputConstraints<n_control_inputs>
      u_constraints_;  // input constraints of system
  const std::array<double, m * n_control_inputs * m * n_control_inputs>
      Ain_;                        // matrix used for rate constraints
  qpOASES::SQProblem qp_problem_;  // qp problem to be solved using qpoases
};

#endif
