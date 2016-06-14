#ifndef DISTRIBUTED_CONTROLLER_H
#define DISTRIBUTED_CONTROLLER_H

#include <Eigen/Eigen>

#include "distributed_solver.h"
#include "input_constraints.h"
#include "observer.h"
#include "constexpr_array.h"

template <class AugLinSys, int p, int m>
class DistributedController {
 private:
  static constexpr int n_delay_states = AugLinSys::n_delay_states;
  static constexpr int n_outputs = AugLinSys::n_outputs;
  static constexpr int n_states = AugLinSys::n_states;
  static constexpr int n_control_inputs = AugLinSys::n_sub_control_inputs;
  static constexpr int n_disturbance_states = AugLinSys::n_disturbance_states;
  static constexpr int n_aug_states = n_delay_states + n_disturbance_states;
  static constexpr int n_obs_states = n_states + n_disturbance_states;
  static constexpr int n_total_states = n_states + n_aug_states;
  static constexpr int n_inputs = AugLinSys::n_inputs;

  static constexpr int n_wsr_max = 10;  // max working set recalculations

 protected:
  using State = Eigen::Matrix<double, n_states, 1>;
  using Input = Eigen::Matrix<double, n_inputs, 1>;
  using ControlInput =
      typename MpcQpSolver<n_total_states, n_outputs, n_control_inputs, p,
                           m>::ControlInput;
  using Output = typename MpcQpSolver<n_total_states, n_outputs,
                                      n_control_inputs, p, m>::Output;
  using AugmentedState =
      typename MpcQpSolver<n_total_states, n_outputs, n_control_inputs, p,
                           m>::AugmentedState;
  using UWeightType = typename MpcQpSolver<n_total_states, n_outputs,
                                           n_control_inputs, p, m>::UWeightType;
  using YWeightType = typename MpcQpSolver<n_total_states, n_outputs,
                                           n_control_inputs, p, m>::YWeightType;
  using QP = typename MpcQpSolver<n_total_states, n_outputs, n_control_inputs,
                                  p, m>::QP;
  using ControlInputPrediction =
      typename MpcQpSolver<n_total_states, n_outputs, n_control_inputs, p,
                           m>::ControlInputPrediction;
  using OutputPrediction =
      typename MpcQpSolver<n_total_states, n_outputs, n_control_inputs, p,
                           m>::OutputPrediction;
  using ControlInputIndex = typename AugLinSys::ControlInputIndex;
  using ObserverMatrix = typename Observer<AugLinSys>::ObserverMatrix;

  AugLinSys auglinsys_;
  Observer<AugLinSys> observer_;
  DistributedSolver<n_total_states, n_outputs, n_control_inputs, p, m>
      qp_solver_;
  State x_;                 // current state of system
  ControlInput u_old_;      // previous optimal input to system
  OutputPrediction y_ref_;  // Reference output
  static constexpr auto n_delay_ =
      typename AugLinSys::DelayType();  // delay states per input

 public:
  /// Constructor
  DistributedController(const AugLinSys& sys,
                        const InputConstraints<n_control_inputs>& constraints,
                        const ObserverMatrix& M);

  /// Set initial output, input and state
  void Initialize(const Output& y_init, const ControlInput& u_init,
                  const Input& full_u_old, const State& x_init,
                  const AugmentedState& dx_init);

  /// Linearize and solve QP for initial input solution
  QP GenerateInitialQP(const Output& y, const Input& full_u_old);

  /// Re-solve QP based on updated inputs from other controllers
  void GetInput(ControlInputPrediction* u_solution, QP* qp,
                const Eigen::VectorXd& du_last);

  /// Set output reference to track
  void SetReferenceOutput(const OutputPrediction& y_ref) { y_ref_ = y_ref; }

  /// Return estimate of current state
  State GetStateEstimate() { return x_; }

  /// Return estimate of current state (no memory allocation performed)
  void GetStateEstimate(double* x_out) {
    for (int i = 0; i < n_states; ++i) x_out[i] = x_(i);
  }
};

/*
 * Get QP solution based on inputs of other controllers
 */
template <class AugLinSys, int p, int m>
void DistributedController<AugLinSys, p, m>::GetInput(
    ControlInputPrediction* u_solution, QP* qp,
    const Eigen::VectorXd& du_last) {
  Eigen::MatrixXd Su_other = auglinsys_.GetSuOther();

  qp_solver_.UpdateAndSolveQP(qp, *u_solution, u_old_, Su_other,
                              du_last.data());
}

#endif
