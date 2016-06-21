#ifndef DISTRIBUTED_CONTROLLER_H
#define DISTRIBUTED_CONTROLLER_H

#include <Eigen/Eigen>

#include "distributed_solver.h"
#include "input_constraints.h"
#include "observer.h"
#include "constexpr_array.h"

template <class AugLinSys, typename StateIndices,
          typename ObserverOutputIndices, typename ControlledOutputIndices,
          int p_in, int m_in>
class DistributedController {
 private:
  static constexpr int n_wsr_max = 10;  // max working set recalculations

 public:
  static constexpr int n_delay_states = AugLinSys::n_delay_states;
  static constexpr int n_observer_outputs = AugLinSys::n_outputs;
  static constexpr int n_states = AugLinSys::n_states;
  static constexpr int n_control_inputs = AugLinSys::n_sub_control_inputs;
  static constexpr int n_disturbance_states = AugLinSys::n_disturbance_states;
  static constexpr int n_aug_states = n_delay_states + n_disturbance_states;
  static constexpr int n_obs_states = n_states + n_disturbance_states;
  static constexpr int n_total_states = n_states + n_aug_states;
  static constexpr int n_inputs = AugLinSys::n_inputs;
  static constexpr int n_full_control_inputs = AugLinSys::n_control_inputs;

  static constexpr int n_controlled_outputs = ControlledOutputIndices::size;
  static constexpr int p = p_in;
  static constexpr int m = m_in;

  using State = Eigen::Matrix<double, n_states, 1>;
  using Input = Eigen::Matrix<double, n_inputs, 1>;
  using ControlInput =
      typename MpcQpSolver<n_total_states, n_controlled_outputs,
                           n_control_inputs, p, m>::ControlInput;
  using FullControlInput = typename AugLinSys::ControlInput;
  using ControlOutput =
      typename MpcQpSolver<n_total_states, n_controlled_outputs,
                           n_control_inputs, p, m>::Output;
  using Output = typename Observer<AugLinSys>::Output;
  using AugmentedState =
      typename MpcQpSolver<n_total_states, n_controlled_outputs,
                           n_control_inputs, p, m>::AugmentedState;
  using UWeightType = typename MpcQpSolver<n_total_states, n_controlled_outputs,
                                           n_control_inputs, p, m>::UWeightType;
  using YWeightType = typename MpcQpSolver<n_total_states, n_controlled_outputs,
                                           n_control_inputs, p, m>::YWeightType;
  using QP = typename MpcQpSolver<n_total_states, n_controlled_outputs,
                                  n_control_inputs, p, m>::QP;
  using ControlInputPrediction =
      typename MpcQpSolver<n_total_states, n_controlled_outputs,
                           n_control_inputs, p, m>::ControlInputPrediction;
  using OutputPrediction =
      typename MpcQpSolver<n_total_states, n_controlled_outputs,
                           n_control_inputs, p, m>::OutputPrediction;
  using ObserverMatrix = typename Observer<AugLinSys>::ObserverMatrix;

  using StateIndexType = StateIndices;
  using ObserverOutputIndexType = ObserverOutputIndices;
  using ControlledOutputIndexType = ControlledOutputIndices;
  using ControlInputIndexType =
      typename AugLinSys::ControlInputIndexType::template IndicesSubArray<
          std::make_integer_sequence<int, n_control_inputs>>;

 protected:
  AugLinSys auglinsys_;
  Observer<AugLinSys> observer_;
  DistributedSolver<n_total_states, n_controlled_outputs, n_control_inputs, p,
                    m> qp_solver_;
  State x_;             // current state of system
  ControlInput u_old_;  // previous optimal input to system
  static constexpr typename AugLinSys::DelayType n_delay_ =
      typename AugLinSys::DelayType();  // delay states per input
  Eigen::MatrixXd su_other_;            // effect of other inputs
  QP qp_;

 public:
  /// Constructor
  DistributedController(const AugLinSys& sys,
                        const InputConstraints<n_control_inputs>& constraints,
                        const ObserverMatrix& M);

  /// Set initial output, input and state
  void Initialize(const State& x_init, const ControlInput& u_init,
                  const Input& full_u_old, const Output& y_init,
                  const AugmentedState& dx_init = AugmentedState::Zero());

  /// Set QP weights
  void SetWeights(const UWeightType& uwt, const YWeightType& ywt) {
    qp_solver_.SetWeights(uwt, ywt);
  }

  /// Set reference output
  void SetOutputReference(const OutputPrediction& y_ref) {
    qp_solver_.SetOutputReference(y_ref);
  }

  /// Linearize QP for initial input solution
  void GenerateInitialQP(const Output& y, const Input& full_u_old);

  /// Re-solve QP based on updated inputs from other controllers
  void GetInput(ControlInputPrediction* u_solution,
                const Eigen::VectorXd& du_last);

  /// Return estimate of current state
  State GetStateEstimate() { return x_; }

  /// Return estimate of current state (no memory allocation performed)
  void GetStateEstimate(double* x_out) {
    for (int i = 0; i < n_states; ++i) x_out[i] = x_(i);
  }
};

// Declaration of static constexpr member
template <class AugLinSys, typename StateIndices,
          typename ObserverOutputIndices, typename ControlledOutputIndices,
          int p, int m>
constexpr typename AugLinSys::DelayType
    // constexpr auto
    DistributedController<AugLinSys, StateIndices, ObserverOutputIndices,
                          ControlledOutputIndices, p, m>::n_delay_;

/*
 * Get QP solution based on inputs of other controllers
 */
template <class AugLinSys, typename StateIndices,
          typename ObserverOutputIndices, typename ControlledOutputIndices,
          int p, int m>
void DistributedController<AugLinSys, StateIndices, ObserverOutputIndices,
                           ControlledOutputIndices, p,
                           m>::GetInput(ControlInputPrediction* u_solution,
                                        const Eigen::VectorXd& du_last) {
  QP qp_new = qp_;
  qp_solver_.UpdateAndSolveQP(&qp_new, u_solution, u_old_, su_other_,
                              du_last.data());
}

#endif
