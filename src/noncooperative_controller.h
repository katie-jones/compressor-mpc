#ifndef NONCOOPERATIVE_CONTROLLER_H
#define NONCOOPERATIVE_CONTROLLER_H

#include <fstream>

#include "controller_interface.h"
#include "input_constraints.h"
#include "prediction.h"
#include "mpc_qp_solver.h"
#include "distributed_solver.h"

template <class System, int n_delay_states, int n_disturbance_states, int p,
          int m, int n_controllers, int n_sub_outputs>
class NonCooperativeController : public ControllerInterface<System, p> {
 protected:
  static constexpr int n_inputs = System::n_inputs;
  static constexpr int n_control_inputs = System::n_control_inputs;
  static constexpr int n_outputs = System::n_outputs;
  static constexpr int n_states = System::n_states;
  static constexpr int n_aug_states = n_delay_states + n_disturbance_states;
  static constexpr int n_obs_states = n_states + n_disturbance_states;

  static constexpr int n_wsr_max = 10;  // max working set recalculations

  using ControllerInterface<System, p>::y_ref_;
  using ControllerInterface<System, p>::u_offset_;
  using ControllerInterface<System, p>::n_delay_;
  using ControllerInterface<System, p>::control_input_index_;

 public:
  /// Number of states in augmented system
  static constexpr int n_total_states = n_aug_states + n_states;

  using SubSolver =
      DistributedSolver<n_total_states, n_sub_outputs,
                        n_control_inputs / n_controllers, p, m, n_controllers>;
  static constexpr int n_sub_control_inputs = SubSolver::n_control_inputs_;

  using State = typename ControllerInterface<System, p>::State;
  using Output = typename ControllerInterface<System, p>::Output;
  using Input = typename ControllerInterface<System, p>::Input;
  using ControlInput = typename ControllerInterface<System, p>::ControlInput;
  using ControlInputIndex =
      typename ControllerInterface<System, p>::ControlInputIndex;
  using OutputPrediction =
      typename ControllerInterface<System, p>::OutputPrediction;
  using AugmentedState =
      typename MpcQpSolver<n_total_states, n_sub_outputs, n_control_inputs, p,
                           m>::AugmentedState;
  using UWeightType = typename MpcQpSolver<n_total_states, n_sub_outputs,
                                           n_control_inputs, p, m>::UWeightType;
  using YWeightType = typename MpcQpSolver<n_total_states, n_sub_outputs,
                                           n_control_inputs, p, m>::YWeightType;
  using QP = typename MpcQpSolver<n_total_states, n_sub_outputs,
                                  n_sub_control_inputs, p, m>::QP;
  using ControlInputPrediction =
      typename MpcQpSolver<n_total_states, n_sub_outputs, n_control_inputs, p,
                           m>::ControlInputPrediction;
  using SubControlInput = typename SubSolver::ControlInput;
  using SubControlInputPrediction = typename SubSolver::ControlInputPrediction;
  using SubOutput = Eigen::Matrix<double, n_sub_outputs, 1>;
  using SubOutputPrediction = Eigen::Matrix<double, p * n_sub_outputs, 1>;

  /// Constructor
  NonCooperativeController(
      const AugmentedLinearizedSystem<System, n_delay_states,
                                      n_disturbance_states>& sys,
      const Observer<System, n_delay_states, n_disturbance_states>& observer,
      const Input& u_offset, const OutputPrediction& y_ref,
      const ControlInputIndex& input_delay,
      const ControlInputIndex& control_input_index,
      const int n_solver_iterations,
      const Eigen::Array<int, n_controllers, n_sub_outputs> sub_output_index,
      const InputConstraints<n_control_inputs>& constraints =
          InputConstraints<n_control_inputs>(),
      const UWeightType& u_weight = UWeightType::Identity(),
      const YWeightType& y_weight = YWeightType::Identity());

  /// Constructor with different YWeights for each sub-controller
  NonCooperativeController(
      const AugmentedLinearizedSystem<System, n_delay_states,
                                      n_disturbance_states>& sys,
      const Observer<System, n_delay_states, n_disturbance_states>& observer,
      const Input& u_offset, const OutputPrediction& y_ref,
      const ControlInputIndex& input_delay,
      const ControlInputIndex& control_input_index,
      const int n_solver_iterations,
      const Eigen::Array<int, n_controllers, n_sub_outputs> sub_output_index,
      const std::array<YWeightType, n_controllers>& y_weights,
      const InputConstraints<n_control_inputs>& constraints =
          InputConstraints<n_control_inputs>(),
      const UWeightType& u_weight = UWeightType::Identity());

  /**
   * Initialize the state, input and optionally state derivative of the system.
   * The QP problem is also initialized so further solutions can be obtained
   * using the hotstart method.
   */
  void SetInitialState(const State& x_init, const Output& y_init,
                       const ControlInput& u_init = ControlInput::Zero());

  /**
   * Compute the next input to apply to the system.
   * Linearizes the system about current state estimate and finds the optimal
   * input value by solving a QP it generates using the MPC formulation.
   */
  virtual const ControlInput GetNextInput(const Output& y) {
    std::ofstream null;
    null.open("/dev/null");
    return GetNextInput(y, null);
  }

  /**
   * Compute the next input to apply to the system.
   * Print the CPU time required for a single sub-controller to perform the
   * optimization in stream given by cpu_time_out.
   */
  const ControlInput GetNextInput(const Output& y, std::ofstream& cpu_time_out);

  /// Output current state estimate
  const AugmentedState GetStateEstimate() const {
    return (AugmentedState() << x_,
            observer_.GetStateEstimate().template tail<n_aug_states>())
        .finished();
  }

 private:
  // Split output into components for each controller
  void SplitOutput(SubOutput y_subs[], const Output& y_full) {
    for (int i = 0; i < n_controllers; i++) {
      for (int j = 0; j < n_sub_outputs; j++) {
        y_subs[i](j) = y_full(sub_output_index_(i, j));
      }
    }
  }

  // Initialize arguments not set in initializer list - called by constructor
  void InitializeArguments(
      const OutputPrediction& y_ref,
      const std::array<YWeightType, n_controllers>& y_weights,
      const InputConstraints<n_control_inputs>& constraints,
      const UWeightType& u_weight);

  AugmentedLinearizedSystem<System, n_delay_states, n_disturbance_states>
      auglinsys_;  // full auglinsys
  Observer<System, n_delay_states, n_disturbance_states>
      observer_;  // observer of entire state
  SubOutputPrediction
      y_sub_refs_[n_controllers];  // reference output for subsystems
  State x_;                        // current augmented state
  ControlInput u_old_;             // previous applied input
  SubControlInputPrediction du_prev_[n_controllers];  // previous QP solution
  std::vector<SubSolver> sub_solvers_;                // solvers of sub QPs
  const int n_solver_iterations_;  // number of QP iterations to solve using
                                   // updated input values
  Eigen::Array<int, n_controllers, n_sub_outputs, Eigen::RowMajor>
      sub_output_index_;  // indices of sub controller outputs
};

#endif
