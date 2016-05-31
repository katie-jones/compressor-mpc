#ifndef COOPERATIVE_CONTROLLER_H
#define COOPERATIVE_CONTROLLER_H

#include "controller_interface.h"
#include "input_constraints.h"
#include "prediction.h"
#include "mpc_qp_solver.h"
#include "distributed_solver.h"

template <class System, int n_delay_states, int n_disturbance_states, int p,
          int m, int n_controllers>
class CooperativeController : public ControllerInterface<System, p> {
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

  using SubSolver = DistributedSolver<n_total_states, n_outputs,
                                      n_control_inputs / n_controllers, p, m,
                                      n_controllers - 1>;
  static constexpr int n_sub_control_inputs = SubSolver::n_control_inputs;

  using State = typename ControllerInterface<System, p>::State;
  using Output = typename ControllerInterface<System, p>::Output;
  using Input = typename ControllerInterface<System, p>::Input;
  using ControlInput = typename ControllerInterface<System, p>::ControlInput;
  using ControlInputIndex =
      typename ControllerInterface<System, p>::ControlInputIndex;
  using OutputPrediction =
      typename ControllerInterface<System, p>::OutputPrediction;
  using AugmentedState =
      typename MpcQpSolver<n_total_states, n_outputs, n_control_inputs, p,
                           m>::AugmentedState;
  using UWeightType = typename MpcQpSolver<n_total_states, n_outputs,
                                           n_control_inputs, p, m>::UWeightType;
  using YWeightType = typename MpcQpSolver<n_total_states, n_outputs,
                                           n_control_inputs, p, m>::YWeightType;
  using QP = typename MpcQpSolver<n_total_states, n_outputs,
                                  n_sub_control_inputs, p, m>::QP;
  using ControlInputPrediction =
      typename MpcQpSolver<n_total_states, n_outputs, n_control_inputs, p,
                           m>::ControlInputPrediction;
  using SubControlInput = typename SubSolver::ControlInput;
  using SubControlInputPrediction = typename SubSolver::ControlInputPrediction;

  /// Constructor
  CooperativeController(
      const AugmentedLinearizedSystem<System, n_delay_states,
                                      n_disturbance_states>& sys,
      const Observer<System, n_delay_states, n_disturbance_states>& observer,
      const OutputPrediction& y_ref, const ControlInputIndex& input_delay,
      const ControlInputIndex& control_input_index, const Input& u_offset,
      const InputConstraints<n_control_inputs>& constraints =
          InputConstraints<n_control_inputs>(),
      const UWeightType& u_weight = UWeightType::Identity(),
      const YWeightType& y_weight = YWeightType::Identity());

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
  virtual const ControlInput GetNextInput(const Output& y);

  /// Output current state estimate
  const AugmentedState GetStateEstimate() const {
    return (AugmentedState() << x_,
            observer_.GetStateEstimate().template tail<n_aug_states>())
        .finished();
  }

 private:
  AugmentedLinearizedSystem<System, n_delay_states, n_disturbance_states>
      auglinsys_;  // full auglinsys
  Observer<System, n_delay_states, n_disturbance_states>
      observer_;                    // observer of entire state
  State x_;                         // current augmented state
  ControlInput u_old_;              // previous applied input
  ControlInputPrediction du_prev_;  // previous QP solution
  std::array<SubSolver, n_controllers> sub_solvers_;  // solvers of sub QPs
};

#endif
