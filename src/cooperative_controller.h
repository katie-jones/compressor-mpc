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
  using QP = typename MpcQpSolver<n_total_states, n_outputs, n_control_inputs,
                                  p, m>::QP;
  using ControlInputPrediction =
      typename MpcQpSolver<n_total_states, n_outputs, n_control_inputs, p,
                           m>::ControlInputPrediction;

  /// Constructor
  CooperativeController(
      const ControlInputIndex& input_delay, const Output& y_init,
      const UWeightType& u_weight = UWeightType().setIdentity(),
      const YWeightType& y_weight = YWeightType().setIdentity(),
      const InputConstraints<n_control_inputs>& constraints =
          InputConstraints<n_control_inputs>(),
      const ControlInput& u_init = ControlInput::Zero());

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
  // Add effect of other inputs to QP
  // void ApplyOtherInputs(QP& qp, const Eigen::VectorXd& du_other,
                        // const Eigen::MatrixXd& Su_other,
                        // const Eigen::MatrixXd& Su) {
    // weight_du_other = Su_other.transpose() * y_weight_ * Su;
    // qp.f += du_other.transpose() * weight_du_other;
  // }

  AugmentedLinearizedSystem<System, n_delay_states, n_disturbance_states>
      auglinsys_;
  Observer<System, n_delay_states, n_disturbance_states> observer_;
  State x_;  // augmented state
  std::array<DistributedSolver<n_total_states, n_outputs,
                               n_control_inputs / n_controllers, p, m,
                               n_controllers - 1>,
             n_controllers> SubSolvers;
};

#endif
