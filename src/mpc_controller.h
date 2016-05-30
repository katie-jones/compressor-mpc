#ifndef MPC_CONTROLLER_H
#define MPC_CONTROLLER_H

#include <Eigen/Eigen>
#include <Eigen/SparseCore>
#include "dynamic_system.h"
#include "aug_lin_sys.h"
#include "observer.h"
#include "mpc_exceptions.h"
#include "qpOASES.hpp"
#include "controller_interface.h"
#include "input_constraints.h"
#include "prediction.h"
#include "mpc_qp_solver.h"

/**
 * MPC Controller for a dynamic system.
 * Template parameters:
 *    System:         class inheriting from DynamicSystem (or implementing the
 * methods found in this class)
 *    n_delay_states: sum of total number of input delay states required by the
 * plant
 *    n_disturbance_states: number of disturbance states present in observer
 *    p: prediction horizon of controller
 *    m: move horizon of controller
 */
template <class System, int n_delay_states, int n_disturbance_states, int p,
          int m>
class MpcController
    : public ControllerInterface<System, p>,
      public MpcQpSolver<System, n_delay_states, n_disturbance_states, p, m> {
 private:
  static constexpr int n_control_inputs = System::n_control_inputs;
  static constexpr int n_inputs = System::n_inputs;
  static constexpr int n_outputs = System::n_outputs;
  static constexpr int n_states = System::n_states;
  static constexpr int n_aug_states = n_delay_states + n_disturbance_states;
  static constexpr int n_obs_states = n_states + n_disturbance_states;
  static constexpr int n_total_states = n_states + n_aug_states;

  static constexpr int n_wsr_max = 10;  // max working set recalculations

  using ControllerInterface<System, p>::y_ref_;
  using ControllerInterface<System, p>::u_offset_;
  using ControllerInterface<System, p>::n_delay_;
  using ControllerInterface<System, p>::control_input_index_;

  using MpcQpSolver<System, n_delay_states, n_disturbance_states, p, m>::u_old_;
  using MpcQpSolver<System, n_delay_states, n_disturbance_states, p,
                    m>::u_weight_;
  using MpcQpSolver<System, n_delay_states, n_disturbance_states, p,
                    m>::y_weight_;
  using MpcQpSolver<System, n_delay_states, n_disturbance_states, p,
                    m>::u_constraints_;
  using MpcQpSolver<System, n_delay_states, n_disturbance_states, p, m>::Ain_;
  using MpcQpSolver<System, n_delay_states, n_disturbance_states, p,
                    m>::qp_problem_;
  using MpcQpSolver<System, n_delay_states, n_disturbance_states, p,
                    m>::auglinsys_;
  using MpcQpSolver<System, n_delay_states, n_disturbance_states, p,
                    m>::observer_;

 public:
  using State = typename ControllerInterface<System, p>::State;
  using Output = typename ControllerInterface<System, p>::Output;
  using Input = typename ControllerInterface<System, p>::Input;
  using ControlInput = typename ControllerInterface<System, p>::ControlInput;
  using ControlInputIndex =
      typename ControllerInterface<System, p>::ControlInputIndex;
  using OutputPrediction =
      typename ControllerInterface<System, p>::OutputPrediction;
  using AugmentedState =
      typename MpcQpSolver<System, n_delay_states, n_disturbance_states, p,
                           m>::AugmentedState;
  using UWeightType =
      typename MpcQpSolver<System, n_delay_states, n_disturbance_states, p,
                           m>::UWeightType;
  using YWeightType =
      typename MpcQpSolver<System, n_delay_states, n_disturbance_states, p,
                           m>::YWeightType;
  using QP = typename MpcQpSolver<System, n_delay_states, n_disturbance_states,
                                  p, m>::QP;

  /// Constructor -- doesn't initialize state or input/output
  MpcController(
      const AugmentedLinearizedSystem<System, n_delay_states,
                                      n_disturbance_states>& sys,
      const Observer<System, n_delay_states, n_disturbance_states>& observer,
      const Input& u_offset, const OutputPrediction& y_ref,
      const ControlInputIndex& input_delay,
      const ControlInputIndex& control_input_index,
      const UWeightType& u_weight = UWeightType().setIdentity(),
      const YWeightType& y_weight = YWeightType().setIdentity(),
      const InputConstraints<n_control_inputs>& constraints =
          InputConstraints<n_control_inputs>());

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

 protected:
  State x_;       // augmented state
  Output y_old_;  // past output
};

#endif
