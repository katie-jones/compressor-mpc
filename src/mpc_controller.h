#ifndef MPC_CONTROLLER_H
#define MPC_CONTROLLER_H

#include <Eigen/Eigen>
#include <Eigen/SparseCore>
#include <fstream>

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
template <class System, typename Delays, int n_disturbance_states, int p, int m>
class MpcController
    : public ControllerInterface<System, Delays, p>,
      public MpcQpSolver<System::n_states + Delays::GetSum() +
                             n_disturbance_states,
                         System::n_outputs, System::n_control_inputs, p, m> {
 private:
  static constexpr int n_delay_states = Delays::GetSum();
  static constexpr int n_control_inputs = System::n_control_inputs;
  static constexpr int n_inputs = System::n_inputs;
  static constexpr int n_outputs = System::n_outputs;
  static constexpr int n_states = System::n_states;
  static constexpr int n_aug_states = n_delay_states + n_disturbance_states;
  static constexpr int n_obs_states = n_states + n_disturbance_states;
  static constexpr int n_total_states = n_states + n_aug_states;

  static constexpr int n_wsr_max = 10;  // max working set recalculations

  using ControllerInterface<System, Delays, p>::u_offset_;
  using ControllerInterface<System, Delays, p>::n_delay_;
  using ControllerInterface<System, Delays, p>::control_input_index_;

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
  using AugLinSys =
      AugmentedLinearizedSystem<System, Delays, n_disturbance_states>;
  using State = typename ControllerInterface<System, Delays, p>::State;
  using Output = typename ControllerInterface<System, Delays, p>::Output;
  using Input = typename ControllerInterface<System, Delays, p>::Input;
  using ControlInput =
      typename ControllerInterface<System, Delays, p>::ControlInput;
  using ControlInputIndex =
      typename ControllerInterface<System, Delays, p>::ControlInputIndex;
  using OutputPrediction =
      typename ControllerInterface<System, Delays, p>::OutputPrediction;
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

  /// Constructor -- doesn't initialize state or input/output
  MpcController(const AugLinSys& sys,
                const Observer<AugLinSys>& observer, const Input& u_offset,
                const OutputPrediction& y_ref,
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
  virtual const ControlInput GetNextInput(const Output& y) {
    std::ofstream null;
    null.open("/dev/null");
    ControlInput u = GetNextInput(y, null);
    null.close();
    return u;
  }

  const ControlInput GetNextInput(const Output& y, std::ofstream& cpu_time_out);

  /// Output current state estimate
  const AugmentedState GetStateEstimate() const {
    return (AugmentedState() << x_,
            observer_.GetStateEstimate().template tail<n_aug_states>())
        .finished();
  }

 protected:
  AugLinSys auglinsys_;
  Observer<AugLinSys> observer_;
  State x_;  // augmented state
  ControlInput u_old_;
};

#endif
