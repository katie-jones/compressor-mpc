#ifndef DISTRIBUTED_CONTROLLER_H
#define DISTRIBUTED_CONTROLLER_H

#include <Eigen/Eigen>
#include <boost/timer/timer.hpp>

#include "constexpr_array.h"
#include "distributed_solver.h"
#include "input_constraints.h"
#include "observer.h"

/**
 * A sub-controller acting as a component in a distributed MPC solution. Manages
 * a QP solver to obtain the next optimal input at each timestep. Updates and
 * re-solves QP as solutions from other sub-controllers are calculated.
 */
template <class AugLinSys, typename StateIndices,
          typename ObserverOutputIndices, typename ControlledOutputIndices,
          int p_in, int m_in>
class DistributedController {
 private:
  static constexpr int n_wsr_max = 10;  // max working set recalculations
  static constexpr bool is_reduced = AugLinSys::n_sub_control_inputs !=
                                     AugLinSys::SystemType::n_control_inputs;

 public:
  /// Number of delayed states
  static constexpr int n_delay_states = AugLinSys::n_delay_states;
  /// Number of total system outputs
  static constexpr int n_observer_outputs = AugLinSys::n_outputs;
  /// Number of system states
  static constexpr int n_states = AugLinSys::n_states;
  /// Number of control inputs for this subcontroller
  static constexpr int n_control_inputs = AugLinSys::n_sub_control_inputs;
  /// Total number of inputs that are delayed
  static constexpr int n_delayed_inputs = AugLinSys::n_delayed_inputs;
  /// Number of disturbance/integrator states
  static constexpr int n_disturbance_states = AugLinSys::n_disturbance_states;
  /// Number of augmented states
  static constexpr int n_aug_states = n_delay_states + n_disturbance_states;
  /// Number of observable states
  static constexpr int n_obs_states = n_states + n_disturbance_states;
  /// Number of total states
  static constexpr int n_total_states = n_states + n_aug_states;
  /// Number of system inputs
  static constexpr int n_inputs = AugLinSys::n_inputs;
  /// Total number of system control inputs (for all subcontrollers)
  static constexpr int n_full_control_inputs = AugLinSys::n_control_inputs;

  /// Number of outputs controlled by this subcontroller
  static constexpr int n_controlled_outputs = ControlledOutputIndices::size;
  /// Prediction horizon of controller
  static constexpr int p = p_in;
  /// Move horizon of controller
  static constexpr int m = m_in;

  /// State of system
  using State = Eigen::Matrix<double, n_states, 1>;
  /// Input to system
  using Input = Eigen::Matrix<double, n_inputs, 1>;
  /// Control input for this subcontroller
  using ControlInput =
      typename MpcQpSolver<n_total_states, n_controlled_outputs,
                           n_control_inputs, p, m>::ControlInput;
  /// Control input for system
  using FullControlInput = typename AugLinSys::ControlInput;
  /// Vector of outputs controlled by this subcontroller 
  using ControlOutput =
      typename MpcQpSolver<n_total_states, n_controlled_outputs,
                           n_control_inputs, p, m>::Output;
  /// Output of system
  using Output = typename Observer<AugLinSys>::Output;
  /// Augmented state of system
  using AugmentedState =
      typename MpcQpSolver<n_total_states, n_controlled_outputs,
                           n_control_inputs, p, m>::AugmentedState;
  /// Matrix of input weights
  using UWeightType = typename MpcQpSolver<n_total_states, n_controlled_outputs,
                                           n_control_inputs, p, m>::UWeightType;
  /// Matrix of output weights
  using YWeightType = typename MpcQpSolver<n_total_states, n_controlled_outputs,
                                           n_control_inputs, p, m>::YWeightType;
  /// QP used in solver
  using QP = typename MpcQpSolver<n_total_states, n_controlled_outputs,
                                  n_control_inputs, p, m>::QP;
  /// Predicted control inputs over move horizon
  using ControlInputPrediction =
      typename MpcQpSolver<n_total_states, n_controlled_outputs,
                           n_control_inputs, p, m>::ControlInputPrediction;
  /// Predicted outputs over prediction horizon
  using OutputPrediction =
      typename MpcQpSolver<n_total_states, n_controlled_outputs,
                           n_control_inputs, p, m>::OutputPrediction;
  /// Observer matrix
  using ObserverMatrix = typename Observer<AugLinSys>::ObserverMatrix;

  /// Indices of current subcontrollers states relative to system states
  using StateIndexType = StateIndices;
  /// Indices of observed outputs relative to system outputs
  using ObserverOutputIndexType = ObserverOutputIndices;
  /// Indices of outputs controlled by this subcontroller relative to system outputs
  using ControlledOutputIndexType = ControlledOutputIndices;
  /// Indices of control inputs controlled by this subcontroller relative to system inputs
  using ControlInputIndexType =
      typename AugLinSys::ControlInputIndexType::template IndicesSubArray<
          std::make_integer_sequence<int, n_control_inputs>>;
  /// Indices of system control inputs relative to system inputs
  using FullControlInputIndexType = typename AugLinSys::ControlInputIndexType;

 protected:
  AugLinSys auglinsys_;
  Observer<AugLinSys> observer_;
  DistributedSolver<n_total_states, n_controlled_outputs, n_control_inputs, p,
                    m>
      qp_solver_;
  State x_;                 // current state of system
  FullControlInput u_old_;  // previous optimal input to system
  static constexpr typename AugLinSys::DelayType n_delay_ =
      typename AugLinSys::DelayType();  // delay states per input
  Eigen::MatrixXd su_other_;            // effect of other inputs
  QP qp_;
  Prediction pred;  // current value of prediction matrices

 public:
  /// Constructor
  DistributedController(const AugLinSys& sys,
                        const InputConstraints<n_control_inputs>& constraints,
                        const ObserverMatrix& M);

  /// Set initial output, input and state
  void Initialize(const State& x_init, const FullControlInput& u_init,
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

  /// Update solution
  void UpdateU(const FullControlInput& du) {
    // Observe system
    observer_.ObserveAPriori(du, u_old_);

    // Update previous input
    u_old_ += du;
  }

  /// Update solution
  void UpdateU(boost::timer::cpu_timer* time_out, const FullControlInput& du) {
    time_out->resume();
    UpdateU(du);
    time_out->stop();
  }

  /// Linearize QP for initial input solution
  void GenerateInitialQP(const Output& y, const Input& full_u_old);

  /// Linearize QP for initial input solution
  void GenerateInitialQP(boost::timer::cpu_timer* time_out, const Output& y,
                         const Input& full_u_old) {
    time_out->resume();
    GenerateInitialQP(y, full_u_old);
    time_out->stop();
  }

  /// Re-solve QP based on updated inputs from other controllers
  void GetInput(ControlInputPrediction* u_solution,
                const Eigen::VectorXd& du_last);

  /// Re-solve QP based on updated inputs from other controllers
  void GetInput(boost::timer::cpu_timer* time_out,
                ControlInputPrediction* u_solution,
                const Eigen::VectorXd& du_last) {
    time_out->resume();
    GetInput(u_solution, du_last);
    time_out->stop();
  }

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
  if (is_reduced) {
    QP qp_new = qp_;
    qp_solver_.UpdateAndSolveQP(&qp_new, u_solution,
                                u_old_.template head<n_control_inputs>(),
                                su_other_, du_last.data());
  } else {
    *u_solution =
        qp_solver_.SolveQP(qp_, u_old_.template head<n_control_inputs>());
  }

  Output y = observer_.GetPreviousOutput();
  OutputPrediction y_pred_u =
      pred.Su * *u_solution + y.template replicate<p, 1>();
}

#endif
