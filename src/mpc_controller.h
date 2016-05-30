#ifndef MPC_CONTROLLER_H
#define MPC_CONTROLLER_H

#include <Eigen/Eigen>
#include <Eigen/SparseCore>
#include "dynamic_system.h"
#include "aug_lin_sys.h"
#include "observer.h"
#include "mpc_exceptions.h"
#include "qpOASES.hpp"

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
class MpcController {
 protected:
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

  /// State of the dynamic system
  typedef Eigen::Matrix<double, n_states, 1> State;

  /// Augmented state including delay and disturbance states
  typedef Eigen::Matrix<double, n_total_states, 1> AugmentedState;

  /// Input to dynamic system
  typedef Eigen::Matrix<double, n_inputs, 1> Input;

  /// Controlled inputs of dynamic system
  typedef Eigen::Matrix<double, n_control_inputs, 1> ControlInput;

  /// Output of dynamic system
  typedef Eigen::Matrix<double, n_outputs, 1> Output;

  /// Output of dynamic system for p time steps
  typedef Eigen::Matrix<double, n_outputs * p, 1> OutputPrediction;

  /// Matrix of input weight terms
  typedef Eigen::Matrix<double, n_control_inputs, n_control_inputs> UWeightType;

  /// Matrix of output weight terms
  typedef Eigen::Matrix<double, n_outputs, n_outputs> YWeightType;

  /// Vector of indices for a ControlInput
  typedef std::array<int, n_control_inputs> ControlInputIndex;

  /// Input constraints controller should respect.
  struct InputConstraints {
    ControlInput lower_bound, upper_bound, lower_rate_bound, upper_rate_bound;
    bool use_rate_constraints;
    InputConstraints();
  };

  /// Constructor -- doesn't initialize state or input/output
  MpcController(
      const AugmentedLinearizedSystem<System, n_delay_states,
                                      n_disturbance_states>& sys,
      const Observer<System, n_delay_states, n_disturbance_states>& observer,
      const ControlInputIndex& input_delay,
      const ControlInputIndex& control_input_index,
      const UWeightType& u_weight = UWeightType().setIdentity(),
      const YWeightType& y_weight = YWeightType().setIdentity(),
      const InputConstraints& constraints = InputConstraints(),
      const Input& u_offset = Input());

  /**
   * Initialize the state, input and optionally state derivative of the system.
   * The QP problem is also initialized so further solutions can be obtained
   * using the hotstart method.
   */
  void SetInitialState(const State& x_init, const Output& y_init,
                       const ControlInput& u_init = ControlInput::Zero(),
                       const AugmentedState& dx_init = AugmentedState::Zero());

  /**
   * Compute the next input to apply to the system.
   * Linearizes the system about current state estimate and finds the optimal
   * input value by solving a QP it generates using the MPC formulation.
   */
  virtual const ControlInput GetNextInput(const Output& y);

  /// Set the reference output trajectory
  void SetReference(const OutputPrediction& y_ref) { y_ref_ = y_ref; }

  /// Get the reference output trajectory
  OutputPrediction GetReference() const { return y_ref_; }

  /// output plant input based on control input and offset
  const Input GetPlantInput(const ControlInput& u_control) const;

  /// output control input based on plant input and offset
  const ControlInput GetControlInput(const Input& u) const;

  /// Output current state estimate
  const AugmentedState GetStateEstimate() const {
    return (AugmentedState() << x_,
            observer_.GetStateEstimate().template tail<n_aug_states>())
        .finished();
  }

 protected:
  // Structure containing QP problem to solve
  struct QP {
    Eigen::Matrix<double, m * n_control_inputs, m * n_control_inputs,
                  Eigen::RowMajor> H;
    Eigen::Matrix<double, m * n_control_inputs, 1> f;
  };

  // generate QP matrices based on linearization
  const QP GenerateQP() const;

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

  Observer<System, n_delay_states, n_disturbance_states> observer_;
  OutputPrediction y_ref_;  // reference trajectory
  State x_;                 // augmented state
  ControlInput u_old_;      // past input
  Output y_old_;            // past output
  const Input u_offset_;    // offset applied to control input
  AugmentedLinearizedSystem<System, n_delay_states, n_disturbance_states>
      auglinsys_;                         // current augmented linearization
  Eigen::Matrix<double, m * n_control_inputs, m * n_control_inputs>
      u_weight_;
  Eigen::SparseMatrix<double> y_weight_;
  const ControlInputIndex n_delay_;       // delay states per input
  const InputConstraints u_constraints_;  // input constraints of system
  const std::array<double, m * n_control_inputs * m * n_control_inputs>
      Ain_;  // matrix used for rate constraints
  // index such that ControlInput[i] -> Input[control_input_index_[i]]
  const ControlInputIndex control_input_index_;
  qpOASES::SQProblem qp_problem_;  // qp problem to be solved using qpoases
};

#endif
