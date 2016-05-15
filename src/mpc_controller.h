#ifndef MPC_CONTROLLER_H
#define MPC_CONTROLLER_H

#include <Eigen/Eigen>
#include "dynamic_system.h"
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

  static constexpr int n_aug_states = n_disturbance_states + n_delay_states;
  static constexpr int n_obs_states = n_states + n_disturbance_states;

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
  typedef Eigen::Matrix<double, p * n_outputs, p * n_outputs> YWeightType;

  /// Observer of dynamic system
  typedef Eigen::Matrix<double, n_obs_states, n_outputs> ObserverMatrix;

  /**
   * Compute the next input to apply to the system.
   * Linearizes the system about current state estimate and finds the optimal
   * input value by solving a QP it generates using the MPC formulation.
   */
  const Input GetNextInput(const Output& y);

  /// Set the reference output trajectory
  void SetReference(const OutputPrediction& y_ref) { y_ref_ = y_ref; }

  /// Get the reference output trajectory
  OutputPrediction GetReference() { return y_ref_; }

 protected:
  // Augmented linearized system
  using AugmentedLinearizedSystem =
      typename DynamicSystem<n_total_states, n_inputs, n_outputs,
                             n_control_inputs>::Linearized;

  // Structure containing prediction matrices of augmented system
  struct Prediction {
    Eigen::Matrix<double, p * n_outputs, n_total_states, Eigen::RowMajor> Sx;
    Eigen::Matrix<double, p * n_outputs, m * n_control_inputs, Eigen::RowMajor>
        Su;
    Eigen::Matrix<double, p * n_outputs, n_states, Eigen::RowMajor> Sf;
  };

  // Structure containing QP problem to solve
  struct QP {
    Eigen::Matrix<double, m * n_control_inputs, m * n_control_inputs,
                  Eigen::RowMajor> H;
    Eigen::Matrix<double, m * n_control_inputs, 1> f;
  };

  // apply observer
  void ObserveAPosteriori(const Output& y_in);

  // Discretize system using runge-kutta 4 method
  static const typename System::Linearized DiscretizeRK4(
      const typename System::Linearized& sys_continuous, const double Ts);

  // update linearized/augmented system
  void LinearizeAndAugment();

  // Generate linearized prediction matrices
  const Prediction GeneratePrediction() const;

  // generate QP matrices based on linearization
  const QP GenerateQP() const;

  // use QPoases to solve QP
  const ControlInput SolveQP(const QP& qp) const;

  // Generate state prediction
  void ObserveAPriori(const ControlInput& u_in);

  // output plant input based on control input and offset
  const Input GetPlantInput(const ControlInput& u_control) const;

  // output control input based on plant input and offset
  const ControlInput GetControlInput(const Input& u) const;

  const System sys_;                     // dynamic system
  OutputPrediction y_ref_;               // reference trajectory
  const double sampling_time_;           // sample time
  AugmentedState x_aug_;                 // augmented state
  AugmentedState dx_aug_;                // differential augmented state
  ControlInput u_old_;                   // past input
  Output y_old_;                         // past output
  Input u_offset_;                       // offset applied to control input
  const ObserverMatrix M_;               // observer matrix used
  AugmentedLinearizedSystem auglinsys_;  // current augmented linearization
  const UWeightType u_weight_;           // input weights
  const YWeightType y_weight_;           // output weights
  const std::array<int, n_control_inputs> n_delay_;  // delay states per input
  const ControlInput lower_bound_, upper_bound_, lower_rate_bound_,
      upper_rate_bound_;  // input constraints
  const Eigen::Matrix<double, 2 * m * n_control_inputs, m * n_control_inputs>
      Ain_;  // matrix used for rate constraints
  // index such that ControlInput[i] -> Input[control_input_index_[i]]
  const std::array<int, n_control_inputs> control_input_index_;
  qpOASES::SQProblem qp_problem_;  // qp problem to be solved using qpoase"
};

#endif
