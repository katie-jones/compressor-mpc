#ifndef AUGMENTED_LINEARIZED_SYSTEM_H
#define AUGMENTED_LINEARIZED_SYSTEM_H

#include <Eigen/Eigen>
#include <Eigen/SparseCore>

#include "null_index_array.h"
#include "prediction.h"

template <class System>
class Observer;

/**
 * Linearized, discretized dynamic (possibly sub-) system with augmented delayed
 * input states and integrator states.
 * This class controls the linearization and generation of augmented system (A,
 * B, C) matrices at each timestep.
 * It also allows linearized prediction matrices predicting the response of the
 * system over a certain prediction horizon to be generated.
 *
 * The augmented A and B matrices are implemented using structs AComposite and
 * BComposite which optimizes linear algebra operators given prior knowledge of
 * their structure.
 *
 * Template parameters:\n
 * - System: class of original dynamic system\n
 * - Delays: ConstexprArray containing delay state info\n
 * - n_disturbance_states_in: number of total disturbance states\n
 * - ControlInputIndices: indices of control inputs of this sub-system relative
 * to those of full system\n
 * - n_sub_control_inputs_in: number of control inputs belonging to this
 * subsystem\n
 */
template <class System, typename Delays, int n_disturbance_states_in,
          typename ControlInputIndices =
              NullIndexArray<System::n_control_inputs>,
          int n_sub_control_inputs_in = -1>
class AugmentedLinearizedSystem {
  friend Observer<
      AugmentedLinearizedSystem<System, Delays, n_disturbance_states_in,
                                ControlInputIndices, n_sub_control_inputs_in>>;
  // tells whether ControlInputIndices is provided or not
  static constexpr bool is_reduced =
      (n_sub_control_inputs_in != System::n_control_inputs);

 public:
  /// Number of delayed states
  static constexpr int n_delay_states = Delays::GetSum();
  /// Number of control inputs
  static constexpr int n_control_inputs = System::n_control_inputs;
  /// Total number of delayed inputs
  static constexpr int n_delayed_inputs = Delays::GetNonzeroEntries();
  /// Number of system inputs
  static constexpr int n_inputs = System::n_inputs;
  /// Number of system outputs
  static constexpr int n_outputs = System::n_outputs;
  /// Number of system states
  static constexpr int n_states = System::n_states;
  /// Number of disturbance/integrator states
  static constexpr int n_disturbance_states = n_disturbance_states_in;
  /// Number of control inputs for this subsystem
  static constexpr int n_sub_control_inputs =
      is_reduced ? n_sub_control_inputs_in : n_control_inputs;
  /// Number of control inputs controlled by other subsystems
  static constexpr int n_other_control_inputs =
      n_control_inputs - n_sub_control_inputs;

  /// Number of augmented states
  static constexpr int n_aug_states = n_disturbance_states + n_delay_states;
  /// Number of observable states
  static constexpr int n_obs_states = n_states + n_disturbance_states;
  /// Number of total states (augmented + real)
  static constexpr int n_total_states = n_aug_states + n_states;

  /// ConstexprArray of Delays to system
  using DelayType = Delays;
  /// ConstexprArray of indices of current subsystem's control inputs
  using ControlInputIndexType = ControlInputIndices;
  /// (Sub-)System being controlled
  using SystemType = System;

  /// State of the dynamic system
  typedef Eigen::Matrix<double, n_states, 1> State;

  /// Input to dynamic system
  typedef typename System::Input Input;

  /// Augmented state including delay and disturbance states
  typedef Eigen::Matrix<double, n_total_states, 1> AugmentedState;

  /// Controlled inputs of dynamic system
  typedef Eigen::Matrix<double, n_control_inputs, 1> ControlInput;

  /// Controlled inputs for given controller
  typedef Eigen::Matrix<double, n_sub_control_inputs, 1> SubControlInput;

  /// Output of dynamic system
  typedef Eigen::Matrix<double, n_outputs, 1> Output;

  /// Vector of indices for a ControlInput
  typedef std::array<int, n_control_inputs> ControlInputIndex;

  /// Constructor
  AugmentedLinearizedSystem(const System& sys, const double sampling_time);

  /// Re-linearize about given state/input
  void Update(const State x, const Input& u);

  /// Generate prediction matrices based on current linearization
  template <typename ControlledOutputIndices>
  const Prediction GeneratePrediction(Eigen::MatrixXd* Su_other, const int p,
                                      const int m) const {
    Prediction pred;
    GeneratePrediction<ControlledOutputIndices>(&pred.Su, &pred.Sx, &pred.Sf,
                                                Su_other, p, m);
    return pred;
  }

  /// Generate prediction matrices for a subset of outputs
  template <typename ControlledOutputIndices>
  void GeneratePrediction(Eigen::MatrixXd* Su, Eigen::MatrixXd* Sx,
                          Eigen::MatrixXd* Sf, Eigen::MatrixXd* Su_other,
                          const int p, const int m) const;

  /// Return current derivative of system
  const State GetDerivative() const { return f; }

  /// Subtract value of input at linearization from delayed states
  static void AdjustFirstDelayedStates(AugmentedState* x,
                                       const ControlInput& u) {
    int index_delayed_inputs = n_obs_states;
    for (int i = 0; i < n_control_inputs; i++) {
      if (n_delay_[i] != 0) {
        (*x)(index_delayed_inputs) -= u[i];
        index_delayed_inputs++;
      }
    }
  }

  /// Subtract value of input at linearization from all delayed states
  static void AdjustAllDelayedStates(AugmentedState* x, const ControlInput& u) {
    int index_delay_states = n_obs_states + n_delayed_inputs;
    int index_delayed_inputs = n_obs_states;
    for (int i = 0; i < n_control_inputs; i++) {
      if (n_delay_[i] != 0) {
        (*x)(index_delayed_inputs) -= u[i];
        for (int j = 1; j < n_delay_[i]; j++) {
          (*x)(index_delay_states + j - 1) -= u[i];
        }
        index_delay_states += n_delay_[i] - 1;
        index_delayed_inputs++;
      }
    }
  }

  /// Subtract value of input at linearization from all delayed states
  static void AdjustAppliedInput(ControlInput* du, const ControlInput& u) {
    for (int i = 0; i < n_control_inputs; i++) {
      if (n_delay_[i] != 0) {
        (*du)(i) += u[i];
      }
    }
  }

 protected:
  /// Augmented A matrix, separated into components for efficiency
  struct AComposite {
    /// Original (non-augmented) A matrix
    Eigen::Matrix<double, n_states, n_states> Aorig;
    /// Component multiplied by last delayed states
    Eigen::Matrix<double, n_states, n_delayed_inputs> Adelay;
    /// Vector containing locations of 1s in augmented part of A matrix
    Eigen::Matrix<int, n_aug_states, 1> Aaug;

    /// Constructor
    AComposite();

    /// Multiply by an augmented state
    inline AugmentedState operator*(const AugmentedState& x) const;

    /// Multiply by only the augmented (delay and disturbance) components of an
    /// augmented state
    AugmentedState TimesAugmentedOnly(
        const Eigen::Matrix<double, n_aug_states, 1> x) const;

    /// Multiply by a C matrix (or one with same dimensions)
    template <int n_sub_outputs>
    void MultiplyC(
        Eigen::Matrix<double, n_sub_outputs, n_total_states>* C) const;
  };

  /// Augmented B matrix, separated into components for efficiency
  struct BComposite {
    /// Non-delayed components of original B-matrix
    Eigen::Matrix<double, n_states, n_control_inputs - n_delayed_inputs> Borig;
    /// Vector containing location of 1s in augmented part of B matrix
    Eigen::Matrix<int, n_control_inputs, 1> Baug;

    /// Constructor
    BComposite();

    /// Multiply B by ControlInput vector
    AugmentedState operator*(const ControlInput& u) const;

    /// Multiply given C matrix by B
    template <int n_sub_outputs>
    Eigen::Matrix<double, n_sub_outputs, System::n_control_inputs> MultiplyC(
        const Eigen::Matrix<double, n_sub_outputs, n_total_states>& C) const;
  };

  // Discretize system using runge-kutta 4 method
  static const typename System::Linearized DiscretizeRK4(
      const typename System::Linearized& sys_continuous, const double Ts);

  AComposite A;  // Augmented A matrix
  BComposite B;  // Augmented B matrix
  Eigen::Matrix<double, n_outputs, n_obs_states>
      C;              // C matrix with delay states removed (all zeros)
  State f;            // Derivative of system at current operating point
  const System sys_;  // Continuous time system to represent
  const double sampling_time_;  // Sampling time of the discretization
  static constexpr Delays n_delay_ = Delays();
};

// Definition of static constexpr member
template <class System, typename Delays, int n_disturbance_states_in,
          typename ControlInputIndices, int n_sub_control_inputs_in>
constexpr Delays AugmentedLinearizedSystem<
    System, Delays, n_disturbance_states_in, ControlInputIndices,
    n_sub_control_inputs_in>::n_delay_;

/*
 * Operator x_out = A*x
 */
template <class System, typename Delays, int n_disturbance_states_in,
          typename ControlInputIndices, int n_sub_control_inputs_in>
inline
    typename AugmentedLinearizedSystem<System, Delays, n_disturbance_states_in,
                                       ControlInputIndices,
                                       n_sub_control_inputs_in>::AugmentedState
        AugmentedLinearizedSystem<System, Delays, n_disturbance_states_in,
                                  ControlInputIndices,
                                  n_sub_control_inputs_in>::AComposite::
        operator*(const AugmentedState& x) const {
  // multiply by augmented part of matrix
  AugmentedState x_out = TimesAugmentedOnly(x.template tail<n_aug_states>());

  // add contribution of Aorig
  x_out.template head<n_states>() += Aorig * x.template head<n_states>();

  return x_out;
}

namespace Eigen {
namespace internal {
template <>
struct scalar_product_traits<double, bool> {
  enum { Defined = 1 };
  typedef double ReturnType;
};

template <>
struct scalar_product_traits<bool, double> {
  enum { Defined = 1 };
  typedef double ReturnType;
};
}
}

#endif
