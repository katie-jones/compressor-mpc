#ifndef TANK_H
#define TANK_H

#include "dynamic_system.h"
#include "constexpr_array.h"

class ParallelCompressors;

/**
 * Tank with with an inlet and outlet valve to be connected to compressors.
 */
class Tank : public DynamicSystem<1, 3, 1, ConstexprArray<>> {
  friend ParallelCompressors;

 public:
  /// Number of system states
  constexpr static int n_states = 1;
  /// Number of system inputs
  constexpr static int n_inputs = 3;
  /// Number of system outputs
  constexpr static int n_outputs = 1;
  /// Number of system control inputs
  constexpr static int n_control_inputs = 0;
  /// Indices of control inputs relative to inputs
  using ControlInputIndex = ConstexprArray<>;

  /// System state
  typedef DynamicSystem<n_states, n_inputs, n_outputs, ControlInputIndex>::State
      State;
  /// System input
  typedef DynamicSystem<n_states, n_inputs, n_outputs, ControlInputIndex>::Input
      Input;
  /// System output
  typedef DynamicSystem<n_states, n_inputs, n_outputs, ControlInputIndex>::Output
      Output;

  /// Parameters determining dynamics of tank.
  struct Params {
    double volume;
    Eigen::Matrix<double, 8, 1> D;
    double m_out_c;
    Params();
  };

  /// Constructor with optional parameters to use
  Tank(Params params = Params()) : params_(params) {}

  /// Get derivative of tank.
  virtual State GetDerivative(const State& x, const Input& u) const;

  /// Linearize system about operating point.
  virtual Linearized GetLinearizedSystem(const State& x, const Input& u) const;

  /// Get output of tank -- same as state
  virtual Output GetOutput(const State& x) const { return x; }

 protected:
  const Params params_;
};

#endif
