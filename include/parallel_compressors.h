#ifndef PARALLEL_COMPRESSORS_H
#define PARALLEL_COMPRESSORS_H

#include "compressor.h"
#include "constexpr_array.h"
#include "dynamic_system.h"
#include "tank.h"

/**
 * Dynamic system of 2 compressors arranged in parallel, connected to a common
 * discharge tank.
 */
class ParallelCompressors
    : public virtual DynamicSystem<11, 9, 4, ConstexprArray<0, 3, 4, 7>> {
 public:
  /// Number of system states
  constexpr static int n_states = 11;
  /// Number of system inputs
  constexpr static int n_inputs = 9;
  /// Number of system outputs
  constexpr static int n_outputs = 4;
  /// Number of system control inputs
  constexpr static int n_control_inputs = 4;
  /// Number of compressors in system
  constexpr static int n_compressors = 2;

  /// Indices of control inputs relative to inputs
  using ControlInputIndex = ConstexprArray<0, 3, 4, 7>;

  /// Compressor with input tank
  using Comp = Compressor<true>;

  /// State of system
  typedef DynamicSystem<n_states, n_inputs, n_outputs, ControlInputIndex>::State
      State;
  /// Input to system
  typedef DynamicSystem<n_states, n_inputs, n_outputs, ControlInputIndex>::Input
      Input;
  /// Output of system
  typedef DynamicSystem<n_states, n_inputs, n_outputs,
                        ControlInputIndex>::Output Output;

  /// Constructor with an array of individual compressors
  ParallelCompressors(const double p_in, const double p_out,
                      const Comp comps[n_compressors], const Tank tank = Tank())
      : p_in_(p_in), p_out_(p_out), tank_(tank) {
    for (int i = 0; i < n_compressors; i++) {
      comps_[i] = comps[i];
    }
  }

  /// Constructor using identical compressors
  ParallelCompressors(const double p_in = 1.0, const double p_out = 1.0,
                      const Comp comp = Comp(), const Tank tank = Tank())
      : p_in_(p_in), p_out_(p_out), tank_(tank) {
    for (int i = 0; i < n_compressors; i++) {
      comps_[i] = comp;
    }
  }

  /// Destructor
  virtual ~ParallelCompressors() {}

  /// Get derivative of compressor system about given operating point.
  virtual State GetDerivative(const State& x, const Input& u) const;

  /// Linearize system about operating point.
  virtual Linearized GetLinearizedSystem(const State& x, const Input& u) const;

  /// Return system output at given state.
  virtual Output GetOutput(const State& x) const;

  /// Return default compressor state.
  static const inline State GetDefaultState() {
    const Comp::State x =
        ((Comp::State() << 0.916, 1.145, 0.152, 440, 0).finished());
    return ((State() << x.replicate(n_compressors, 1), 1.12).finished());
  }

  /// Return default compressor input.
  static const inline Input GetDefaultInput() {
    Eigen::Array<double, 4, 1> u;
    u << 0.304, 0.43, 1.0, 0;

    return ((Input() << u.replicate(n_compressors, 1), 0.7).finished());
  }

 protected:
  Comp comps_[n_compressors];
  Tank tank_;
  const double p_in_;
  const double p_out_;
  constexpr static int n_comp_inputs = 4;
  inline Comp::Input GetCompressorInput(Input u_in, int i, State x) const {
    Comp::Input u;
    u << u_in.segment<n_comp_inputs>(i * n_comp_inputs), p_in_, x.tail<1>();
    return u;
  }
};

#endif
