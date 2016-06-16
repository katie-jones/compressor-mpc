#ifndef TANK_H
#define TANK_H

#include "dynamic_system.h"
#include "global.h"
#include "constexpr_array.h"

class ParallelCompressors;

class Tank : public DynamicSystem<1, 3, 1, ConstexprArray<>> {
  friend ParallelCompressors;

 public:
  constexpr static int n_states = 1;
  constexpr static int n_inputs = 3;
  constexpr static int n_outputs = 1;
  constexpr static int n_control_inputs = 0;
  using ControlInputIndex = ConstexprArray<>;

  typedef DynamicSystem<n_states, n_inputs, n_outputs, ControlInputIndex>::State
      State;
  typedef DynamicSystem<n_states, n_inputs, n_outputs, ControlInputIndex>::Input
      Input;
  typedef DynamicSystem<n_states, n_inputs, n_outputs, ControlInputIndex>::Output
      Output;

  /// Parameters determining dynamics of tank.
  struct Params {
    double volume;
    Vec<8> D;
    double m_out_c;
    Params();
  };

  /// Optionally give parameters to use
  Tank(Params params = Params()) : params_(params) {}

  /// Get derivative of tank.
  virtual State GetDerivative(const State x, const Input u) const;

  /// Linearize system about operating point.
  virtual Linearized GetLinearizedSystem(const State x, const Input u) const;

  /// Get output of tank -- same as state
  virtual Output GetOutput(const State x) const { return x; }

 protected:
  const Params params_;
};

#endif
