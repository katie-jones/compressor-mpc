#ifndef TANK_H
#define TANK_H

#include "dynamic_system.h"
#include "global.h"

class Tank : public DynamicSystem<1, 3, 1, 0> {
 public:
  constexpr static int n_states = 1;
  constexpr static int n_inputs = 3;
  constexpr static int n_outputs = 1;
  constexpr static int n_control_inputs = 0;

  typedef DynamicSystem<n_states, n_inputs, n_outputs, n_control_inputs>::State
      State;
  typedef DynamicSystem<n_states, n_inputs, n_outputs, n_control_inputs>::Input
      Input;
  typedef DynamicSystem<n_states, n_inputs, n_outputs, n_control_inputs>::Output
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
