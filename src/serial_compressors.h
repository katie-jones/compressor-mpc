#ifndef SERIAL_COMPRESSORS_H
#define SERIAL_COMPRESSORS_H

#include "dynamic_system.h"
#include "compressor.h"
#include "constexpr_array.h"
#include "valve_eqs.h"

// TODO: make n_compressors a template argument
class SerialCompressors
    : public virtual DynamicSystem<10, 8, 4, ConstexprArray<0, 3, 4, 7>> {
 public:
  constexpr static int n_states = 10;
  constexpr static int n_inputs = 8;
  constexpr static int n_outputs = 4;
  constexpr static int n_control_inputs = 4;
  constexpr static int n_compressors = 2;
  constexpr static int n_follower_compressors = n_compressors - 1;

  using ControlInputIndex = ConstexprArray<0, 3, 4, 7>;

  /// Compressor with input tank
  using FirstComp = Compressor<true>;
  using FollowerComp = Compressor<false>;
  using Comp = CompressorBase;

  typedef DynamicSystem<n_states, n_inputs, n_outputs, ControlInputIndex>::State
      State;
  typedef DynamicSystem<n_states, n_inputs, n_outputs, ControlInputIndex>::Input
      Input;
  typedef DynamicSystem<n_states, n_inputs, n_outputs,
                        ControlInputIndex>::Output Output;

  SerialCompressors(const double p_in, const double p_out,
                    const FirstComp& first_comp, const FollowerComp comps[])
      : p_in_(p_in), p_out_(p_out), first_comp_(first_comp) {
    p_comps_[0] = &first_comp_;
    for (int i = 0; i < n_follower_compressors; i++) {
      comps_[i] = comps[i];
      p_comps_[i + 1] = &comps_[i];
    }
  }

  SerialCompressors(const double p_in = 1.0, const double p_out = 1.0,
                    const FirstComp& first_comp = FirstComp(),
                    const FollowerComp& comp = FollowerComp())
      : p_in_(p_in), p_out_(p_out), first_comp_(first_comp) {
    p_comps_[0] = &first_comp_;
    for (int i = 0; i < n_follower_compressors; i++) {
      comps_[i] = comp;
      p_comps_[i] = &comps_[i];
    }
  }

  virtual ~SerialCompressors() {}

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
  const FirstComp first_comp_;
  FollowerComp comps_[n_compressors - 1];
  const Comp* p_comps_[n_compressors];
  const double p_in_;
  const double p_out_;
  constexpr static int n_comp_inputs = 4;

  inline Comp::Input GetCompressorInput(Input u_in, int i, State x) const {
    Comp::Input u;
    double p_out;

    // Get output pressure depending on compressor number
    if (i == n_compressors - 1) {
      p_out = p_out_;
    } else {
      p_out = x((i + 1) * Comp::n_states);
    }

    // Input with no entry for p_in/m_in
    u << u_in.template segment<n_comp_inputs>(i * n_comp_inputs), -1, p_out;

    if (i == 1) {
      // Set p_in
      u(n_comp_inputs) = p_in_;
    } else {
      // set m_in
      u(n_comp_inputs) = ValveEqs::CalculateValveMassFlow(
          x((i - 1) * Comp::n_states + 1), x(i * Comp::n_states),
          u((i - 1) * Comp::n_states + 2), p_comps_[i - 1]->params_.D,
          p_comps_[i - 1]->params_.m_out_c);
    }

    return u;
  }
};

#endif
