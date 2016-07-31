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
      p_comps_[i + 1] = &comps_[i];
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
  // TODO: generalize for 3+ compressors
  static const inline State GetDefaultState() {
    const Comp::State x =
        ((Comp::State() << 0.916, 1.145, 0.152, 440, 0).finished());
    return ((State() << 0.867, 1.03, 0.176, 395, 0, 0.999, 1.19, 0.176, 395, 0)
                .finished());
  }

  /// Return default compressor input.
  // TODO: generalize for 3+ compressors
  static const inline Input GetDefaultInput() {
    return ((Input() << 0.304, 0.405, 1, 0, 0.304, -1, 0.393, 0).finished());
  }

 protected:
  const FirstComp first_comp_;
  FollowerComp comps_[n_compressors - 1];
  const Comp* p_comps_[n_compressors];
  const double p_in_;
  const double p_out_;
  constexpr static int n_comp_inputs = 4;

  Comp::Input GetCompressorInput(const Input& u_in, const int i,
                                        const State& x) const {
    Comp::Input u;
    double p_out,p_in;

    // Get output pressure depending on compressor number
    if (i==0){
      p_in = p_in_;
    } else {
      p_in = -1;
    }

    if (i == n_compressors - 1) {
      p_out = p_out_;
    } else {
      p_out = x((i + 1) * Comp::n_states);
    }

    // Input with no entry for p_in/m_in
    u << u_in.template segment<n_comp_inputs>(i * n_comp_inputs), p_in, p_out;

    return u;
  }
};

#endif
