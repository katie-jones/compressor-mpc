#ifndef COMPRESSOR_H
#define COMPRESSOR_H

#include "dynamic_system.h"
#include "global.h"

class Compressor : public virtual DynamicSystem<5, 6, 2, 2> {
 public:
  static constexpr int n_states = 5;
  static constexpr int n_inputs = 6;
  static constexpr int n_outputs = 2;
  static constexpr int n_control_inputs = 2;

  typedef DynamicSystem<n_states, n_inputs, n_outputs, n_control_inputs>::State
      State;
  typedef DynamicSystem<n_states, n_inputs, n_outputs, n_control_inputs>::Input
      Input;
  typedef DynamicSystem<n_states, n_inputs, n_outputs, n_control_inputs>::Output
      Output;

  /// Coefficients describing dynamics of compressor.
  struct Coefficients {
    double J, tau_r, m_in_c, m_out_c, torque_drive_c, delta_bar, n_bar;
    Vec<8> C, D;
    Vec<12> A;
    Vec<2> m_rec_ss_c, SD_c;
    Vec<3> T_ss_c;
    Coefficients();
  };

  /// Characteristics of fluid flow in compressor.
  struct FlowConstants {
    double a, V1, V2, AdivL;
    FlowConstants();
  };

  /// Optionally specify pressures, coefficients and flow constants.
  Compressor(Coefficients coeffs = Coefficients(),
             FlowConstants flow_constants = FlowConstants())
      : coeffs_(coeffs), flow_constants_(flow_constants) {}

  virtual ~Compressor() {}

  /**
   * Get derivative and mass flow of compressor about given operating point.
   * Return derivative of compressor about state x, input u and pressure at the
   * outlet pout. Also return mass flow through the compressor in variable
   * m_out.
   */
  State GetDerivative(const State x, const Input u, double &m_out) const;

  /**
   * Get derivative of compressor about given operating point.
   * Define a local variable to give as m_out argument to GetDerivative.
   */
  virtual inline State GetDerivative(const State x, const Input u) const {
    double m_out = 0;
    return GetDerivative(x, u, m_out);
  }

  /// Return default compressor state.
  static const inline State GetDefaultState() {
    return ((State() << 0.898, 1.126, 0.15, 440, 0).finished());
  }

  /// Return default compressor input.
  static const inline Input GetDefaultInput() {
    return ((Input() << 0.304, 0.405, 0.393, 0, 1.0, 1.0).finished());
  }

  /// Linearize system about operating point.
  virtual Linearized GetLinearizedSystem(const State x, const Input u) const;

  /// Return system output at given state.
  virtual Output GetOutput(const State x) const;

 protected:
  Coefficients coeffs_;  // parameters determining dynamics of compressor
  FlowConstants flow_constants_;  // flow parameters
};

#endif
