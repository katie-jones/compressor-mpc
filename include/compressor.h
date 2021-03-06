#ifndef COMPRESSOR_H
#define COMPRESSOR_H

#include "constexpr_array.h"
#include "dynamic_system.h"

class ParallelCompressors;
class SerialCompressors;

template <bool has_input_tank>
class Compressor;

/**
 * Virtual base class of compressor, combining elements common to compressors
 * with and without an input tank. Inherits from the DynamicSystem class.
 * Defines compressor parameters and provides default values.
 */
class CompressorBase
    : public virtual DynamicSystem<5, 6, 2, ConstexprArray<0, 3>> {
  friend ParallelCompressors;
  friend SerialCompressors;

 public:
  static constexpr int n_states = 5;
  static constexpr int n_inputs = 6;
  static constexpr int n_outputs = 2;
  static constexpr int n_control_inputs = 2;
  using ControlInputIndex = ConstexprArray<0, 3>;

  typedef DynamicSystem<n_states, n_inputs, n_outputs, ControlInputIndex>::State
      State;
  typedef DynamicSystem<n_states, n_inputs, n_outputs, ControlInputIndex>::Input
      Input;
  typedef DynamicSystem<n_states, n_inputs, n_outputs,
                        ControlInputIndex>::Output Output;

  /// Parameters describing dynamics of compressor.
  struct Parameters {
    double J, tau_r, m_in_c, m_out_c, torque_drive_c, delta_bar, n_bar;
    double V1, V2, AdivL, SD_multiplier;
    Eigen::Matrix<double, 8, 1> C, D;
    Eigen::Matrix<double, 12, 1> A;
    Eigen::Matrix<double, 2, 1> m_rec_ss_c, SD_c;
    Eigen::Matrix<double, 3, 1> T_ss_c;
    Parameters();
  };

  /// Optionally specify pressures, coefficients and flow constants.
  CompressorBase(const Parameters& params = Parameters()) : params_(params) {}

  /// Copy constructor
  CompressorBase(const CompressorBase& x) : params_(x.params_) {}

  /// Casting constructor
  template <bool has_input_tank>
  CompressorBase(const Compressor<has_input_tank>& x) : params_(x.params_) {}

  /// Equals operator
  CompressorBase& operator=(const CompressorBase& x) {
    params_ = x.params_;
    return *this;
  }

  /// Equals operator
  template <bool has_input_tank>
  CompressorBase& operator=(const Compressor<has_input_tank>& x) {
    return *this;
  }

  /// Destructor
  virtual ~CompressorBase() {}

  /**
   * Get derivative and mass flow of compressor about given operating point.
   * Return derivative of compressor about state x, input u and pressure at the
   * outlet pout. Also return mass flow through the compressor in variable
   * m_out.
   */
  virtual State GetDerivative(double* m_out, const State& x,
                              const Input& u) const = 0;

  /**
   * Get derivative of compressor about given operating point.
   * Define a local variable to give as m_out argument to GetDerivative.
   */
  virtual inline State GetDerivative(const State& x, const Input& u) const {
    double m_out = 0;
    return this->GetDerivative(&m_out, x, u);
  }

  /// Return default compressor state.
  static const inline State GetDefaultState() {
    return ((State() << 0.898, 1.126, 0.15, 439.5, 0).finished());
  }

  /// Return default compressor input.
  static const inline Input GetDefaultInput() {
    return ((Input() << 0.304, 0.405, 0.393, 0, 1.0, 1.0).finished());
  }

  /// Linearize system about operating point.
  virtual Linearized GetLinearizedSystem(const State& x,
                                         const Input& u) const = 0;

  /// Linearize system about operating point.
  virtual Linearized GetLinearizedSystem(double* m_out, const State& x,
                                         const Input& u) const = 0;

  /// Return system output at given state.
  virtual Output GetOutput(const State& x) const;

 protected:
  Parameters params_;
};

/**
 * Implementation of CompressorBase with input tank or without.
 * Implements all required virtual methods to be dynamic system.
 */
template <bool has_input_tank>
class Compressor : public CompressorBase {
 public:
  /// Type of this object
  using type = Compressor<has_input_tank>;

  using CompressorBase::GetOutput;
  using CompressorBase::GetDerivative;
  using CompressorBase::GetDefaultInput;
  using CompressorBase::GetDefaultState;

  /// Optionally specify pressures, coefficients and flow constants.
  Compressor(const Parameters& params = Parameters())
      : CompressorBase(params) {}

  /// Copy constructor
  Compressor(const CompressorBase& x) : CompressorBase(x) {}

  /// Equals operator
  Compressor& operator=(const Compressor& x) {
    params_ = x.params_;
    return *this;
  }

  /// Equals operator
  Compressor& operator=(const CompressorBase& x) { return *this; }

  /// Destructor
  virtual ~Compressor() {}

  /**
   * Get derivative and mass flow of compressor about given operating point.
   * Return derivative of compressor about state x, input u and pressure at the
   * outlet pout. Also return mass flow through the compressor in variable
   * m_out.
   */
  virtual State GetDerivative(double* m_out, const State& x,
                              const Input& u) const;

  /// Linearize system about operating point.
  virtual Linearized GetLinearizedSystem(const State& x, const Input& u) const {
    double m_out;
    return GetLinearizedSystem(&m_out, x, u);
  }

  /// Linearize system about operating point.
  virtual Linearized GetLinearizedSystem(double* m_out, const State& x,
                                         const Input& u) const;
};

#endif
