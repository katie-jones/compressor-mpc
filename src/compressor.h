#ifndef CHILD_COMPRESSOR_H
#define CHILD_COMPRESSOR_H

#include <Eigen/Eigen>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/algebra/vector_space_algebra.hpp>

#include "global_constants.h"

// Forward declaration for friend class
class ParallelCompressors;

/**
 * Minimal implementation of single compressor.
 * Only contains parameters as member variables. Derivative can be calculated
 * but input and state must always be specified.
 */
class Comp {
 public:
  constexpr static int n_states = 5;   ///< Number of states in model.
  constexpr static int n_inputs = 4;   ///< Number of inputs to model.
  constexpr static int n_outputs = 2;  ///< Number of outputs from model.

  /// Type describing state of compressor.
  typedef Eigen::Array<double, n_states, 1> CompressorState;

  /// Type describing inputs to compressor.
  typedef Eigen::Array<double, n_inputs, 1> CompressorInput;

  /// Type describing outputs of compressor.
  typedef Vec<n_outputs> CompressorOutput;

  /// Coefficients describing dynamics of compressor.
  struct Coefficients {
    double J, tau_r, m_in_c, m_out_c, torque_drive_c;
    Vec<8> C, D;
    Vec<12> A;
    Vec<2> m_rec_ss_c, SD_c;
    Vec<3> T_ss_c;
    Coefficients(const Coefficients &x);
    Coefficients();
  };

  /// Characteristics of fluid flow in compressor.
  struct FlowConstants {
    double a, Pin, Pout, V1, V2, AdivL;
    FlowConstants();
    FlowConstants(const FlowConstants &x);
  };

  /// Initialize compressor with coefficients and flow parameters.
  Comp(Coefficients coeffs = Coefficients(),
       FlowConstants flow_constants = FlowConstants())
      : coeffs(coeffs), flow_constants(flow_constants) {}

  /// Copy constructor.
  Comp(const Comp &x) : coeffs(x.coeffs), flow_constants(x.flow_constants) {}

  /**
   * Get derivative and mass flow of compressor about given operating point.
   * Return derivative of compressor about state x, input u and pressure at the
   * outlet pout. Also return mass flow through the compressor in variable
   * m_out.
   */
  CompressorState GetDerivative(const CompressorState x,
                                const CompressorInput u, const double pout,
                                double &m_out) const;

  /**
   * Get derivative of compressor about given operating point.
   * Define a local variable to give as m_out argument to GetDerivative.
   */
  inline CompressorState GetDerivative(const CompressorState x,
                                       const CompressorInput u,
                                       const double pout) const {
    double m_out = 0;
    return GetDerivative(x, u, pout, m_out);
  }

  /// Return output values of compressor for given state.
  CompressorOutput GetOutput(const CompressorState x) const;

  /// Return default compressor state.
  const static inline CompressorState GetDefaultState() {
    return ((CompressorState() << 0.898, 1.126, 0.15, 440, 0).finished());
  }

  /// Return default compressor input.
  const static inline CompressorInput GetDefaultInput() {
    return ((CompressorInput() << 0.304, 0.405, 0.393, 0).finished());
  }

 protected:
  Coefficients coeffs; // parameters determining dynamics of compressor
  FlowConstants flow_constants; // flow parameters

  friend ParallelCompressors;
};

#endif
