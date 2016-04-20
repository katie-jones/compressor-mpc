#ifndef CHILD_COMPRESSOR_H
#define CHILD_COMPRESSOR_H

#include <Eigen/Eigen>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/algebra/vector_space_algebra.hpp>

#include "global_constants.h"

// Forward declaration for friend class
class ParallelCompressors;

/*
 * Minimal implementation of single compressor
 * Only contains parameters as member variables. Input and state must always be specified.
 */
class Comp {
 public:
  constexpr static int n_states = 5;
  constexpr static int n_inputs = 4;
  constexpr static int n_outputs = 2;

  typedef Eigen::Array<double, n_states, 1> CompressorState;
  typedef Eigen::Array<double, n_inputs, 1> CompressorInput;
  typedef Vec<n_outputs> CompressorOutput;

  struct Coefficients {
    double J, tau_r, m_in_c, m_out_c, torque_drive_c;
    Vec<8> C, D;
    Vec<12> A;
    Vec<2> m_rec_ss_c, SD_c;
    Vec<3> T_ss_c;
    Coefficients(const Coefficients &x);
    Coefficients();
  };

  struct FlowConstants {
    double a, Pin, Pout, V1, V2, AdivL;
    FlowConstants();
    FlowConstants(const FlowConstants &x);
  };

  Comp(Coefficients coeffs = Coefficients(),
       FlowConstants flow_constants = FlowConstants())
      : coeffs(coeffs), flow_constants(flow_constants) {}

  Comp(const Comp &x) : coeffs(x.coeffs), flow_constants(x.flow_constants) {}

  CompressorState GetDerivative(const CompressorState x,
                                const CompressorInput u, const double pout,
                                double &m_out) const;
  inline CompressorState GetDerivative(const CompressorState x,
                                       const CompressorInput u,
                                       const double pout) const {
    double m_out = 0;
    return GetDerivative(x, u, pout, m_out);
  }
  
  CompressorOutput GetOutput(const CompressorState x) const;

  const static inline CompressorState GetDefaultState() {
    return ((CompressorState() << 0.898, 1.126, 0.15, 440, 0).finished());
  }
  const static inline CompressorInput GetDefaultInput() {
    return ((CompressorInput() << 0.304, 0.405, 0.393, 0).finished());
  }


 protected:
  Coefficients coeffs;
  FlowConstants flow_constants;

  friend ParallelCompressors;
};


#endif
