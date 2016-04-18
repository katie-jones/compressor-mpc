#ifndef COMPRESSOR_H
#define COMPRESSOR_H

#include <Eigen/Eigen>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/algebra/vector_space_algebra.hpp>
#include <iostream>

#include "global_constants.h"

/*
 * Class containing single compressor
 * Overloaded () operator to give state derivative
 */
class Compressor {
 public:
  constexpr static int n_states = 5;
  constexpr static int n_inputs = 4;
  constexpr static int n_outputs = 2;

  typedef Eigen::Array<double, n_states, 1> CompressorState;
  typedef Eigen::Array<double, n_inputs, 1> CompressorInput;
  typedef Vec<n_outputs> CompressorOutput;
  typedef void (*IntegrationCallbackPtr)(const CompressorState, const double);

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

  const static CompressorInput default_input;
  const static CompressorState default_initial_state;
  constexpr static double default_pout = 1.0;

  Compressor(CompressorState x = default_initial_state,
             CompressorInput u = default_input, double pout = default_pout,
             Coefficients coeffs = Coefficients(),
             FlowConstants flow_constants = FlowConstants())
      : x(x), u(u), pout(pout), coeffs(coeffs), flow_constants(flow_constants) {}

  void operator()(const CompressorState &x_in, CompressorState &dxdt,
                  const double /* t */) const {
    dxdt = GetDerivative(x_in, u, true);
  }

  CompressorState GetDerivative(const CompressorState x,
                                const CompressorInput u, const bool flag) const;

  CompressorOutput GetOutput() const;

  double GetMassFlowOut(const CompressorState x_in) const;
  double GetMassFlowOut() const;

  friend void IntegrateCompressor(Compressor comp, const double t0,
                                  const double tf, const double dt,
                                  IntegrationCallbackPtr callback,
                                  const double rel_error = 1e-6,
                                  const double abs_error = 1e-6) {
    ControlledStepper stepper =
        make_controlled(rel_error, abs_error, Compressor::Dopri5Stepper());
    integrate_const(stepper, comp, comp.x, t0, tf, dt, callback);
  }

  const Coefficients coeffs;
  const FlowConstants flow_constants;
  CompressorInput u;
  CompressorState x;
  double pout;

 private:
  typedef boost::numeric::odeint::runge_kutta_dopri5<
      CompressorState, double, CompressorState, double,
      boost::numeric::odeint::vector_space_algebra> Dopri5Stepper;
  typedef boost::numeric::odeint::controlled_runge_kutta<Dopri5Stepper>
      ControlledStepper;
};

// Define norm of Eigen::Array
namespace boost {
namespace numeric {
namespace odeint {
template <>
struct vector_space_norm_inf<Compressor::CompressorState> {
  typedef double result_type;
  double operator()(Compressor::CompressorState x) const {
    double absval = 0;
    for (int i = 0; i < Compressor::n_states; i++) absval += x[i] * x[i];
    return sqrt(absval);
  }
};
}
}
}


#endif

