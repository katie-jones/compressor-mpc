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
  const static int n_states = 5;
  const static int n_inputs = 5;

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

  typedef Eigen::Array<double, n_states, 1> CompressorState;
  typedef Eigen::Array<double, n_inputs, 1> CompressorInput;
  typedef void (*IntegrationCallbackPtr)(const CompressorState, const double);

  CompressorInput u;
  CompressorState x;
  const Coefficients coeffs;
  const FlowConstants flow_constants;
  Compressor(CompressorState x, CompressorInput u,
             Coefficients coeffs = Coefficients(),
             FlowConstants flow_constants = FlowConstants())
      : x(x), u(u), coeffs(coeffs), flow_constants(flow_constants) {}

  void operator()(const CompressorState &x_in, CompressorState &dxdt,
                  const double /* t */) const {
  
    dxdt = GetDerivative(x_in, u, true);
  }

  CompressorState GetDerivative(const CompressorState x,
                                const CompressorInput u, const bool flag) const;

  const static CompressorInput default_input;
  const static CompressorState default_initial_state;

  friend void IntegrateCompressor(Compressor comp, const double t0,
                                  const double tf, const double dt,
                                  IntegrationCallbackPtr callback,
                                  const double rel_error = 1e-6,
                                  const double abs_error = 1e-6);

 private:
  typedef boost::numeric::odeint::runge_kutta_dopri5<
      CompressorState, double, CompressorState, double,
      boost::numeric::odeint::vector_space_algebra> Dopri5Stepper;
  typedef boost::numeric::odeint::controlled_runge_kutta<Dopri5Stepper>
      ControlledStepper;
};


#endif

