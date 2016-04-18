#ifndef COMPRESSOR_SIMULATION_H
#define COMPRESSOR_SIMULATION_H

#include <Eigen/Eigen>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/algebra/vector_space_algebra.hpp>
#include <iostream>

#include "global_constants.h"
#include "compressor.h"

/*
 * Class containing single compressor
 * Overloaded () operator to give state derivative
 */
class CompressorSimulation : public Comp {
 public:
  typedef void (*IntegrationCallbackPtr)(const CompressorState, const double);

  constexpr static double default_pout = 1.0;

  CompressorSimulation(CompressorState x = GetDefaultState(), CompressorInput u = GetDefaultInput(),
                       double pout = default_pout,
                       Coefficients coeffs = Coefficients(),
                       FlowConstants flow_constants = FlowConstants())
      : x(x), u(u), pout(pout), Comp(coeffs, flow_constants) {}

  void operator()(const CompressorState &x_in, CompressorState &dxdt,
                  const double /* t */) const {
    dxdt = GetDerivative(x_in, u, pout);
  }

  inline CompressorOutput GetOutput() { return Comp::GetOutput(x); }

  static inline CompressorState GetDefaultState() {
    return ((CompressorState() << 0.898, 1.126, 0.15, 440, 0).finished());
  }
  static inline CompressorInput GetDefaultInput() {
    return ((CompressorInput() << 0.304, 0.405, 0.393, 0).finished());
  }

  friend void IntegrateCompressor(CompressorSimulation comp, const double t0,
                                  const double tf, const double dt,
                                  IntegrationCallbackPtr callback,
                                  const double rel_error = 1e-6,
                                  const double abs_error = 1e-6) {
    ControlledStepper stepper =
        make_controlled(rel_error, abs_error, Dopri5Stepper());
    integrate_const(stepper, comp, comp.x, t0, tf, dt, callback);
  }

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
struct vector_space_norm_inf<Comp::CompressorState> {
  typedef double result_type;
  double operator()(Comp::CompressorState x) const {
    double absval = 0;
    for (int i = 0; i < Comp::n_states; i++) absval += x[i] * x[i];
    return sqrt(absval);
  }
};
}
}
}

#endif
