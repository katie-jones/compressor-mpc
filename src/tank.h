#ifndef TANK_H
#define TANK_H

#include <Eigen/Eigen>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/algebra/vector_space_algebra.hpp>
#include <iostream>

#include "global_constants.h"

/*
 * Class containing compressor tank
 * Overloaded () operator to give state derivative
 */
class Tank {
 public:
  typedef double TankInput;
  typedef double TankState;

  const static int n_states = 1;
  const static int n_inputs = 1;

  struct Params {
    double pout;
    double volume;
    Vec<8> D;
    double m_out_c;
    Params();
    Params(const Params &x);
  };

  void operator()(const TankState &x_in, TankState &dxdt,
                  const double /* t */) const {
  
    dxdt = GetDerivative(x_in, u);
  }

  TankInput u;
  TankState x;
  Params params;

  Tank(TankState xinit, TankInput uinit) : x(xinit), u(uinit) {}
  Tank() : x(default_initial_state), u(default_input) {}

  TankState GetDerivative(const TankState x,
                                const TankInput u) const;

  const static TankInput default_input;
  const static TankState default_initial_state;

 private:
  typedef boost::numeric::odeint::runge_kutta_dopri5<double> Dopri5Stepper;
  typedef boost::numeric::odeint::controlled_runge_kutta<Dopri5Stepper>
      ControlledStepper;
};


#endif

