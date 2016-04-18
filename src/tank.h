#ifndef TANK_H
#define TANK_H

#include <Eigen/Eigen>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/algebra/vector_space_algebra.hpp>
#include <iostream>

/*
 * Class containing compressor tank
 * Overloaded () operator to give state derivative
 */
class Tank {
 public:
  const static int n_states = 1;
  const static int n_inputs = 1;

  void operator()(const double &x_in, double &dxdt,
                  const double /* t */) const {
  
    dxdt = GetDerivative(x_in, u);
  }

  typedef double TankInput;
  typedef double TankState;

  TankInput u;
  TankState x;

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

