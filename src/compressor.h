#ifndef COMPRESSOR_H
#define COMPRESSOR_H

#include <Eigen/Eigen>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/algebra/vector_space_algebra.hpp>
#include <iostream>

#include "const_sim.h"
#include "defs.h"

using namespace std;
using namespace boost::numeric::odeint;
using namespace Eigen;

/*
 * Class containing single compressor
 * Overloaded () operator to give state derivative
 */
class compressor {
  comp_input u;
  double t;

 public:
  const bool flag;
  compressor(comp_input u_in, bool flag_in = true) : u(u_in), flag(flag_in) {
    t = 0.0;
  }

  void operator()(const comp_state &x, comp_state &dxdt, const double /* t */);

  void reset(double time) { t = time; }
  // void do_step(
  // friend void integrate(compressor comp, comp_state xinit);
};

// friend void integrate(compressor comp, comp_state xinit) {
// runge_kutta4<comp_state> stepper;
// comp_input u;
// integrate_const(stepper, comp, xinit, 0.0, 10.0, 0.05);
// }

#endif
