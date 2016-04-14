#include "compressor.h"
#include <iostream>
#include <boost/numeric/odeint.hpp>
#include "const_sim.h"

using namespace std;

void write_states(const comp_state &x, const double t) {
  cout << "time: " << t << endl
       << "states: " << x[0] << '\t' << x[1] << '\t' << x[2] << '\t' << x[3]
       << '\t' << x[4] << endl;
}

int main(int argc, char *argv[]) {
  const comp_input uinit(
      (comp_input() << 0.304, 0.405, 0.393, 0, 0).finished());
  comp_state xinit((comp_state() << 0.898, 1.126, 0.15, 439.5, 0.0).finished());
  const double T0 = 0.0;
  const double Tf = 1.0;

  compressor comp(uinit);

  typedef runge_kutta_dopri5<comp_state, double, comp_state, double,
                             vector_space_algebra> dopri5_type;
  typedef controlled_runge_kutta<dopri5_type> ctrl_stepper;

  ctrl_stepper stpr = make_controlled(1e-6, 1e-6, dopri5_type());
  integrate_const(stpr, comp, xinit, T0, Tf, Ts, write_states);

  return 0;
}
