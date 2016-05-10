#include "tank.h"

#include "global.h"
#include "valve_eqs.h"

using namespace ValveEqs;

constexpr double pi = 3.14159265358979323846;
constexpr double speed_sound = 340.;

Tank::State Tank::GetDerivative(const State x, const Input u) const {
  State dxdt;

  const double p_d = x(0);
  const double u_tank = u(0);
  const double p_out = u(1);
  const double mass_flow_in = u(2);

  const double m_out =
      CalculateValveMassFlow(p_d, p_out, u_tank, params_.D, params_.m_out_c);

  dxdt << speed_sound *speed_sound / params_.volume *(mass_flow_in - m_out) *
              1e-5;
  return dxdt;
}

Tank::Linearized Tank::GetLinearizedSystem(const State x, const Input u) const {
  Linearized linsys;

  const double p_d = x(0);
  const double u_tank = u(0);
  const double p_out = u(1);
  const double mass_flow_in = u(2);

  linsys.A << CalculateValveDerivative(p_d, p_out, u_tank, params_.D,
                                       params_.volume);
  linsys.C << 1;
  linsys.f = GetDerivative(x, u);
  return linsys;
}

Tank::Params::Params() {
  volume = 20 * pi * (0.60 / 2) * (0.60 / 2) * 2 +
           pi * (0.08 / 2) * (0.08 / 2) * 5.940;
  D = (Vec<8>() << -0.0083454, -0.0094965, 0.16826, -0.032215, -0.61199,
       0.94175, -0.48522, 0.10369).finished();
  m_out_c = 0.017;
}
