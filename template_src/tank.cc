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

  linsys.A << CalculateValveDerivative(p_d, p_out, u_tank, params_.D);
  linsys.C << 1;
  return linsys;

  // Att = -speed_sound * speed_sound / params_.volume * 1e-5 *
        // (100/2. / sqrt(abs(p_d * 100 - p_out * 100))) *
        // (D2_t(1) * udt ^ 3 + D2_t(2) * udt ^ 2 + D2_t(3) * udt + D2_t(4));
}

Tank::Params::Params() {
  volume = 20 * pi * (0.60 / 2) * (0.60 / 2) * 2 +
           pi * (0.08 / 2) * (0.08 / 2) * 5.940;
  D = (Vec<8>() << -0.0083454, -0.0094965, 0.16826, -0.032215, -0.61199,
       0.94175, -0.48522, 0.10369).finished();
  m_out_c = 0.017;
}
