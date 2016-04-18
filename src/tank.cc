#include "tank.h"

#include <cmath>
#include <boost/math/special_functions/sign.hpp>

using namespace Eigen;

Tank::TankState Tank::GetDerivative(const TankState x, const TankInput u,
                                    const double mass_flow_compressors) const {
  const double p_d = x(0);
  const double u_tank = u(0);
  const double dp_sqrt2 =
      10 * sqrt(abs(p_d - params.pout)) * boost::math::sign(p_d - params.pout);

  const Vec<8> M5((Vec<8>() << dp_sqrt2 * pow(u_tank, 3),
                   dp_sqrt2 * u_tank * u_tank, dp_sqrt2 * u_tank, dp_sqrt2,
                   pow(u_tank, 3), u_tank * u_tank, u_tank, 1).finished());

  const double m_out = params.D.transpose() * M5 + params.m_out_c;

  TankState dxdt = TankState(speed_sound * speed_sound / params.volume *
                   (mass_flow_compressors - m_out));
  return dxdt;
}

Tank::Params::Params() {
  pout = 1;
  volume = 20 * pi * (0.60 / 2) * (0.60 / 2) * 2 +
           pi * (0.08 / 2) * (0.08 / 2) * 5.940;
  D = (Vec<8>() << -0.0083454, -0.0094965, 0.16826, -0.032215, -0.61199,
       0.94175, -0.48522, 0.10369).finished();
  m_out_c = 0.017;
}

Tank::Params::Params(const Params &x) {
  pout = x.pout;
  volume = x.volume;
  D = x.D;
  m_out_c = x.m_out_c;
}
