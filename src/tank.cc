#include "tank.h"

using namespace Eigen;

const double pi = 3.14159265358979323846;

const Tank::TankInput Tank::default_input = 0.7;
const Tank::TankState Tank::default_initial_state = 1.12;

Tank::TankState Tank::GetDerivative(const TankState x,
                                    const TankInput u) const {
  TankState dxdt = x;
  return dxdt;
}

Tank::Params::Params() {
  pout = 1;
  volume = 20 * pi * (0.60 / 2) * (0.60 / 2) * 2 +
           pi * (0.08 / 2) * (0.08 / 2) * 5.940;
  D = (Eigen::Matrix<double, 8, 1>() << -0.0083454, -0.0094965, 0.16826, -0.032215,
       -0.61199, 0.94175, -0.48522, 0.10369).finished();
  m_out_c = 0.017;
}

Tank::Params::Params(const Params &x) {
  pout = x.pout;
  volume = x.volume;
  D = x.D;
  m_out_c = x.m_out_c;
}
