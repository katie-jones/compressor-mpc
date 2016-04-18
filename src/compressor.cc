#include "compressor.h"

#include <cmath>
#include <boost/math/special_functions/sign.hpp>

using namespace std;
using namespace boost::numeric::odeint;

typedef boost::numeric::odeint::runge_kutta_dopri5<
    Compressor::CompressorState, double, Compressor::CompressorState, double,
    boost::numeric::odeint::vector_space_algebra> Dopri5Type;
typedef boost::numeric::odeint::controlled_runge_kutta<Dopri5Type>
    ControlledStepper;

const Compressor::CompressorInput Compressor::default_input =
    ((CompressorInput() << 0.304, 0.405, 0.393, 0).finished());

const Compressor::CompressorState Compressor::default_initial_state =
    ((CompressorState() << 0.898, 1.126, 0.15, 440, 0).finished());

double Compressor::GetMassFlowOut() const { return GetMassFlowOut(x); }

double Compressor::GetMassFlowOut(
    const Compressor::CompressorState x_in) const {
  const double p2 = x_in(1);
  const double dp_sqrt2 =
      10 * sqrt(abs(p2 - pout)) * boost::math::sign(p2 - pout);

  const Vec<8> M5((Vec<8>() << dp_sqrt2 * pow(u(2), 3), dp_sqrt2 * u(2) * u(2),
                   dp_sqrt2 * u(2), dp_sqrt2, pow(u(2), 3), u(2) * u(2), u(2),
                   1).finished());

  return coeffs.D.transpose() * M5 + coeffs.m_out_c;
}

Compressor::CompressorState Compressor::GetDerivative(
    const CompressorState x_in, const CompressorInput u,
    const bool flag) const {
  CompressorState dxdt;

  const double p1 = x_in(0);
  const double p2 = x_in(1);
  const double mc = x_in(2);
  const double wc = x_in(3);
  const double mr = x_in(4);

  const double td = u(0) * coeffs.torque_drive_c / wc;
  const double u_in = u(1);
  const double u_out = u(2);
  const double u_rec = u(3);
  const double pout = u(4);

  const double dp_sqrt = 10 * sqrt(abs(flow_constants.Pin - p1)) *
                         boost::math::sign(flow_constants.Pin - p1);

  const Vec<8> M3((Vec<8>() << dp_sqrt * pow(u_in, 3), dp_sqrt * u_in * u_in,
                   dp_sqrt * u_in, dp_sqrt, pow(u_in, 3), u_in * u_in, u_in,
                   1).finished());

  const double m_in = coeffs.C.transpose() * M3 + coeffs.m_in_c;

  double m_rec_ss;
  if (flag)
    m_rec_ss =
        coeffs.m_rec_ss_c.dot(
            (Vec<2>() << sqrt(p2 * 1e5 - p1 * 1e5) * u_rec, 1).finished()) *
        (u_rec > 1e-2);
  else
    m_rec_ss = coeffs.m_rec_ss_c(0) * sqrt(p2 * 1e5 - p1 * 1e5) * u_rec;

  const double m_out = GetMassFlowOut(x_in);

  const double mc2 = mc * mc;
  const double mc3 = mc * mc2;
  const double wc2 = wc * wc;
  const Vec<12> M((Vec<12>() << wc2 * mc3, wc2 * mc2, wc2 * mc, wc2, wc * mc3,
                   wc * mc2, wc * mc, wc, mc3, mc2, mc, 1).finished());

  const double p_ratio = coeffs.A.transpose() * M;

  const double T_ss_model =
      coeffs.T_ss_c(0) + coeffs.T_ss_c(1) * mc + coeffs.T_ss_c(2);

  dxdt(0) = flow_constants.a * flow_constants.a / flow_constants.V1 *
            (m_in + mr - mc) * 1e-5;
  dxdt(1) = flow_constants.a * flow_constants.a / flow_constants.V2 *
            (mc - mr - m_out) * 1e-5;
  dxdt(2) = flow_constants.AdivL * (p_ratio * p1 - p2) * 1e5;
  dxdt(3) = (td - T_ss_model) / coeffs.J;
  dxdt(4) = coeffs.tau_r * (m_rec_ss - mr);

  return dxdt;
}

Compressor::CompressorOutput Compressor::GetOutput() const {
  const double p1 = x(0);
  const double p2 = x(1);
  const double mass_flow = x(2);
  const double surge_distance =
      -(p2 / p1) / coeffs.SD_c(0) - coeffs.SD_c(1) / coeffs.SD_c(0) - mass_flow;
  Compressor::CompressorOutput y;
  y << p2, surge_distance;
  return y;
}

Compressor::Coefficients::Coefficients(const Compressor::Coefficients &x) {
  J = x.J;
  tau_r = x.tau_r;
  m_in_c = x.m_in_c;
  m_out_c = x.m_out_c;
  torque_drive_c = x.torque_drive_c;
  C = x.C;
  D = x.D;
  A = x.A;
  m_rec_ss_c = x.m_rec_ss_c;
  SD_c = x.SD_c;
  T_ss_c = x.T_ss_c;
}

Compressor::Coefficients::Coefficients() {
  J = (0.4 + 0.2070) * 0.4;
  tau_r = 1 / 0.5 + 1;
  A = (Vec<12>() << 0.000299749505193654, -0.000171254191089237,
       3.57321648097597e-05, -9.1783572200945e-07, -0.252701086129365,
       0.136885752773673, -0.02642368327081, 0.00161012740365743,
       54.8046725371143, -29.9550791497765, 5.27827499839098,
       0.693826282579158).finished();

  C = (Vec<8>() << -0.423884232813775, 0.626400271518973, -0.0995040168384753,
       0.0201535563630318, -0.490814924104294, 0.843580880467905,
       -0.423103455111209, 0.0386841406482887).finished();

  D = (Vec<8>() << -0.0083454, -0.0094965, 0.16826, -0.032215, -0.61199,
       0.94175, -0.48522, 0.10369).finished();

  m_in_c = 0.0051;
  m_rec_ss_c = (Vec<2>() << 0.0047, 0.0263).finished();

  m_out_c = 0.017;

  T_ss_c = (Vec<3>() << 2.5543945754982, 47.4222669576423, 0.6218).finished();

  SD_c = (Vec<2>() << 5.55, 0.66).finished();

  torque_drive_c = 15000;
}

Compressor::FlowConstants::FlowConstants() {
  a = 340;
  Pin = 1;
  Pout = 1;
  V1 = 2 * pi * (0.60 / 2) * (0.60 / 2) * 2 +
       pi * (0.08 / 2) * (0.08 / 2) * 8.191;
  V2 = 0.5 * V1;

  AdivL = pi * (0.08 / 2) * (0.08 / 2) / 3 * 0.1;
}

Compressor::FlowConstants::FlowConstants(const FlowConstants &x) {
  a = x.a;
  Pin = x.Pin;
  Pout = x.Pout;
  V1 = x.V1;
  V2 = x.V2;
  AdivL = x.AdivL;
}

