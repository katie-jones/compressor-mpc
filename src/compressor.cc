#include "child_compressor.h"

Comp::CompressorState Comp::GetDerivative(const CompressorState x,
                                          const CompressorInput u,
                                          const double pout,
                                          double &m_out) {
  CompressorState dxdt;

  const double p1 = x(0);
  const double p2 = x(1);
  const double mc = x(2);
  const double wc = x(3);
  const double mr = x(4);

  const double td = u(0) * coeffs.torque_drive_c / wc;
  const double u_input = u(1);
  const double u_out = u(2);
  const double u_rec = u(3);

  const double dp_sqrt = 10 * sqrt(abs(flow_constants.Pin - p1)) *
                         boost::math::sign(flow_constants.Pin - p1);

  const Vec<8> M3((Vec<8>() << dp_sqrt * pow(u_input, 3),
                   dp_sqrt * u_input * u_input, dp_sqrt * u_input, dp_sqrt,
                   pow(u_input, 3), u_input * u_input, u_input, 1).finished());

  const double m_in = coeffs.C.transpose() * M3 + coeffs.m_in_c;

  double m_rec_ss =
      coeffs.m_rec_ss_c.dot(
          (Vec<2>() << sqrt(p2 * 1e5 - p1 * 1e5) * u_rec, 1).finished()) *
      (u_rec > 1e-2);

  const double dp_sqrt2 =
      10 * sqrt(abs(p2 - pout)) * boost::math::sign(p2 - pout);

  const Vec<8> M5((Vec<8>() << dp_sqrt2 * pow(u(2), 3), dp_sqrt2 * u(2) * u(2),
                   dp_sqrt2 * u(2), dp_sqrt2, pow(u(2), 3), u(2) * u(2), u(2),
                   1).finished());

  m_out = coeffs.D.transpose() * M5 + coeffs.m_out_c;

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

