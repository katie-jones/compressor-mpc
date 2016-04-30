#include "compressor.h"

#include <iostream>
#include <cmath>
#include <boost/math/special_functions/sign.hpp>

using namespace Eigen;

// Calculate the mass flow across a valve
// p_in: upstream pressure
// p_out: downstream pressure
// u_valve: valve setting
// C: coefficients of valve map
// m_offset: mass flow offset
double CalculateValveMassFlow(double p_in, double p_out, double u_valve,
                              Vec<8> C, double m_offset);
// Calculate partial derivative wrt pressure across a valve
// p_in: upstream pressure
// p_out: downstream pressure
// u_valve: valve setting
// C: coefficients of valve map
double CalculateValveDerivative(double p_in, double p_out, double u_valve,
                                Vec<8> C);

Compressor::State Compressor::GetDerivative(const State x, const Input u,
                                            double &m_out) const {
  State dxdt;

  const double p1 = x(0);
  const double p2 = x(1);
  const double mc = x(2);
  const double wc = x(3);
  const double mr = x(4);

  const double td = u(0) * coeffs_.torque_drive_c / wc;
  const double u_input = u(1);
  const double u_out = u(2);
  const double u_rec = u(3);
  const double p_in = u(4);
  const double p_out = u(5);

  const double m_in =
      CalculateValveMassFlow(p_in, p1, u_input, coeffs_.C, coeffs_.m_in_c);
  m_out = CalculateValveMassFlow(p2, p_out, u_out, coeffs_.D, coeffs_.m_out_c);

  double m_rec_ss =
      coeffs_.m_rec_ss_c.dot(
          (Vec<2>() << sqrt(p2 * 1e5 - p1 * 1e5) * u_rec, 1).finished()) *
      (u_rec > 1e-2);

  const double mc2 = mc * mc;
  const double mc3 = mc * mc2;
  const double wc2 = wc * wc;
  const Vec<12> M((Vec<12>() << wc2 * mc3, wc2 * mc2, wc2 * mc, wc2, wc * mc3,
                   wc * mc2, wc * mc, wc, mc3, mc2, mc, 1).finished());

  const double p_ratio = coeffs_.A.transpose() * M;

  const double T_ss_model =
      coeffs_.T_ss_c(0) + coeffs_.T_ss_c(1) * mc + coeffs_.T_ss_c(2);

  dxdt(0) = flow_constants_.a * flow_constants_.a / flow_constants_.V1 *
            (m_in + mr - mc) * 1e-5;
  dxdt(1) = flow_constants_.a * flow_constants_.a / flow_constants_.V2 *
            (mc - mr - m_out) * 1e-5;
  dxdt(2) = flow_constants_.AdivL * (p_ratio * p1 - p2) * 1e5;
  dxdt(3) = (td - T_ss_model) / coeffs_.J;
  dxdt(4) = coeffs_.tau_r * (m_rec_ss - mr);

  return dxdt;
}

Compressor::Output Compressor::GetOutput(const State x) const {
  const double p1 = x(0);
  const double p2 = x(1);
  const double mass_flow = x(2);
  const double surge_distance = -(p2 / p1) / coeffs_.SD_c(0) +
                                coeffs_.SD_c(1) / coeffs_.SD_c(0) + mass_flow;
  Compressor::Output y;
  y << p2, surge_distance;
  return y;
}

Compressor::Linearized Compressor::GetLinearizedSystem(const State x,
                                                       const Input u) const {
  Linearized linsys;

  const double p1 = x(0);
  const double p2 = x(1);
  const double mc = x(2);
  const double wc = x(3);
  const double mr = x(4);

  const double td_in = u(0);
  const double u_input = u(1);
  const double u_out = u(2);
  const double u_rec = u(3);
  const double p_in = u(4);
  const double p_out = u(5);

  // Partial derivatives of p1
  linsys.A.row(0) << -speed_sound * speed_sound / flow_constants_.V1 * 1e-5 *
                         CalculateValveDerivative(p_in, p1, u_input, coeffs_.C),
      0, -speed_sound * speed_sound / flow_constants_.V1 * 1e-5, 0,
      speed_sound * speed_sound / flow_constants_.V1 * 1e-5;

  // Partial derivatives of p2
  linsys.A.row(1) << 0,
      -speed_sound * speed_sound / flow_constants_.V2 * 1e-5 *
          CalculateValveDerivative(p2, p_out, u_out, coeffs_.D),
      speed_sound * speed_sound / flow_constants_.V2 * 1e-5, 0,
      -speed_sound * speed_sound / flow_constants_.V2 * 1e-5;

  // Intermediate calculations
  const double wc2 = wc * wc;
  const double wc3 = wc * wc2;
  const double mc2 = mc * mc;
  const double mc3 = mc * mc2;

  Vec<12> M, dM_dmcomp, dM_dwcomp;

  dM_dmcomp << 3 * wc2 *mc2, 2 * wc2 *mc, wc2, 0, 3 * wc *mc2, 2 * wc *mc, wc,
      0, 3 * mc2, 2 * mc, 1, 0;
  dM_dwcomp << 2 * wc *mc3, 2 * wc *mc2, 2 * wc *mc, 2 * wc, mc3, mc2, mc, 1, 0,
      0, 0, 0;

  M << wc2 *mc3, wc2 *mc2, wc2 *mc, wc2, wc *mc3, wc *mc2, wc *mc, wc, mc3, mc2,
      mc, 1;

  double p_ratio = coeffs_.A.dot(M);

  // Partial derivatives of qc
  linsys.A.row(2) << flow_constants_.AdivL * (p_ratio * 1e5),
      -flow_constants_.AdivL * 1e5,
      flow_constants_.AdivL * (p1 * 1e5) * coeffs_.A.dot(dM_dmcomp),
      flow_constants_.AdivL * (p1 * 1e5) * coeffs_.A.dot(dM_dwcomp), 0;

  // Partial derivatives of omegac
  linsys.A.row(3) << 0, 0, -1.0 / coeffs_.J * coeffs_.T_ss_c(0),
      -1.0 / coeffs_.J * td_in * coeffs_.torque_drive_c / wc2, 0;

  // Partial derivatives of qr
  linsys.A.row(4) << -coeffs_.tau_r * (coeffs_.m_rec_ss_c(0) * 1 / 2 * u_rec /
                                       sqrt(p2 * 1e5 - p1 * 1e5) * 1e5),
      coeffs_.tau_r * (coeffs_.m_rec_ss_c(0) * 1 / 2 * u_rec /
                       sqrt(p2 * 1e5 - p1 * 1e5) * 1e5),
      0, 0, -coeffs_.tau_r;

  // B matrix
  // Approximate deadzone using exponentials
  constexpr double x0 = 1e-2;
  double dmr_ur =
      coeffs_.tau_r * coeffs_.m_rec_ss_c(0) * sqrt(p2 * 1e5 - p1 * 2e5);
  if (u_rec < 2 * x0) {
    double a;
    if (u_rec >= x0) {
      a = coeffs_.delta_bar +
          (1 - coeffs_.delta_bar) * exp(coeffs_.n_bar * (u_rec - x0));
    } else {
      a = 2 - (1 - coeffs_.delta_bar) * exp(-coeffs_.n_bar * u_rec);
    }
    dmr_ur = a * dmr_ur;
  }

  linsys.B << 0, 0, 0, 0, 0, 0, 1.0 / coeffs_.J *coeffs_.torque_drive_c / wc, 0,
      0, dmr_ur;

  linsys.C << 0, 1, 0, 0, 0, 100 * p2 / (coeffs_.SD_c(0) * p1 * p1),
      -100. / (coeffs_.SD_c(0) * p1), 100, 0, 0;

  return linsys;
}

double CalculateValveMassFlow(double p_in, double p_out, double u_valve,
                              Vec<8> C, double m_offset) {
  const double dp_sqrt =
      10 * sqrt(abs(p_in - p_out)) * boost::math::sign(p_in - p_out);

  const Vec<8> M3((Vec<8>() << dp_sqrt * pow(u_valve, 3),
                   dp_sqrt * u_valve * u_valve, dp_sqrt * u_valve, dp_sqrt,
                   pow(u_valve, 3), u_valve * u_valve, u_valve, 1).finished());

  return C.dot(M3) + m_offset;
}

double CalculateValveDerivative(double p_in, double p_out, double u_valve,
                                Vec<8> C) {
  Vec<4> M;
  M << u_valve *u_valve *u_valve, u_valve *u_valve, u_valve, 1;
  return (-boost::math::sign(p_in - p_out) / 2. * 100 /
          sqrt(abs(p_in * 100 - p_out * 100))) *
         M.dot(C.head<4>());
}

Compressor::FlowConstants::FlowConstants() {
  a = 340;
  V1 = 2 * pi * (0.60 / 2) * (0.60 / 2) * 2 +
       pi * (0.08 / 2) * (0.08 / 2) * 8.191;
  V2 = 0.5 * V1;

  AdivL = pi * (0.08 / 2) * (0.08 / 2) / 3 * 0.1;
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

  delta_bar = 0.1;

  n_bar = 1e2;
}
