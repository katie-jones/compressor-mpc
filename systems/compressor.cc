#include "compressor.h"

#include <cmath>
#include <iostream>

#include "valve_eqs.h"

constexpr double pi = 3.14159265358979323846;
constexpr double speed_sound = 340.;

using namespace Eigen;
using namespace ValveEqs;

template <bool has_input_tank>
auto Compressor<has_input_tank>::GetDerivative(double *m_out, const State &x,
                                               const Input &u) const -> State {
  State dxdt;

  const double p1 = x(0);
  const double p2 = x(1);
  const double mc = x(2);
  const double wc = x(3);
  const double mr = x(4);

  const double td = u(0) * params_.torque_drive_c / wc;
  const double u_input = u(1);
  const double u_out = u(2);
  const double u_rec = u(3);
  const double p_out = u(5);

  double m_in;
  if (has_input_tank) {
    m_in = CalculateValveMassFlow(u(4), p1, u_input, params_.C, params_.m_in_c);
  } else {
    m_in = u(4);
  }

  *m_out = CalculateValveMassFlow(p2, p_out, u_out, params_.D, params_.m_out_c);

  double m_rec_ss =
      params_.m_rec_ss_c.dot(
          (Eigen::Matrix<double, 2, 1>() << sqrt(p2 * 1e5 - p1 * 1e5) * u_rec,
           1).finished()) *
      (u_rec > 1e-2);

  const double mc2 = mc * mc;
  const double mc3 = mc * mc2;
  const double wc2 = wc * wc;
  const Eigen::Matrix<double, 12, 1> M(
      (Eigen::Matrix<double, 12, 1>() << wc2 * mc3, wc2 * mc2, wc2 * mc, wc2,
       wc * mc3, wc * mc2, wc * mc, wc, mc3, mc2, mc,
       1).finished());

  const double p_ratio = params_.A.transpose() * M;

  const double T_ss_model =
      params_.T_ss_c(0) + params_.T_ss_c(1) * mc + params_.T_ss_c(2);

  dxdt(0) = speed_sound * speed_sound / params_.V1 * (m_in + mr - mc) * 1e-5;
  dxdt(1) = speed_sound * speed_sound / params_.V2 * (mc - mr - *m_out) * 1e-5;
  dxdt(2) = params_.AdivL * (p_ratio * p1 - p2) * 1e5;
  dxdt(3) = (td - T_ss_model) / params_.J;
  dxdt(4) = params_.tau_r * (m_rec_ss - mr);

  return dxdt;
}

auto CompressorBase::GetOutput(const State &x) const -> Output {
  const double p1 = x(0);
  const double p2 = x(1);
  const double mass_flow = x(2);
  const double surge_distance =
      params_.SD_multiplier * (-(p2 / p1) / params_.SD_c(0) +
                               params_.SD_c(1) / params_.SD_c(0) + mass_flow);
  Output y;
  y << p2, surge_distance;
  return y;
}

template <bool has_input_tank>
auto Compressor<has_input_tank>::GetLinearizedSystem(double *m_out,
                                                     const State &x,
                                                     const Input &u) const
    -> Linearized {
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
  const double p_out = u(5);

  // Partial derivatives of p1
  linsys.A.row(0) << -1, 0, -speed_sound * speed_sound / params_.V1 * 1e-5, 0,
      speed_sound * speed_sound / params_.V1 * 1e-5;

  // if there's an input tank, calculate the derivative
  if (has_input_tank) {
    linsys.A(0, 0) =
        -CalculateValveDerivative(u(4), p1, u_input, params_.C, params_.V1);
  }

  // Partial derivatives of p2
  linsys.A.row(1) << 0,
      -CalculateValveDerivative(p2, p_out, u_out, params_.D, params_.V2),
      speed_sound * speed_sound / params_.V2 * 1e-5, 0,
      -speed_sound * speed_sound / params_.V2 * 1e-5;

  // Intermediate calculations
  const double wc2 = wc * wc;
  const double wc3 = wc * wc2;
  const double mc2 = mc * mc;
  const double mc3 = mc * mc2;

  Eigen::Matrix<double, 12, 1> M, dM_dmcomp, dM_dwcomp;

  dM_dmcomp << 3 * wc2 * mc2, 2 * wc2 * mc, wc2, 0, 3 * wc * mc2, 2 * wc * mc,
      wc, 0, 3 * mc2, 2 * mc, 1, 0;
  dM_dwcomp << 2 * wc * mc3, 2 * wc * mc2, 2 * wc * mc, 2 * wc, mc3, mc2, mc, 1,
      0, 0, 0, 0;

  M << wc2 * mc3, wc2 * mc2, wc2 * mc, wc2, wc * mc3, wc * mc2, wc * mc, wc,
      mc3, mc2, mc, 1;

  double p_ratio = params_.A.dot(M);

  // Partial derivatives of qc
  linsys.A.row(2) << params_.AdivL * (p_ratio * 1e5), -params_.AdivL * 1e5,
      params_.AdivL * (p1 * 1e5) * params_.A.dot(dM_dmcomp),
      params_.AdivL * (p1 * 1e5) * params_.A.dot(dM_dwcomp), 0;

  // Partial derivatives of omegac
  linsys.A.row(3) << 0, 0, -1.0 / params_.J * params_.T_ss_c(1),
      -1.0 / params_.J * td_in * params_.torque_drive_c / wc2, 0;

  // Partial derivatives of qr
  linsys.A.row(4) << -params_.tau_r * (params_.m_rec_ss_c(0) * 1 / 2 * u_rec /
                                       sqrt(p2 * 1e5 - p1 * 1e5) * 1e5),
      params_.tau_r * (params_.m_rec_ss_c(0) * 1 / 2 * u_rec /
                       sqrt(p2 * 1e5 - p1 * 1e5) * 1e5),
      0, 0, -params_.tau_r;

  // B matrix
  // Approximate deadzone using exponentials
  double dmr_ur =
      params_.tau_r * params_.m_rec_ss_c(0) * sqrt(p2 * 1e5 - p1 * 1e5);

  constexpr double x0 = 1e-2;
  if (u_rec < 2 * x0) {
    double a;
    if (u_rec >= x0) {
      a = params_.delta_bar +
          (1 - params_.delta_bar) * exp(params_.n_bar * (u_rec - x0));
    } else {
      a = 2 - (1 - params_.delta_bar) * exp(-params_.n_bar * u_rec);
    }
    dmr_ur = a * dmr_ur;
  }

  linsys.B << 0, 0, 0, 0, 0, 0, 1.0 / params_.J * params_.torque_drive_c / wc,
      0, 0, dmr_ur;

  linsys.C << 0, 1, 0, 0, 0, 100 * p2 / (params_.SD_c(0) * p1 * p1),
      -100. / (params_.SD_c(0) * p1), 100, 0, 0;

  linsys.f = GetDerivative(m_out, x, u);

  return linsys;
}

CompressorBase::Parameters::Parameters() {
  J = (0.4 + 0.2070) * 0.4;
  tau_r = 1 / 0.5;
  A = (Eigen::Matrix<double, 12, 1>() << 0.000299749505193654,
       -0.000171254191089237, 3.57321648097597e-05, -9.1783572200945e-07,
       -0.252701086129365, 0.136885752773673, -0.02642368327081,
       0.00161012740365743, 54.8046725371143, -29.9550791497765,
       5.27827499839098, 0.693826282579158)
          .finished();

  C = (Eigen::Matrix<double, 8, 1>() << -0.423884232813775, 0.626400271518973,
       -0.0995040168384753, 0.0201535563630318, -0.490814924104294,
       0.843580880467905, -0.423103455111209, 0.0386841406482887)
          .finished();

  D = (Eigen::Matrix<double, 8, 1>() << -0.0083454, -0.0094965, 0.16826,
       -0.032215, -0.61199, 0.94175, -0.48522, 0.10369)
          .finished();

  m_in_c = 0.0051;
  m_rec_ss_c = (Eigen::Matrix<double, 2, 1>() << 0.0047, 0.0263).finished();

  m_out_c = 0.017;

  T_ss_c = (Eigen::Matrix<double, 3, 1>() << 2.5543945754982, 47.4222669576423,
            0.6218)
               .finished();

  SD_c = (Eigen::Matrix<double, 2, 1>() << 5.55, 0.66).finished();
  SD_multiplier = 100;

  torque_drive_c = 15000;

  delta_bar = 0.1;

  n_bar = 1e2;

  V1 = 2 * pi * (0.60 / 2.0) * (0.60 / 2.0) * 2.0 +
       pi * (0.08 / 2.0) * (0.08 / 2.0) * 8.191;

  V2 = pi * (0.60 / 2.0) * (0.60 / 2.0) * 2.0 +
       pi * (0.08 / 2.0) * (0.08 / 2.0) * 5.940;

  AdivL = pi * (0.08 / 2) * (0.08 / 2) / 3 * 0.1;
}

#include "compressor_list.h"
