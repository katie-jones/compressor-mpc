#include <cmath>
#include <boost/math/special_functions/sign.hpp>

#include "compressor.h"
#include "comp_coeffs.h"
#include "const_flow.h"

using namespace std;
using namespace boost::numeric::odeint;
using namespace Eigen;
using namespace comp;
using namespace boost;

void compressor::operator()(const comp_state &x, comp_state &dxdt,
                            const double /* t */) {
  const double td = u(0) * coeff::torque_drive_c / x(3);
  const double u_in = u(1);
  const double u_out = u(2);
  const double u_rec = u(3);
  flow::Pout = 1;

  const double p1 = x(0);
  const double p2 = x(1);
  const double mc = x(2);
  const double wc = x(3);
  const double mr = x(4);

  const double dp_sqrt =
      10 * sqrt(abs(flow::Pin - p1)) * math::sign(flow::Pin - p1);

  const Vec<8> M3((Vec<8>() << dp_sqrt * pow(u_in, 3), dp_sqrt * u_in * u_in,
                   dp_sqrt * u_in, dp_sqrt, pow(u_in, 3), u_in * u_in, u_in,
                   1).finished());

  const double m_in = coeff::C.transpose() * M3 + coeff::m_in_c;

  double m_rec_ss;
  if (flag)
    m_rec_ss =
        coeff::m_rec_ss_c.dot(
            (Vec<2>() << sqrt(p2 * 1e5 - p1 * 1e5) * u_rec, 1).finished()) *
        (u_rec > 1e-2);
  else
    m_rec_ss = coeff::m_rec_ss_c(0) * sqrt(p2 * 1e5 - p1 * 1e5) * u_rec;

  const double dp_sqrt2 =
      10 * sqrt(abs(p2 - flow::Pout)) * math::sign(p2 - flow::Pout);

  const Vec<8> M5((Vec<8>() << dp_sqrt2 * pow(u_out, 3),
                   dp_sqrt2 * u_out * u_out, dp_sqrt2 * u_out, dp_sqrt2,
                   pow(u_out, 3), u_out * u_out, u_out, 1).finished());

  const double m_out = coeff::D.transpose() * M5 + coeff::m_out_c;

  const double mc2 = mc * mc;
  const double mc3 = mc * mc2;
  const double wc2 = wc * wc;
  const Vec<12> M((Vec<12>() << wc2 * mc3, wc2 * mc2, wc2 * mc, wc2, wc * mc3,
                   wc * mc2, wc * mc, wc, mc3, mc2, mc, 1).finished());

  const double p_ratio = coeff::A.transpose() * M;

  const double T_ss_model =
      coeff::T_ss_c(0) + coeff::T_ss_c(1) * mc + coeff::T_ss_c(2);

  dxdt(0) = flow::a * flow::a / flow::V1 * (m_in + mr - mc) * 1e-5;
  dxdt(1) = flow::a * flow::a / flow::V2 * (mc - mr - m_out) * 1e-5;
  dxdt(2) = flow::AdivL * (p_ratio * p1 - p2) * 1e5;
  dxdt(3) = (td - T_ss_model) / coeff::J;
  dxdt(4) = coeff::tau_r * (m_rec_ss - mr);

}


