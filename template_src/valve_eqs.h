#ifndef VALVE_EQS_H
#define VALVE_EQS_H

namespace ValveEqs {
// Calculate partial derivative wrt pressure across a valve
// p_in: upstream pressure
// p_out: downstream pressure
// u_valve: valve setting
// C: coefficients of valve map
inline double CalculateValveDerivative(double p_in, double p_out, double u_valve,
                                Vec<8> C) {
  Vec<4> M;
  M << u_valve *u_valve *u_valve, u_valve *u_valve, u_valve, 1;
  return (-boost::math::sign(p_in - p_out) / 2. * 100 /
          sqrt(std::abs(p_in * 100 - p_out * 100))) *
         M.dot(C.head<4>());
}

// Calculate the mass flow across a valve
// p_in: upstream pressure
// p_out: downstream pressure
// u_valve: valve setting
// C: coefficients of valve map
// m_offset: mass flow offset
inline double CalculateValveMassFlow(double p_in, double p_out,
                                             double u_valve, Vec<8> C,
                                             double m_offset) {
  const double dp_sqrt =
      10 * std::sqrt(std::abs(p_in - p_out)) * boost::math::sign(p_in - p_out);

  const Vec<8> M3((Vec<8>() << dp_sqrt * u_valve * u_valve * u_valve,
                   dp_sqrt * u_valve * u_valve, dp_sqrt * u_valve, dp_sqrt,
                   u_valve * u_valve * u_valve, u_valve * u_valve, u_valve,
                   1).finished());

  return C.dot(M3) + m_offset;
}
}

#endif
