#include "parallel_compressors.h"

#include "valve_eqs.h"

typedef Compressor Cp;

using namespace ValveEqs;

ParallelCompressors::State ParallelCompressors::GetDerivative(
    const State x, const Input u) const {
  State dxdt;
  double mass_flow;
  double mass_flow_total = 0;

  for (int i = 0; i < n_compressors; i++) {
    dxdt.segment<Cp::n_states>(i * Cp::n_states) =
        comps_[i].GetDerivative(x.segment<Cp::n_states>(i * Cp::n_states),
                                GetCompressorInput(u, i, x), mass_flow);
    mass_flow_total += mass_flow;
  }
  dxdt.tail<Tank::n_states>() = tank_.GetDerivative(
      x.tail<Tank::n_states>(),
      (Tank::Input() << u(n_inputs - 1), p_out_, mass_flow_total).finished());
  return dxdt;
}

ParallelCompressors::Linearized ParallelCompressors::GetLinearizedSystem(
    const State x, const Input u) const {
  Linearized linsys;
  double mass_flow, mass_flow_total = 0;
  double p_compressor;

  linsys.A.setZero();
  linsys.B.setZero();
  linsys.C.setZero();
  linsys.f.setZero();

  Cp::Linearized comps_linsys[n_compressors];
  for (int i = 0; i < n_compressors; i++) {
    // linearized compressor
    comps_linsys[i] = comps_[i].GetLinearizedSystem(
        x.segment<Cp::n_states>(i * Cp::n_states), GetCompressorInput(u, i, x));

    // A part of compressor i
    linsys.A.block<Cp::n_states, Cp::n_states>(
        i * Cp::n_states, i * Cp::n_states) = comps_linsys[i].A;

    // B part of compressor i
    linsys.B.block<Cp::n_states, Cp::n_control_inputs>(
        i * Cp::n_states, i * Cp::n_control_inputs) = comps_linsys[i].B;

    // C part of compressor i
    // linsys.C.block<Cp::n_outputs, Cp::n_states>(
    // i * Cp::n_outputs, i * Cp::n_states) = comps_linsys[i].C;
    linsys.C.block<1, Cp::n_states>(i, i * Cp::n_states) =
        comps_linsys[i].C.row(1);  // surge distances

    // mass flow out of compressor i
    mass_flow_total += CalculateValveMassFlow(
        x(i * Cp::n_states + 1), x(n_states - 1), u(i * n_comp_inputs + 2),
        comps_[i].params_.D, comps_[i].params_.m_out_c);

    // effect of compressor i on tank
    linsys.A.block<Tank::n_states, Cp::n_states>(n_compressors * Cp::n_states,
                                                 i * Cp::n_states)
        << 0,
        CalculateValveDerivative(x(i * Cp::n_states + 1), x(n_states - 1),
                                 u(i * n_comp_inputs + 2), comps_[i].params_.D,
                                 tank_.params_.volume),
        0, 0, 0;

    // effect of tank on compressor i
    linsys.A.block<Cp::n_states, Tank::n_states>(i * Cp::n_states,
                                                 n_compressors * Cp::n_states)
        << 0,
        CalculateValveDerivative(x(i * Cp::n_states + 1), x(n_states - 1),
                                 u(i * n_comp_inputs + 2), comps_[i].params_.D,
                                 comps_[i].params_.V2),
        0, 0, 0;
  }

  Tank::Linearized tank_linsys = tank_.GetLinearizedSystem(
      x.tail<Tank::n_states>(),
      (Tank::Input() << u(n_inputs - 1), p_out_, mass_flow_total).finished());

  linsys.A.block<Tank::n_states, Tank::n_states>(n_compressors * Cp::n_states,
                                                 n_compressors * Cp::n_states) =
      tank_linsys.A;

  // Delta pressure
  linsys.C.row(n_compressors) << comps_linsys[0].C.row(0),
      -comps_linsys[1].C.row(0),
      Eigen::Matrix<double, Tank::n_outputs, Tank::n_states>::Zero();

  // Tank pressure
  linsys.C.block<Tank::n_outputs, Tank::n_states>(
      n_outputs - Tank::n_outputs - 1, n_compressors * Cp::n_states) =
      tank_linsys.C;

  linsys.f = GetDerivative(x, u);

  return linsys;
}

ParallelCompressors::Output ParallelCompressors::GetOutput(
    const State x) const {
  Output y;
  Compressor::Output y_comp;
  Eigen::Matrix<double, n_compressors, 1> pressures, surge_distances;

  for (int i = 0; i < n_compressors; i++) {
    y_comp = comps_[i].GetOutput(
        x.segment<Compressor::n_states>(i * Compressor::n_states));
    // y.segment<Compressor::n_outputs>(i * Compressor::n_outputs) =
    // comps_[i].GetOutput(
    // x.segment<Compressor::n_states>(i * Compressor::n_states));
    pressures[i] = y_comp(0);
    surge_distances[i] = y_comp[1];
  }
  y.template head<n_compressors>() = surge_distances;
  y(n_compressors) = pressures(0) - pressures(1);
  y.tail<1>() = tank_.GetOutput(x.tail<1>());  // tank pressure
  return y;
}
