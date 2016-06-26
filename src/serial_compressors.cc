#include "serial_compressors.h"

#include "valve_eqs.h"

using namespace ValveEqs;
using Cp = CompressorBase;

SerialCompressors::State SerialCompressors::GetDerivative(
    const State& x, const Input& u) const {
  State dxdt;
  Cp::Input u_sub;
  double m_out = -1;

  for (int i = 0; i < n_compressors; i++) {
    u_sub = GetCompressorInput(u, i, x);
    if (i > 0) {
      u_sub(n_comp_inputs) = m_out;
    }

    dxdt.template segment<Cp::n_states>(i * Cp::n_states) =
        p_comps_[i]->GetDerivative(
            &m_out, x.template segment<Cp::n_states>(i * Cp::n_states), u_sub);
  }

  return dxdt;
}

SerialCompressors::Linearized SerialCompressors::GetLinearizedSystem(
    const State& x, const Input& u) const {
  Linearized linsys;
  double m_out = 0;
  double p_compressor;
  Cp::Input u_sub;

  linsys.A.setZero();
  linsys.B.setZero();
  linsys.C.setZero();
  linsys.f.setZero();

  Cp::Linearized comp_linsys;

  // First compressor
  comp_linsys = first_comp_.GetLinearizedSystem(
      &m_out, x.template head<Cp::n_states>(), GetCompressorInput(u, 0, x));
  linsys.A.topLeftCorner<Cp::n_states, Cp::n_states>() = comp_linsys.A;
  linsys.B.topLeftCorner<Cp::n_states, Cp::n_control_inputs>() = comp_linsys.B;
  linsys.C.topLeftCorner<Cp::n_outputs, Cp::n_states>() = comp_linsys.C;

  // Effect of first on second compressor
  linsys.A(1 * Cp::n_states, 1) = CalculateValveDerivative(
      x(1), x(Cp::n_states), u(2), first_comp_.params_.D, comps_[0].params_.V1);

  // Following compressors
  for (int i = 1; i < n_compressors; i++) {
    // Get input
    u_sub = GetCompressorInput(u, i, x);
    u_sub(n_comp_inputs) = m_out;

    // Linearized compressor model
    comp_linsys = p_comps_[i]->GetLinearizedSystem(
        &m_out, x.template segment<Cp::n_states>(i * Cp::n_states),
        GetCompressorInput(u, i, x));

    // A part of compressor i
    linsys.A.template block<Cp::n_states, Cp::n_states>(
        i * Cp::n_states, i * Cp::n_states) = comp_linsys.A;

    // Fix dp1/p1 in A
    linsys.A(i * Cp::n_states, i * Cp::n_states) =
        -CalculateValveDerivative(
            x((i - 1) * Cp::n_states + 1), x(i * Cp::n_states),
            u((i - 1) * n_comp_inputs + 2), p_comps_[i - 1]->params_.D,
            p_comps_[i]->params_.V1);

    // Cross terms
    // Effect of current p1 on previous p2
    linsys.A((i - 1) * Cp::n_states + 1, i * Cp::n_states) =
        CalculateValveDerivative(
            x((i - 1) * Cp::n_states + 1), x(i * Cp::n_states),
            u((i - 1) * n_comp_inputs + 2), p_comps_[i - 1]->params_.D,
            p_comps_[i - 1]->params_.V2);

    // Effect of current p2 on next p1
    if (i < n_follower_compressors) {
      linsys.A((i + 1) * Cp::n_states, i * Cp::n_states + 1) =
          CalculateValveDerivative(
              x(i * Cp::n_states + 1), x((i + 1) * Cp::n_states),
              u(i * n_comp_inputs + 2), p_comps_[i]->params_.D,
              p_comps_[i + 1]->params_.V1);
    }

    // B part of compressor i
    linsys.B.block<Cp::n_states, Cp::n_control_inputs>(
        i * Cp::n_states, i * Cp::n_control_inputs) = comp_linsys.B;

    // C part of compressor i
    linsys.C.block<1, Cp::n_states>(i, i * Cp::n_states) =
        comp_linsys.C.row(1);  // surge distances

    // Compressor derivative
    linsys.f = GetDerivative(x, u);
  }

  return linsys;
}

SerialCompressors::Output SerialCompressors::GetOutput(const State& x) const {
  Output y;

  for (int i = 0; i < n_compressors; i++) {
    y.template segment<Cp::n_outputs>(i * Cp::n_outputs) =
        p_comps_[i]->GetOutput(
            x.template segment<Cp::n_states>(i * Cp::n_states));
  }
  return y;
}
