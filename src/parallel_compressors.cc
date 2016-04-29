#include "parallel_compressors.h"

ParallelCompressors::SysState ParallelCompressors::GetDerivative(
    const SysState x, const SysInput u) const {
  SysState dxdt;
  double mass_flow = 0;
  double mass_flow_total = 0;

  // Calculate derivative for each compressor
  for (int i = 0; i < n_compressors; i++) {
    dxdt.segment<n_comp_states>(i * n_comp_states) = comps[i].GetDerivative(
        // states of current compressor
        x.segment<n_comp_states>(i * n_comp_states),
        // inputs of current compressor
        u.segment<n_comp_inputs>(i * n_comp_inputs),
        // pressure in output tank
        x(x.size() - n_tank_states),
        // return mass flow through compressor in this variable
        mass_flow);

    // Update total mass flow through compressors
    mass_flow_total += mass_flow;
  }

  // Get derivative of output tank
  dxdt.tail<n_tank_states>() = tank.GetDerivative(
      x.tail<n_tank_states>(), u.tail<n_tank_inputs>(), mass_flow_total);

  return dxdt;
}
