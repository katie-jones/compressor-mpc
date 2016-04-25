#include "parallel_compressors.h"

ParallelCompressors::SysState ParallelCompressors::GetDerivative(
    const SysState x, const SysInput u) const {
  SysState dxdt;
  Tank::TankState ptank = x.tail<n_tank_states>();
  Comp::CompressorState comp_deriv;
  Comp::CompressorState xcomp;
  Comp::CompressorInput ucomp;
  double mass_flow = 0;
  double mass_flow_total = 0;

  for (int i = 0; i < n_compressors; i++) {
    xcomp = x.segment<n_comp_states>(i * n_comp_states);
    ucomp = u.segment<n_comp_inputs>(i * n_comp_inputs);
    comp_deriv = comps[i].GetDerivative(xcomp, ucomp, ptank(0), mass_flow);
    std::cout << "mass flow: " << mass_flow << std::endl;
    dxdt.segment<n_comp_states>(i * n_comp_states) = comp_deriv;
    mass_flow_total += mass_flow;
  }
  dxdt.tail<n_tank_states>() =
      tank.GetDerivative(ptank, u.tail<n_tank_inputs>(), mass_flow_total);

  return dxdt;
}

