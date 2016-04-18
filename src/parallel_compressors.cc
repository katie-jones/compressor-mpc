#include "parallel_compressors.h"


// ParallelCompressors::SysState ParallelCompressors::GetDerivative(const
// SysState x, const SysInput u) const {
// Compressor::CompressorState comp_derivs[n_compressors];
// for (int i=0; i<n_compressors; i++) {
// comps[i].u = u.block<n_comp_inputs,1>(i*n_comp_inputs,0);
// comps[i].pout = x(n_states-1);
// comps[i].x = x.block<n_comp_states,1>(i*n_comp_states,0);
// comp_derivs[i] = comps[i].GetDerivative();
// }

// }
