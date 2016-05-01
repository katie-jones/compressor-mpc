#ifndef SIMULATION_COMPRESSOR_H
#define SIMULATION_COMPRESSOR_H

#include "compressor.h"
#include "simulation_system.h"

class SimulationCompressor : public SimulationSystem<5, 6, 2, 2>,
                             public Compressor {
 public:
  SimulationCompressor(Compressor::State x_in = GetDefaultState(),
                       Compressor::Input u_in = GetDefaultInput())
      : SimulationSystem<n_states, n_inputs, n_outputs, n_control_inputs>(
            x_in, u_in) {}
};

// Define vector_space_norm_inf for the state used in order for odeint to work
namespace boost {
namespace numeric {
namespace odeint {
template <>
struct vector_space_norm_inf<Compressor::State> {
  typedef double result_type;
  double operator()(Compressor::State x) const {
    double absval = 0;

    for (int i = 0; i < SimulationCompressor::n_states; i++)
      absval += x[i] * x[i];
    return sqrt(absval);
  }
};
}
}
}
#endif
