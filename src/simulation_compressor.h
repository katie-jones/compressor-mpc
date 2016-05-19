#ifndef SIMULATION_COMPRESSOR_H
#define SIMULATION_COMPRESSOR_H

#include "compressor.h"
#include "simulation_system.h"

class SimulationCompressor : public SimulationSystem<5, 6, 2, 2>,
                             public Compressor {
  typedef SimulationSystem<5, 6, 2, 2>::ControlInputIndex ControlInputIndex;
  static const ControlInputIndex GetDefaultInputIndex() { return {0, 3}; }

 public:
  SimulationCompressor(
      Compressor::Input u_offset = GetDefaultInput(),
      Compressor::State x_in = GetDefaultState(),
      Compressor::ControlInput u_init = Compressor::ControlInput::Zero(),
      Compressor::Parameters params = Compressor::Parameters())
      : SimulationSystem<n_states, n_inputs, n_outputs, n_control_inputs>(
            u_offset, GetDefaultInputIndex(), x_in, u_init),
        Compressor(params) {}
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
