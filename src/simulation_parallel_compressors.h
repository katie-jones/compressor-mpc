#ifndef SIMULATION_PARALLEL_COMPRESSORS_H
#define SIMULATION_PARALLEL_COMPRESSORS_H

#include "parallel_compressors.h"
#include "simulation_system.h"

class SimulationParallelCompressors : public SimulationSystem<11, 9, 5, 4>,
                                      public ParallelCompressors {
  typedef SimulationSystem<11, 9, 5, 4>::ControlInputIndex ControlInputIndex;
  static const ControlInputIndex GetDefaultInputIndex() { return {0, 3, 5, 8}; }

 public:
  SimulationParallelCompressors(
      ParallelCompressors::Input u_offset = GetDefaultInput(),
      ParallelCompressors::State x_in = GetDefaultState(),
      ParallelCompressors::ControlInput u_in =
          ParallelCompressors::ControlInput::Zero())
      : SimulationSystem<n_states, n_inputs, n_outputs, n_control_inputs>(
            u_offset, GetDefaultInputIndex(), x_in, u_in) {}
};

// Define vector_space_norm_inf for the state used in order for odeint to work
namespace boost {
namespace numeric {
namespace odeint {
template <>
struct vector_space_norm_inf<ParallelCompressors::State> {
  typedef double result_type;
  double operator()(ParallelCompressors::State x) const {
    double absval = 0;

    for (int i = 0; i < SimulationParallelCompressors::n_states; i++)
      absval += x[i] * x[i];
    return sqrt(absval);
  }
};
}
}
}

#endif
