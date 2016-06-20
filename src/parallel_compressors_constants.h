#ifndef PARALLEL_COMPRESSORS_CONSTANTS_H
#define PARALLEL_COMPRESSORS_CONSTANTS_H

#include "parallel_compressors.h"
#include "constexpr_array.h"
#include "null_index_array.h"

#define AUGMENTEDSYSTEM1                                           \
  AugmentedLinearizedSystem<                                       \
      ParallelCompressors, PARALLEL_COMPRESSORS_CONSTANTS::Delays, \
      PARALLEL_COMPRESSORS_CONSTANTS::n_disturbance_states,        \
      PARALLEL_COMPRESSORS_CONSTANTS::ControlInputIndices1,        \
      PARALLEL_COMPRESSORS_CONSTANTS::n_sub_control_inputs>

#define AUGMENTEDSYSTEM2                                           \
  AugmentedLinearizedSystem<                                       \
      ParallelCompressors, PARALLEL_COMPRESSORS_CONSTANTS::Delays, \
      PARALLEL_COMPRESSORS_CONSTANTS::n_disturbance_states,        \
      PARALLEL_COMPRESSORS_CONSTANTS::ControlInputIndices2,        \
      PARALLEL_COMPRESSORS_CONSTANTS::n_sub_control_inputs>

#define OBSERVER1 Observer<AUGMENTEDSYSTEM1>
#define OBSERVER2 Observer<AUGMENTEDSYSTEM2>

#define CONTROLLER1                                                   \
  DistributedController<                                              \
      AUGMENTEDSYSTEM1, PARALLEL_COMPRESSORS_CONSTANTS::StateIndices, \
      PARALLEL_COMPRESSORS_CONSTANTS::ObserverOutputIndices,          \
      PARALLEL_COMPRESSORS_CONSTANTS::ControlledOutputIndices,        \
      PARALLEL_COMPRESSORS_CONSTANTS::p, PARALLEL_COMPRESSORS_CONSTANTS::m>

#define CONTROLLER2                                                   \
  DistributedController<                                              \
      AUGMENTEDSYSTEM2, PARALLEL_COMPRESSORS_CONSTANTS::StateIndices, \
      PARALLEL_COMPRESSORS_CONSTANTS::ObserverOutputIndices,          \
      PARALLEL_COMPRESSORS_CONSTANTS::ControlledOutputIndices,        \
      PARALLEL_COMPRESSORS_CONSTANTS::p, PARALLEL_COMPRESSORS_CONSTANTS::m>

namespace PARALLEL_COMPRESSORS_CONSTANTS {

constexpr int n_delay_states = 80;
constexpr int n_disturbance_states = 4;
constexpr int p = 100;
constexpr int m = 2;
constexpr int n_controllers = 2;
constexpr int n_sub_outputs = 3;
constexpr int n_sub_control_inputs = 2;
constexpr int n_total_states =
    ParallelCompressors::n_states + n_delay_states + n_disturbance_states;

using Delays = ConstexprArray<0, 40, 0, 40>;
using OutputIndices = ConstexprArray<0, 1, 3>;
using StateIndices = NullIndexArray<ParallelCompressors::n_states +
                                    n_delay_states + n_disturbance_states>;
using ObserverOutputIndices = NullIndexArray<ParallelCompressors::n_outputs>;
using InputIndices = ConstexprArray<0, 3, 4, 7>;
using ControlledOutputIndices = ConstexprArray<0, 1, 3>;

using OutputIndexList = ConstexprArrayList<OutputIndices, OutputIndices>;
using StateIndexList = ConstexprArrayList<StateIndices, StateIndices>;

using ControlInputIndices1 = ConstexprArray<0, 1, 2, 3>;
using ControlInputIndices2 = ConstexprArray<2, 3, 0, 1>;
}

#endif
