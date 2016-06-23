#ifndef PARALLEL_COMPRESSORS_CONSTANTS_H
#define PARALLEL_COMPRESSORS_CONSTANTS_H

#include "constexpr_array.h"
#include "null_index_array.h"
#include "parallel_compressors.h"

#define AUGMENTEDSYSTEM_DIST1                                      \
  AugmentedLinearizedSystem<                                       \
      ParallelCompressors, PARALLEL_COMPRESSORS_CONSTANTS::Delays, \
      PARALLEL_COMPRESSORS_CONSTANTS::n_disturbance_states,        \
      PARALLEL_COMPRESSORS_CONSTANTS::ControlInputIndices1,        \
      PARALLEL_COMPRESSORS_CONSTANTS::n_sub_control_inputs>

#define AUGMENTEDSYSTEM_DIST2                                      \
  AugmentedLinearizedSystem<                                       \
      ParallelCompressors, PARALLEL_COMPRESSORS_CONSTANTS::Delays, \
      PARALLEL_COMPRESSORS_CONSTANTS::n_disturbance_states,        \
      PARALLEL_COMPRESSORS_CONSTANTS::ControlInputIndices2,        \
      PARALLEL_COMPRESSORS_CONSTANTS::n_sub_control_inputs>

#define AUGMENTEDSYSTEM_CENTRALIZED                                \
  AugmentedLinearizedSystem<                                       \
      ParallelCompressors, PARALLEL_COMPRESSORS_CONSTANTS::Delays, \
      PARALLEL_COMPRESSORS_CONSTANTS::n_disturbance_states,        \
      PARALLEL_COMPRESSORS_CONSTANTS::ControlInputIndices1,        \
      ParallelCompressors::n_control_inputs>

#define OBSERVER_DIST1 Observer<AUGMENTEDSYSTEM_DIST1>
#define OBSERVER_DIST2 Observer<AUGMENTEDSYSTEM_DIST2>
#define OBSERVER_CENTRALIZED Observer<AUGMENTEDSYSTEM_CENTRALIZED>

#define CONTROLLER_COOP1                                                   \
  DistributedController<                                                   \
      AUGMENTEDSYSTEM_DIST1, PARALLEL_COMPRESSORS_CONSTANTS::StateIndices, \
      PARALLEL_COMPRESSORS_CONSTANTS::ObserverOutputIndices,               \
      PARALLEL_COMPRESSORS_CONSTANTS::ControlledOutputIndices,             \
      PARALLEL_COMPRESSORS_CONSTANTS::p, PARALLEL_COMPRESSORS_CONSTANTS::m>

#define CONTROLLER_COOP2                                                   \
  DistributedController<                                                   \
      AUGMENTEDSYSTEM_DIST2, PARALLEL_COMPRESSORS_CONSTANTS::StateIndices, \
      PARALLEL_COMPRESSORS_CONSTANTS::ObserverOutputIndices,               \
      PARALLEL_COMPRESSORS_CONSTANTS::ControlledOutputIndices,             \
      PARALLEL_COMPRESSORS_CONSTANTS::p, PARALLEL_COMPRESSORS_CONSTANTS::m>

#define CONTROLLER_NONCOOP1                                                   \
  DistributedController<                                                   \
      AUGMENTEDSYSTEM_DIST1, PARALLEL_COMPRESSORS_CONSTANTS::StateIndices, \
      PARALLEL_COMPRESSORS_CONSTANTS::ObserverOutputIndices,               \
      PARALLEL_COMPRESSORS_CONSTANTS::NCControlledOutputIndices1,             \
      PARALLEL_COMPRESSORS_CONSTANTS::p, PARALLEL_COMPRESSORS_CONSTANTS::m>

#define CONTROLLER_NONCOOP2                                                   \
  DistributedController<                                                   \
      AUGMENTEDSYSTEM_DIST2, PARALLEL_COMPRESSORS_CONSTANTS::StateIndices, \
      PARALLEL_COMPRESSORS_CONSTANTS::ObserverOutputIndices,               \
      PARALLEL_COMPRESSORS_CONSTANTS::NCControlledOutputIndices2,             \
      PARALLEL_COMPRESSORS_CONSTANTS::p, PARALLEL_COMPRESSORS_CONSTANTS::m>

#define CONTROLLER_CENTRALIZED                                 \
  DistributedController<                                       \
      AUGMENTEDSYSTEM_CENTRALIZED,                             \
      PARALLEL_COMPRESSORS_CONSTANTS::StateIndices,            \
      PARALLEL_COMPRESSORS_CONSTANTS::ObserverOutputIndices,   \
      PARALLEL_COMPRESSORS_CONSTANTS::ControlledOutputIndices, \
      PARALLEL_COMPRESSORS_CONSTANTS::p, PARALLEL_COMPRESSORS_CONSTANTS::m>

namespace PARALLEL_COMPRESSORS_CONSTANTS {

constexpr int n_delay_states = 80;
constexpr int n_disturbance_states = 4;
constexpr int p = 100;
constexpr int m = 2;
constexpr int n_controllers = 2;
constexpr int n_sub_outputs = 3;
constexpr int n_sub_outputs_nc = 2;
constexpr int n_sub_control_inputs = 2;
constexpr int n_total_states =
    ParallelCompressors::n_states + n_delay_states + n_disturbance_states;

using Delays = ConstexprArray<0, 40, 0, 40>;
using StateIndices = NullIndexArray<ParallelCompressors::n_states +
                                    n_delay_states + n_disturbance_states>;
using ObserverOutputIndices = NullIndexArray<ParallelCompressors::n_outputs>;
using InputIndices = ConstexprArray<0, 3, 4, 7>;

using ControlledOutputIndices = ConstexprArray<0, 1, 3>;
using NCControlledOutputIndices1 = ConstexprArray<0, 3>;
using NCControlledOutputIndices2 = ConstexprArray<1, 3>;

using ControlInputIndices1 = ConstexprArray<0, 1, 2, 3>;
using ControlInputIndices2 = ConstexprArray<2, 3, 0, 1>;
}

#endif
