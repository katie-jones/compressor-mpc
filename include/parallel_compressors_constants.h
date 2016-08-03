#ifndef PARALLEL_COMPRESSORS_CONSTANTS_H
#define PARALLEL_COMPRESSORS_CONSTANTS_H

#include "constexpr_array.h"
#include "null_index_array.h"
#include "parallel_compressors.h"

#define PARALLEL_AUGSYS_DIST1                                      \
  AugmentedLinearizedSystem<                                       \
      ParallelCompressors, PARALLEL_COMPRESSORS_CONSTANTS::Delays, \
      PARALLEL_COMPRESSORS_CONSTANTS::n_disturbance_states,        \
      PARALLEL_COMPRESSORS_CONSTANTS::ControlInputIndices1,        \
      PARALLEL_COMPRESSORS_CONSTANTS::n_sub_control_inputs>

#define PARALLEL_AUGSYS_DIST2                                      \
  AugmentedLinearizedSystem<                                       \
      ParallelCompressors, PARALLEL_COMPRESSORS_CONSTANTS::Delays, \
      PARALLEL_COMPRESSORS_CONSTANTS::n_disturbance_states,        \
      PARALLEL_COMPRESSORS_CONSTANTS::ControlInputIndices2,        \
      PARALLEL_COMPRESSORS_CONSTANTS::n_sub_control_inputs>

#define PARALLEL_AUGSYS_CENT                                       \
  AugmentedLinearizedSystem<                                       \
      ParallelCompressors, PARALLEL_COMPRESSORS_CONSTANTS::Delays, \
      PARALLEL_COMPRESSORS_CONSTANTS::n_disturbance_states,        \
      PARALLEL_COMPRESSORS_CONSTANTS::ControlInputIndices1,        \
      ParallelCompressors::n_control_inputs>

#define PARALLEL_OBS_DIST1 Observer<PARALLEL_AUGSYS_DIST1>
#define PARALLEL_OBS_DIST2 Observer<PARALLEL_AUGSYS_DIST2>
#define PARALLEL_OBS_CENT Observer<PARALLEL_AUGSYS_CENT>

#define PARALLEL_CTRL_COOP1                                                \
  DistributedController<                                                   \
      PARALLEL_AUGSYS_DIST1, PARALLEL_COMPRESSORS_CONSTANTS::StateIndices, \
      PARALLEL_COMPRESSORS_CONSTANTS::ObserverOutputIndices,               \
      PARALLEL_COMPRESSORS_CONSTANTS::ControlledOutputIndices,             \
      PARALLEL_COMPRESSORS_CONSTANTS::p, PARALLEL_COMPRESSORS_CONSTANTS::m>

#define PARALLEL_CTRL_COOP2                                                \
  DistributedController<                                                   \
      PARALLEL_AUGSYS_DIST2, PARALLEL_COMPRESSORS_CONSTANTS::StateIndices, \
      PARALLEL_COMPRESSORS_CONSTANTS::ObserverOutputIndices,               \
      PARALLEL_COMPRESSORS_CONSTANTS::ControlledOutputIndices,             \
      PARALLEL_COMPRESSORS_CONSTANTS::p, PARALLEL_COMPRESSORS_CONSTANTS::m>

#define PARALLEL_CTRL_NONCOOP1                                             \
  DistributedController<                                                   \
      PARALLEL_AUGSYS_DIST1, PARALLEL_COMPRESSORS_CONSTANTS::StateIndices, \
      PARALLEL_COMPRESSORS_CONSTANTS::ObserverOutputIndices,               \
      PARALLEL_COMPRESSORS_CONSTANTS::NCControlledOutputIndices1,          \
      PARALLEL_COMPRESSORS_CONSTANTS::p, PARALLEL_COMPRESSORS_CONSTANTS::m>

#define PARALLEL_CTRL_NONCOOP2                                             \
  DistributedController<                                                   \
      PARALLEL_AUGSYS_DIST2, PARALLEL_COMPRESSORS_CONSTANTS::StateIndices, \
      PARALLEL_COMPRESSORS_CONSTANTS::ObserverOutputIndices,               \
      PARALLEL_COMPRESSORS_CONSTANTS::NCControlledOutputIndices2,          \
      PARALLEL_COMPRESSORS_CONSTANTS::p, PARALLEL_COMPRESSORS_CONSTANTS::m>

#define PARALLEL_CTRL_CENT                                                \
  DistributedController<                                                  \
      PARALLEL_AUGSYS_CENT, PARALLEL_COMPRESSORS_CONSTANTS::StateIndices, \
      PARALLEL_COMPRESSORS_CONSTANTS::ObserverOutputIndices,              \
      PARALLEL_COMPRESSORS_CONSTANTS::ControlledOutputIndices,            \
      PARALLEL_COMPRESSORS_CONSTANTS::p, PARALLEL_COMPRESSORS_CONSTANTS::m>

namespace PARALLEL_COMPRESSORS_CONSTANTS {

using Delays = ConstexprArray<0, 40, 0, 40>;
constexpr int n_delay_states = Delays::GetSum();

constexpr int n_disturbance_states = 4;
constexpr int p = 100;
constexpr int m = 2;
constexpr int n_controllers = 2;
constexpr int n_sub_outputs = 3;
constexpr int n_sub_outputs_nc = 2;
constexpr int n_sub_control_inputs = 2;
constexpr int n_total_states =
    ParallelCompressors::n_states + n_delay_states + n_disturbance_states;

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
