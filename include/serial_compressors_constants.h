#ifndef SERIAL_COMPRESSORS_CONSTANTS_H
#define SERIAL_COMPRESSORS_CONSTANTS_H

#include "constexpr_array.h"
#include "null_index_array.h"
#include "serial_compressors.h"

#define SERIAL_AUGSYS_DIST1                                    \
  AugmentedLinearizedSystem<                                   \
      SerialCompressors, SERIAL_COMPRESSORS_CONSTANTS::Delays, \
      SERIAL_COMPRESSORS_CONSTANTS::n_disturbance_states,      \
      SERIAL_COMPRESSORS_CONSTANTS::ControlInputIndices1,      \
      SERIAL_COMPRESSORS_CONSTANTS::n_sub_control_inputs>

#define SERIAL_AUGSYS_DIST2                                    \
  AugmentedLinearizedSystem<                                   \
      SerialCompressors, SERIAL_COMPRESSORS_CONSTANTS::Delays, \
      SERIAL_COMPRESSORS_CONSTANTS::n_disturbance_states,      \
      SERIAL_COMPRESSORS_CONSTANTS::ControlInputIndices2,      \
      SERIAL_COMPRESSORS_CONSTANTS::n_sub_control_inputs>

#define SERIAL_AUGSYS_CENT                              \
  AugmentedLinearizedSystem<                                   \
      SerialCompressors, SERIAL_COMPRESSORS_CONSTANTS::Delays, \
      SERIAL_COMPRESSORS_CONSTANTS::n_disturbance_states,      \
      SERIAL_COMPRESSORS_CONSTANTS::ControlInputIndices1,      \
      SerialCompressors::n_control_inputs>

#define SERIAL_OBS_DIST1 Observer<SERIAL_AUGSYS_DIST1>
#define SERIAL_OBS_DIST2 Observer<SERIAL_AUGSYS_DIST2>
#define SERIAL_OBS_CENT Observer<SERIAL_AUGSYS_CENT>

#define SERIAL_CTRL_COOP1                                              \
  DistributedController<                                               \
      SERIAL_AUGSYS_DIST1, SERIAL_COMPRESSORS_CONSTANTS::StateIndices, \
      SERIAL_COMPRESSORS_CONSTANTS::ObserverOutputIndices,             \
      SERIAL_COMPRESSORS_CONSTANTS::ControlledOutputIndices,           \
      SERIAL_COMPRESSORS_CONSTANTS::p, SERIAL_COMPRESSORS_CONSTANTS::m>

#define SERIAL_CTRL_COOP2                                              \
  DistributedController<                                               \
      SERIAL_AUGSYS_DIST2, SERIAL_COMPRESSORS_CONSTANTS::StateIndices, \
      SERIAL_COMPRESSORS_CONSTANTS::ObserverOutputIndices,             \
      SERIAL_COMPRESSORS_CONSTANTS::ControlledOutputIndices,           \
      SERIAL_COMPRESSORS_CONSTANTS::p, SERIAL_COMPRESSORS_CONSTANTS::m>

#define SERIAL_CTRL_NONCOOP_OLD1                                           \
  DistributedController<                                               \
      SERIAL_AUGSYS_DIST1, SERIAL_COMPRESSORS_CONSTANTS::StateIndices, \
      SERIAL_COMPRESSORS_CONSTANTS::ObserverOutputIndices,             \
      SERIAL_COMPRESSORS_CONSTANTS::OldNCControlledOutputIndices1,        \
      SERIAL_COMPRESSORS_CONSTANTS::p, SERIAL_COMPRESSORS_CONSTANTS::m>

#define SERIAL_CTRL_NONCOOP_OLD2                                           \
  DistributedController<                                               \
      SERIAL_AUGSYS_DIST2, SERIAL_COMPRESSORS_CONSTANTS::StateIndices, \
      SERIAL_COMPRESSORS_CONSTANTS::ObserverOutputIndices,             \
      SERIAL_COMPRESSORS_CONSTANTS::OldNCControlledOutputIndices2,        \
      SERIAL_COMPRESSORS_CONSTANTS::p, SERIAL_COMPRESSORS_CONSTANTS::m>

#define SERIAL_CTRL_NONCOOP1                                           \
  DistributedController<                                               \
      SERIAL_AUGSYS_DIST1, SERIAL_COMPRESSORS_CONSTANTS::StateIndices, \
      SERIAL_COMPRESSORS_CONSTANTS::ObserverOutputIndices,             \
      SERIAL_COMPRESSORS_CONSTANTS::NCControlledOutputIndices1,        \
      SERIAL_COMPRESSORS_CONSTANTS::p, SERIAL_COMPRESSORS_CONSTANTS::m>

#define SERIAL_CTRL_NONCOOP2                                           \
  DistributedController<                                               \
      SERIAL_AUGSYS_DIST2, SERIAL_COMPRESSORS_CONSTANTS::StateIndices, \
      SERIAL_COMPRESSORS_CONSTANTS::ObserverOutputIndices,             \
      SERIAL_COMPRESSORS_CONSTANTS::NCControlledOutputIndices2,        \
      SERIAL_COMPRESSORS_CONSTANTS::p, SERIAL_COMPRESSORS_CONSTANTS::m>

#define SERIAL_CTRL_CENT                                              \
  DistributedController<                                              \
      SERIAL_AUGSYS_CENT, SERIAL_COMPRESSORS_CONSTANTS::StateIndices, \
      SERIAL_COMPRESSORS_CONSTANTS::ObserverOutputIndices,            \
      SERIAL_COMPRESSORS_CONSTANTS::ControlledOutputIndices,          \
      SERIAL_COMPRESSORS_CONSTANTS::p, SERIAL_COMPRESSORS_CONSTANTS::m>

namespace SERIAL_COMPRESSORS_CONSTANTS {

using Delays = ConstexprArray<0, 40, 0, 40>;
constexpr int n_delay_states = Delays::GetSum();

constexpr int n_disturbance_states = 4;
constexpr int p = 100;
constexpr int m = 2;
constexpr int n_controllers = 2;
constexpr int n_sub_outputs = 4;
constexpr int n_sub_outputs_nc = 3;
constexpr int n_sub_control_inputs = 2;
constexpr int n_total_states =
    SerialCompressors::n_states + n_delay_states + n_disturbance_states;

using StateIndices = NullIndexArray<SerialCompressors::n_states +
                                    n_delay_states + n_disturbance_states>;
using ObserverOutputIndices = NullIndexArray<SerialCompressors::n_outputs>;
using InputIndices = ConstexprArray<0, 3, 4, 7>;

using ControlledOutputIndices = NullIndexArray<SerialCompressors::n_outputs>;
using OldNCControlledOutputIndices1 = ConstexprArray<0, 1, 2>;
using OldNCControlledOutputIndices2 = ConstexprArray<2, 3, 1>;
using NCControlledOutputIndices1 = ConstexprArray<0, 1>;
using NCControlledOutputIndices2 = ConstexprArray<2, 3>;

using ControlInputIndices1 = ConstexprArray<0, 1, 2, 3>;
using ControlInputIndices2 = ConstexprArray<2, 3, 0, 1>;
}

#endif
