#ifndef COMMON_VARIABLES_H
#define COMMON_VARIABLES_H

#include <boost/timer/timer.hpp>
#include <fstream>
#include <iostream>
#include "aug_lin_sys.h"
#include "constexpr_array.h"
#include "distributed_controller.h"
#include "input_constraints.h"
#include "nerve_center.h"
#include "null_index_array.h"
#include "observer.h"
#include "read_files.h"
#include "simulation_system.h"

//------------------------------------------------------------------------------
//-------------------------- PARALLEL COMPRESSORS ------------------------------
//------------------------------------------------------------------------------
#ifdef SYSTEM_TYPE_PARALLEL
#include "parallel_compressors.h"
#include "parallel_compressors_constants.h"

using namespace PARALLEL_COMPRESSORS_CONSTANTS;
using SimSystem = SimulationSystem<ParallelCompressors, Delays, InputIndices>;

//------------------------------ CENTRALIZED -----------------------------------
#ifdef CONTROLLER_TYPE_CENTRALIZED
using AugmentedSystem = PARALLEL_AUGSYS_CENT;
using Obsv = PARALLEL_OBS_CENT;
using Controller = PARALLEL_CTRL_CENT;
using NvCtr = NerveCenter<ParallelCompressors, n_total_states, Controller>;

#else   // ifdef CONTROLLER_TYPE_CENTRALIZED
//------------------------------ DISTRIBUTED -----------------------------------
using Obsv1 = PARALLEL_OBS_DIST1;
using Obsv2 = PARALLEL_OBS_DIST2;
using Obsv = Obsv1;

using AugmentedSystem1 = PARALLEL_AUGSYS_DIST1;
using AugmentedSystem2 = PARALLEL_AUGSYS_DIST2;
#endif  // ifdef CONTROLLER_TYPE_CENTRALIZED

//------------------------------ COOPERATIVE -----------------------------------
#ifdef CONTROLLER_TYPE_COOP
using NvCtr = NerveCenter<ParallelCompressors, n_total_states,
                          PARALLEL_CTRL_COOP1, PARALLEL_CTRL_COOP2>;

using Controller1 = PARALLEL_CTRL_COOP1;
using Controller2 = PARALLEL_CTRL_COOP2;

#endif  // ifdef CONTROLLER_TYPE_COOP
//---------------------------- NON-COOPERATIVE ---------------------------------
#ifdef CONTROLLER_TYPE_NCOOP
using NvCtr = NerveCenter<ParallelCompressors, n_total_states,
                          PARALLEL_CTRL_COOP1, PARALLEL_CTRL_COOP2>;

using Controller1 = PARALLEL_CTRL_COOP1;
using Controller2 = PARALLEL_CTRL_COOP2;
#endif  // ifdef CONTROLLER_TYPE_NCOOP

#endif  // ifdef SYSTEM_TYPE_PARALLEL

//------------------------------------------------------------------------------
//--------------------------- SERIAL COMPRESSORS -------------------------------
//------------------------------------------------------------------------------
#ifdef SYSTEM_TYPE_SERIAL
#include "serial_compressors.h"
#include "serial_compressors_constants.h"

using namespace SERIAL_COMPRESSORS_CONSTANTS;
using SimSystem = SimulationSystem<SerialCompressors, Delays, InputIndices>;

//------------------------------ CENTRALIZED -----------------------------------
#ifdef CONTROLLER_TYPE_CENTRALIZED
using AugmentedSystem = SERIAL_AUGSYS_CENT;
using Obsv = SERIAL_OBS_CENT;
using Controller = SERIAL_CTRL_CENT;
using NvCtr = NerveCenter<SerialCompressors, n_total_states, Controller>;

#else
//------------------------------ DISTRIBUTED -----------------------------------
using Obsv1 = SERIAL_OBS_DIST1;
using Obsv2 = SERIAL_OBS_DIST2;
using Obsv = Obsv1;

using AugmentedSystem1 = SERIAL_AUGSYS_DIST1;
using AugmentedSystem2 = SERIAL_AUGSYS_DIST2;
#endif  // ifdef CONTROLLER_TYPE_CENTRALIZED

//------------------------------ COOPERATIVE -----------------------------------
#ifdef CONTROLLER_TYPE_COOP
using NvCtr = NerveCenter<SerialCompressors, n_total_states, SERIAL_CTRL_COOP1,
                          SERIAL_CTRL_COOP2>;

using Controller1 = SERIAL_CTRL_COOP1;
using Controller2 = SERIAL_CTRL_COOP2;
#endif  // ifdef CONTROLLER_TYPE_COOP

//---------------------------- NON-COOPERATIVE ---------------------------------
#ifdef CONTROLLER_TYPE_NCOOP
using NvCtr = NerveCenter<SerialCompressors, n_total_states,
                          SERIAL_CTRL_NONCOOP1, SERIAL_CTRL_NONCOOP2>;

using Controller1 = SERIAL_CTRL_NONCOOP1;
using Controller2 = SERIAL_CTRL_NONCOOP2;

#endif  // ifdef CONTROLLER_TYPE_NCOOP

#ifdef CONTROLLER_TYPE_NCOOP_UNSTABLE
// using NvCtr =
    // NerveCenter<SerialCompressors, n_total_states,
                // SERIAL_CTRL_NONCOOP1, SERIAL_CTRL_NONCOOP2>;

// using Controller1 = SERIAL_CTRL_NONCOOP1;
// using Controller2 = SERIAL_CTRL_NONCOOP2;

#endif  // ifdef CONTROLLER_TYPE_NCOOP_UNSTABLE

#endif  // ifdef SYSTEM_TYPE_SERIAL

using namespace ReadFiles;

#endif
