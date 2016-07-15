#ifndef MPC_QP_SOLVER_LIST_H
#define MPC_QP_SOLVER_LIST_H

#include "mpc_qp_solver.h"
#include "parallel_compressors_constants.h"
#include "serial_compressors_constants.h"

// Parallel
template class MpcQpSolver<PARALLEL_COMPRESSORS_CONSTANTS::n_total_states,
                           PARALLEL_COMPRESSORS_CONSTANTS::n_sub_outputs,
                           ParallelCompressors::n_control_inputs,
                           PARALLEL_COMPRESSORS_CONSTANTS::p,
                           PARALLEL_COMPRESSORS_CONSTANTS::m>;

template class MpcQpSolver<PARALLEL_COMPRESSORS_CONSTANTS::n_total_states,
                           PARALLEL_COMPRESSORS_CONSTANTS::n_sub_outputs,
                           PARALLEL_COMPRESSORS_CONSTANTS::n_sub_control_inputs,
                           PARALLEL_COMPRESSORS_CONSTANTS::p,
                           PARALLEL_COMPRESSORS_CONSTANTS::m>;

template class MpcQpSolver<PARALLEL_COMPRESSORS_CONSTANTS::n_total_states,
                           PARALLEL_COMPRESSORS_CONSTANTS::n_sub_outputs_nc,
                           PARALLEL_COMPRESSORS_CONSTANTS::n_sub_control_inputs,
                           PARALLEL_COMPRESSORS_CONSTANTS::p,
                           PARALLEL_COMPRESSORS_CONSTANTS::m>;

// Serial
template class MpcQpSolver<SERIAL_COMPRESSORS_CONSTANTS::n_total_states,
                           SERIAL_COMPRESSORS_CONSTANTS::n_sub_outputs,
                           SerialCompressors::n_control_inputs,
                           SERIAL_COMPRESSORS_CONSTANTS::p,
                           SERIAL_COMPRESSORS_CONSTANTS::m>;

template class MpcQpSolver<SERIAL_COMPRESSORS_CONSTANTS::n_total_states,
                           SERIAL_COMPRESSORS_CONSTANTS::n_sub_outputs,
                           SERIAL_COMPRESSORS_CONSTANTS::n_sub_control_inputs,
                           SERIAL_COMPRESSORS_CONSTANTS::p,
                           SERIAL_COMPRESSORS_CONSTANTS::m>;

template class MpcQpSolver<SERIAL_COMPRESSORS_CONSTANTS::n_total_states,
                           SERIAL_COMPRESSORS_CONSTANTS::n_sub_outputs_nc,
                           SERIAL_COMPRESSORS_CONSTANTS::n_sub_control_inputs,
                           SERIAL_COMPRESSORS_CONSTANTS::p,
                           SERIAL_COMPRESSORS_CONSTANTS::m>;

template class MpcQpSolver<SERIAL_COMPRESSORS_CONSTANTS::n_total_states, 2,
                           SERIAL_COMPRESSORS_CONSTANTS::n_sub_control_inputs,
                           SERIAL_COMPRESSORS_CONSTANTS::p,
                           SERIAL_COMPRESSORS_CONSTANTS::m>;
#endif
