#ifndef MPC_QP_SOLVER_LIST_H
#define MPC_QP_SOLVER_LIST_H

#include "mpc_qp_solver.h"
#include "parallel_compressors_constants.h"

template class MpcQpSolver<
    PARALLEL_COMPRESSORS_CONSTANTS::n_total_states,
    ParallelCompressors::n_outputs, ParallelCompressors::n_control_inputs,
    PARALLEL_COMPRESSORS_CONSTANTS::p, PARALLEL_COMPRESSORS_CONSTANTS::m>;

template class MpcQpSolver<PARALLEL_COMPRESSORS_CONSTANTS::n_total_states,
                           PARALLEL_COMPRESSORS_CONSTANTS::n_sub_outputs,
                           PARALLEL_COMPRESSORS_CONSTANTS::n_sub_control_inputs,
                           PARALLEL_COMPRESSORS_CONSTANTS::p,
                           PARALLEL_COMPRESSORS_CONSTANTS::m>;

#endif
