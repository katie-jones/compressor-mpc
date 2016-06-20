#ifndef DISTRIBUTED_CONTROLLER_LIST_H
#define DISTRIBUTED_CONTROLLER_LIST_H

#include "distributed_controller.h"
#include "parallel_compressors.h"
#include "constexpr_array.h"
#include "null_index_array.h"
#include "aug_lin_sys.h"
#include "parallel_compressors_constants.h"

extern template class AUGMENTEDSYSTEM1;
extern template class AUGMENTEDSYSTEM2;

template int CONTROLLER1::InitializeFull<
    ParallelCompressors::n_states, ParallelCompressors::n_outputs,
    PARALLEL_COMPRESSORS_CONSTANTS::n_total_states>(
    const Eigen::Matrix<double, ParallelCompressors::n_states, 1>& x_init_full,
    const FullControlInput& u_init_full, const Input& full_u_old,
    const Eigen::Matrix<double, ParallelCompressors::n_outputs, 1>& y_init_full,
    const Eigen::Matrix<double, PARALLEL_COMPRESSORS_CONSTANTS::n_total_states,
                        1> dx_init_full);

template int CONTROLLER2::InitializeFull<
    ParallelCompressors::n_states, ParallelCompressors::n_outputs,
    PARALLEL_COMPRESSORS_CONSTANTS::n_total_states>(
    const Eigen::Matrix<double, ParallelCompressors::n_states, 1>& x_init_full,
    const FullControlInput& u_init_full, const Input& full_u_old,
    const Eigen::Matrix<double, ParallelCompressors::n_outputs, 1>& y_init_full,
    const Eigen::Matrix<double, PARALLEL_COMPRESSORS_CONSTANTS::n_total_states,
                        1> dx_init_full);

template class CONTROLLER1; 
template class CONTROLLER2;


#endif
