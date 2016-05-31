#ifndef MPC_CONTROLLER_LIST_H
#define MPC_CONTROLLER_LIST_H

#include "parallel_compressors.h"
#include "aug_lin_sys.h"
#include "observer.h"
#include "mpc_qp_solver.h"

// Single compressor
// template class MpcController<Compressor, 40, 2, 100, 2>;

// template class MpcController<Compressor, 5, 2, 100, 2>;
// template class MpcController<Compressor, 40, 2, 5, 2>;
// template class MpcController<Compressor, 5, 2, 5, 2>;
// template class MpcController<Compressor, 2, 2, 5, 2>;

// Parallel compressors
extern template class AugmentedLinearizedSystem<ParallelCompressors, 80, 4>;
extern template class Observer<ParallelCompressors, 80, 4>;
template class MpcQpSolver<95, 4, 4, 100, 2>;
template class MpcController<ParallelCompressors, 80, 4, 100, 2>;

#endif
