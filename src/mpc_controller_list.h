#ifndef MPC_CONTROLLER_LIST_H
#define MPC_CONTROLLER_LIST_H
// #include "compressor.h"
#include "parallel_compressors.h"

// Single compressor
// template class MpcController<Compressor, 40, 2, 100, 2>;

// extern template class AugmentedLinearizedSystem<Compressor, 5, 2>;
extern template class AugmentedLinearizedSystem<ParallelCompressors, 80, 4>;

// template class MpcController<Compressor, 5, 2, 100, 2>;
// template class MpcController<Compressor, 40, 2, 5, 2>;
// template class MpcController<Compressor, 5, 2, 5, 2>;
// template class MpcController<Compressor, 2, 2, 5, 2>;

// Parallel compressors
template class MpcController<ParallelCompressors, 80, 4, 100, 2>;

#endif
