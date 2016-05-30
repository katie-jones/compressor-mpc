#ifndef AUG_LIN_SYS_LIST_H
#define AUG_LIN_SYS_LIST_H
// #include "compressor.h"
#include "parallel_compressors.h"

// Single compressor
// template class AugmentedLinearizedSystem<Compressor, 5, 2>;

// Parallel compressors
template class AugmentedLinearizedSystem<ParallelCompressors, 80, 4>;

#endif
