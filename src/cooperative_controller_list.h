#ifndef COOPERATIVE_CONTROLLER_LIST_H
#define COOPERATIVE_CONTROLLER_LIST_H
#include "parallel_compressors.h"
#include "cooperative_controller.h"

extern template class AugmentedLinearizedSystem<ParallelCompressors, 80, 4>;
extern template class Observer<ParallelCompressors, 80, 4>;
extern template class MpcQpSolver<95, 4, 2, 100, 2>;
template class CooperativeController<ParallelCompressors, 80, 4, 100, 2,
                                            2>;


#endif 
