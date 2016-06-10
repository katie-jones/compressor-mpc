#ifndef NONCOOPERATIVE_CONTROLLER_LIST_H
#define NONCOOPERATIVE_CONTROLLER_LIST_H
#include "parallel_compressors.h"
#include "noncooperative_controller.h"
#include "constexpr_array.h"

extern template class AugmentedLinearizedSystem<
    ParallelCompressors, ConstexprArray<0, 40, 0, 40>, 4>;
extern template class Observer<ParallelCompressors,
                               ConstexprArray<0, 40, 0, 40>, 4>;

extern template class MpcQpSolver<95, 2, 2, 100, 2>;
extern template class MpcQpSolver<95, 3, 2, 100, 2>;

template class NonCooperativeController<
    ParallelCompressors, ConstexprArray<0, 40, 0, 40>, 4, 100, 2, 2, 2>;

template class NonCooperativeController<
    ParallelCompressors, ConstexprArray<0, 40, 0, 40>, 4, 100, 2, 2, 3>;

#endif
