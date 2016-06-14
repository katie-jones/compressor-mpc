#ifndef OBSERVER_LIST_H
#define OBSERVER_LIST_H
#include "parallel_compressors.h"
#include "constexpr_array.h"
#include "aug_lin_sys.h"

extern template class AugmentedLinearizedSystem<
    ParallelCompressors, ConstexprArray<0, 40, 0, 40>, 4>;

template class Observer<
    AugmentedLinearizedSystem<ParallelCompressors, ConstexprArray<0, 40, 0, 40>,
                              4>>;

#endif
