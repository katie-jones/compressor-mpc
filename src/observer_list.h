#ifndef OBSERVER_LIST_H
#define OBSERVER_LIST_H
#include "parallel_compressors.h"
#include "constexpr_array.h"
#include "aug_lin_sys.h"

template class Observer<
    AugmentedLinearizedSystem<ParallelCompressors, ConstexprArray<0, 40, 0, 40>,
                              4, ConstexprArray<0, 1, 2, 3>, 2>>;
template class Observer<
    AugmentedLinearizedSystem<ParallelCompressors, ConstexprArray<0, 40, 0, 40>,
                              4, ConstexprArray<2, 3, 0, 1>, 2>>;

#endif
