#ifndef OBSERVER_LIST_H
#define OBSERVER_LIST_H
#include "parallel_compressors.h"
#include "constexpr_array.h"

template class ConstexprArray<0,40,0,40>;
template class Observer<ParallelCompressors, ConstexprArray<0,40,0,40>, 4>;

#endif
