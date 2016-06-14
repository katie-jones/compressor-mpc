#ifndef DISTRIBUTED_CONTROLLER_LIST_H
#define DISTRIBUTED_CONTROLLER_LIST_H

#include "distributed_controller.h"
#include "parallel_compressors.h"
#include "constexpr_array.h"
#include "aug_lin_sys.h"

extern template class AugmentedLinearizedSystem<
    ParallelCompressors, ConstexprArray<0, 40, 0, 40>, 4>;
template class DistributedController<
    AugmentedLinearizedSystem<ParallelCompressors, ConstexprArray<0, 40, 0, 40>,
                              4>,
    100, 2>;

#endif
