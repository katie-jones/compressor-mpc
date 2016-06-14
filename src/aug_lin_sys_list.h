#ifndef AUG_LIN_SYS_LIST_H
#define AUG_LIN_SYS_LIST_H
#include "parallel_compressors.h"
#include "prediction.h"
#include "constexpr_array.h"

// Parallel compressors
template class AugmentedLinearizedSystem<ParallelCompressors,
                                         ConstexprArray<0, 40, 0, 40>, 4,
                                         ConstexprArray<0, 1, 2, 3>, 2>;
template class AugmentedLinearizedSystem<ParallelCompressors,
                                         ConstexprArray<0, 40, 0, 40>, 4,
                                         ConstexprArray<2, 3, 0, 1>, 2>;

#endif
