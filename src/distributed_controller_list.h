#ifndef DISTRIBUTED_CONTROLLER_LIST_H
#define DISTRIBUTED_CONTROLLER_LIST_H

#include "distributed_controller.h"
#include "parallel_compressors.h"

template class DistributedController<ParallelCompressors, 40, 2, 100, 2, 2>;


#endif
