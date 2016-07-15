#ifndef DISTRIBUTED_CONTROLLER_LIST_H
#define DISTRIBUTED_CONTROLLER_LIST_H

#include "distributed_controller.h"
#include "parallel_compressors.h"
#include "serial_compressors.h"
#include "constexpr_array.h"
#include "null_index_array.h"
#include "aug_lin_sys.h"
#include "parallel_compressors_constants.h"
#include "serial_compressors_constants.h"

// Parallel
template class PARALLEL_CTRL_COOP1;
template class PARALLEL_CTRL_COOP2;

template class PARALLEL_CTRL_NONCOOP1;
template class PARALLEL_CTRL_NONCOOP2;

template class PARALLEL_CTRL_CENT;

// Serial
template class SERIAL_CTRL_COOP1;
template class SERIAL_CTRL_COOP2;

template class SERIAL_CTRL_NONCOOP1;
template class SERIAL_CTRL_NONCOOP2;

template class SERIAL_CTRL_CENT;

#endif
