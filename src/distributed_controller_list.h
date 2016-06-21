#ifndef DISTRIBUTED_CONTROLLER_LIST_H
#define DISTRIBUTED_CONTROLLER_LIST_H

#include "distributed_controller.h"
#include "parallel_compressors.h"
#include "constexpr_array.h"
#include "null_index_array.h"
#include "aug_lin_sys.h"
#include "parallel_compressors_constants.h"

extern template class AUGMENTEDSYSTEM1;
extern template class AUGMENTEDSYSTEM2;

template class CONTROLLER1; 
template class CONTROLLER2;


#endif
