#ifndef MPC_CONTROLLER_LIST_H
#define MPC_CONTROLLER_LIST_H
#include "augmented_system.h"
#include "compressor.h"

template class MpcController<AugmentedSystem<Compressor, 2, 3>, 5, 2>;

#endif
