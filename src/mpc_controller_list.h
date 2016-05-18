#ifndef MPC_CONTROLLER_LIST_H
#define MPC_CONTROLLER_LIST_H
#include "compressor.h"

template class MpcController<Compressor, 40, 2, 100, 2>;
template class MpcController<Compressor, 40, 2, 5, 2>;
template class MpcController<Compressor, 5, 2, 5, 2>;
template class MpcController<Compressor, 2, 2, 5, 2>;

#endif
