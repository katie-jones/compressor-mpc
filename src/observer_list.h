#ifndef OBSERVER_LIST_H
#define OBSERVER_LIST_H
#include "parallel_compressors.h"
#include "serial_compressors.h"
#include "constexpr_array.h"
#include "aug_lin_sys.h"
#include "parallel_compressors_constants.h"
#include "serial_compressors_constants.h"

template class OBSERVER_DIST1;
template class OBSERVER_DIST2;
template class OBSERVER_CENTRALIZED;

template class SERIAL_OBS_DIST1;
template class SERIAL_OBS_DIST2;
template class SERIAL_OBS_CENT;

extern template class SERIAL_AUGSYS_DIST1;
extern template class SERIAL_AUGSYS_DIST2;
extern template class SERIAL_AUGSYS_CENT;


#endif
