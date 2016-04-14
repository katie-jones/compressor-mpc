#ifndef CONST_SIM_H
#define CONST_SIM_H
#include "defs.h"
const double Ts = 0.05;
const comp_input uoff1((comp_input() << 0.304, 0.43, 1, 0, 0).finished());
const comp_input uoff2 = uoff1;

const double ud = 0.7;

const int ysize = 2;

#endif
