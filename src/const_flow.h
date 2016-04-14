#ifndef COMP_CONST_FLOW_H
#define COMP_CONST_FLOW_H

namespace comp {
namespace flow {
const double pi = 3.14159265358979323846;
const double a = 340.;
const double Pin = 1;
double Pout = -1;

const double V1 =
    2 * pi * (0.60 / 2) * (0.60 / 2) * 2 + pi * (0.08 / 2) * (0.08 / 2) * 8.191;
const double V2 = 0.5 * V1;

const double AdivL = pi * (0.08 / 2) * (0.08 / 2) / 3 * 0.1;
}
}

#endif

