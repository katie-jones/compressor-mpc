#include "tank.h"

const Tank::TankInput Tank::default_input = 0.7;
const Tank::TankState Tank::default_initial_state = 1.12;

Tank::TankState Tank::GetDerivative(const TankState x,
                                const TankInput u) const {
  TankState dxdt = x;
  return dxdt;
}

