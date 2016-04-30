#include <iostream>
#include "compressor.h"

int main(void) {
  Compressor x;

  Compressor::State xin = Compressor::GetDefaultState();
  Compressor::Linearized linsys = x.GetLinearizedSystem(xin, Compressor::GetDefaultInput());

  std::cout << xin(0) << std::endl;
  std::cout << linsys.A(0,0) << std::endl;
  std::cout << linsys.B(3,1) << std::endl;

  return 0;
}
