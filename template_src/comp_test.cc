#include <iostream>
#include "compressor.h"
#include "simulation_compressor.h"

void Callback(Compressor::State x, double t) {
  std::cout << "t: " << t << "\tx: ";
  for (int i=0; i<5; i++) 
    std::cout << x[i] << "\t";
  std::cout << std::endl;
}

int main(void) {
  Compressor x;

  Compressor::State xin = Compressor::GetDefaultState();
  Compressor::Linearized linsys = x.GetLinearizedSystem(xin, Compressor::GetDefaultInput());

  std::cout << xin(0) << std::endl;
  std::cout << linsys.A(0,0) << std::endl;
  std::cout << linsys.B(3,1) << std::endl;

  SimulationCompressor y;
  y.Integrate(0, 1, 0.05, &Callback);

  return 0;
}
