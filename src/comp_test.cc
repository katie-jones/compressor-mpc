#include <iostream>
#include <fstream>
#include "compressor.h"
#include "simulation_compressor.h"
#include "tank.h"
#include "parallel_compressors.h"
#include "simulation_parallel_compressors.h"
#include "augmented_system.h"

extern template class AugmentedSystem<Compressor, 2, 3>;
using AugSys = AugmentedSystem<Compressor, 2, 3>;

void Callback(Compressor::State x, double t) {
  std::cout << "t: " << t << "\tx: ";
  for (int i = 0; i < 5; i++) std::cout << x[i] << "\t";
  std::cout << std::endl;
}

void CallbackSys(ParallelCompressors::State x, double t) {
  std::cout << "t: " << t << "\tx: ";
  for (int i = 0; i < ParallelCompressors::n_states; i++)
    std::cout << x[i] << "\t";
  std::cout << std::endl;
}

int main(void) {
  Compressor x;

  Compressor::State xin = Compressor::GetDefaultState();
  Compressor::Linearized linsys =
      x.GetLinearizedSystem(xin, Compressor::GetDefaultInput());

  std::cout << xin(0) << std::endl;
  std::cout << linsys.A(0, 0) << std::endl;
  std::cout << linsys.B(3, 1) << std::endl;

  SimulationCompressor y;
  y.Integrate(0, 1, 0.05, &Callback);

  Tank z;
  Tank::State xt;
  xt << 1.12;
  Tank::Input yt;
  yt << 0.7, 1, 0.3;
  Tank::Linearized zlin = z.GetLinearizedSystem(xt, yt);
  // std::cout << zlin.A(0) << std::endl;

  ParallelCompressors parcomp = ParallelCompressors();
  ParallelCompressors::Linearized plin =
      parcomp.GetLinearizedSystem(ParallelCompressors::GetDefaultState(),
                                  ParallelCompressors::GetDefaultInput());
  std::cout << plin.A(10, 1) << std::endl;

  std::ifstream finaltime("finaltime.txt");
  double tf;
  finaltime >> tf;
  SimulationParallelCompressors compsys = SimulationParallelCompressors();
  compsys.Integrate(0, tf, 0.05, CallbackSys);

  std::array<int, 2> n_delay = {0, 3};
  std::array<int, 2> control_input_index = {0, 3};

  AugSys augsys =
      AugSys(x, 0.05, Compressor::GetDefaultState(), AugSys::ObserverMatrix(),
             n_delay, control_input_index);

  // AugSys::AugmentedLinearizedSystem auglin =
  // augsys.LinearizeAndAugment(linsys);
  std::cout << linsys << std::endl;
  // std::cout << auglin << std::endl;

  return 0;
}
