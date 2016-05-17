#include <iostream>
#include <fstream>
#include "compressor.h"
#include "simulation_compressor.h"
#include "tank.h"
#include "parallel_compressors.h"
#include "simulation_parallel_compressors.h"
#include "mpc_controller.h"
#include "print_matrix.h"

extern template class MpcController<Compressor, 2, 2, 5, 2>;
using Controller = MpcController<Compressor, 2, 2, 5, 2>;

namespace Control {
constexpr int n_delay_states = 2;
constexpr int n_disturbance_states = 2;
constexpr int p = 5;
constexpr int m = 2;
}

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

  // Compressor::State xin = Compressor::GetDefaultState();
  // Compressor::Linearized linsys =
  // x.GetLinearizedSystem(xin, Compressor::GetDefaultInput());

  // std::cout << xin(0) << std::endl;
  // std::cout << linsys.A(0, 0) << std::endl;
  // std::cout << linsys.B(3, 1) << std::endl;

  // SimulationCompressor y;
  // y.Integrate(0, 1, 0.05, &Callback);

  // Tank z;
  // Tank::State xt;
  // xt << 1.12;
  // Tank::Input yt;
  // yt << 0.7, 1, 0.3;
  // Tank::Linearized zlin = z.GetLinearizedSystem(xt, yt);
  // // std::cout << zlin.A(0) << std::endl;

  // ParallelCompressors parcomp = ParallelCompressors();
  // ParallelCompressors::Linearized plin =
  // parcomp.GetLinearizedSystem(ParallelCompressors::GetDefaultState(),
  // ParallelCompressors::GetDefaultInput());
  // std::cout << plin.A(10, 1) << std::endl;

  // std::ifstream finaltime("finaltime.txt");
  // double tf;
  // finaltime >> tf;
  // SimulationParallelCompressors compsys = SimulationParallelCompressors();
  // compsys.Integrate(0, tf, 0.05, CallbackSys);

  // std::array<int, 2> n_delay = {0, 3};
  // std::array<int, 2> control_input_index = {0, 3};

  // std::array<double, 2*2*2*3*3> A = MpcController<Compressor, 4, 2, 20,
  // 3>::GetConstraintMatrix();

  // std::cout << "A = " << std::endl;
  // for (int i=0; i<2*3*2; i++) {
  // for (int j=0; j<2*3; j++) {
  // std::cout << A[j+6*i] << "\t";
  // }
  // std::cout << std::endl;
  // }

  // Controller ctrl(augsys);
  // Controller::Prediction pred = ctrl.GetPredictionMatrices();
  const Controller::ObserverMatrix M =
      (Controller::ObserverMatrix()
           << Eigen::Matrix<double, x.n_states,
                            x.n_control_inputs>::Zero(),
       Eigen::Matrix<double, Control::n_delay_states,
                     Control::n_delay_states>::Identity()).finished();
  const Controller::ControlInputIndex delay = {0, 2};
  const Controller::ControlInputIndex index = {0, 3};

  Controller::InputConstraints constraints;
  constraints.lower_bound << -0.3, 0;
  constraints.upper_bound << 0.3, 1;
  constraints.lower_rate_bound << -0.1, -0.1;
  constraints.upper_rate_bound << 0.1, 1;

  Controller ctrl(x, M, 0.05, delay, index);
  ctrl.SetInitialState(x.GetDefaultState(), Controller::ControlInput::Zero());
  Compressor::Input u = ctrl.GetNextInput(Compressor::Output::Zero());
  print_matrix(std::cout, u.matrix(), "u");

  return 0;
}
