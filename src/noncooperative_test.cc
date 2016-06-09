#include <iostream>
#include <fstream>
#include <boost/timer/timer.hpp>
#include "parallel_compressors.h"
#include "simulation_system.h"
#include "aug_lin_sys.h"
#include "observer.h"
#include "print_matrix.h"
#include "input_constraints.h"
#include "noncooperative_controller.h"
#include "prediction.h"

namespace Control {
constexpr int n_delay_states = 80;
constexpr int n_disturbance_states = 4;
constexpr int p = 100;
constexpr int m = 2;
constexpr int n_controllers = 2;
constexpr int n_sub_outputs = 2;
}

extern template class AugmentedLinearizedSystem<ParallelCompressors, 80, 4>;
extern template class Observer<ParallelCompressors, 80, 4>;
extern template class NonCooperativeController<ParallelCompressors, 80, 4, 100,
                                               2, 2, 2>;
extern template class MpcQpSolver<95, 2, 2, 100, 2>;

extern template const Prediction
AugmentedLinearizedSystem<ParallelCompressors, 80, 4>::GenerateSubPrediction<2>(
    const int p, const int m, const int* output_index) const;

using AugmentedSystem =
    AugmentedLinearizedSystem<ParallelCompressors, Control::n_delay_states,
                              Control::n_disturbance_states>;

using Obsv = Observer<ParallelCompressors, Control::n_delay_states,
                      Control::n_disturbance_states>;

using Controller = NonCooperativeController<
    ParallelCompressors, Control::n_delay_states, Control::n_disturbance_states,
    Control::p, Control::m, Control::n_controllers, Control::n_sub_outputs>;

using SimSystem =
    SimulationSystem<ParallelCompressors, Control::n_delay_states>;

SimSystem *p_sim_compressor;
ParallelCompressors *p_compressor;
Controller *p_controller;
std::ofstream output_file;
std::ofstream cpu_times_file;

void Callback(ParallelCompressors::State x, double t) {
  output_file << t << std::endl;
  output_file << x.transpose() << std::endl;

  Controller::Output y = p_compressor->GetOutput(x);

  output_file << y.transpose() << std::endl;

  // Get and apply next input
  Controller::ControlInput u =
      p_controller->GetNextInput(p_compressor->GetOutput(x), cpu_times_file);
  p_sim_compressor->SetInput(u);

  output_file << u.transpose() << std::endl
              << std::endl;
}

int main(void) {
  output_file.open("noncoop_output.dat");
  cpu_times_file.open("noncoop_cpu_times.dat");

  ParallelCompressors compressor;
  p_compressor = &compressor;

  ParallelCompressors::Input u_default = ParallelCompressors::GetDefaultInput();
  ParallelCompressors::State x_init = ParallelCompressors::GetDefaultState();

  // index of controlled inputs
  const Controller::ControlInputIndex index = {0, 3, 4, 7};
  const Controller::ControlInputIndex delay = {0, Control::n_delay_states / 2,
                                               0, Control::n_delay_states / 2};

  // index of outputs per subcontroller
  const Eigen::Matrix<int, Control::n_controllers, Control::n_sub_outputs>
      output_index =
          (Eigen::Matrix<int, Control::n_controllers, Control::n_sub_outputs>()
               << 0,
           3, 1, 3).finished();

  SimSystem sim_comp(p_compressor, u_default, index, delay, x_init);
  p_sim_compressor = &sim_comp;

  const double sampling_time = 0.05;

  const Obsv::ObserverMatrix M =
      (Obsv::ObserverMatrix() << Eigen::Matrix<double, compressor.n_states,
                                               compressor.n_outputs>::Zero(),
       Eigen::Matrix<double, Control::n_disturbance_states,
                     compressor.n_outputs>::Identity()).finished();

  Controller::UWeightType uwt = Controller::UWeightType::Zero();
  Controller::YWeightType ywt = Controller::YWeightType::Zero();

  std::ifstream weight_file;
  weight_file.open("uweight_noncoop");
  for (int i = 0; i < uwt.rows(); i++) {
    weight_file >> uwt(i, i);
  }
  weight_file.close();

  weight_file.open("yweight_noncoop");
  for (int i = 0; i < ywt.rows(); i++) {
    weight_file >> ywt(i, i);
  }
  weight_file.close();

  const Controller::Input offset = u_default;

  const Controller::OutputPrediction y_ref =
      (Controller::Output() << 4.5, 4.5, 0, 1.12)
          .finished()
          .replicate<Control::p, 1>();

  // Input constraints
  InputConstraints<ParallelCompressors::n_control_inputs> constraints;
  constraints.lower_bound << -0.3, 0, -0.3, 0;
  constraints.upper_bound << 0.3, 1, 0.3, 1;
  // constraints.lower_bound << 0, 0, 0, 0;
  // constraints.upper_bound << 0, 0, 0, 0;
  constraints.lower_rate_bound << -0.1, -0.1, -0.1, -0.1;
  constraints.upper_rate_bound << 0.1, 1, 0.1, 1;
  constraints.use_rate_constraints = true;

  // Setup controller
  AugmentedSystem sys(compressor, sampling_time, delay);
  Obsv observer(M, compressor.GetOutput(x_init));

  const int n_solver_iterations = 3;

  Controller ctrl(sys, observer, u_default, y_ref, delay, index,
                  n_solver_iterations, output_index, constraints, uwt, ywt);
  p_controller = &ctrl;

  ctrl.SetInitialState(x_init, compressor.GetOutput(x_init));

  // Integrate system and time it
  boost::timer::cpu_timer integrate_timer;
  sim_comp.Integrate(0, 50, sampling_time, &Callback);
  integrate_timer.stop();

  // Apply disturbance
  ParallelCompressors::Input u_disturbance = u_default;
  u_disturbance(8) -= 0.3;

  sim_comp.SetOffset(u_disturbance);
  integrate_timer.resume();
  sim_comp.Integrate(50 + sampling_time, 500, sampling_time, &Callback);
  integrate_timer.stop();
  boost::timer::cpu_times int_elapsed = integrate_timer.elapsed();
  boost::timer::nanosecond_type elapsed_ns(int_elapsed.system +
                                           int_elapsed.user);
  output_file.close();
  cpu_times_file.close();
  std::cout << "CPU time: " << elapsed_ns << std::endl;
  std::cout << "Wall time: " << int_elapsed.wall << std::endl
            << std::endl;

  return 0;
}
