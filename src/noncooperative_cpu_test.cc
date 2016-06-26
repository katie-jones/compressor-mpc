#include <boost/timer/timer.hpp>
#include <sstream>
#include <fstream>
#include <iostream>
#include "aug_lin_sys.h"
#include "constexpr_array.h"
#include "distributed_controller.h"
#include "input_constraints.h"
#include "nerve_center.h"
#include "null_index_array.h"
#include "observer.h"
#include "parallel_compressors.h"
#include "parallel_compressors_constants.h"
#include "simulation_system.h"

using namespace PARALLEL_COMPRESSORS_CONSTANTS;

using SimSystem = SimulationSystem<ParallelCompressors, Delays, InputIndices>;

using NvCtr = NerveCenter<ParallelCompressors, n_total_states,
                          CONTROLLER_NONCOOP1, CONTROLLER_NONCOOP2>;

using AugmentedSystem1 = AUGMENTEDSYSTEM_DIST1;
using AugmentedSystem2 = AUGMENTEDSYSTEM_DIST2;

using Obsv1 = OBSERVER_DIST1;
using Obsv2 = OBSERVER_DIST2;

using Controller1 = CONTROLLER_NONCOOP1;
using Controller2 = CONTROLLER_NONCOOP2;

SimSystem *p_sim_compressor;
ParallelCompressors *p_compressor;
NvCtr *p_controller;
std::ofstream output_file;
std::ofstream cpu_times_file;

boost::timer::cpu_timer timer;

void Callback(ParallelCompressors::State x, double t) {
  output_file << t << std::endl;
  output_file << x.transpose() << std::endl;

  ParallelCompressors::Output y = p_compressor->GetOutput(x);

  output_file << y.transpose() << std::endl;

  // Get and apply next input
  timer.resume();
  NvCtr::ControlInput u =
      p_controller->GetNextInput(p_compressor->GetOutput(x));
  timer.stop();

  boost::timer::cpu_times elapsed = timer.elapsed();
  boost::timer::nanosecond_type elapsed_ns(elapsed.system + elapsed.user);

  cpu_times_file << elapsed_ns/2.0 << std::endl;

  p_sim_compressor->SetInput(u);

  output_file << u.transpose() << std::endl
              << std::endl;
}

int main(int argc, char **argv) {
  timer.stop();

  int n_solver_iterations;

  if (argc < 2) {
    std::cout << "Number of solver iterations: ";
    std::cin >> n_solver_iterations;
  } else {
    std::istringstream ss(argv[1]);
    if (!(ss >> n_solver_iterations))
      std::cerr << "Invalid number " << argv[1] << '\n';
  }

  std::cout << "Running non-cooperative simulation using "
            << n_solver_iterations << " solver iterations... ";

  // Time entire simulation
  boost::timer::cpu_timer simulation_timer;

  boost::timer::cpu_times time_offset = timer.elapsed();
  boost::timer::nanosecond_type offset_ns(time_offset.system + time_offset.user);

  output_file.open("noncoop_output" + std::to_string(n_solver_iterations) +
                   ".dat");
  cpu_times_file.open("noncoop_cpu_times" +
                      std::to_string(n_solver_iterations) + ".dat");

  cpu_times_file << offset_ns << std::endl;

  ParallelCompressors compressor;
  p_compressor = &compressor;

  ParallelCompressors::Input u_default = ParallelCompressors::GetDefaultInput();
  ParallelCompressors::State x_init = ParallelCompressors::GetDefaultState();

  SimSystem sim_comp(p_compressor, u_default, x_init);
  p_sim_compressor = &sim_comp;

  const double sampling_time = 0.05;

  const Obsv1::ObserverMatrix M =
      (Obsv1::ObserverMatrix() << Eigen::Matrix<double, compressor.n_states,
                                                compressor.n_outputs>::Zero(),
       Eigen::Matrix<double, n_disturbance_states,
                     compressor.n_outputs>::Identity()).finished();

  NvCtr::UWeightType uwt = NvCtr::UWeightType::Zero();
  NvCtr::YWeightType ywt = NvCtr::YWeightType::Zero();

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

  const AugmentedSystem1::Input u_offset = u_default;

  const NvCtr::OutputPrediction y_ref =
      (NvCtr::Output() << 4.5, 4.5, 0, 1.12).finished().replicate<p, 1>();

  // Input constraints
  InputConstraints<n_sub_control_inputs> constraints;
  constraints.lower_bound << -0.3, 0;
  constraints.upper_bound << 0.3, 1;
  constraints.lower_rate_bound << -0.1, -0.1;
  constraints.upper_rate_bound << 0.1, 1;
  constraints.use_rate_constraints = true;

  // Setup controller
  AugmentedSystem1 sys1(compressor, sampling_time);
  AugmentedSystem2 sys2(compressor, sampling_time);
  Controller1 ctrl1(sys1, constraints, M);
  Controller2 ctrl2(sys2, constraints, M);

  // Create a nerve center
  std::tuple<Controller1, Controller2> ctrl_tuple(ctrl1, ctrl2);
  NvCtr nerve_center(compressor, ctrl_tuple, n_solver_iterations);
  p_controller = &nerve_center;

  // Test functions
  nerve_center.SetWeights(uwt, ywt);
  nerve_center.SetOutputReference(y_ref);

  nerve_center.Initialize(compressor.GetDefaultState(),
                          AugmentedSystem1::ControlInput::Zero(),
                          compressor.GetDefaultInput(),
                          compressor.GetOutput(compressor.GetDefaultState()));

  // Integrate system
  sim_comp.Integrate(0, 50, sampling_time, &Callback);

  // Apply disturbance
  ParallelCompressors::Input u_disturbance = u_default;
  u_disturbance(8) -= 0.3;

  sim_comp.SetOffset(u_disturbance);

  sim_comp.Integrate(50 + sampling_time, 1000, sampling_time, &Callback);

  output_file.close();
  cpu_times_file.close();

  boost::timer::cpu_times simulation_cpu_time = simulation_timer.elapsed();
  boost::timer::nanosecond_type simulation_ns(simulation_cpu_time.system + simulation_cpu_time.user);

  std::cout << "Finished." << std::endl
            << "Total time required:\t"
            << static_cast<double>(simulation_ns) / 1.0e6 << " ms." << std::endl
            << std::endl;

  return 0;
}
