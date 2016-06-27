#include <boost/timer/timer.hpp>
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

constexpr int n_solver_iterations = 1;

using namespace PARALLEL_COMPRESSORS_CONSTANTS;

using SimSystem = SimulationSystem<ParallelCompressors, Delays, InputIndices>;

using AugmentedSystem = AUGMENTEDSYSTEM_CENTRALIZED;

using Obsv = OBSERVER_CENTRALIZED;

using Controller = CONTROLLER_CENTRALIZED;

using NvCtr = NerveCenter<ParallelCompressors, n_total_states, Controller>;

SimSystem *p_sim_compressor;
ParallelCompressors *p_compressor;
NvCtr *p_controller;
std::ofstream output_file;

void Callback(ParallelCompressors::State x, double t) {
  output_file << t << std::endl;
  output_file << x.transpose() << std::endl;

  ParallelCompressors::Output y = p_compressor->GetOutput(x);

  output_file << y.transpose() << std::endl;

  // Get and apply next input
  NvCtr::ControlInput u =
      p_controller->GetNextInput(p_compressor->GetOutput(x));

  p_sim_compressor->SetInput(u);

  output_file << u.transpose() << std::endl
              << std::endl;
}

int main(void) {

  std::cout << "Running centralized simulation... ";
  std::cout.flush();

  output_file.open("parallel/output/cent_output.dat");

  // Time entire simulation
  boost::timer::cpu_timer simulation_timer;

  ParallelCompressors compressor;
  p_compressor = &compressor;

  ParallelCompressors::Input u_default = ParallelCompressors::GetDefaultInput();
  ParallelCompressors::State x_init = ParallelCompressors::GetDefaultState();

  SimSystem sim_comp(p_compressor, u_default, x_init);
  p_sim_compressor = &sim_comp;

  const double sampling_time = 0.05;

  const Obsv::ObserverMatrix M =
      (Obsv::ObserverMatrix() << Eigen::Matrix<double, compressor.n_states,
                                               compressor.n_outputs>::Zero(),
       Eigen::Matrix<double, n_disturbance_states,
                     compressor.n_outputs>::Identity()).finished();

  NvCtr::UWeightType uwt = NvCtr::UWeightType::Zero();
  NvCtr::YWeightType ywt = NvCtr::YWeightType::Zero();

  std::ifstream read_file;
  read_file.open("parallel/uweight_cent");
  for (int i = 0; i < uwt.rows(); i++) {
    if (!(read_file >> uwt(i, i))) {
      std::cerr << "Error reading input weight from file parallel/uweight_cent"
                << std::endl;
      return -1;
    }
  }
  read_file.close();

  read_file.open("parallel/yweight_cent");
  for (int i = 0; i < ywt.rows(); i++) {
    if (!(read_file >> ywt(i, i))) {
      std::cerr << "Error reading output weight from file parallel/yweight_cent"
                << std::endl;
      return -1;
    }
  }
  read_file.close();

  const AugmentedSystem::Input u_offset = u_default;

  // Read in reference output
  ParallelCompressors::Output y_ref_sub;
  read_file.open("parallel/yref");
  for (int i = 0; i < y_ref_sub.size(); i++) {
    if (!(read_file >> y_ref_sub(i))) {
      std::cerr << "Error reading reference output from file parallel/yref"
                << std::endl;
      return -1;
    }
  }
  read_file.close();

  const NvCtr::OutputPrediction y_ref = y_ref_sub.replicate<Controller::p, 1>();

  // Input constraints
  InputConstraints<Controller::n_control_inputs> constraints;
  constraints.use_rate_constraints = true;
  read_file.open("parallel/constraints");

  for (int i = 0; i < Controller::n_control_inputs; i++) {
    if (!(read_file >> constraints.lower_bound(i))) {
      std::cerr
          << "Error reading input constraints from file parallel/constraints"
          << std::endl;
      return -1;
    }
  }
  for (int i = 0; i < Controller::n_control_inputs; i++) {
    if (!(read_file >> constraints.upper_bound(i))) {
      std::cerr
          << "Error reading input constraints from file parallel/constraints"
          << std::endl;
      return -1;
    }
  }
  for (int i = 0; i < Controller::n_control_inputs; i++) {
    if (!(read_file >> constraints.lower_rate_bound(i))) {
      std::cerr
          << "Error reading input constraints from file parallel/constraints"
          << std::endl;
      return -1;
    }
  }
  for (int i = 0; i < Controller::n_control_inputs; i++) {
    if (!(read_file >> constraints.upper_rate_bound(i))) {
      std::cerr
          << "Error reading input constraints from file parallel/constraints"
          << std::endl;
      return -1;
    }
  }

  // Setup controller
  AugmentedSystem sys(compressor, sampling_time);
  Controller ctrl(sys, constraints, M);

  // Create a nerve center
  std::tuple<Controller> ctrl_tuple(ctrl);
  NvCtr nerve_center(compressor, ctrl_tuple, n_solver_iterations);
  p_controller = &nerve_center;

  // Test functions
  nerve_center.SetWeights(uwt, ywt);
  nerve_center.SetOutputReference(y_ref);

  nerve_center.Initialize(compressor.GetDefaultState(),
                          AugmentedSystem::ControlInput::Zero(),
                          compressor.GetDefaultInput(),
                          compressor.GetOutput(compressor.GetDefaultState()));

  // Integrate system
  sim_comp.Integrate(0, 50, sampling_time, &Callback);

  // Apply disturbance
  ParallelCompressors::Input u_disturbance = u_default;
  u_disturbance(8) -= 0.3;

  sim_comp.SetOffset(u_disturbance);

  sim_comp.Integrate(50 + sampling_time, 500, sampling_time, &Callback);

  output_file.close();

  boost::timer::cpu_times simulation_cpu_time = simulation_timer.elapsed();
  boost::timer::nanosecond_type simulation_ns(simulation_cpu_time.system +
                                              simulation_cpu_time.user);

  std::cout << "Finished." << std::endl
            << "Total time required:\t"
            << static_cast<double>(simulation_ns) / 1.0e6 << " ms." << std::endl
            << std::endl;

  std::ofstream info_file;
  info_file.open("parallel/output/cent_info.dat");
  info_file << uwt << std::endl << ywt << std::endl << y_ref;
  info_file.close();

  return 0;
}
