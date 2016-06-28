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
#include "read_files.h"

constexpr int n_solver_iterations = 1;

using namespace PARALLEL_COMPRESSORS_CONSTANTS;
using namespace ReadFiles;

using SimSystem = SimulationSystem<ParallelCompressors, Delays, InputIndices>;

using AugmentedSystem = AUGMENTEDSYSTEM_CENTRALIZED;

using Obsv = OBSERVER_CENTRALIZED;

using Controller = CONTROLLER_CENTRALIZED;

using NvCtr = NerveCenter<ParallelCompressors, n_total_states, Controller>;

SimSystem *p_sim_compressor;
ParallelCompressors *p_compressor;
NvCtr *p_controller;
std::ofstream cpu_times_file;

boost::timer::cpu_timer timer;

void Callback(ParallelCompressors::State x, double t) {
  ParallelCompressors::Output y = p_compressor->GetOutput(x);

  // Get and apply next input
  timer.resume();
  NvCtr::ControlInput u =
      p_controller->GetNextInput(p_compressor->GetOutput(x));
  timer.stop();

  boost::timer::cpu_times elapsed = timer.elapsed();
  boost::timer::nanosecond_type elapsed_ns(elapsed.system + elapsed.user);

  cpu_times_file << elapsed_ns << std::endl;

  p_sim_compressor->SetInput(u);
}

int main(void) {
  timer.stop();

  const std::string folder_name = "parallel/";

  const std::string constraints_fname = folder_name + "constraints";
  const std::string cpu_times_fname = folder_name + "output/cent_cpu_times.dat";
  const std::string yref_fname = folder_name + "yref";
  const std::string ywt_fname = folder_name + "yweight_cent";
  const std::string uwt_fname = folder_name + "uweight_cent";


  boost::timer::cpu_times time_offset = timer.elapsed();
  boost::timer::nanosecond_type offset_ns(time_offset.system +
                                          time_offset.user);

  std::cout << "Running centralized simulation... ";
  std::cout.flush();

  cpu_times_file.open(cpu_times_fname);

  cpu_times_file << offset_ns << std::endl;

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

  const AugmentedSystem::Input u_offset = u_default;

  // Weights
  NvCtr::UWeightType uwt = NvCtr::UWeightType::Zero();
  NvCtr::YWeightType ywt = NvCtr::YWeightType::Zero();

  if (!(ReadDataFromFile(uwt.data(), uwt.rows(), uwt_fname,
                         uwt.rows() + 1))) {
    return -1;
  }
  if (!(ReadDataFromFile(ywt.data(), ywt.rows(), ywt_fname,
                         ywt.rows() + 1))) {
    return -1;
  }

  // Read in reference output
  ParallelCompressors::Output y_ref_sub;
  if (!(ReadDataFromFile(y_ref_sub.data(), y_ref_sub.size(),
                         yref_fname))) {
    return -1;
  }
  const NvCtr::OutputPrediction y_ref = y_ref_sub.replicate<Controller::p, 1>();

  // Input constraints
  InputConstraints<Controller::n_control_inputs> constraints;
  constraints.use_rate_constraints = true;
  if (!(ReadConstraintsFromFile(&constraints, constraints_fname))) {
    return -1;
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

  cpu_times_file.close();

  boost::timer::cpu_times simulation_cpu_time = simulation_timer.elapsed();
  boost::timer::nanosecond_type simulation_ns(simulation_cpu_time.system +
                                              simulation_cpu_time.user);

  std::cout << "Finished." << std::endl
            << "Total time required:\t"
            << static_cast<double>(simulation_ns) / 1.0e6 << " ms." << std::endl
            << std::endl;

  return 0;
}
