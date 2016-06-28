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
#include "read_files.h"
#include "parallel_compressors.h"
#include "parallel_compressors_constants.h"
#include "simulation_system.h"

using namespace PARALLEL_COMPRESSORS_CONSTANTS;
using namespace ReadFiles;

using SimSystem = SimulationSystem<ParallelCompressors, Delays, InputIndices>;

using NvCtr = NerveCenter<ParallelCompressors, n_total_states, CONTROLLER_COOP1,
                          CONTROLLER_COOP2>;

using AugmentedSystem1 = AUGMENTEDSYSTEM_DIST1;
using AugmentedSystem2 = AUGMENTEDSYSTEM_DIST2;

using Obsv1 = OBSERVER_DIST1;
using Obsv2 = OBSERVER_DIST2;

using Controller1 = CONTROLLER_COOP1;
using Controller2 = CONTROLLER_COOP2;

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

int main(int argc, char **argv) {
  int n_solver_iterations;

  if (argc < 2) {
    std::cout << "Number of solver iterations: ";
    std::cin >> n_solver_iterations;
  } else {
    std::istringstream ss(argv[1]);
    if (!(ss >> n_solver_iterations))
      std::cerr << "Invalid number " << argv[1] << '\n';
  }

  const std::string folder_name = "parallel/";

  const std::string constraints_fname = folder_name + "dist_constraints";
  const std::string output_fname = folder_name + "output/coop_output" +
                                      std::to_string(n_solver_iterations) +
                                      ".dat";
  const std::string info_fname = folder_name + "output/coop_info" +
                                      std::to_string(n_solver_iterations) +
                                      ".dat";
  const std::string yref_fname = folder_name + "yref";
  const std::string ywt_fname = folder_name + "yweight_coop";
  const std::string uwt_fname = folder_name + "uweight_coop";


  std::cout << "Running cooperative simulation using " << n_solver_iterations
            << " solver iterations... " << std::endl;
  std::cout.flush();

  // Time entire simulation
  boost::timer::cpu_timer simulation_timer;

  output_file.open(output_fname);

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

  const AugmentedSystem1::Input u_offset = u_default;

  // Weights
  NvCtr::UWeightType uwt = NvCtr::UWeightType::Zero();
  NvCtr::YWeightType ywt = NvCtr::YWeightType::Zero();

  if (!(ReadDataFromFile(uwt.data(), uwt.rows(), uwt_fname, uwt.rows() + 1))) {
    return -1;
  }
  if (!(ReadDataFromFile(ywt.data(), ywt.rows(), ywt_fname, ywt.rows() + 1))) {
    return -1;
  }

  // Read in reference output
  ParallelCompressors::Output y_ref_sub;
  if (!(ReadDataFromFile(y_ref_sub.data(), y_ref_sub.size(), yref_fname))) {
    return -1;
  }
  const NvCtr::OutputPrediction y_ref = y_ref_sub.replicate<Controller1::p, 1>();

  // Input constraints
  InputConstraints<Controller1::n_control_inputs> constraints;
  constraints.use_rate_constraints = true;
  if (!(ReadConstraintsFromFile(&constraints, constraints_fname))) {
    return -1;
  }

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
  info_file.open(info_fname);
  info_file << uwt << std::endl << ywt << std::endl << y_ref;
  info_file.close();
  return 0;
}
