#undef CONTROLLER_TYPE_COOP
#undef CONTROLLER_TYPE_CENTRALIZED
#undef SYSTEM_TYPE_SERIAL

#define CONTROLLER_TYPE_NCOOP
#define SYSTEM_TYPE_PARALLEL

#include "common-variables.h"

SimSystem *p_sim_compressor;
ParallelCompressors *p_compressor;
NvCtr *p_controller;
std::ofstream timing_file;

boost::timer::nanosecond_type time_initialize[2];
boost::timer::nanosecond_type time_solve[2];
boost::timer::nanosecond_type time_central;

void Callback(ParallelCompressors::State x, double t) {
  // Get and apply next input
  NvCtr::ControlInput u =
      p_controller
          ->GetNextInputWithTiming<NvCtr::TimerType::SPLIT_INITIAL_SOLVE_TIMES>(
              p_compressor->GetOutput(x), &time_central, time_initialize,
              time_solve);

  p_sim_compressor->SetInput(u);
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
  const std::string timing_fname = folder_name + "timing_ncoop" + std::to_string(n_solver_iterations) + ".dat";
  const std::string yref_fname = folder_name + "yref";
  const std::string ywt_fname = folder_name + "yweight_ncoop";
  const std::string uwt_fname = folder_name + "uweight_ncoop";

  timing_file.open(timing_fname, std::fstream::out | std::fstream::app);

  std::cout << "Running parallel non-cooperative simulation... ";
  std::cout.flush();

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
  const NvCtr::OutputPrediction y_ref =
      y_ref_sub.replicate<Controller1::p, 1>();

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
  u_disturbance(6) -= 0.1;

  sim_comp.SetOffset(u_disturbance);

  sim_comp.Integrate(50 + sampling_time, 500, sampling_time, &Callback);

  std::cout << "Finished." << std::endl;

  timing_file << time_central << std::endl;
  timing_file << time_initialize[0] << std::endl;
  timing_file << time_solve[0] << std::endl;
  timing_file << time_initialize[1] << std::endl;
  timing_file << time_solve[1] << std::endl;

  timing_file.close();

  return 0;
}
