// Test of serial centralized controller
#define CONTROLLER_TYPE_CENTRALIZED
#define SYSTEM_TYPE_SERIAL

#include "common-variables.h"

constexpr int n_solver_iterations = 1;

SimSystem *p_sim_compressor;
SerialCompressors *p_compressor;
NvCtr *p_controller;
std::ofstream timing_file;

boost::timer::nanosecond_type time_initialize;
boost::timer::nanosecond_type time_solve;
boost::timer::nanosecond_type time_central;

boost::timer::cpu_timer timer;

void Callback(SerialCompressors::State x, double t) {
  // Get and apply next input
  NvCtr::ControlInput u =
      p_controller
          ->GetNextInputWithTiming<NvCtr::TimerType::SPLIT_INITIAL_SOLVE_TIMES>(
              p_compressor->GetOutput(x), &time_central, &time_initialize,
              &time_solve);

  p_sim_compressor->SetInput(u);
}

int main(void) {
  timer.stop();

  const std::string folder_name = "serial/";

  const std::string constraints_fname = folder_name + "constraints";
  const std::string timing_fname =
      folder_name + "timing_centralized.dat";
  const std::string cpu_times_fname = folder_name + "output/cent_cpu_times.dat";
  const std::string yref_fname = folder_name + "yref";
  const std::string ywt_fname = folder_name + "yweight_cent";
  const std::string uwt_fname = folder_name + "uweight_cent";

  timing_file.open(timing_fname, std::fstream::out | std::fstream::app);

  std::cout << "Running serial centralized simulation... ";
  std::cout.flush();

  SerialCompressors compressor;
  p_compressor = &compressor;

  SerialCompressors::Input u_default = SerialCompressors::GetDefaultInput();
  SerialCompressors::State x_init = SerialCompressors::GetDefaultState();

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

  if (!(ReadDataFromFile(uwt.data(), uwt.rows(), uwt_fname, uwt.rows() + 1))) {
    return -1;
  }
  if (!(ReadDataFromFile(ywt.data(), ywt.rows(), ywt_fname, ywt.rows() + 1))) {
    return -1;
  }

  // Read in reference output
  SerialCompressors::Output y_ref_sub;
  if (!(ReadDataFromFile(y_ref_sub.data(), y_ref_sub.size(), yref_fname))) {
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
  SerialCompressors::Input u_disturbance = u_default;
  u_disturbance(6) -= 0.1;

  sim_comp.SetOffset(u_disturbance);

  sim_comp.Integrate(50 + sampling_time, 500, sampling_time, &Callback);

  std::cout << "Finished." << std::endl;

  timing_file << time_central << std::endl;
  timing_file << time_initialize << std::endl;
  timing_file << time_solve << std::endl;

  timing_file.close();

  return 0;
}
