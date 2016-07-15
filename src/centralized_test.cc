#define CONTROLLER_TYPE_CENTRALIZED
#define SYSTEM_TYPE_PARALLEL

#include "common-variables.h"

constexpr int n_solver_iterations = 1;

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
  const std::string folder_name = "parallel/";

  const std::string constraints_fname = folder_name + "constraints";
  const std::string output_fname = folder_name + "output/cent_output.dat";
  const std::string info_fname = folder_name + "output/cent_info.dat";
  const std::string yref_fname = folder_name + "yref";
  const std::string ywt_fname = folder_name + "yweight_cent";
  const std::string uwt_fname = folder_name + "uweight_cent";

  std::cout << "Running centralized simulation... ";
  std::cout.flush();

  output_file.open(output_fname);

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
  info_file << uwt << std::endl
            << ywt << std::endl
            << y_ref;
  info_file.close();

  return 0;
}
