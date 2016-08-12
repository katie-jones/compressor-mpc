#define CONTROLLER_TYPE_CENTRALIZED
#define SYSTEM_TYPE_PARALLEL
#define EIGEN_DONT_PARALLELIZE 1

#include "common-variables.h"

int n_solver_iterations;

SimSystem *p_sim_compressor;
ParallelCompressors *p_compressor;
NvCtr *p_controller;
std::ofstream cpu_times_file;

boost::timer::nanosecond_type time_total;

void Callback(ParallelCompressors::State x, double t) {
  ParallelCompressors::Output y = p_compressor->GetOutput(x);

  // Get and apply next input
  NvCtr::ControlInput u =
      p_controller->GetNextInputWithTiming<NvCtr::TOTAL_TIME>(
          p_compressor->GetOutput(x), &time_total);

  p_sim_compressor->SetInput(u);
}

int main(int argc, char **argv) {
  time_total = 0;
  if (argc < 2) {
    n_solver_iterations = 1;
  } else {
    std::istringstream ss(argv[1]);
    if (!(ss >> n_solver_iterations)) {
      std::cerr << "Invalid number " << argv[1] << std::endl;
      return 1;
    }
    if (n_solver_iterations < 0 || n_solver_iterations > 1) {
      std::cout << "Warning: Number of solver iterations for centralized "
                   "controller should be 0 or 1. Using 1 solver iteration."
                << std::endl;
      n_solver_iterations = 1;
    }
  }

  const std::string folder_name = "parallel/";

  const std::string constraints_fname = folder_name + "constraints";
  const std::string cpu_times_fname = folder_name + "output/cent_cpu_times.dat";
  const std::string yref_fname = folder_name + "yref";
  const std::string ywt_fname = folder_name + "yweight_cent";
  const std::string uwt_fname = folder_name + "uweight_cent";
  const std::string disturbances_fname = folder_name + "disturbances_cent";

  std::cout << "Running centralized simulation... ";
  std::cout.flush();

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
  NvCtr nerve_center(ctrl_tuple, n_solver_iterations);
  p_controller = &nerve_center;

  // Test functions
  nerve_center.SetWeights(uwt, ywt);
  nerve_center.SetOutputReference(y_ref);

  nerve_center.Initialize(compressor.GetDefaultState(),
                          AugmentedSystem::ControlInput::Zero(),
                          compressor.GetDefaultInput(),
                          compressor.GetOutput(compressor.GetDefaultState()));

  // Initialize disturbance
  ParallelCompressors::Input u_disturbance;
  double t_past = -sampling_time;
  double t_next;

  std::ifstream disturbances_file;
  disturbances_file.open(disturbances_fname);

  // Read inputs and times from file
  while (ReadDataFromStream(u_disturbance.data(), disturbances_file,
                            u_disturbance.size())) {
    if (!(ReadDataFromStream(&t_next, disturbances_file, 1))) {
      std::cerr << "Simulation time could not be read." << std::endl;
      break;
    }
    std::cout << "Simulating from time " << t_past << " to time " << t_next
              << " with offset:" << std::endl;
    for (int i = 0; i < u_disturbance.size(); i++) {
      std::cout << u_disturbance(i) << "\t";
    }
    std::cout << std::endl;

    sim_comp.SetOffset(u_default + u_disturbance);

    sim_comp.Integrate(t_past + sampling_time, t_next, sampling_time,
                       &Callback);

    t_past = t_next;
  }

  disturbances_file.close();
  boost::timer::cpu_times simulation_cpu_time = simulation_timer.elapsed();
  boost::timer::nanosecond_type simulation_ns(simulation_cpu_time.system +
                                              simulation_cpu_time.user);

  cpu_times_file.open(cpu_times_fname);
  cpu_times_file << time_total << std::endl;
  cpu_times_file.close();

  std::cout << "Finished." << std::endl
            << "Total time required:\t"
            << static_cast<double>(simulation_ns) / 1.0e6 << " ms." << std::endl
            << std::endl;

  return 0;
}
