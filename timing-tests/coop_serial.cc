#define CONTROLLER_TYPE_COOP
#define SYSTEM_TYPE_SERIAL

#include "common-variables.h"

SimSystem *p_sim_compressor;
SerialCompressors *p_compressor;
NvCtr *p_controller;
std::ofstream cpu_times_file;
boost::timer::nanosecond_type time_out;
int n_solver_iterations;
constexpr int n_max_solver_iterations = 9;

void Callback(SerialCompressors::State x, double t) {
  SerialCompressors::Output y = p_compressor->GetOutput(x);

  // Get and apply next input
  NvCtr::ControlInput u = p_controller->GetNextInputWithTiming(
      p_compressor->GetOutput(x), n_solver_iterations, &time_out);
  cpu_times_file << time_out << std::endl;

  p_sim_compressor->SetInput(u);
}

int main(int argc, char **argv) {
  if (argc < 2) {
    std::cout << "Number of solver iterations for timing: ";
    std::cin >> n_solver_iterations;
  } else {
    std::istringstream ss(argv[1]);
    if (!(ss >> n_solver_iterations))
      std::cerr << "Invalid number " << argv[1] << '\n';
  }

  const std::string folder_name = "serial/";

  const std::string constraints_fname = folder_name + "dist_constraints";
  const std::string cpu_times_fname = folder_name + "output/coop_cpu_times" +
                                      std::to_string(n_solver_iterations) +
                                      ".dat";
  const std::string yref_fname = folder_name + "yref";
  const std::string ywt_fname = folder_name + "yweight_coop";
  const std::string uwt_fname = folder_name + "uweight_coop";
  const std::string disturbances_fname = folder_name + "disturbances_coop";
  const std::string xinit_fname = folder_name + "xinit_coop";
  const std::string uinit_fname = folder_name + "uinit_coop";

  cpu_times_file.open(cpu_times_fname);

  std::cout << "Running serial cooperative simulation using "
            << n_max_solver_iterations << " solver iterations (timing "
            << n_solver_iterations << " iterations)... " << std::endl;

  // Time entire simulation
  boost::timer::cpu_timer simulation_timer;

  SerialCompressors compressor;
  p_compressor = &compressor;

  SerialCompressors::Input u_default;
  SerialCompressors::State x_init;

  // Read initial state
  if (!(ReadDataFromFile(u_default.data(), u_default.size(), uinit_fname))) {
    std::cerr << "Initial inputs could not be read" << std::endl;
    return 1;
  }

  u_default += compressor.GetDefaultInput();

  if (!(ReadDataFromFile(x_init.data(), x_init.size(), xinit_fname))) {
    std::cerr << "Initial state could not be read" << std::endl;
    return 1;
  }

  SimSystem sim_comp(p_compressor, u_default, x_init);
  p_sim_compressor = &sim_comp;

  const double sampling_time = 0.05;

  const Obsv1::ObserverMatrix M =
      (Obsv1::ObserverMatrix() << Eigen::Matrix<double, compressor.n_states,
                                                compressor.n_outputs>::Zero(),
       Eigen::Matrix<double, n_disturbance_states,
                     compressor.n_outputs>::Identity())
          .finished();

  // Weights
  NvCtr::UWeightType uwt = NvCtr::UWeightType::Zero();
  NvCtr::YWeightType ywt = NvCtr::YWeightType::Zero();

  std::ifstream yweight_file;
  yweight_file.open(ywt_fname);

  if (!(ReadDataFromFile(uwt.data(), uwt.rows(), uwt_fname, uwt.rows() + 1))) {
    std::cerr << "U weights could not be read" << std::endl;
    return -1;
  }
  if (!(ReadDataFromStream(ywt.data(), yweight_file, ywt.rows(),
                           ywt.rows() + 1))) {
    std::cerr << "Y1 weights could not be read" << std::endl;
    return -1;
  }

  yweight_file.close();

  // Read in reference output
  SerialCompressors::Output y_ref_sub;
  if (!(ReadDataFromFile(y_ref_sub.data(), y_ref_sub.size(), yref_fname))) {
    std::cerr << "Reference signal could not be read" << std::endl;
    return -1;
  }
  const NvCtr::OutputPrediction y_ref =
      y_ref_sub.replicate<Controller1::p, 1>();

  // Input constraints
  InputConstraints<n_sub_control_inputs> constraints;
  constraints.use_rate_constraints = true;
  if (!(ReadConstraintsFromFile(&constraints, constraints_fname))) {
    std::cerr << "Input constraints could not be read" << std::endl;
    return -1;
  }

  // Setup controller
  AugmentedSystem1 sys1(compressor, sampling_time);
  AugmentedSystem2 sys2(compressor, sampling_time);
  Controller1 ctrl1(sys1, constraints, M);
  Controller2 ctrl2(sys2, constraints, M);

  // Create a nerve center
  std::tuple<Controller1, Controller2> ctrl_tuple(ctrl1, ctrl2);
  NvCtr nerve_center(ctrl_tuple, n_max_solver_iterations);
  p_controller = &nerve_center;

  // Set up controllers
  nerve_center.SetWeights(uwt, ywt);
  nerve_center.SetOutputReference(y_ref);

  nerve_center.Initialize(x_init, AugmentedSystem1::ControlInput::Zero(),
                          u_default, compressor.GetOutput(x_init));

  // Initialize disturbance
  SerialCompressors::Input u_disturbance = SerialCompressors::Input::Zero();
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
