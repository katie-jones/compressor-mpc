#define CONTROLLER_TYPE_CENTRALIZED
#define SYSTEM_TYPE_SERIAL

#include "common-variables.h"

constexpr int n_solver_iterations = 0;

SimSystem *p_sim_compressor;
SerialCompressors *p_compressor;
NvCtr *p_controller;
std::ofstream output_file;
std::ofstream cpu_times_file;

boost::timer::nanosecond_type time_initialize;
boost::timer::nanosecond_type time_solve;
boost::timer::nanosecond_type time_central;

boost::timer::cpu_timer timer;

void Callback(SerialCompressors::State x, double t) {
  output_file << t << std::endl;
  output_file << x.transpose() << std::endl;

  SerialCompressors::Output y = p_compressor->GetOutput(x);

  output_file << y.transpose() << std::endl;

  // Get and apply next input
  timer.resume();
  NvCtr::ControlInput u =
      p_controller
          ->GetNextInputWithTiming<NvCtr::TimerType::SPLIT_INITIAL_SOLVE_TIMES>(
              p_compressor->GetOutput(x), &time_central, &time_initialize,
              &time_solve);
  timer.stop();

  boost::timer::cpu_times elapsed = timer.elapsed();
  boost::timer::nanosecond_type elapsed_ns(elapsed.system + elapsed.user);

  cpu_times_file << elapsed_ns << std::endl;

  p_sim_compressor->SetInput(u);

  output_file << u.transpose() << std::endl << std::endl;
}

int main(void) {
  timer.stop();

  boost::timer::cpu_times time_offset = timer.elapsed();
  boost::timer::nanosecond_type offset_ns(time_offset.system +
                                          time_offset.user);
  const std::string folder_name = "serial/";

  const std::string constraints_fname = folder_name + "constraints";
  const std::string output_fname = folder_name + "output/cent_output0.dat";
  const std::string info_fname = folder_name + "output/cent_info0.dat";
  const std::string cpu_times_fname = folder_name + "output/cent_cpu_times0.dat";
  const std::string yref_fname = folder_name + "yref";
  const std::string ywt_fname = folder_name + "yweight_cent";
  const std::string uwt_fname = folder_name + "uweight_cent";
  const std::string disturbances_fname = folder_name + "disturbances_cent";

  std::cout << "Running serial centralized simulation... ";
  std::cout.flush();

  cpu_times_file.open(cpu_times_fname);

  cpu_times_file << offset_ns << std::endl;
  output_file.open(output_fname);

  // Time entire simulation
  boost::timer::cpu_timer simulation_timer;

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
                     compressor.n_outputs>::Identity())
          .finished();

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
  SerialCompressors::Input u_disturbance;
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
  output_file.close();
  cpu_times_file.close();

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

  info_file.open("centralized_serial_times.dat");

  info_file << "Central time: " << time_central << std::endl;
  info_file << "Initialize time: " << time_initialize << std::endl;
  info_file << "Solve time: " << time_solve << std::endl;

  info_file.close();

  return 0;
}
