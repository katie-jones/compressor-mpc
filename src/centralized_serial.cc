// Test of serial centralized controller

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
#include "serial_compressors.h"
#include "serial_compressors_constants.h"
#include "simulation_system.h"

constexpr int n_solver_iterations = 1;

using namespace SERIAL_COMPRESSORS_CONSTANTS;

using SimSystem = SimulationSystem<SerialCompressors, Delays, InputIndices>;

using AugmentedSystem = SERIAL_AUGSYS_CENT;

using Obsv = SERIAL_OBS_CENT;
using Controller = SERIAL_CTRL_CENT;

using NvCtr = NerveCenter<SerialCompressors, n_total_states, Controller>;

SimSystem *p_sim_compressor;
SerialCompressors *p_compressor;
NvCtr *p_controller;
std::ofstream output_file;
std::ofstream cpu_times_file;

boost::timer::cpu_timer timer;

void Callback(SerialCompressors::State x, double t) {
  output_file << t << std::endl;
  output_file << x.transpose() << std::endl;

  SerialCompressors::Output y = p_compressor->GetOutput(x);

  output_file << y.transpose() << std::endl;

  // Get and apply next input
  timer.resume();
  NvCtr::ControlInput u =
      p_controller->GetNextInput(p_compressor->GetOutput(x));
  timer.stop();

  boost::timer::cpu_times elapsed = timer.elapsed();
  boost::timer::nanosecond_type elapsed_ns(elapsed.system + elapsed.user);

  cpu_times_file << elapsed_ns << std::endl;

  p_sim_compressor->SetInput(u);

  output_file << u.transpose() << std::endl
              << std::endl;
}

int main(void) {
  timer.stop();

  boost::timer::cpu_times time_offset = timer.elapsed();
  boost::timer::nanosecond_type offset_ns(time_offset.system +
                                          time_offset.user);

  std::cout << "Running serial centralized simulation... ";
  std::cout.flush();

  cpu_times_file.open("serial/output/cent_cpu_times.dat");

  cpu_times_file << offset_ns << std::endl;
  output_file.open("serial/output/cent_output.dat");

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
                     compressor.n_outputs>::Identity()).finished();

  NvCtr::UWeightType uwt = NvCtr::UWeightType::Zero();
  NvCtr::YWeightType ywt = NvCtr::YWeightType::Zero();

  std::ifstream read_file;
  read_file.open("serial/uweight_cent");
  for (int i = 0; i < uwt.rows(); i++) {
    if (!(read_file >> uwt(i, i))) {
      std::cerr << "Error reading input weight from file serial/uweight_cent"
                << std::endl;
      return -1;
    }
  }
  read_file.close();

  read_file.open("serial/yweight_cent");
  for (int i = 0; i < ywt.rows(); i++) {
    if (!(read_file >> ywt(i, i))) {
      std::cerr << "Error reading output weight from file serial/yweight_cent"
                << std::endl;
      return -1;
    }
  }
  read_file.close();

  const AugmentedSystem::Input u_offset = u_default;

  // Read in reference output
  SerialCompressors::Output y_ref_sub;
  read_file.open("serial/yref");
  for (int i = 0; i < y_ref_sub.size(); i++) {
    if (!(read_file >> y_ref_sub(i))) {
      std::cerr << "Error reading reference output from file serial/yref_cent"
                << std::endl;
      return -1;
    }
  }
  read_file.close();

  const NvCtr::OutputPrediction y_ref = y_ref_sub.replicate<Controller::p, 1>();

  // Input constraints
  InputConstraints<SerialCompressors::n_control_inputs> constraints;
  constraints.use_rate_constraints = true;
  read_file.open("serial/constraints");

  for (int i=0; i<SerialCompressors::n_control_inputs; i++) {
    if (!(read_file >> constraints.lower_bound(i))) {
      std::cerr << "Error reading input constraints from file serial/constraints" << std::endl;
      return -1;
    }
  }
  for (int i=0; i<SerialCompressors::n_control_inputs; i++) {
    if (!(read_file >> constraints.upper_bound(i))) {
      std::cerr << "Error reading input constraints from file serial/constraints" << std::endl;
      return -1;
    }
  }
  for (int i=0; i<SerialCompressors::n_control_inputs; i++) {
    if (!(read_file >> constraints.lower_rate_bound(i))) {
      std::cerr << "Error reading input constraints from file serial/constraints" << std::endl;
      return -1;
    }
  }
  for (int i=0; i<SerialCompressors::n_control_inputs; i++) {
    if (!(read_file >> constraints.upper_rate_bound(i))) {
      std::cerr << "Error reading input constraints from file serial/constraints" << std::endl;
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
  SerialCompressors::Input u_disturbance = u_default;
  u_disturbance(6) -= 0.1;

  sim_comp.SetOffset(u_disturbance);

  sim_comp.Integrate(50 + sampling_time, 500, sampling_time, &Callback);

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
  info_file.open("serial/output/cent_info.dat");
  info_file << uwt << std::endl
            << ywt << std::endl
            << y_ref;
  info_file.close();

  return 0;
}
