#include <iostream>
#include <fstream>
#include "compressor.h"
#include "simulation_compressor.h"
#include "tank.h"
#include "parallel_compressors.h"
#include "simulation_parallel_compressors.h"
#include "mpc_controller.h"
#include "print_matrix.h"

namespace Control {
constexpr int n_delay_states = 40;
constexpr int n_disturbance_states = 2;
constexpr int p = 100;
constexpr int m = 2;
}

extern template class MpcController<Compressor, Control::n_delay_states,
                                    Control::n_disturbance_states, Control::p,
                                    Control::m>;
using Controller =
    MpcController<Compressor, Control::n_delay_states,
                  Control::n_disturbance_states, Control::p, Control::m>;

SimulationCompressor *p_sim_compressor;
Compressor *p_compressor;
Controller *p_controller;
std::ofstream output_file;

void Callback(Compressor::State x, double t) {
  output_file << "t: " << t << "\tx: ";
  for (int i = 0; i < 5; i++) output_file << x[i] << "\t";
  output_file << std::endl;

  Controller::Output yref =
      p_controller->GetReference().template head<Compressor::n_outputs>();

  Controller::Output y = p_compressor->GetOutput(x);
  
  // Get and apply next input
  Controller::Input u = p_controller->GetNextInput(p_compressor->GetOutput(x));
  p_sim_compressor->SetInput(u);
  print_matrix(std::cout, u, "u");

  Controller::ControlInput uctrl = p_controller->GetControlInput(u);
}

int main(void) {
  output_file.open("output.txt");
  Compressor compressor;
  SimulationCompressor sim_comp;
  p_compressor = &compressor;
  p_sim_compressor = &sim_comp;

  Compressor::Input u_default = Compressor::GetDefaultInput();

  const double sampling_time = 0.05;

  const Controller::ObserverMatrix M =
      (Controller::ObserverMatrix()
           << Eigen::Matrix<double, compressor.n_states,
                            compressor.n_control_inputs>::Zero(),
       Eigen::Matrix<double, Control::n_disturbance_states,
                     Control::n_disturbance_states>::Identity()).finished();

  const Controller::ControlInputIndex delay = {0, Control::n_delay_states};
  const Controller::ControlInputIndex index = {0, 3};

  const Controller::UWeightType uwt =
      (Controller::UWeightType() << 100, 0, 0, 1e4).finished();
  const Controller::YWeightType ywt =
      (Controller::YWeightType() << 0.1, 0, 0, 1).finished();
  const Controller::Input offset = Compressor::GetDefaultInput();
  const Controller::OutputPrediction y_ref =
      compressor.GetOutput(Compressor::GetDefaultState())
          .replicate<Control::p, 1>();
  print_matrix(std::cout, y_ref, "yref");

  Controller::InputConstraints constraints;
  constraints.lower_bound << -0.3, 0;
  constraints.upper_bound << 0.3, 1;
  // constraints.lower_bound << 0, 0;
  // constraints.upper_bound << 0, 0;
  constraints.lower_rate_bound << -0.1, -0.1;
  constraints.upper_rate_bound << 0.1, 1;
  constraints.use_rate_constraints = true;

  Controller ctrl(compressor, M, 0.05, delay, index, uwt, ywt, constraints,
                  offset);
  p_controller = &ctrl;

  ctrl.SetReference(y_ref);
  ctrl.SetInitialState(compressor.GetDefaultState(),
                       Controller::ControlInput::Zero());

  sim_comp.Integrate(0, 2, sampling_time, &Callback);
  // Compressor::Input u_disturbance = u_default;
  // u_disturbance(1) -= 0.1;
  // sim_comp.SetInput(u_disturbance);
  // sim_comp.Integrate(10+sampling_time, 100, sampling_time, &Callback);

  output_file.close();

  return 0;
}
