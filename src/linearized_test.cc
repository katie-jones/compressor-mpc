#include <iostream>
#include <fstream>
#include <boost/timer/timer.hpp>
#include "compressor.h"
#include "simulation_compressor.h"
#include "tank.h"
#include "parallel_compressors.h"
#include "simulation_parallel_compressors.h"
#include "mpc_controller.h"
#include "print_matrix.h"

namespace Control {
constexpr int n_delay_states = 5;
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
  // Update state
  // p_sim_compressor->SetState(x);

  output_file << t << std::endl;
  for (int i = 0; i < Compressor::n_states; i++) output_file << x[i] << "\t";
  output_file << std::endl;

  Controller::Output yref =
      p_controller->GetReference().template head<Compressor::n_outputs>();
  Controller::Output y = p_compressor->GetOutput(x);

  for (int i = 0; i < Compressor::n_outputs; i++)
    output_file << yref(i) - y(i) << "\t";
  output_file << std::endl;

  // Get and apply next input
  Controller::ControlInput u =
      p_controller->GetNextInput(p_compressor->GetOutput(x));
  p_sim_compressor->SetInput(u);

  for (int i = 0; i < Compressor::n_control_inputs; i++)
    output_file << u(i) << "\t";
  output_file << std::endl
              << std::endl;
}

int main(void) {
  output_file.open("output.txt");

  Compressor compressor;
  SimulationCompressor sim_comp;
  p_compressor = &compressor;
  p_sim_compressor = &sim_comp;

  Compressor::Input u_default = Compressor::GetDefaultInput();
  Compressor::State x_init = (Compressor::State() << 0.898978, 1.12501,
                              0.151226, 440.679, 0).finished();

  sim_comp.SetOffset(u_default);
  sim_comp.SetState(x_init);

  const double sampling_time = 0.05;

  const Controller::ObserverMatrix M =
      (Controller::ObserverMatrix()
           << Eigen::Matrix<double, compressor.n_states,
                            compressor.n_control_inputs>::Zero(),
       Eigen::Matrix<double, Control::n_disturbance_states,
                     Control::n_disturbance_states>::Identity()).finished();

  const Controller::ControlInputIndex delay = {0, Control::n_delay_states};

  // index of controlled states
  const Controller::ControlInputIndex index = {0, 3};

  const Controller::UWeightType uwt =
      (Controller::UWeightType() << 100, 0, 0, 1e4).finished();
  const Controller::YWeightType ywt =
      (Controller::YWeightType() << 0.1, 0, 0, 1).finished();

  const Controller::Input offset = u_default;

  const Controller::OutputPrediction y_ref =
      compressor.GetOutput(x_init).replicate<Control::p, 1>();

  // Input constraints
  Controller::InputConstraints constraints;
  constraints.lower_bound << -0.3, 0;
  constraints.upper_bound << 0.3, 1;
  constraints.lower_rate_bound << -0.1, -0.1;
  constraints.upper_rate_bound << 0.1, 1;
  constraints.use_rate_constraints = true;

  // Setup controller
  Controller ctrl(compressor, M, sampling_time, delay, index, uwt, ywt,
                  constraints, offset);
  p_controller = &ctrl;

  ctrl.SetReference(y_ref);
  ctrl.SetInitialState(x_init, Controller::ControlInput::Zero());

  // Integrate system and time it
  boost::timer::cpu_timer integrate_timer;
  sim_comp.Integrate(0, 50, sampling_time, &Callback);
  const boost::timer::cpu_times int_elapsed = integrate_timer.elapsed();
  const boost::timer::nanosecond_type elapsed_ns(int_elapsed.system +
                                                 int_elapsed.user);
  std::cout << "CPU time: " << elapsed_ns << std::endl;

  // Apply disturbance
  Compressor::Input u_disturbance = u_default;
  u_disturbance(2) -= 0.1;
  sim_comp.SetOffset(u_disturbance);
  sim_comp.Integrate(50 + sampling_time, 500, sampling_time, &Callback);

  output_file.close();

  return 0;
}
