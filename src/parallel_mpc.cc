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
constexpr int n_delay_states = 80;
constexpr int n_disturbance_states = 4;
constexpr int p = 100;
constexpr int m = 2;
}

extern template class MpcController<
    ParallelCompressors, Control::n_delay_states, Control::n_disturbance_states,
    Control::p, Control::m>;
using Controller =
    MpcController<ParallelCompressors, Control::n_delay_states,
                  Control::n_disturbance_states, Control::p, Control::m>;

SimulationParallelCompressors *p_sim_compressor;
ParallelCompressors *p_compressor;
Controller *p_controller;
std::ofstream output_file;

void Callback(ParallelCompressors::State x, double t) {
  // Update state

  output_file << t << std::endl;
  for (int i = 0; i < ParallelCompressors::n_states; i++)
    output_file << x[i] << "\t";
  output_file << std::endl;

  Controller::Output yref =
      p_controller->GetReference()
          .template head<ParallelCompressors::n_outputs>();
  Controller::Output y = p_compressor->GetOutput(x);

  for (int i = 0; i < ParallelCompressors::n_outputs; i++)
    output_file << yref(i) - y(i) << "\t";
  output_file << std::endl;

  // Get and apply next input
  // Controller::ControlInput u =
      // p_controller->GetNextInput(p_compressor->GetOutput(x));
  // p_sim_compressor->SetInput(u);

  // for (int i = 0; i < ParallelCompressors::n_control_inputs; i++)
    // output_file << u(i) << "\t";
  // output_file << std::endl
              // << std::endl;
}

int main(void) {
  output_file.open("output.txt");

  ParallelCompressors compressor;
  SimulationParallelCompressors sim_comp;
  p_compressor = &compressor;
  p_sim_compressor = &sim_comp;

  ParallelCompressors::Input u_default = ParallelCompressors::GetDefaultInput();
  ParallelCompressors::State x_init =
      (ParallelCompressors::State() << 0.915654, 1.14501, 0.151568, 439.989, 0,
       0.915654, 1.14501, 0.151568, 439.989, 0, 1.12016).finished();

  sim_comp.SetOffset(u_default);
  sim_comp.SetState(x_init);

  const double sampling_time = 0.05;

  const Controller::ObserverMatrix M =
      (Controller::ObserverMatrix()
           << Eigen::Matrix<double, compressor.n_states,
                            compressor.n_outputs>::Zero(),
       Eigen::Matrix<double, Control::n_disturbance_states,
                     compressor.n_outputs>::Identity()).finished();

  const Controller::ControlInputIndex delay = {0, Control::n_delay_states / 2,
                                               0, Control::n_delay_states / 2};

  // index of controlled states
  const Controller::ControlInputIndex index = {0, 3, 4, 7};

  const Controller::UWeightType uwt =
      (Controller::UWeightType() << 2e4, 0, 0, 0, 0, 2e5, 0, 0, 0, 0, 2e4, 0, 0,
       0, 0, 2e5).finished();
  // (Controller::YWeightType() << 0.1, 0, -0.1, 0, 0, 0, 1, 0, 0, 0, -0.1, 0,
  // 0.1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 5e2).finished();
  Controller::YWeightType ywt = Controller::YWeightType::Zero();
  ywt.diagonal() << 1, 1, 0.1, 5e2;
      // (Eigen::Diagonal<double, compressor.n_outputs>() << 1, 1, 0.1, 5e2)
          // .finished();

  const Controller::Input offset = u_default;

  const Controller::OutputPrediction y_ref =
      compressor.GetOutput(x_init).replicate<Control::p, 1>();

  // Input constraints
  Controller::InputConstraints constraints;
  constraints.lower_bound << -0.3, 0, -0.3, 0;
  constraints.upper_bound << 0.3, 1, 0.3, 1;
  constraints.lower_rate_bound << -0.1, -0.1, -0.1, -0.1;
  constraints.upper_rate_bound << 0.1, 1, 0.1, 1;
  constraints.use_rate_constraints = true;

  // Setup controller
  Controller ctrl(compressor, M, sampling_time, delay, index, uwt, ywt,
                  constraints, offset);
  p_controller = &ctrl;

  ctrl.SetReference(y_ref);
  ctrl.SetInitialState(x_init, Controller::ControlInput::Zero());

  // Integrate system and time it
  boost::timer::cpu_timer integrate_timer;
  sim_comp.Integrate(0, 200, sampling_time, &Callback);
  const boost::timer::cpu_times int_elapsed = integrate_timer.elapsed();
  const boost::timer::nanosecond_type elapsed_ns(int_elapsed.system +
                                                 int_elapsed.user);

  // Apply disturbance
  // ParallelCompressors::Input u_disturbance = u_default;
  // u_disturbance(2) -= 0.1;
  // sim_comp.SetOffset(u_disturbance);
  // sim_comp.Integrate(50 + sampling_time, 500, sampling_time, &Callback);

  output_file.close();
  std::cout << "CPU time: " << elapsed_ns << std::endl;

  return 0;
}
