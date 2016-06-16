#include <boost/timer/timer.hpp>
#include <fstream>
#include <iostream>
#include "aug_lin_sys.h"
#include "input_constraints.h"
#include "noncooperative_controller.h"
#include "observer.h"
#include "parallel_compressors.h"
#include "print_matrix.h"
#include "simulation_system.h"
#include "constexpr_array.h"
#include "distributed_controller.h"
#include "nerve_center.h"

namespace Control {
constexpr int n_delay_states = 80;
constexpr int n_disturbance_states = 4;
constexpr int p = 100;
constexpr int m = 2;
constexpr int n_controllers = 2;
constexpr int n_sub_outputs = 3;
constexpr int n_sub_control_inputs = 2;

using Delays = ConstexprArray<0, 40, 0, 40>;
using OutputIndices = ConstexprArray<0, 1, 3>;
using StateIndices = ConstexprArray<0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10>;
using ControlInputIndices1 = ConstexprArray<0, 1, 2, 3>;
using ControlInputIndices2 = ConstexprArray<2, 3, 0, 1>;
using InputIndices = ConstexprArray<0, 3, 4, 7>;
using ControlledOutputIndices = ConstexprArray<0, 1, 3>;

using OutputIndexList = ConstexprArrayList<OutputIndices, OutputIndices>;
using StateIndexList = ConstexprArrayList<StateIndices, StateIndices>;
}

using namespace Control;

using AugmentedSystem1 =
    AugmentedLinearizedSystem<ParallelCompressors, Delays, n_disturbance_states,
                              ControlInputIndices1, n_sub_control_inputs>;
using AugmentedSystem2 =
    AugmentedLinearizedSystem<ParallelCompressors, Delays, n_disturbance_states,
                              ControlInputIndices2, n_sub_control_inputs>;

using Obsv1 = Observer<AugmentedSystem1>;
using Obsv2 = Observer<AugmentedSystem2>;

using Controller1 =
    DistributedController<AugmentedSystem1, ControlledOutputIndices, p, m>;
using Controller2 =
    DistributedController<AugmentedSystem2, ControlledOutputIndices, p, m>;

using SimSystem = SimulationSystem<ParallelCompressors, Delays, InputIndices>;

using NvCtr = NerveCenter<ParallelCompressors, StateIndexList, OutputIndexList,
                          Controller1, Controller2>;

extern template class AugmentedLinearizedSystem<
    ParallelCompressors, Delays, n_disturbance_states, ControlInputIndices1,
    n_sub_control_inputs>;
extern template class AugmentedLinearizedSystem<
    ParallelCompressors, Delays, n_disturbance_states, ControlInputIndices2,
    n_sub_control_inputs>;

extern template class Observer<AugmentedSystem1>;
extern template class Observer<AugmentedSystem2>;

extern template class DistributedController<AugmentedSystem1,
                                            ControlledOutputIndices, p, m>;
extern template class DistributedController<AugmentedSystem2,
                                            ControlledOutputIndices, p, m>;

SimSystem *p_sim_compressor;
ParallelCompressors *p_compressor;
// Controller *p_controller;
std::ofstream output_file;
std::ofstream cpu_times_file;

void Callback(ParallelCompressors::State x, double t) {
  output_file << t << std::endl;
  output_file << x.transpose() << std::endl;

  ParallelCompressors::Output y = p_compressor->GetOutput(x);

  output_file << y.transpose() << std::endl;

  // Get and apply next input
  // Controller::ControlInput u =
  // p_controller->GetNextInput(p_compressor->GetOutput(x), cpu_times_file);
  // p_sim_compressor->SetInput(u);

  // output_file << u.transpose() << std::endl
  // << std::endl;
}

int main(void) {
  output_file.open("coop_output.dat");
  cpu_times_file.open("coop_cpu_times.dat");

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
       Eigen::Matrix<double, Control::n_disturbance_states,
                     compressor.n_outputs>::Identity()).finished();

  // Controller::UWeightType uwt = Controller::UWeightType::Zero();
  // Controller::YWeightType ywt = Controller::YWeightType::Zero();

  // std::ifstream weight_file;
  // weight_file.open("uweight_coop");
  // for (int i = 0; i < uwt.rows(); i++) {
  // weight_file >> uwt(i, i);
  // }
  // weight_file.close();

  // weight_file.open("yweight_coop");
  // for (int i = 0; i < ywt.rows(); i++) {
  // weight_file >> ywt(i, i);
  // }
  // weight_file.close();

  const AugmentedSystem1::Input u_offset = u_default;

  const Controller1::OutputPrediction y_ref =
      (Controller1::ControlOutput() << 4.5, 4.5, 1.12)
          .finished()
          .replicate<p, 1>();

  // Input constraints
  InputConstraints<n_sub_control_inputs> constraints;
  constraints.lower_bound << -0.3, 0;
  constraints.upper_bound << 0.3, 1;
  constraints.lower_rate_bound << -0.1, -0.1;
  constraints.upper_rate_bound << 0.1, 1;
  constraints.use_rate_constraints = true;

  // Setup controller
  AugmentedSystem1 sys1(compressor, sampling_time);
  AugmentedSystem2 sys2(compressor, sampling_time);
  Controller1 ctrl1(sys1, constraints, M);
  Controller2 ctrl2(sys2, constraints, M);

  // Test functions
  // ctrl1.Initialize(compressor.GetOutput(x_init),
  // AugmentedSystem1::SubControlInput::Zero(), u_offset, x_init);
  // ctrl2.Initialize(compressor.GetOutput(x_init),
  // AugmentedSystem2::SubControlInput::Zero(), u_offset, x_init);
  // ctrl1.SetOutputReference(y_ref);
  // ctrl2.SetOutputReference(y_ref);

  // Generate a QP
  // MpcQpSolver<
  // n_disturbance_states + n_delay_states + ParallelCompressors::n_states,
  // ControlledOutputIndices::size, n_sub_control_inputs, p, m>::QP qp =
  // ctrl1.GenerateInitialQP(compressor.GetOutput(x_init), u_offset);

  // Test QP solver
  // Controller1::ControlInputPrediction u_sol;
  // for (int i = 0; i < 20; i++) {
  // ctrl1.GetInput(&u_sol, &qp,
  // Eigen::Matrix<double, n_sub_control_inputs * m, 1>::Zero());

  // std::cout << u_sol << std::endl;
  // }

  // Create a nerve center
  std::tuple<Controller1, Controller2> ctrl_tuple(ctrl1, ctrl2);
  NvCtr nerve_center(compressor, ctrl_tuple);

  output_file.close();
  cpu_times_file.close();
  return 0;
}
