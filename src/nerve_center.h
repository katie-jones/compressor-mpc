#ifndef NERVE_CENTER_H
#define NERVE_CENTER_H

// #include <tuple>
// #include <utility>

#include "distributed_controller.h"

template <typename System, typename StateIndices, typename OutputIndices,
          typename... SubControllers>
class NerveCenter {
 protected:
  static constexpr int n_states = System::n_states;
  static constexpr int n_inputs = System::n_inputs;
  static constexpr int n_outputs = System::n_outputs;
  static constexpr int n_control_inputs = System::n_control_inputs;
  static constexpr int n_controllers = sizeof...(SubControllers);

  using State = Eigen::Matrix<double, System::n_states, 1>;
  using Input = Eigen::Matrix<double, System::n_inputs, 1>;
  using Output = Eigen::Matrix<double, System::n_outputs, 1>;
  using ControlInput = Eigen::Matrix<double, System::n_control_inputs, 1>;

  // Sub controllers
  // DistributedControllerBase* sub_controllers_[n_controllers];
  std::tuple<SubControllers...> sub_controllers_;

  // Current state estimate
  State x_;

  // Previous applied input
  ControlInput u_old_;

  // Previous output from plant
  Output y_old_;

  // Dynamic system to control
  System sys_;

 public:
  NerveCenter(const System& sys, std::tuple<SubControllers...>& controllers)
      : sys_(sys), sub_controllers_(controllers) {}

  void Initialize(const State& x_init, const ControlInput& u_init,
                  const Output& y_init) {}

 private:
  template <int... Ints>
  InitializeHelper(const State& x_init, const ControlInput& u_init,
                   const Output& y_init,
                   const std::integer_sequence<int, Ints...>);
};

#endif
