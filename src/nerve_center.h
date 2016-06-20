#ifndef NERVE_CENTER_H
#define NERVE_CENTER_H

#include <tuple>

#include "distributed_controller.h"
#include "controller_interface.h"

template <typename System, int n_total_states, typename... SubControllers>
class NerveCenter : public ControllerInterface<System> {
 protected:
  static constexpr int n_states = System::n_states;
  static constexpr int n_inputs = System::n_inputs;
  static constexpr int n_outputs = System::n_outputs;
  static constexpr int n_control_inputs = System::n_control_inputs;
  static constexpr int n_controllers = sizeof...(SubControllers);

  using State = Eigen::Matrix<double, System::n_states, 1>;
  using AugmentedState = Eigen::Matrix<double, n_total_states, 1>;
  using Input = Eigen::Matrix<double, System::n_inputs, 1>;
  using Output = Eigen::Matrix<double, System::n_outputs, 1>;
  using ControlInput = Eigen::Matrix<double, System::n_control_inputs, 1>;

  using ControllerInterface<System>::u_offset_;

  // Sub controllers
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
      : ControllerInterface<System>(Input::Zero()),
        sys_(sys),
        sub_controllers_(controllers) {}

  /// Initialize all sub controllers based on given initial conditions
  void Initialize(const State& x_init, const ControlInput& u_init,
                  const Input& u_init_full, const Output& y_init,
                  const AugmentedState& dx_init = AugmentedState::Zero()) {
    using expander = int[];
    expander{(std::get<SubControllers>(sub_controllers_)
                  .template InitializeFull<n_states, n_outputs, n_total_states>(
                      x_init, u_init, u_init_full, y_init, dx_init))...};
    u_offset_ = u_init_full;
  }

  // TODO: Write this function
  virtual const ControlInput GetNextInput(const Output& y) {
    return ControlInput::Zero();
  }
};

#endif
