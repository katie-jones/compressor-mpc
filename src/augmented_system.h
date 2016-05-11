#ifndef AUGMENTED_SYSTEM_H
#define AUGMENTED_SYSTEM_H

#include <Eigen/Eigen>

#include "dynamic_system.h"

template <class System, int n_disturbance_states, int n_delay_states>
class AugmentedSystem {
 protected:
  static constexpr int n_control_inputs = System::n_control_inputs;
  static constexpr int n_inputs = System::n_inputs;
  static constexpr int n_outputs = System::n_outputs;
  static constexpr int n_states = System::n_states;

 public:
  typedef Eigen::Matrix<double, n_states + n_disturbance_states, n_outputs>
      ObserverMatrix;

  using AugmentedLinearizedSystem =
      typename DynamicSystem<n_states + n_disturbance_states + n_delay_states,
                             n_inputs, n_outputs, n_control_inputs>::Linearized;

 protected:
  typedef typename System::Linearized LinearizedSystem;
  typedef typename System::State State;
  typedef typename System::Output Output;
  typedef typename System::Input Input;
  typedef Eigen::Matrix<double, n_control_inputs, 1> ControlInput;
  typedef Eigen::Matrix<double,
                        n_states + n_disturbance_states + n_delay_states,
                        1> AugmentedState;

  System sys_;
  const double Ts_;
  AugmentedState x_aug_;
  AugmentedState dx_aug_;
  ControlInput u_old_;
  Input u_offset_;
  Output y_old_;
  const ObserverMatrix M_;
  AugmentedLinearizedSystem auglinsys_;
  std::array<int, n_control_inputs> n_delay_;
  std::array<int, n_control_inputs> control_input_index_;

  LinearizedSystem DiscretizeRK4(LinearizedSystem &sys_continuous);
  AugmentedLinearizedSystem LinearizeAndAugment(
      LinearizedSystem &sys_continuous);

  ControlInput GetControlInput(Input u, Input offset = Input::Zero()) {
    ControlInput u_control;
    for (int i = 0; i < n_control_inputs; i++) {
      u_control(i) =
          u(control_input_index_[i]) - offset(control_input_index_[i]);
    }
    return u_control;
  }

  Input GetPlantInput(const ControlInput &u_control,
                      const Input &offset = Input::Zero()) {
    Input u = offset;
    for (int i = 0; i < n_control_inputs; i++) {
      u(control_input_index_[i]) += u_control(i);
    }
    return u;
  }

 public:
  AugmentedSystem(const System &sys, const double Ts, const State &x_init,
                  const ObserverMatrix &M,
                  std::array<int, n_control_inputs> n_delay,
                  std::array<int, n_control_inputs> control_input_index,
                  const Input &u_init = Input(),
                  const Output &y_init = Output(),
                  const AugmentedState &dx_init = AugmentedState())
      : sys_(sys),
        Ts_(Ts),
        x_aug_((AugmentedState() << x_init,
                Eigen::MatrixXd::Zero(n_delay_states + n_disturbance_states, 1))
                   .finished()),
        dx_aug_(dx_init),
        y_old_(y_init),
        M_(M),
        n_delay_(n_delay),
        control_input_index_(control_input_index) {
    u_old_ = GetControlInput(u_init);
  }

  AugmentedSystem(const System &sys, const double Ts,
                  const AugmentedState &x_init, const ObserverMatrix &M,
                  std::array<int, n_control_inputs> n_delay,
                  std::array<int, n_control_inputs> control_input_index,
                  const Input &u_init = Input(),
                  const Output &y_init = Output(),
                  const AugmentedState &dx_init = AugmentedState())
      : sys_(sys),
        Ts_(Ts),
        x_aug_(x_init),
        dx_aug_(dx_init),
        y_old_(y_init),
        M_(M),
        n_delay_(n_delay),
        control_input_index_(control_input_index) {
    u_old_ = GetControlInput(u_init);
  }

  void ObserveAPriori(ControlInput &u_new);
  void ObserveAPosteriori(Output &y_new);
};

#endif
