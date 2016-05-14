#ifndef AUGMENTED_SYSTEM_H
#define AUGMENTED_SYSTEM_H

#include <Eigen/Eigen>

#include "dynamic_system.h"

template <class System, int n_disturbance_states, int n_delay_states>
class AugmentedSystem {
 public:
  static constexpr int n_control_inputs = System::n_control_inputs;
  static constexpr int n_inputs = System::n_inputs;
  static constexpr int n_outputs = System::n_outputs;
  static constexpr int n_states = System::n_states;

  static constexpr int n_aug_states = n_disturbance_states + n_delay_states;
  static constexpr int n_obs_states = n_states + n_disturbance_states;
  static constexpr int n_total_states = n_aug_states + n_states;

  typedef typename System::Linearized LinearizedSystem;
  typedef typename System::State State;
  typedef typename System::Output Output;
  typedef typename System::Input Input;
  typedef Eigen::Matrix<double, n_control_inputs, 1> ControlInput;
  typedef Eigen::Matrix<double, n_total_states, 1> AugmentedState;
  using AugmentedLinearizedSystem =
      typename DynamicSystem<n_total_states, n_inputs, n_outputs,
                             n_control_inputs>::Linearized;

 public:
  typedef Eigen::Matrix<double, n_obs_states, n_outputs> ObserverMatrix;

 private:
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

  LinearizedSystem DiscretizeRK4(const LinearizedSystem &sys_continuous);
  AugmentedLinearizedSystem LinearizeAndAugment(
      const LinearizedSystem &sys_continuous);

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
                  const std::array<int, n_control_inputs> &n_delay,
                  const std::array<int, n_control_inputs> &control_input_index,
                  const Input &u_init = Input::Zero(),
                  const Input &u_offset = Input::Zero(),
                  const Output &y_init = Output::Zero(),
                  const AugmentedState &dx_init = AugmentedState::Zero())
      : AugmentedSystem(
            sys, Ts,
            (AugmentedState() << x_init,
             Eigen::Matrix<double, n_aug_states, 1>::Zero()).finished(),
            M, n_delay, control_input_index, u_init, u_offset, y_init,
            dx_init) {
    auglinsys_ =
        LinearizeAndAugment((sys_.GetLinearizedSystem(x_init, u_init)));
  }

  AugmentedSystem(const System &sys, const double Ts,
                  const AugmentedState &x_init, const ObserverMatrix &M,
                  const std::array<int, n_control_inputs> &n_delay,
                  const std::array<int, n_control_inputs> &control_input_index,
                  const Input &u_init = Input::Zero(),
                  const Input &u_offset = Input::Zero(),
                  const Output &y_init = Output::Zero(),
                  const AugmentedState &dx_init = AugmentedState::Zero())
      : sys_(sys),
        Ts_(Ts),
        x_aug_(x_init),
        dx_aug_(dx_init),
        u_offset_(u_offset),
        y_old_(y_init),
        M_(M),
        n_delay_(n_delay),
        control_input_index_(control_input_index) {
    u_old_ = GetControlInput(u_init);
    auglinsys_ = LinearizeAndAugment(
        sys_.GetLinearizedSystem(x_init.template head<n_states>(), u_init));
  }

  void ObserveAPriori(ControlInput &u_new);
  void ObserveAPosteriori(Output &y_new);
  AugmentedLinearizedSystem GetLinearization() const { return auglinsys_; }
  inline AugmentedState GetDeltaState() const { return dx_aug_; }
};

#endif
