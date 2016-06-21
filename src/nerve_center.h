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

  using expander = int[];

 public:
  using State = Eigen::Matrix<double, System::n_states, 1>;
  using AugmentedState = Eigen::Matrix<double, n_total_states, 1>;
  using Input = Eigen::Matrix<double, System::n_inputs, 1>;
  using Output = Eigen::Matrix<double, System::n_outputs, 1>;
  using ControlInput = Eigen::Matrix<double, System::n_control_inputs, 1>;

 protected:
  using ControllerInterface<System>::u_offset_;

  static constexpr int GetMaxP() {
    int p = 0;
    int dummy[]{(p = (p > SubControllers::p) ? p : SubControllers::p)...};
    return p;
  };

  static constexpr int p_max = GetMaxP();

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
  using UWeightType = Eigen::Matrix<double, n_control_inputs, n_control_inputs>;
  using YWeightType = Eigen::Matrix<double, n_outputs, n_outputs>;
  using OutputPrediction = Eigen::Matrix<double, p_max * n_outputs, 1>;

  NerveCenter(const System& sys, std::tuple<SubControllers...>& controllers)
      : ControllerInterface<System>(Input::Zero()),
        sys_(sys),
        sub_controllers_(controllers) {}

  /// Initialize all sub controllers based on given initial conditions
  void Initialize(const State& x_init, const ControlInput& u_init,
                  const Input& u_init_full, const Output& y_init,
                  const AugmentedState& dx_init = AugmentedState::Zero()) {
    expander{InitializeHelper(&std::get<SubControllers>(sub_controllers_),
                              x_init, u_init, u_init_full, y_init, dx_init)...};
    u_offset_ = u_init_full;
  }

  /// Set input/output weights for all controllers
  void SetWeights(const UWeightType& uwt, const YWeightType& ywt) {
    expander{SetWeightsHelper(&std::get<SubControllers>(sub_controllers_), uwt,
                              ywt)...};
  }

  /// Set reference output
  void SetOutputReference(const OutputPrediction& y_ref) {
    expander{SetOutputReferenceHelper(
        &std::get<SubControllers>(sub_controllers_), y_ref)...};
  }

  // TODO: Write this function
  virtual const ControlInput GetNextInput(const Output& y) {
    return ControlInput::Zero();
  }

 private:
  // Initialize each controller
  template <typename T>
  int InitializeHelper(T* controller, const State& x_init,
                       const ControlInput& u_init, const Input& u_init_full,
                       const Output& y_init, const AugmentedState& dx_init) {
    typename T::State x_init_sub;
    typename T::ControlInput u_init_sub;
    typename T::Output y_init_sub;
    typename T::AugmentedState dx_init_sub;

    T::StateIndexType::template IndicesSubArray<std::make_integer_sequence<
        int, T::n_states>>::GetSubArray(x_init_sub.data(), x_init.data());
    T::ControlInputIndexType::GetSubArray(u_init_sub.data(), u_init.data());
    T::ObserverOutputIndexType::GetSubArray(y_init_sub.data(), y_init.data());
    T::StateIndexType::GetSubArray(dx_init_sub.data(), dx_init.data());

    controller->Initialize(x_init_sub, u_init_sub, u_init_full, y_init_sub,
                           dx_init_sub);
  }

  // Set weights for each controller
  template <typename T>
  int SetWeightsHelper(T* controller, const UWeightType& uwt,
                       const YWeightType& ywt) {
    typename T::UWeightType uwt_sub;
    typename T::YWeightType ywt_sub;

    T::ControlInputIndexType::template GetSubMatrix<n_control_inputs>(
        uwt_sub.data(), uwt.data());
    T::ControlledOutputIndexType::template GetSubMatrix<n_outputs>(
        ywt_sub.data(), ywt.data());

    controller->SetWeights(uwt_sub, ywt_sub);
  }

  // Set reference output for each controller
  template <typename T>
  int SetOutputReferenceHelper(T* controller, const OutputPrediction& y_ref) {
    Eigen::Matrix<double, T::p * T::n_controlled_outputs, 1> y_ref_sub;

    const double* orig_data = y_ref.data();
    double* new_data = y_ref_sub.data();

    // Loop through each prediction number and take correct outputs
    for (auto i = 0; i < T::p; i++) {
      T::ControlledOutputIndexType::GetSubArray(new_data, orig_data);
      new_data += T::n_controlled_outputs;
      orig_data += n_outputs;
    }

    controller->SetOutputReference(y_ref_sub);
  }
};

#endif
