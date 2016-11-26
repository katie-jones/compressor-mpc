#ifndef NERVE_CENTER_H
#define NERVE_CENTER_H

#include <tuple>

#include "controller_interface.h"
#include "distributed_controller.h"

/**
 * Nerve center and communication hub for distributed MPC controller.
 * Coordinates exchange of information between the plant and controllers, as
 * well as between the different sub-controllers.
 * Template parameters:\n
 * - System: class of dynamic system to be controlled\n
 * - n_total_states: total number of (augmented + real) states\n
 * - SubControllers: parameter pack of sub-controllers that make up the\n
 * distributed MPC controller\n
 */
template <typename System, int n_total_states, typename... SubControllers>
class NerveCenter : public ControllerInterface<System> {
 protected:
  static constexpr int n_states = System::n_states;
  static constexpr int n_inputs = System::n_inputs;
  static constexpr int n_outputs = System::n_outputs;
  static constexpr int n_control_inputs = System::n_control_inputs;
  static constexpr int n_controllers = sizeof...(SubControllers);

  using expander = int[];

  static constexpr int GetPredictionControlInputs() {
    int u_total = 0;
    int dummy[]{
        (u_total += SubControllers::m * SubControllers::n_control_inputs)...};
    return u_total;
  }

  static constexpr int n_prediction_control_inputs =
      GetPredictionControlInputs();

 public:
  /// State of system
  using State = Eigen::Matrix<double, System::n_states, 1>;
  /// Augmented state of system
  using AugmentedState = Eigen::Matrix<double, n_total_states, 1>;
  /// Input to system
  using Input = Eigen::Matrix<double, System::n_inputs, 1>;
  /// Output of system
  using Output = Eigen::Matrix<double, System::n_outputs, 1>;
  /// Control input to system
  using ControlInput = Eigen::Matrix<double, System::n_control_inputs, 1>;
  /// Prediction of control input over move horizon
  using ControlInputPrediction =
      Eigen::Matrix<double, n_prediction_control_inputs, 1>;

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

  // Previous applied control input
  ControlInput u_old_;

  // Previous control input predictions
  ControlInputPrediction du_old_;

  // Number of solver iterations
  const int n_solver_iterations_;

 public:
  /// Matrix of weights on inputs
  using UWeightType = Eigen::Matrix<double, n_control_inputs, n_control_inputs>;
  /// Matrix of weights on outputs
  using YWeightType = Eigen::Matrix<double, n_outputs, n_outputs>;
  /// Tuple of weights on outputs for each sub-controller
  using SubYWeightType = std::tuple<typename SubControllers::YWeightType...>;
  /// Prediction of output over maximum prediction horizon
  using OutputPrediction = Eigen::Matrix<double, p_max * n_outputs, 1>;

  /// Constructor
  NerveCenter(std::tuple<SubControllers...>& controllers,
              const int n_solver_iterations)
      : ControllerInterface<System>(Input::Zero()),
        n_solver_iterations_(n_solver_iterations),
        u_old_(ControlInput::Zero()),
        du_old_(ControlInputPrediction::Zero()),
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

  /// Set input/output weights individually for each controller
  void SetWeights(const UWeightType& uwt, const SubYWeightType& ywts) {
    expander{SetWeightsSubHelper(&std::get<SubControllers>(sub_controllers_),
                                 uwt, ywts)...};
  }

  /// Set reference output
  void SetOutputReference(const OutputPrediction& y_ref) {
    expander{SetOutputReferenceHelper(
        &std::get<SubControllers>(sub_controllers_), y_ref)...};
  }

  /// Solve QPs and get next input to apply
  virtual const ControlInput GetNextInput(const Output& y) {
    return GetNextInputWithTiming(y);
  }

  /// Solve QPs and get next input to apply.
  /// If n_timing_iterations is given, the CPU time required to generate the
  /// QPs, solve them n_timing_iterations times and return the solution is
  /// measured and returned in time_out. (0 <= n_timing_iterations <
  /// n_solver_iterations_)
  const ControlInput GetNextInputWithTiming(
      const Output& y, const int n_timing_iterations = -1,
      boost::timer::nanosecond_type* time_out = NULL) {
    // Start timer
    boost::timer::cpu_timer timer;

    const Input& u_full_old = System::GetPlantInput(u_old_, u_offset_);
    // Initialize all QPs
    expander{InitializeQPHelper(&std::get<SubControllers>(sub_controllers_), y,
                                u_full_old)...};

    // Solve QPs n_solver_iterations_ times
    ControlInputPrediction du_prev = du_old_;
    ControlInputPrediction du_new = ControlInputPrediction::Zero();

    int prediction_index = 0;
    for (int i = 0; i < n_solver_iterations_; i++) {
      if (i == n_timing_iterations) timer.stop();

      expander{SolveQPHelper(&std::get<SubControllers>(sub_controllers_),
                             &du_new, &prediction_index, du_prev)...};

      du_prev = du_new;
      prediction_index = 0;
    }
    if (timer.is_stopped()) timer.resume();

    // Update du_old_, u_old_
    du_old_ = du_prev;

    ControlInput du = -u_old_;  // difference between curr and prev solution

    int input_index = 0;
    expander{UpdateUOld<SubControllers>(&prediction_index, &input_index)...};

    // Send solutions to subcontrollers
    du += u_old_;

    expander{SendUHelper(&std::get<SubControllers>(sub_controllers_), du)...};

    // Assign total central time
    timer.stop();
    if (time_out) {
      boost::timer::cpu_times elapsed = timer.elapsed();
      *time_out = elapsed.wall;
    }

    return u_old_;
  }

 private:
  // Initialize each controller
  template <typename T>
  int InitializeHelper(T* controller, const State& x_init,
                       const ControlInput& u_init, const Input& u_init_full,
                       const Output& y_init, const AugmentedState& dx_init) {
    typename T::State x_init_sub;
    typename T::FullControlInput u_init_sub;
    typename T::Output y_init_sub;
    typename T::AugmentedState dx_init_sub;

    T::StateIndexType::template IndicesSubArray<std::make_integer_sequence<
        int, T::n_states> >::GetSubArray(x_init_sub.data(), x_init.data());
    T::FullControlInputIndexType::GetSubArray(u_init_sub.data(), u_init.data());
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

  // Set weights for each controller
  template <typename T>
  int SetWeightsSubHelper(T* controller, const UWeightType& uwt,
                          const SubYWeightType& ywts) {
    typename T::YWeightType ywt_sub = std::get<GetControllerNumber<T>()>(ywts);
    typename T::UWeightType uwt_sub;

    T::ControlInputIndexType::template GetSubMatrix<n_control_inputs>(
        uwt_sub.data(), uwt.data());

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

  // Initialize QP for each controller
  template <typename T>
  int InitializeQPHelper(T* controller, const Output& y,
                         const Input& u_full_old) {
    typename T::Output y_sub;
    T::ObserverOutputIndexType::GetSubArray(y_sub.data(), y.data());
    controller->GenerateInitialQP(y_sub, u_full_old);
  }

  // Time InitializeQPHelper
  template <typename T>
  int TimeInitializeQPHelper(boost::timer::nanosecond_type*& time_out,
                             T* controller, const Output& y,
                             const Input& u_full_old) {
    boost::timer::cpu_timer timer;

    InitializeQPHelper(controller, y, u_full_old);

    boost::timer::cpu_times elapsed = timer.elapsed();
    *time_out += elapsed.wall;
    time_out++;
  }

  // Solve QP for each controller
  template <typename T>
  int SolveQPHelper(T* controller, ControlInputPrediction* du_new,
                    int* prediction_index,
                    const ControlInputPrediction& du_old) {
    // Take previous solution from other controllers (not this one)
    Eigen::Matrix<double,
                  n_prediction_control_inputs - T::m * T::n_control_inputs, 1>
        du_old_sub;
    du_old_sub << du_old.head(*prediction_index),
        du_old.tail(n_prediction_control_inputs - *prediction_index -
                    T::m * T::n_control_inputs);

    typename T::ControlInputPrediction du_new_sub;

    // Get next control input from controller
    controller->GetInput(&du_new_sub, du_old_sub);

    // Copy new values and update index
    du_new->template segment<T::m * T::n_control_inputs>(*prediction_index) =
        du_new_sub;
    *prediction_index += T::m * T::n_control_inputs;
  }

  // Time SolveQPHelper
  template <typename T>
  int TimeSolveQPHelper(boost::timer::nanosecond_type*& time_out, T* controller,
                        ControlInputPrediction* du_new, int* prediction_index,
                        const ControlInputPrediction& du_old) {
    boost::timer::cpu_timer timer;

    SolveQPHelper(controller, du_new, prediction_index, du_old);

    boost::timer::cpu_times elapsed = timer.elapsed();
    *time_out += elapsed.wall;
    time_out++;
  }

  // Add new solution (in du_old_) to u_old_
  template <typename T>
  int UpdateUOld(int* prediction_index, int* input_index) {
    u_old_.template segment<T::n_control_inputs>(*input_index) +=
        du_old_.template segment<T::n_control_inputs>(*prediction_index);
    *prediction_index += T::m * T::n_control_inputs;
    *input_index += T::n_control_inputs;
  }

  // Send new solution to subcontrollers
  template <typename T>
  int SendUHelper(T* controller, const ControlInput& du) {
    ControlInput du_reordered;
    du_reordered.setZero();
    T::ControlInputIndexType::GetSubArray(du_reordered.data(), du.data());
    controller->UpdateU(du_reordered);
  }

  // Get the index of controller of typename T
  template <typename T>
  static constexpr int GetControllerNumber(const T* controller = NULL) {
    int ind = -1;
    bool is_right_type[]{(std::is_same<T, SubControllers>::value)...};

    for (auto i = 0; i < n_controllers; i++) {
      if (is_right_type[i]) {
        ind = i;
      }
    }
    return ind;
  }
};

#endif
