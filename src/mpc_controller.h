#ifndef MPC_CONTROLLER_H
#define MPC_CONTROLLER_H

#include <Eigen/Eigen>

template <class AugmentedSystem, int p, int m>
class MpcController {
 protected:
  static constexpr int n_control_inputs = AugmentedSystem::n_control_inputs;
  static constexpr int n_inputs = AugmentedSystem::n_inputs;
  static constexpr int n_outputs = AugmentedSystem::n_outputs;
  static constexpr int n_states = AugmentedSystem::n_total_states;

 public:
  typedef Eigen::Matrix<double, n_states, 1> State;
  typedef Eigen::Matrix<double, n_inputs, 1> Input;
  typedef Eigen::Matrix<double, n_control_inputs, 1> ControlInput;
  typedef Eigen::Matrix<double, n_outputs, 1> Output;
  typedef Eigen::Matrix<double, n_outputs * p, 1> OutputPrediction;
  typedef Eigen::Matrix<double, n_control_inputs, n_control_inputs> UWeightType;
  typedef Eigen::Matrix<double, n_outputs, n_outputs> YWeightType;

  MpcController(AugmentedSystem augsys, Output y_ref,
                UWeightType u_weight = UWeightType::Identity(),
                YWeightType y_weight = YWeightType::Identity())
      : MpcController(augsys, static_cast<OutputPrediction>(
                                  y_ref.template replicate<p, 1>()),
                      u_weight, y_weight) {}

  MpcController(AugmentedSystem augsys,
                OutputPrediction y_ref = OutputPrediction::Zero(),
                UWeightType u_weight = UWeightType::Identity(),
                YWeightType y_weight = YWeightType::Identity())
      : augsys_(augsys),
        y_ref_(y_ref),
        u_weight_(u_weight),
        y_weight_(y_weight) {}

  Input GetNextInput(Output y);

  void SetReference(OutputPrediction y_ref) { y_ref_ = y_ref; }
  OutputPrediction GetReference() { return y_ref_; }

  struct Prediction {
    Eigen::Matrix<double, p * n_outputs, n_states, Eigen::RowMajor> Sx;
    Eigen::Matrix<double, p * n_outputs, m * n_control_inputs, Eigen::RowMajor>
        Su;
    Eigen::Matrix<double, p * n_outputs, n_states, Eigen::RowMajor> Sf;
  };

  Prediction GetPredictionMatrices() const;

 protected:
  AugmentedSystem augsys_;
  OutputPrediction y_ref_;

  const UWeightType u_weight_;
  const YWeightType y_weight_;
  const Eigen::Matrix<double, 2 * m * n_control_inputs, m * n_control_inputs>
      Ain_;

  void UpdateStateEstimate(Output y_in) { augsys_.ObserveAPosteriori(y_in); }

  Input SolveQP() const;

  void MakeStatePrediction(ControlInput u_in) { augsys_.ObserveAPriori(u_in); }
};

#endif
