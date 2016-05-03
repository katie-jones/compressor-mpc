#ifndef MPC_CONTROLLER_H
#define MPC_CONTROLLER_H

#include <Eigen/Eigen>

template <int n_states, int n_inputs, int n_outputs, int n_control_inputs>
class MpcController {
 public:
  typedef Eigen::Matrix<double, n_states, 1> State;
  typedef Eigen::Matrix<double, n_inputs, 1> PlantInput;
  typedef Eigen::Matrix<double, n_control_inputs, 1> Input;
  typedef Eigen::Matrix<double, n_outputs, 1> Output;

 protected:
  DynamicSystem<n_states, n_inputs, n_outputs, n_control_inputs> *p_sys_;
  int p_;
  int m_;
  State xbar_;
  State dxbar_;
  Input upast_;
  Output y_old_;
  Output y_ref_;

  PlantInput u_offset_;
  int u_control_index_[n_control_inputs];

  const Eigen::Matrix<double, n_states, n_outputs> M_;
  const Eigen::Matrix<double, n_control_inputs, n_control_inputs> u_weight_;
  const Eigen::Matrix<double, n_outputs, n_outputs> y_weight_;
  const static Eigen::Matrix<double, 2 * m * n_control_inputs,
                             m * n_control_inputs> Ain_;

  struct Prediction {
    Eigen::Matrix<double, p * n_outputs, n_states, Eigen::RowMajor> Sx;
    Eigen::Matrix<double, p * n_outputs, n_control_inputs, Eigen::RowMajor> Su;
    Eigen::Matrix<double, p * n_outputs, n_states, Eigen::RowMajor> Sf;
  };

  void UpdateStateEstimate(Output y_in);

  Prediction GetPredictionMatrices() const;

  PlantInput SolveQP() const;

  void MakeStatePrediction(Input u_in);
};

#endif
