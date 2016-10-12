#ifndef PREDICTION_H
#define PREDICTION_H

#include <Eigen/Eigen>

/**
 * Linearized prediction of dynamic system.
 */
struct Prediction {
  /// Prediction matrix for control inputs of a sub-controller
  Eigen::MatrixXd Su;
  /// Prediction matrix for initial augmented state
  Eigen::MatrixXd Sx;
  /// Prediction matrix for derivative of system at linearization point
  Eigen::MatrixXd Sf;
  /// Prediction matrix for control inputs from other sub-controllers
  Eigen::MatrixXd Su_other;

  /// Empty constructor
  Prediction() : Su(), Sx(), Sf(), Su_other() {}

  /// Constructor with prediction matrices given
  Prediction(const Eigen::MatrixXd& Su_in, const Eigen::MatrixXd& Sx_in,
             const Eigen::MatrixXd& Sf_in,
             const Eigen::MatrixXd& Su_other_in = Eigen::MatrixXd())
      : Su(Su_in), Sx(Sx_in), Sf(Sf_in), Su_other(Su_other_in) {}
};

#endif
