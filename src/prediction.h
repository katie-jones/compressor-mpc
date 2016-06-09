#ifndef PREDICTION_H
#define PREDICTION_H

#include <Eigen/Eigen>

/// Linearized prediction of dynamic system
struct Prediction {
  Eigen::MatrixXd Su, Sx, Sf, Su_other;
  Prediction() : Su(), Sx(), Sf(), Su_other() {}
  Prediction(const Eigen::MatrixXd& Su_in, const Eigen::MatrixXd& Sx_in,
             const Eigen::MatrixXd& Sf_in,
             const Eigen::MatrixXd& Su_other_in = Eigen::MatrixXd())
      : Su(Su_in), Sx(Sx_in), Sf(Sf_in), Su_other(Su_other_in) {}
};

#endif
