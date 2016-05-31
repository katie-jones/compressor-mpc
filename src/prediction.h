#ifndef PREDICTION_H
#define PREDICTION_H

#include <Eigen/Eigen>

/// Linearized prediction of dynamic system
struct Prediction {
  Eigen::MatrixXd Su, Sx, Sf;
  Prediction() : Su(), Sx(), Sf() {}
  Prediction(const Eigen::MatrixXd& Su_in, const Eigen::MatrixXd& Sx_in,
             const Eigen::MatrixXd& Sf_in)
      : Su(Su_in), Sx(Sx_in), Sf(Sf_in) {}
};

#endif
