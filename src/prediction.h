#ifndef PREDICTION_H
#define PREDICTION_H

#include <Eigen/Eigen>

/// Linearized prediction of dynamic system
struct Prediction {
  Eigen::MatrixXd Su, Sx, Sf;
};

#endif
