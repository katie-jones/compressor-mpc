#ifndef AUG_LIN_SYS_LIST_H
#define AUG_LIN_SYS_LIST_H
#include "parallel_compressors.h"
#include "prediction.h"
#include "constexpr_array.h"
#include "parallel_compressors_constants.h"

// templates for GeneratePrediction function
template const Prediction AUGMENTEDSYSTEM1::GeneratePrediction<ConstexprArray<0, 1, 3>>(
    Eigen::MatrixXd* Su_other, const int p, const int m) const;

template void AUGMENTEDSYSTEM1::GeneratePrediction<ConstexprArray<0, 1, 3>>(
    Eigen::MatrixXd* Su, Eigen::MatrixXd* Sx, Eigen::MatrixXd* Sf,
    Eigen::MatrixXd* Su_other, const int p, const int m) const;

template const Prediction AUGMENTEDSYSTEM2::GeneratePrediction<ConstexprArray<0, 1, 3>>(
    Eigen::MatrixXd* Su_other, const int p, const int m) const;

template void AUGMENTEDSYSTEM2::GeneratePrediction<ConstexprArray<0, 1, 3>>(
    Eigen::MatrixXd* Su, Eigen::MatrixXd* Sx, Eigen::MatrixXd* Sf,
    Eigen::MatrixXd* Su_other, const int p, const int m) const;

// class template explicit instantiation
template class AUGMENTEDSYSTEM1;
template class AUGMENTEDSYSTEM2;

#endif
