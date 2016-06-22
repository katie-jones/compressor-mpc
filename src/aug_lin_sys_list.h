#ifndef AUG_LIN_SYS_LIST_H
#define AUG_LIN_SYS_LIST_H
#include "constexpr_array.h"
#include "parallel_compressors.h"
#include "parallel_compressors_constants.h"
#include "prediction.h"

// templates for GeneratePrediction function
template const Prediction
AUGMENTEDSYSTEM_COOP1::GeneratePrediction<ConstexprArray<0, 1, 3>>(
    Eigen::MatrixXd* Su_other, const int p, const int m) const;

template void AUGMENTEDSYSTEM_COOP1::GeneratePrediction<
    ConstexprArray<0, 1, 3>>(Eigen::MatrixXd* Su, Eigen::MatrixXd* Sx,
                             Eigen::MatrixXd* Sf, Eigen::MatrixXd* Su_other,
                             const int p, const int m) const;

template const Prediction
AUGMENTEDSYSTEM_COOP2::GeneratePrediction<ConstexprArray<0, 1, 3>>(
    Eigen::MatrixXd* Su_other, const int p, const int m) const;

template void AUGMENTEDSYSTEM_COOP2::GeneratePrediction<
    ConstexprArray<0, 1, 3>>(Eigen::MatrixXd* Su, Eigen::MatrixXd* Sx,
                             Eigen::MatrixXd* Sf, Eigen::MatrixXd* Su_other,
                             const int p, const int m) const;

template const Prediction
AUGMENTEDSYSTEM_CENTRALIZED::GeneratePrediction<ConstexprArray<0, 1, 3>>(
    Eigen::MatrixXd* Su_other, const int p, const int m) const;

template void AUGMENTEDSYSTEM_CENTRALIZED::GeneratePrediction<
    ConstexprArray<0, 1, 3>>(Eigen::MatrixXd* Su, Eigen::MatrixXd* Sx,
                             Eigen::MatrixXd* Sf, Eigen::MatrixXd* Su_other,
                             const int p, const int m) const;

// class template explicit instantiation
template class AUGMENTEDSYSTEM_COOP1;
template class AUGMENTEDSYSTEM_COOP2;
template class AUGMENTEDSYSTEM_CENTRALIZED;

#endif
