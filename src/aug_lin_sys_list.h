#ifndef AUG_LIN_SYS_LIST_H
#define AUG_LIN_SYS_LIST_H
#include "parallel_compressors.h"
#include "prediction.h"
#include "constexpr_array.h"

using System1 =
    AugmentedLinearizedSystem<ParallelCompressors, ConstexprArray<0, 40, 0, 40>,
                              4, ConstexprArray<0, 1, 2, 3>, 2>;
using System2 =
    AugmentedLinearizedSystem<ParallelCompressors, ConstexprArray<0, 40, 0, 40>,
                              4, ConstexprArray<2, 3, 0, 1>, 2>;

// templates for GeneratePrediction function
template const Prediction System1::GeneratePrediction<ConstexprArray<0, 1, 3>>(
    Eigen::MatrixXd* Su_other, const int p, const int m) const;

template void System1::GeneratePrediction<ConstexprArray<0, 1, 3>>(
    Eigen::MatrixXd* Su, Eigen::MatrixXd* Sx, Eigen::MatrixXd* Sf,
    Eigen::MatrixXd* Su_other, const int p, const int m) const;

template const Prediction System2::GeneratePrediction<ConstexprArray<0, 1, 3>>(
    Eigen::MatrixXd* Su_other, const int p, const int m) const;

template void System2::GeneratePrediction<ConstexprArray<0, 1, 3>>(
    Eigen::MatrixXd* Su, Eigen::MatrixXd* Sx, Eigen::MatrixXd* Sf,
    Eigen::MatrixXd* Su_other, const int p, const int m) const;

// class template explicit instantiation
template class AugmentedLinearizedSystem<ParallelCompressors,
                                         ConstexprArray<0, 40, 0, 40>, 4,
                                         ConstexprArray<0, 1, 2, 3>, 2>;
template class AugmentedLinearizedSystem<ParallelCompressors,
                                         ConstexprArray<0, 40, 0, 40>, 4,
                                         ConstexprArray<2, 3, 0, 1>, 2>;

#endif
