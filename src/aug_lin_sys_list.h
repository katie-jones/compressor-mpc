#ifndef AUG_LIN_SYS_LIST_H
#define AUG_LIN_SYS_LIST_H
#include "parallel_compressors.h"
#include "prediction.h"

template const Prediction
AugmentedLinearizedSystem<ParallelCompressors, 80, 4>::GenerateSubPrediction<2>(
    const int p, const int m, const int* output_index) const;

template const Prediction
AugmentedLinearizedSystem<ParallelCompressors, 80, 4>::GenerateSubPrediction<3>(
    const int p, const int m, const int* output_index) const;

template void
AugmentedLinearizedSystem<ParallelCompressors, 80, 4>::GenerateSubPrediction<2>(
    Eigen::MatrixXd* Su, Eigen::MatrixXd* Sx, Eigen::MatrixXd* Sf, const int p,
    const int m, const int* output_index) const;

template const Prediction
AugmentedLinearizedSystem<ParallelCompressors, 80, 4>::GenerateSubPrediction<4>(
    const int p, const int m, const int* output_index) const;

template void
AugmentedLinearizedSystem<ParallelCompressors, 80, 4>::GenerateSubPrediction<4>(
    Eigen::MatrixXd* Su, Eigen::MatrixXd* Sx, Eigen::MatrixXd* Sf, const int p,
    const int m, const int* output_index) const;

// Parallel compressors
template class AugmentedLinearizedSystem<ParallelCompressors, 80, 4>;

#endif
