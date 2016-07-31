#ifndef AUG_LIN_SYS_LIST_H
#define AUG_LIN_SYS_LIST_H
#include "constexpr_array.h"
#include "parallel_compressors.h"
#include "parallel_compressors_constants.h"
#include "serial_compressors.h"
#include "serial_compressors_constants.h"
#include "prediction.h"

// -----------------------------------------------------------------------------
// ------------------------- PARALLEL ------------------------------------------
// -----------------------------------------------------------------------------

// templates for GeneratePrediction function
// Cooperative
template const Prediction PARALLEL_AUGSYS_DIST1::GeneratePrediction<
    PARALLEL_COMPRESSORS_CONSTANTS::ControlledOutputIndices>(
    Eigen::MatrixXd* Su_other, const int p, const int m) const;

template void PARALLEL_AUGSYS_DIST1::GeneratePrediction<
    PARALLEL_COMPRESSORS_CONSTANTS::ControlledOutputIndices>(
    Eigen::MatrixXd* Su, Eigen::MatrixXd* Sx, Eigen::MatrixXd* Sf,
    Eigen::MatrixXd* Su_other, const int p, const int m) const;

template const Prediction PARALLEL_AUGSYS_DIST2::GeneratePrediction<
    PARALLEL_COMPRESSORS_CONSTANTS::ControlledOutputIndices>(
    Eigen::MatrixXd* Su_other, const int p, const int m) const;

template void PARALLEL_AUGSYS_DIST2::GeneratePrediction<
    PARALLEL_COMPRESSORS_CONSTANTS::ControlledOutputIndices>(
    Eigen::MatrixXd* Su, Eigen::MatrixXd* Sx, Eigen::MatrixXd* Sf,
    Eigen::MatrixXd* Su_other, const int p, const int m) const;

// Noncooperative
template const Prediction PARALLEL_AUGSYS_DIST1::GeneratePrediction<
    PARALLEL_COMPRESSORS_CONSTANTS::NCControlledOutputIndices1>(
    Eigen::MatrixXd* Su_other, const int p, const int m) const;

template void PARALLEL_AUGSYS_DIST1::GeneratePrediction<
    PARALLEL_COMPRESSORS_CONSTANTS::NCControlledOutputIndices1>(
    Eigen::MatrixXd* Su, Eigen::MatrixXd* Sx, Eigen::MatrixXd* Sf,
    Eigen::MatrixXd* Su_other, const int p, const int m) const;

template const Prediction PARALLEL_AUGSYS_DIST2::GeneratePrediction<
    PARALLEL_COMPRESSORS_CONSTANTS::NCControlledOutputIndices2>(
    Eigen::MatrixXd* Su_other, const int p, const int m) const;

template void PARALLEL_AUGSYS_DIST2::GeneratePrediction<
    PARALLEL_COMPRESSORS_CONSTANTS::NCControlledOutputIndices2>(
    Eigen::MatrixXd* Su, Eigen::MatrixXd* Sx, Eigen::MatrixXd* Sf,
    Eigen::MatrixXd* Su_other, const int p, const int m) const;

// Centralized
template const Prediction PARALLEL_AUGSYS_CENT::GeneratePrediction<
    PARALLEL_COMPRESSORS_CONSTANTS::ControlledOutputIndices>(
    Eigen::MatrixXd* Su_other, const int p, const int m) const;

template void PARALLEL_AUGSYS_CENT::GeneratePrediction<
    PARALLEL_COMPRESSORS_CONSTANTS::ControlledOutputIndices>(
    Eigen::MatrixXd* Su, Eigen::MatrixXd* Sx, Eigen::MatrixXd* Sf,
    Eigen::MatrixXd* Su_other, const int p, const int m) const;

// class template explicit instantiation
template class PARALLEL_AUGSYS_DIST1;
template class PARALLEL_AUGSYS_DIST2;
template class PARALLEL_AUGSYS_CENT;

// -----------------------------------------------------------------------------
// --------------------------- SERIAL ------------------------------------------
// -----------------------------------------------------------------------------

// templates for GeneratePrediction function
// Cooperative
template const Prediction SERIAL_AUGSYS_DIST1::GeneratePrediction<
    SERIAL_COMPRESSORS_CONSTANTS::ControlledOutputIndices>(
    Eigen::MatrixXd* Su_other, const int p, const int m) const;

template void SERIAL_AUGSYS_DIST1::GeneratePrediction<
    SERIAL_COMPRESSORS_CONSTANTS::ControlledOutputIndices>(
    Eigen::MatrixXd* Su, Eigen::MatrixXd* Sx, Eigen::MatrixXd* Sf,
    Eigen::MatrixXd* Su_other, const int p, const int m) const;

template const Prediction SERIAL_AUGSYS_DIST2::GeneratePrediction<
    SERIAL_COMPRESSORS_CONSTANTS::ControlledOutputIndices>(
    Eigen::MatrixXd* Su_other, const int p, const int m) const;

template void SERIAL_AUGSYS_DIST2::GeneratePrediction<
    SERIAL_COMPRESSORS_CONSTANTS::ControlledOutputIndices>(
    Eigen::MatrixXd* Su, Eigen::MatrixXd* Sx, Eigen::MatrixXd* Sf,
    Eigen::MatrixXd* Su_other, const int p, const int m) const;

// Noncooperative
template const Prediction SERIAL_AUGSYS_DIST1::GeneratePrediction<
    SERIAL_COMPRESSORS_CONSTANTS::NCControlledOutputIndices1>(
    Eigen::MatrixXd* Su_other, const int p, const int m) const;

template void SERIAL_AUGSYS_DIST1::GeneratePrediction<
    SERIAL_COMPRESSORS_CONSTANTS::NCControlledOutputIndices1>(
    Eigen::MatrixXd* Su, Eigen::MatrixXd* Sx, Eigen::MatrixXd* Sf,
    Eigen::MatrixXd* Su_other, const int p, const int m) const;

template const Prediction SERIAL_AUGSYS_DIST2::GeneratePrediction<
    SERIAL_COMPRESSORS_CONSTANTS::NCControlledOutputIndices2>(
    Eigen::MatrixXd* Su_other, const int p, const int m) const;

template void SERIAL_AUGSYS_DIST2::GeneratePrediction<
    SERIAL_COMPRESSORS_CONSTANTS::NCControlledOutputIndices2>(
    Eigen::MatrixXd* Su, Eigen::MatrixXd* Sx, Eigen::MatrixXd* Sf,
    Eigen::MatrixXd* Su_other, const int p, const int m) const;

// Noncooperative unstable
template const Prediction SERIAL_AUGSYS_DIST1::GeneratePrediction<
    SERIAL_COMPRESSORS_CONSTANTS::OldNCControlledOutputIndices1>(
    Eigen::MatrixXd* Su_other, const int p, const int m) const;

template void SERIAL_AUGSYS_DIST1::GeneratePrediction<
    SERIAL_COMPRESSORS_CONSTANTS::OldNCControlledOutputIndices1>(
    Eigen::MatrixXd* Su, Eigen::MatrixXd* Sx, Eigen::MatrixXd* Sf,
    Eigen::MatrixXd* Su_other, const int p, const int m) const;

template const Prediction SERIAL_AUGSYS_DIST2::GeneratePrediction<
    SERIAL_COMPRESSORS_CONSTANTS::OldNCControlledOutputIndices2>(
    Eigen::MatrixXd* Su_other, const int p, const int m) const;

template void SERIAL_AUGSYS_DIST2::GeneratePrediction<
    SERIAL_COMPRESSORS_CONSTANTS::OldNCControlledOutputIndices2>(
    Eigen::MatrixXd* Su, Eigen::MatrixXd* Sx, Eigen::MatrixXd* Sf,
    Eigen::MatrixXd* Su_other, const int p, const int m) const;

// Centralized
template const Prediction SERIAL_AUGSYS_CENT::GeneratePrediction<
    SERIAL_COMPRESSORS_CONSTANTS::ControlledOutputIndices>(
    Eigen::MatrixXd* Su_other, const int p, const int m) const;

template void SERIAL_AUGSYS_CENT::GeneratePrediction<
    SERIAL_COMPRESSORS_CONSTANTS::ControlledOutputIndices>(
    Eigen::MatrixXd* Su, Eigen::MatrixXd* Sx, Eigen::MatrixXd* Sf,
    Eigen::MatrixXd* Su_other, const int p, const int m) const;

// class template explicit instantiation
template class SERIAL_AUGSYS_DIST1;
template class SERIAL_AUGSYS_DIST2;
template class SERIAL_AUGSYS_CENT;

#endif
