#ifndef PRINT_MATRIX_H
#define PRINT_MATRIX_H

#include <iostream>
#include <string>
#include <Eigen/Eigen>
#include <iomanip>

template <typename Derived>
static std::ostream& print_matrix(std::ostream& os,
                                  const Eigen::MatrixBase<Derived>& x,
                                  const std::string name) {
  if (!name.empty()) os << name << " = " << std::endl;
  for (int i = 0; i < x.rows(); i++) {
    for (int j = 0; j < x.cols(); j++) {
      os << std::setprecision(3) << std::setw(10) << x(i, j) << " ";
    }
    os << std::endl;
  }
  os << std::endl;
  return os;
}

static std::ostream& print_pointer(std::ostream& os, const double* p_data,
                                   const int m, const int n,
                                   const std::string name) {
  os << name << " = " << std::endl;
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++)
      os << std::setprecision(3) << std::setw(10) << p_data[j + i * n] << " ";
    os << std::endl;
  }
  os << std::endl;
  return os;
}
#endif
