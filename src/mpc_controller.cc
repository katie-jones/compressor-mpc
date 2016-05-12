#include "mpc_controller.h"
#include "print_matrix.h"
#include <iostream>

template <class AugmentedSystem, int p, int m>
typename MpcController<AugmentedSystem, p, m>::Prediction
MpcController<AugmentedSystem, p, m>::GetPredictionMatrices() const {
  Prediction pred;
  typename AugmentedSystem::AugmentedLinearizedSystem linsys =
      augsys_.GetLinearization();
  Eigen::Matrix<double, n_outputs, n_states> c_times_a = linsys.C;
  Eigen::Matrix<double, n_outputs, n_control_inputs> to_add;
  int ind_col, ind_row;

  pred.Su.setZero();
  pred.Sx.setZero();
  pred.Sf.setZero();

  pred.Sf.template topRows<n_outputs>() = linsys.C;

  for (int i = 0; i < p; i++) {
    if (i > 0) {
      pred.Sf.template block<n_outputs, n_states>(i * n_outputs, 0) =
          pred.Sf.template block<n_outputs, n_states>((i - 1) * n_outputs, 0) +
          c_times_a;
    }

    to_add = c_times_a * linsys.B;
    for (int j = 0; j < p - i; j++) {
      ind_row = i + j;
      if (j < m)
        ind_col = j;
      else
        ind_col = m - 1;
      pred.Su.template block<n_outputs, n_control_inputs>(
          ind_row * n_outputs, ind_col * n_control_inputs) += to_add;
    }
    c_times_a *= linsys.A;
    pred.Sx.template block<n_outputs, n_states>(i * n_outputs, 0) = c_times_a;
  }
  print_matrix(std::cout, linsys.A, "A");
  print_matrix(std::cout, linsys.B, "B");
  print_matrix(std::cout, linsys.C, "C");
  print_matrix(std::cout, linsys.f, "f");

  print_matrix(std::cout, pred.Su, "Su");
  print_matrix(std::cout, pred.Sf, "Sf");
  print_matrix(std::cout, pred.Sx, "Sx");
}

#include "mpc_controller_list.h"
