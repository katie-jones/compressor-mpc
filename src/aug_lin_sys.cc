#include "mpc_controller.h"

/*
 * Default values for linearization
 */
template <class System, int n_delay_states, int n_disturbance_states, int p,
          int m>
MpcController<System, n_delay_states, n_disturbance_states, p, m>::
    AugmentedLinearizedSystemTwo::AugmentedLinearizedSystemTwo(
        const ControlInputIndex& n_delay_in)
    : A(AComposite(n_delay_in)),
      B(Eigen::MatrixXd::Zero(n_total_states, n_control_inputs)),
      C(Eigen::MatrixXd::Zero(n_outputs, n_obs_states)),
      f(System::State::Zero()),
      n_delay(n_delay_in) {
  C.template rightCols<n_disturbance_states>() =
      Eigen::Matrix<double, n_outputs, n_disturbance_states>::Identity();

  int index_delay_states = n_obs_states;
  for (int i = 0; i < n_control_inputs; i++) {
    if (n_delay[i] != 0) {
      index_delay_states += n_delay[i];
      B(index_delay_states - 1, i) = 1;
    }
  }
}

template <class System, int n_delay_states, int n_disturbance_states, int p,
          int m>
MpcController<System, n_delay_states, n_disturbance_states, p, m>::
    AugmentedLinearizedSystemTwo::AComposite::AComposite(
        const ControlInputIndex& n_delay)
    : Aaug(Eigen::SparseMatrix<bool>(n_aug_states, n_aug_states)),
      Aorig(Eigen::Matrix<double, n_states, n_states>::Zero()),
      Adelay(Eigen::Matrix<double, n_states, n_delay_states>::Zero()) {
  Eigen::Matrix<double, n_aug_states, 1> reserve_a;
  reserve_a.setConstant(1);

  Aaug.reserve(reserve_a);

  for (int i = 0; i < n_disturbance_states; i++) {
    Aaug.insert(i, i) = 1;
  }

  cumulative_delay[0] = 0;
  for (int i = 0; i < n_control_inputs; i++) {
    if (i > 0) {
      cumulative_delay[i] = cumulative_delay[i - 1];
    }
    if (n_delay[i] != 0) {
      const int size_block = n_delay[i] - 1;
      for (int j = 0; j < size_block; j++) {
        Aaug.insert(n_disturbance_states + cumulative_delay[i] + j,
                    n_disturbance_states + cumulative_delay[i] + j + 1) = 1;
      }
      Aaug.insert(n_disturbance_states + cumulative_delay[i] + size_block, 0) =
          0;
      cumulative_delay[i] += n_delay[i];
    }
  }
}

template <class System, int n_delay_states, int n_disturbance_states, int p,
          int m>
inline typename MpcController<System, n_delay_states, n_disturbance_states, p,
                              m>::AugmentedLinearizedSystemTwo::Ctype&
    MpcController<System, n_delay_states, n_disturbance_states, p,
                  m>::AugmentedLinearizedSystemTwo::Ctype::
    operator*=(const AComposite& a) {
  // const Ctype c = *this;
  // this->template leftCols<n_states>() =
  // c.template leftCols<n_states>() * a.Aorig;

  // this->template rightCols<n_aug_states>() =
  // c.template rightCols<n_aug_states>() * a.Aaug;

  // this->template rightCols<n_delay_states>() =
  // c.template leftCols<n_states>() * a.Adelay;

  this->template leftCols<n_states>() *= a.Aorig;
  this->template rightCols<n_aug_states>() *= a.Aaug;
  this->template rightCols<n_delay_states>() += this->template leftCols<n_states>() * a.Adelay;

  // for (int i = 0; i < n_control_inputs; i++) {
    // this->col(n_obs_states + a.cumulative_delay[i]) +=
        // this->template leftCols<n_states>() * a.Adelay.col(i);
  // }

  return *this;
}

// template <class System, int n_delay_states, int n_disturbance_states, int p,
// int m>
// inline Eigen::Matrix<double, System::n_outputs, System::n_control_inputs>
// MpcController<System, n_delay_states, n_disturbance_states, p,
// m>::AugmentedLinearizedSystemTwo::Ctype::
// operator*(const Eigen::Matrix<double, n_total_states, n_control_inputs,
// Eigen::RowMajor>& b) const {
// // const Eigen::Matrix<double, n_outputs, n_total_states> c = *this;
// return (*this) * b;
// }

template <class System, int n_delay_states, int n_disturbance_states, int p,
          int m>
inline typename MpcController<System, n_delay_states, n_disturbance_states, p,
                              m>::AugmentedState
    MpcController<System, n_delay_states, n_disturbance_states, p,
                  m>::AugmentedLinearizedSystemTwo::AComposite::
    operator*(const AugmentedState& x) const {
  AugmentedState x_out;
  x_out.template tail<n_aug_states>() = Aaug * x.template tail<n_aug_states>();
  x_out.template head<n_states>() = Aorig * x.template head<n_states>() +
                                    Adelay * x.template tail<n_delay_states>();

  // for (int i = 0; i < n_control_inputs; i++) {
  // x_out.template head<n_states>() +=
  // x(n_obs_states + cumulative_delay[i]) * Adelay.col(i);
  // }
  return x_out;
}

template <class System, int n_delay_states, int n_disturbance_states, int p,
          int m>
void MpcController<System, n_delay_states, n_disturbance_states, p, m>::
    AugmentedLinearizedSystemTwo::Update(
        const typename System::Linearized& sys_discrete) {
  A.Aorig = sys_discrete.A;
  int index_delay_states = 0;
  for (int i = 0; i < n_control_inputs; i++) {
    if (n_delay[i] == 0) {
      B.template block<n_states, 1>(0, i) = sys_discrete.B.col(i);
    } else {
      A.Adelay.col(index_delay_states) = sys_discrete.B.col(i);
      index_delay_states += n_delay[i];
    }
  }

  C.template leftCols<n_states>() = sys_discrete.C;
  f = sys_discrete.f;
}

/*
 * Generate prediction matrices
 */
// template <class System, int n_delay_states, int n_disturbance_states, int p,
// int m>
// inline typename MpcController<System, n_delay_states, n_disturbance_states,
// p,
// m>::Prediction
// MpcController<System, n_delay_states, n_disturbance_states, p, m>::
// AugmentedLinearizedSystemTwo::GetPrediction() const {
// }

// A.reserve(reserve_a);

// A.setZero();
// B.setZero();
// C.setZero();
// f.setZero();

// // Disturbance states
// A.template block<n_disturbance_states, n_disturbance_states>(
// n_states, n_states) = Eigen::Matrix<double, n_disturbance_states,
// n_disturbance_states>::Identity();
// C.template block<n_outputs, n_disturbance_states>(0, n_states) =
// Eigen::Matrix<double, n_outputs, n_disturbance_states>::Identity();

// // Delay states
// int index_delay_states = n_obs_states;
// for (int i = 0; i < n_control_inputs; i++) {
// if (n_delay[i] != 0) {
// const int size_block = n_delay[i] - 1;
// A.block(index_delay_states, index_delay_states + 1, size_block,
// size_block) = Eigen::MatrixXd::Identity(size_block, size_block);
// index_delay_states += n_delay[i];
// B(index_delay_states - 1, i) = 1;
// }
// }
// }
#include "mpc_controller_list.h"
