#include "mpc_controller.h"

/*
 * Default values for linearization
 */
template <class System, int n_delay_states, int n_disturbance_states, int p,
          int m>
MpcController<System, n_delay_states, n_disturbance_states, p, m>::
    AugmentedLinearizedSystemTwo::AugmentedLinearizedSystemTwo(
        const ControlInputIndex& n_delay)
    : A(AComposite(n_delay)),
      B(Eigen::MatrixXd::Zero(n_total_states, n_control_inputs)),
      C(Eigen::MatrixXd::Zero(n_outputs, n_obs_states)),
      f(System::State::Zero()) {
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
        const ControlInputIndex& n_delay_in)
    : Aaug(Eigen::SparseMatrix<bool>(n_aug_states, n_aug_states)),
      Aorig(Eigen::Matrix<double, n_states, n_states>::Zero()),
      Adelay(Eigen::Matrix<double, n_states, n_control_inputs>::Zero()),
      n_delay(n_delay_in) {
  Eigen::Matrix<double, n_aug_states, 1> reserve_a;
  reserve_a.setConstant(1);

  Aaug.reserve(reserve_a);

  for (int i = 0; i < n_disturbance_states; i++) {
    Aaug.insert(i, i) = 1;
  }

  int index_delay_states = n_disturbance_states;
  for (int i = 0; i < n_control_inputs; i++) {
    if (n_delay[i] != 0) {
      const int size_block = n_delay[i] - 1;
      for (int j = 0; j < size_block; j++) {
        Aaug.insert(index_delay_states + j, index_delay_states + j + 1) = 1;
      }
      Aaug.insert(index_delay_states + size_block, 0) = 0;
      index_delay_states += n_delay[i];
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
  this->template leftCols<n_states>() *= a.Aorig;
  this->template rightCols<n_aug_states>() *= a.Aaug;

  int index_delay_states = 0;
  for (int i = 0; i < n_control_inputs; i++) {
    if (a.n_delay[i] != 0) {
      this->col(n_obs_states + index_delay_states) +=
          this->template leftCols<n_states>() * a.Adelay.col(i);
      // a.Adelay.col(index_delay_states);
      index_delay_states += a.n_delay[i];
    }
  }

  return *this;
}

template <class System, int n_delay_states, int n_disturbance_states, int p,
          int m>
inline typename MpcController<System, n_delay_states, n_disturbance_states, p,
                              m>::AugmentedState
    MpcController<System, n_delay_states, n_disturbance_states, p,
                  m>::AugmentedLinearizedSystemTwo::AComposite::
    operator*(const AugmentedState& x) const {
  AugmentedState x_out;
  x_out.template tail<n_aug_states>() = Aaug * x.template tail<n_aug_states>();
  x_out.template head<n_states>() = Aorig * x.template head<n_states>();

  int index_delay_states = 0;
  for (int i = 0; i < n_control_inputs; i++) {
    if (n_delay[i] != 0) {
      x_out.template head<n_states>() +=
          Adelay.col(i) * x(n_obs_states + index_delay_states);
      index_delay_states += n_delay[i];
    }
  }
  return x_out;
}

template <class System, int n_delay_states, int n_disturbance_states, int p,
          int m>
void MpcController<System, n_delay_states, n_disturbance_states, p, m>::
    AugmentedLinearizedSystemTwo::Update(
        const typename System::Linearized& sys_discrete) {
  A.Aorig = sys_discrete.A;
  for (int i = 0; i < n_control_inputs; i++) {
    if (A.n_delay[i] == 0) {
      B.template block<n_states, 1>(0, i) = sys_discrete.B.col(i);
    } else {
      A.Adelay.col(i) = sys_discrete.B.col(i);
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

#include "mpc_controller_list.h"
