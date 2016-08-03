#include "aug_lin_sys.h"
#ifndef NDEBUG
#include <iostream>
#endif

/*
 * Default values for linearization
 */
template <class System, typename Delays, int n_disturbance_states_in,
          typename ControlInputIndices, int n_sub_control_inputs_in>
AugmentedLinearizedSystem<System, Delays, n_disturbance_states_in,
                          ControlInputIndices, n_sub_control_inputs_in>::
    AugmentedLinearizedSystem(const System& sys, const double sampling_time)
    : sys_(sys),
      sampling_time_(sampling_time),
      A(AComposite()),
      B(BComposite()),
      C(Eigen::MatrixXd::Zero(n_outputs, n_obs_states)),
      f(System::State::Zero()) {
  C.template rightCols<n_disturbance_states>() =
      Eigen::Matrix<double, n_outputs, n_disturbance_states>::Identity();
}

/*
 * Initialize AComposite
 */
template <class System, typename Delays, int n_disturbance_states_in,
          typename ControlInputIndices, int n_sub_control_inputs_in>
AugmentedLinearizedSystem<System, Delays, n_disturbance_states_in,
                          ControlInputIndices,
                          n_sub_control_inputs_in>::AComposite::AComposite()
    : Aaug(Eigen::SparseMatrix<bool>(n_aug_states, n_aug_states)),
      Aorig(Eigen::Matrix<double, n_states, n_states>::Zero()),
      Adelay(Eigen::Matrix<double, n_states, n_delayed_inputs>::Zero()) {
  Eigen::Matrix<double, n_aug_states, 1> reserve_a;
  reserve_a.setConstant(1);

  for (int i = 0; i < n_delayed_inputs; i++) {
    reserve_a[n_disturbance_states + i] = 0;
  }

  Aaug.reserve(reserve_a);

  // Make identity matrix for disturbance states
  for (int i = 0; i < n_disturbance_states; i++) {
    Aaug.insert(i, i) = 1;
  }

  int index_delay_states = n_disturbance_states + n_delayed_inputs;
  int index_delayed_inputs = n_disturbance_states;
  for (int i = 0; i < n_control_inputs; i++) {
    if (n_delay_[i] != 0) {
      const int size_block = n_delay_[i] - 1;
      Aaug.insert(index_delayed_inputs, index_delay_states) = 1;
      for (int j = 1; j < size_block; j++) {
        Aaug.insert(index_delay_states + j - 1, index_delay_states + j) = 1;
      }
      index_delay_states += n_delay_[i] - 1;
      index_delayed_inputs++;
    }
  }
}

/*
 * Operator C *= A
 */
template <class System, typename Delays, int n_disturbance_states_in,
          typename ControlInputIndices, int n_sub_control_inputs_in>
template <int n_sub_outputs>
void AugmentedLinearizedSystem<System, Delays, n_disturbance_states_in,
                               ControlInputIndices, n_sub_control_inputs_in>::
    AComposite::MultiplyC(
        Eigen::Matrix<double, n_sub_outputs, n_total_states>* C) const {
  const Eigen::Matrix<double, n_sub_outputs, n_states> temp =
      C->template leftCols<n_states>();

  C->template leftCols<n_states>() *= Aorig;
  C->template rightCols<n_aug_states>() *= Aaug;
  C->template block<n_sub_outputs, n_delayed_inputs>(0, n_obs_states) +=
      temp * Adelay;
}

/*
 * Operator output = C * B
 */
template <class System, typename Delays, int n_disturbance_states_in,
          typename ControlInputIndices, int n_sub_control_inputs_in>
template <int n_sub_outputs>
Eigen::Matrix<double, n_sub_outputs, System::n_control_inputs>
AugmentedLinearizedSystem<System, Delays, n_disturbance_states_in,
                          ControlInputIndices, n_sub_control_inputs_in>::
    BComposite::MultiplyC(
        const Eigen::Matrix<double, n_sub_outputs, n_total_states>& C) const {
  Eigen::Matrix<double, n_sub_outputs, n_control_inputs> output =
      C.template leftCols<n_states>() * Borig;
  output += C.template rightCols<n_delay_states>() * Baug;

  return output;
}

/*
 * Multiply A by augmented part of x
 */
template <class System, typename Delays, int n_disturbance_states_in,
          typename ControlInputIndices, int n_sub_control_inputs_in>
typename AugmentedLinearizedSystem<System, Delays, n_disturbance_states_in,
                                   ControlInputIndices,
                                   n_sub_control_inputs_in>::AugmentedState
AugmentedLinearizedSystem<System, Delays, n_disturbance_states_in,
                          ControlInputIndices, n_sub_control_inputs_in>::
    AComposite::TimesAugmentedOnly(
        const Eigen::Matrix<double, n_aug_states, 1> x) const {
  // Augmented part
  AugmentedState x_out = AugmentedState::Zero();
  x_out.template tail<n_aug_states>() = Aaug * x.template tail<n_aug_states>();
  x_out.template head<n_states>() +=
      Adelay * x.template segment<n_delayed_inputs>(n_disturbance_states);

  return x_out;
}

/*
 * Update augmented system with new linearization values
 */
template <class System, typename Delays, int n_disturbance_states_in,
          typename ControlInputIndices, int n_sub_control_inputs_in>
void AugmentedLinearizedSystem<
    System, Delays, n_disturbance_states_in, ControlInputIndices,
    n_sub_control_inputs_in>::Update(const State x, const Input& u) {
  typename System::Linearized sys_continuous = sys_.GetLinearizedSystem(x, u);

  // linearize and discretize dynamic system
  typename System::Linearized sys_discrete =
      DiscretizeRK4(sys_continuous, sampling_time_);

  A.Aorig = sys_discrete.A;
  int index;
  int index_delayed_inputs = 0;
  for (int i = 0; i < n_control_inputs; i++) {
    // if system has re-ordered control inputs, re-order columns of A and B
    if (is_reduced)
      index = ControlInputIndices::GetEntry(i);
    else
      index = i;
    if (n_delay_[i] == 0) {
      B.Borig.col(index) = sys_discrete.B.col(i);
    } else {
      A.Adelay.col(index_delayed_inputs) = sys_discrete.B.col(index);
      index_delayed_inputs++;
    }
  }

  C.template leftCols<n_states>() = sys_discrete.C;
  f = sys_discrete.f;
}

/*
 * Initialize BComposite
 */
template <class System, typename Delays, int n_disturbance_states_in,
          typename ControlInputIndices, int n_sub_control_inputs_in>
AugmentedLinearizedSystem<System, Delays, n_disturbance_states_in,
                          ControlInputIndices,
                          n_sub_control_inputs_in>::BComposite::BComposite()
    : Baug(Eigen::SparseMatrix<bool>(n_delay_states, n_control_inputs)),
      Borig(Eigen::Matrix<double, n_states, n_control_inputs>::Zero()) {
  Baug.reserve(n_delayed_inputs);

  int index_delay_states = n_delayed_inputs;
  for (int i = 0; i < n_control_inputs; i++) {
    if (n_delay_[i] != 0) {
      index_delay_states += n_delay_[i] - 1;
      Baug.insert(index_delay_states - 1, i) = 1;
    }
  }
}

/*
 * Operator x_out = B*u
 */
template <class System, typename Delays, int n_disturbance_states_in,
          typename ControlInputIndices, int n_sub_control_inputs_in>
typename AugmentedLinearizedSystem<System, Delays, n_disturbance_states_in,
                                   ControlInputIndices,
                                   n_sub_control_inputs_in>::AugmentedState
    AugmentedLinearizedSystem<System, Delays, n_disturbance_states_in,
                              ControlInputIndices,
                              n_sub_control_inputs_in>::BComposite::
    operator*(const ControlInput& u) const {
  AugmentedState x_out = AugmentedState::Zero();
  x_out.template head<n_states>() = Borig * u;
  x_out.template segment<n_disturbance_states>(n_states).setZero();
  x_out.template tail<n_delay_states>() = Baug * u;

  return x_out;
}

/*
 * Discretize using Runge-Kutta 4
 */
template <class System, typename Delays, int n_disturbance_states_in,
          typename ControlInputIndices, int n_sub_control_inputs_in>
const typename System::Linearized AugmentedLinearizedSystem<
    System, Delays, n_disturbance_states_in, ControlInputIndices,
    n_sub_control_inputs_in>::DiscretizeRK4(const typename System::Linearized&
                                                sys_continuous,
                                            const double Ts) {
  typename System::Linearized sys_discrete;
  Eigen::Matrix<double, n_states, n_states> A2, A3, Acommon;

  A2 = sys_continuous.A * sys_continuous.A;
  A3 = A2 * sys_continuous.A;
  Acommon = Ts * Eigen::Matrix<double, n_states, n_states>::Identity() +
            Ts * Ts / 2.0 * sys_continuous.A + Ts * Ts * Ts / 6.0 * A2 +
            Ts * Ts * Ts * Ts / 24.0 * A3;

  sys_discrete.A = Eigen::Matrix<double, n_states, n_states>::Identity() +
                   Acommon * sys_continuous.A;
  sys_discrete.B = Acommon * sys_continuous.B;
  sys_discrete.C = sys_continuous.C;
  sys_discrete.f = Acommon * sys_continuous.f;

  return sys_discrete;
}

/*
 * Generate linearized prediction matrices
 */
template <class System, typename Delays, int n_disturbance_states_in,
          typename ControlInputIndices, int n_sub_control_inputs_in>
template <typename ControlledOutputIndices>
void AugmentedLinearizedSystem<
    System, Delays, n_disturbance_states_in, ControlInputIndices,
    n_sub_control_inputs_in>::GeneratePrediction(Eigen::MatrixXd* Su,
                                                 Eigen::MatrixXd* Sx,
                                                 Eigen::MatrixXd* Sf,
                                                 Eigen::MatrixXd* Su_other,
                                                 const int p,
                                                 const int m) const {
  // number of outputs being controlled
  constexpr int n_controlled_outputs = ControlledOutputIndices::size;

  // Initialize c_times_a to C (with only desired outputs)
  Eigen::Matrix<double, n_controlled_outputs, n_total_states> c_times_a;
  c_times_a.template rightCols<n_delay_states>().setZero();

  for (int i = 0; i < ControlledOutputIndices::size; i++) {
    c_times_a.template block<1, n_obs_states>(i, 0) =
        C.row(ControlledOutputIndices::GetEntry(i));
  }

  // Initialize prediction to zero
  Su->resize(p * n_controlled_outputs, m * n_sub_control_inputs);
  Sx->resize(p * n_controlled_outputs, n_aug_states);
  Sf->resize(p * n_controlled_outputs, n_states);
  if (is_reduced)
    Su_other->resize(p * n_controlled_outputs, m * n_other_control_inputs);

  Su->setZero();
  Sx->setZero();
  Sf->setZero();
  if (is_reduced) Su_other->setZero();

  Eigen::Matrix<double, n_controlled_outputs, n_control_inputs> to_add;
  int ind_col, ind_row;

  Sf->template topRows<n_controlled_outputs>() =
      c_times_a.template leftCols<n_states>();

  for (int i = 0; i < p; i++) {
    if (i > 0) {
      // Sf is additive
      Sf->template block<n_controlled_outputs, n_states>(
          i * n_controlled_outputs, 0) =
          Sf->template block<n_controlled_outputs, n_states>(
              (i - 1) * n_controlled_outputs, 0) +
          c_times_a.template leftCols<n_states>();
    }

    to_add = B.MultiplyC(c_times_a);
    for (int j = 0; j < p - i; j++) {
      ind_row = i + j;
      if (j < m)
        ind_col = j;
      else
        ind_col = m - 1;

      // Split Su matrix: first own inputs, then other controllers'
      Su->template block<n_controlled_outputs, n_sub_control_inputs>(
          ind_row * n_controlled_outputs, ind_col * n_sub_control_inputs) +=
          to_add.template leftCols<n_sub_control_inputs>();

      if (is_reduced)
        Su_other->template block<n_controlled_outputs, n_other_control_inputs>(
            ind_row * n_controlled_outputs, ind_col * n_other_control_inputs) +=
            to_add.template rightCols<n_other_control_inputs>();
    }
    A.MultiplyC(&c_times_a);
    Sx->template block<n_controlled_outputs, n_aug_states>(
        i * n_controlled_outputs, 0) =
        c_times_a.template rightCols<n_aug_states>();
  }
}

#include "aug_lin_sys_list.h"
