#include "aug_lin_sys.h"

/*
 * Default values for linearization
 */
template <class System, int n_delay_states, int n_disturbance_states>
AugmentedLinearizedSystem<System, n_delay_states, n_disturbance_states>::
    AugmentedLinearizedSystem(const System& sys, const double sampling_time,
                              const ControlInputIndex& n_delay)
    : sys_(sys),
      sampling_time_(sampling_time),
      A(AComposite(n_delay)),
      B(BComposite(n_delay)),
      C(Eigen::MatrixXd::Zero(n_outputs, n_obs_states)),
      f(System::State::Zero()) {
  C.template rightCols<n_disturbance_states>() =
      Eigen::Matrix<double, n_outputs, n_disturbance_states>::Identity();
}

/*
 * Initialize AComposite
 */
template <class System, int n_delay_states, int n_disturbance_states>
AugmentedLinearizedSystem<System, n_delay_states, n_disturbance_states>::
    AComposite::AComposite(const ControlInputIndex& n_delay_in)
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

/*
 * Operator C *= A
 */
template <class System, int n_delay_states, int n_disturbance_states>
inline typename AugmentedLinearizedSystem<System, n_delay_states,
                                          n_disturbance_states>::Ctype&
    AugmentedLinearizedSystem<System, n_delay_states,
                              n_disturbance_states>::Ctype::
    operator*=(const AComposite& a) {
  this->template leftCols<n_states>() *= a.Aorig;
  this->template rightCols<n_aug_states>() *= a.Aaug;

  int index_delay_states = 0;
  for (int i = 0; i < n_control_inputs; i++) {
    if (a.n_delay[i] != 0) {
      this->col(n_obs_states + index_delay_states) +=
          this->template leftCols<n_states>() * a.Adelay.col(i);
      index_delay_states += a.n_delay[i];
    }
  }

  return *this;
}

/*
 * Operator output = C * B
 */
template <class System, int n_delay_states, int n_disturbance_states>
inline Eigen::Matrix<double, System::n_outputs, System::n_control_inputs>
    AugmentedLinearizedSystem<System, n_delay_states,
                              n_disturbance_states>::Ctype::
    operator*(const BComposite& b) {
  Eigen::Matrix<double, n_outputs, n_control_inputs> output =
      this->template leftCols<n_states>() * b.Borig;
  output += this->template rightCols<n_delay_states>() * b.Baug;

  return output;
}

/*
 * Operator x_out = A*x
 */
template <class System, int n_delay_states, int n_disturbance_states>
inline typename AugmentedLinearizedSystem<System, n_delay_states,
                                          n_disturbance_states>::AugmentedState
    AugmentedLinearizedSystem<System, n_delay_states,
                              n_disturbance_states>::AComposite::
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

/*
 * Update augmented system with new linearization values
 */
template <class System, int n_delay_states, int n_disturbance_states>
void AugmentedLinearizedSystem<System, n_delay_states,
                               n_disturbance_states>::Update(const State x,
                                                             const Input& u) {
  typename System::Linearized sys_continuous = sys_.GetLinearizedSystem(x, u);

  // linearize and discretize dynamic system
  typename System::Linearized sys_discrete =
      DiscretizeRK4(sys_continuous, sampling_time_);

  A.Aorig = sys_discrete.A;
  for (int i = 0; i < n_control_inputs; i++) {
    if (A.n_delay[i] == 0) {
      B.Borig.col(i) = sys_discrete.B.col(i);
    } else {
      A.Adelay.col(i) = sys_discrete.B.col(i);
    }
  }

  C.template leftCols<n_states>() = sys_discrete.C;
  f = sys_discrete.f;
}

/*
 * Initialize BComposite
 */
template <class System, int n_delay_states, int n_disturbance_states>
AugmentedLinearizedSystem<System, n_delay_states, n_disturbance_states>::
    BComposite::BComposite(const ControlInputIndex& n_delay)
    : Baug(Eigen::SparseMatrix<bool>(n_delay_states, n_control_inputs)),
      Borig(Eigen::Matrix<double, n_states, n_control_inputs>::Zero()) {
  Eigen::Matrix<double, n_control_inputs, 1> reserve_b;
  reserve_b.setConstant(1);
  Baug.reserve(reserve_b);

  int index_delay_states = 0;
  for (int i = 0; i < n_control_inputs; i++) {
    if (n_delay[i] != 0) {
      index_delay_states += n_delay[i];
      Baug.insert(index_delay_states - 1, i) = 1;
    }
  }
}

/*
 * Operator x_out = B*u
 */
template <class System, int n_delay_states, int n_disturbance_states>
typename AugmentedLinearizedSystem<System, n_delay_states,
                                   n_disturbance_states>::AugmentedState
    AugmentedLinearizedSystem<System, n_delay_states,
                              n_disturbance_states>::BComposite::
    operator*(const ControlInput& u) const {
  AugmentedState x_out;
  x_out.template tail<n_delay_states>() = Baug * u;
  x_out.template head<n_states>() = Borig * u;

  return x_out;
}

/*
 * Discretize using Runge-Kutta 4
 */
template <class System, int n_delay_states, int n_disturbance_states>
const typename System::Linearized
AugmentedLinearizedSystem<System, n_delay_states, n_disturbance_states>::
    DiscretizeRK4(const typename System::Linearized& sys_continuous,
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
template <class System, int n_delay_states, int n_disturbance_states>
const typename AugmentedLinearizedSystem<System, n_delay_states,
                                         n_disturbance_states>::Prediction
AugmentedLinearizedSystem<System, n_delay_states,
                          n_disturbance_states>::GeneratePrediction(const int p, const int m) const {
  Prediction pred;

  // pred.Sx = Eigen::MatrixXd::Zero(p * n_outputs, n_total_states);
  pred.Sx = Eigen::MatrixXd::Zero(p*n_outputs, n_aug_states);
  pred.Sf = Eigen::MatrixXd::Zero(p * n_outputs, n_states);
  pred.Su = Eigen::MatrixXd::Zero(p * n_outputs, m * n_control_inputs);

  Ctype c_times_a;
  c_times_a.template leftCols<n_obs_states>() = C;
  c_times_a.template rightCols<n_delay_states>().setZero();

  Eigen::Matrix<double, n_outputs, n_control_inputs> to_add;
  int ind_col, ind_row;

  pred.Sf.template topRows<n_outputs>() =
      C.template leftCols<n_states>();

  for (int i = 0; i < p; i++) {
    if (i > 0) {
      pred.Sf.template block<n_outputs, n_states>(i * n_outputs, 0) =
          pred.Sf.template block<n_outputs, n_states>((i - 1) * n_outputs, 0) +
          c_times_a.template leftCols<n_states>();
    }

    to_add = c_times_a * B;
    for (int j = 0; j < p - i; j++) {
      ind_row = i + j;
      if (j < m)
        ind_col = j;
      else
        ind_col = m - 1;
      pred.Su.template block<n_outputs, n_control_inputs>(
          ind_row * n_outputs, ind_col * n_control_inputs) += to_add;
    }
    c_times_a *= A;
    pred.Sx.template block<n_outputs, n_aug_states>(i * n_outputs, 0) =
        c_times_a.rightCols<n_aug_states>();
  }

  return pred;
}

#include "aug_lin_sys_list.h"
