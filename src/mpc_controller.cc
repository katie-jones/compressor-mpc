#include "mpc_controller.h"
#include "print_matrix.h"
#include <iostream>

/*
 * Observe system given new output values from plant
 */
template <class System, int n_delay_states, int n_disturbance_states, int p,
          int m>
void MpcController<System, n_delay_states, n_disturbance_states, p,
                   m>::ObserveAPosteriori(const Output& y_in) {
  // apply observer to non-delay states
  dx_aug_.template head<n_obs_states>() =
      dx_aug_.template head<n_obs_states>() +
      M_ *
          (y_in - y_old_ -
           (Output)(auglinsys_.C.template leftCols<n_obs_states>() *
                    (dx_aug_.template head<n_obs_states>()))).matrix();

  // add change in normal states
  x_aug_.template head<n_states>() =
      x_aug_.template head<n_states>() + dx_aug_.template head<n_states>();

  // take augmented states (delay/disturbance) directly from dx_aug_
  x_aug_.template tail<n_aug_states>() = dx_aug_.template tail<n_aug_states>();

  // store input value
  y_old_ = y_in;
}

/*
 * Discretize linearized system using runge-kutta 4
 */
template <class System, int n_delay_states, int n_disturbance_states, int p,
          int m>
const typename System::Linearized MpcController<
    System, n_delay_states, n_disturbance_states, p,
    m>::DiscretizeRK4(const typename System::Linearized& sys_continuous,
                      const double Ts) {
  typename System::Linearized sys_discrete;
  Eigen::Matrix<double, n_states, n_states> A2, A3, Acommon;

  A2 = sys_continuous.A * sys_continuous.A;
  A3 = A2 * sys_continuous.A;
  Acommon = Ts * Eigen::Matrix<double, n_states, n_states>::Identity() +
            Ts * Ts / 2.0 * sys_continuous.A + Ts * Ts * Ts / 6.0 * A2 +
            Ts * Ts * Ts / 24.0 * A3;

  sys_discrete.A = Eigen::Matrix<double, n_states, n_states>::Identity() +
                   Acommon * sys_continuous.A;
  sys_discrete.B = Acommon * sys_continuous.B;
  sys_discrete.C = sys_continuous.C;
  sys_discrete.f = Acommon * sys_continuous.f;
  return sys_discrete;
}

/*
 * Linearize and augment the system about current state estimate
 */
template <class System, int n_delay_states, int n_disturbance_states, int p,
          int m>
void MpcController<System, n_delay_states, n_disturbance_states, p,
                   m>::LinearizeAndAugment() {
  // linearize and discretize dynamic system
  typename System::Linearized sys_discrete =
      DiscretizeRK4(sys_.GetLinearizedSystem(x_aug_.template head<n_states>(),
                                             GetPlantInput(u_old_)),
                    sampling_time_);

  AugmentedLinearizedSystem sys_out;
  sys_out.A.setZero();
  sys_out.B.setZero();
  sys_out.C.setZero();
  sys_out.f.setZero();

  sys_out.A.template topLeftCorner<n_states, n_states>() = sys_discrete.A;
  sys_out.A.template block<n_disturbance_states, n_disturbance_states>(
      n_states, n_states) = Eigen::Matrix<double, n_disturbance_states,
                                          n_disturbance_states>::Identity();

  int index_delay_states = n_obs_states;

  for (int i = 0; i < n_control_inputs; i++) {
    if (n_delay_[i] == 0) {
      sys_out.B.template block<n_states, 1>(0, i) = sys_discrete.B.col(i);
    } else {
      sys_out.A.template block<n_states, 1>(0, index_delay_states) =
          sys_discrete.B.col(i);
      const int size_block = n_delay_[i] - 1;
      sys_out.A.template block(index_delay_states, index_delay_states + 1,
                               size_block, size_block) =
          Eigen::MatrixXd::Identity(size_block, size_block);
      index_delay_states += n_delay_[i];
      sys_out.B(index_delay_states - 1, i) = 1;
    }
  }

  sys_out.C.template leftCols<n_states>() = sys_discrete.C;
  sys_out.C.template block<n_outputs, n_disturbance_states>(0, n_states) =
      Eigen::Matrix<double, n_outputs, n_disturbance_states>::Identity();

  sys_out.f.template head<n_states>() = sys_discrete.f;

  auglinsys_ = sys_out;
  print_matrix(std::cout, auglinsys_.A, "A");
  print_matrix(std::cout, auglinsys_.B, "B");
  print_matrix(std::cout, auglinsys_.C, "C");
  print_matrix(std::cout, auglinsys_.f, "f");
}

/*
 * Generate QP matrices from MPC formulation
 */
template <class System, int n_delay_states, int n_disturbance_states, int p,
          int m>
const typename MpcController<System, n_delay_states, n_disturbance_states, p,
                             m>::QP
MpcController<System, n_delay_states, n_disturbance_states, p, m>::GenerateQP()
    const {
  QP qp;
  const Prediction pred = GeneratePrediction();
  const OutputPrediction dy_ref = y_ref_ - y_old_.template replicate<p, 1>();

  AugmentedState delta_x0 = dx_aug_;
  static_cast<State>(delta_x0.template head<n_states>()).setZero();

  int index_delay_states = n_obs_states;
  for (int i = 0; i < n_control_inputs; i++) {
    if (n_delay_[i] != 0) {
      for (int j = 0; j < n_delay_[i]; j++) {
        delta_x0(index_delay_states + j) -= u_old_[i];
      }
      index_delay_states += n_delay_[i];
    }
  }

  qp.H = pred.Su.transpose() * y_weight_ * pred.Su;

  qp.f = auglinsys_.f.template head<n_states>().transpose() *
         pred.Sf.transpose() * y_weight_ * pred.Su;
  qp.f -= dy_ref.transpose() * y_weight_ * pred.Su;
  qp.f += delta_x0.transpose() * pred.Sx.transpose() * y_weight_ * pred.Su;

  return qp;
}

/*
 * Generate linearized prediction matrices
 */
template <class System, int n_delay_states, int n_disturbance_states, int p,
          int m>
const typename MpcController<System, n_delay_states, n_disturbance_states, p,
                             m>::Prediction
MpcController<System, n_delay_states, n_disturbance_states, p,
              m>::GeneratePrediction() const {
  Prediction pred;
  Eigen::Matrix<double, n_outputs, n_total_states> c_times_a = auglinsys_.C;
  Eigen::Matrix<double, n_outputs, n_control_inputs> to_add;
  int ind_col, ind_row;

  pred.Su.setZero();
  pred.Sx.setZero();
  pred.Sf.setZero();

  pred.Sf.template topRows<n_outputs>() =
      auglinsys_.C.template leftCols<n_states>();

  for (int i = 0; i < p; i++) {
    if (i > 0) {
      pred.Sf.template block<n_outputs, n_states>(i * n_outputs, 0) =
          pred.Sf.template block<n_outputs, n_states>((i - 1) * n_outputs, 0) +
          c_times_a.template leftCols<n_states>();
    }

    to_add = c_times_a * auglinsys_.B;
    for (int j = 0; j < p - i; j++) {
      ind_row = i + j;
      if (j < m)
        ind_col = j;
      else
        ind_col = m - 1;
      pred.Su.template block<n_outputs, n_control_inputs>(
          ind_row * n_outputs, ind_col * n_control_inputs) += to_add;
    }
    c_times_a *= auglinsys_.A;
    pred.Sx.template block<n_outputs, n_total_states>(i * n_outputs, 0) =
        c_times_a;
  }

  print_matrix(std::cout, pred.Su, "Su");
  print_matrix(std::cout, pred.Sf, "Sf");
  print_matrix(std::cout, pred.Sx, "Sx");

  return pred;
}

/*
 * Solve QP and obtain next input to apply
 */
template <class System, int n_delay_states, int n_disturbance_states, int p,
          int m>
const typename MpcController<System, n_delay_states, n_disturbance_states, p,
                             m>::ControlInput
MpcController<System, n_delay_states, n_disturbance_states, p, m>::SolveQP(
    const QP& qp) const {
  return ControlInput();
}

/*
 * Predict next state of system based on input value to be applied
 */
template <class System, int n_delay_states, int n_disturbance_states, int p,
          int m>
void MpcController<System, n_delay_states, n_disturbance_states, p,
                   m>::ObserveAPriori(const ControlInput& u_in) {}

/*
 * Calculate a plant input based on control input and offset
 */
template <class System, int n_delay_states, int n_disturbance_states, int p,
          int m>
const typename MpcController<System, n_delay_states, n_disturbance_states, p,
                             m>::Input
MpcController<System, n_delay_states, n_disturbance_states, p,
              m>::GetPlantInput(const ControlInput& u_control) const {
  Input u = u_offset_;
  for (int i = 0; i < n_control_inputs; i++) {
    u(control_input_index_[i]) += u_control(i);
  }
  return u;
}

/*
 * Calculate a control input based on plant input and offset
 */
template <class System, int n_delay_states, int n_disturbance_states, int p,
          int m>
const typename MpcController<System, n_delay_states, n_disturbance_states, p,
                             m>::ControlInput
MpcController<System, n_delay_states, n_disturbance_states, p,
              m>::GetControlInput(const Input& u) const {
  ControlInput u_control;
  for (int i = 0; i < n_control_inputs; i++) {
    u_control(i) =
        u(control_input_index_[i]) - u_offset_(control_input_index_[i]);
  }
  return u_control;
}

#include "mpc_controller_list.h"
