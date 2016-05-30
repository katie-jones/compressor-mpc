#include "mpc_controller.h"

/*
 * Observe system given new output values from plant
 */
template <class AugmentedLinearizedSystem, int p, int m>
void MpcController<AugmentedLinearizedSystem, p, m>::ObserveAPosteriori(
    const Output& y_in) {
  // apply observer to non-delay states
  dx_aug_.template head<n_obs_states>() =
      dx_aug_.template head<n_obs_states>() +
      M_ *
          (y_in - y_old_ -
           (Output)(auglinsys_.C * (dx_aug_.template head<n_obs_states>())))
              .matrix();

  // add change in normal states
  x_ = x_ + dx_aug_.template head<n_states>();

  // store input value
  y_old_ = y_in;
}

/*
 * Generate QP matrices from MPC formulation
 */
template <class AugmentedLinearizedSystem, int p, int m>
const typename MpcController<AugmentedLinearizedSystem, p, m>::QP
MpcController<AugmentedLinearizedSystem, p, m>::GenerateQP() const {
  QP qp;
  const Prediction pred = GeneratePrediction();
  const OutputPrediction dy_ref = y_ref_ - y_old_.template replicate<p, 1>();

  AugmentedState delta_x0 = dx_aug_;
  delta_x0.template head<n_states>().setZero();

  int index_delay_states = n_obs_states;
  for (int i = 0; i < n_control_inputs; i++) {
    if (n_delay_[i] != 0) {
      for (int j = 0; j < n_delay_[i]; j++) {
        delta_x0(index_delay_states + j) -= u_old_[i];
      }
      index_delay_states += n_delay_[i];
    }
  }

  Eigen::SparseMatrix<double> y_weight_full(p * n_outputs, p * n_outputs);
  const Eigen::Matrix<double, p * n_outputs, 1> reserve_values =
      Eigen::Matrix<double, p * n_outputs, 1>::Constant(n_outputs);
  y_weight_full.reserve(reserve_values);
  for (int i = 0; i < p; i++) {
    for (int j = 0; j < n_outputs * n_outputs; j++) {
      y_weight_full.insert(i * n_outputs + j % n_outputs,
                           i * n_outputs + j / n_outputs) =
          y_weight_(j % n_outputs, j / n_outputs);
    }
  }

  Eigen::Matrix<double, m * n_control_inputs, m * n_control_inputs>
      u_weight_full;
  u_weight_full.setZero();
  for (int i = 0; i < m; i++) {
    u_weight_full.template block<n_control_inputs, n_control_inputs>(
        i * n_control_inputs, i * n_control_inputs) = u_weight_;
  }

  qp.H = pred.Su.transpose() * y_weight_full * pred.Su + u_weight_full;

  qp.f =
      auglinsys_.f.transpose() * pred.Sf.transpose() * y_weight_full * pred.Su -
      dy_ref.transpose() * y_weight_full * pred.Su +
      delta_x0.transpose() * pred.Sx.transpose() * y_weight_full * pred.Su;

  return qp;
}

/*
 * Generate linearized prediction matrices
 */
template <class AugmentedLinearizedSystem, int p, int m>
const typename MpcController<AugmentedLinearizedSystem, p, m>::Prediction
MpcController<AugmentedLinearizedSystem, p, m>::GeneratePrediction() const {
  Prediction pred;

  pred.Sx = Eigen::MatrixXd::Zero(p * n_outputs, n_total_states);
  pred.Sf = Eigen::MatrixXd::Zero(p * n_outputs, n_states);
  pred.Su = Eigen::MatrixXd::Zero(p * n_outputs, m * n_control_inputs);

  typename AugmentedLinearizedSystem::Ctype c_times_a;
  c_times_a.template leftCols<n_obs_states>() = auglinsys_.C;
  c_times_a.template rightCols<n_delay_states>().setZero();

  Eigen::Matrix<double, n_outputs, n_control_inputs> to_add;
  int ind_col, ind_row;

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

  return pred;
}

/*
 * Solve QP and obtain next input to apply
 */
template <class AugmentedLinearizedSystem, int p, int m>
const typename MpcController<AugmentedLinearizedSystem, p, m>::ControlInput
MpcController<AugmentedLinearizedSystem, p, m>::SolveQP(const QP& qp) {
  int n_wsr = n_wsr_max;

  // Replicate constraints for number of moves
  const Eigen::Matrix<double, m * n_control_inputs, 1> lb =
      (u_constraints_.lower_bound - u_old_).template replicate<m, 1>();
  const Eigen::Matrix<double, m * n_control_inputs, 1> ub =
      (u_constraints_.upper_bound - u_old_).template replicate<m, 1>();
  const Eigen::Matrix<double, m * n_control_inputs, 1> lbA =
      u_constraints_.lower_rate_bound.template replicate<m, 1>();
  const Eigen::Matrix<double, m * n_control_inputs, 1> ubA =
      u_constraints_.upper_rate_bound.template replicate<m, 1>();

  qpOASES::returnValue status =
      qp_problem_.hotstart(qp.H.data(), qp.f.data(), Ain_.data(), lb.data(),
                           ub.data(), lbA.data(), ubA.data(), n_wsr, NULL);

  if (status != qpOASES::SUCCESSFUL_RETURN) {
    // QP not solved, return zeros
    return ControlInput::Zero();
  }

  ControlInput u_solution;
  qp_problem_.getPrimalSolution(u_solution.data());
  return u_solution;
}

/*
 * Predict next state of system based on input value to be applied
 */
template <class AugmentedLinearizedSystem, int p, int m>
void MpcController<AugmentedLinearizedSystem, p, m>::ObserveAPriori(
    const ControlInput& du_in) {
  ControlInput du;
  AugmentedState dx = dx_aug_;

  dx.template head<n_states>().setZero();
  int index_delay_states = n_obs_states;

  for (int i = 0; i < n_control_inputs; i++) {
    if (n_delay_[i] == 0) {
      du(i) = du_in(i);
    } else {
      du(i) = u_old_(i) + du_in(i);
      dx(index_delay_states) -= u_old_(i);
      index_delay_states += n_delay_[i];
    }
  }

  dx_aug_ = auglinsys_.B * du + auglinsys_.A * dx;
  dx_aug_.template head<n_states>() += auglinsys_.f;
}

/*
 * Calculate a plant input based on control input and offset
 */
template <class AugmentedLinearizedSystem, int p, int m>
const typename MpcController<AugmentedLinearizedSystem, p, m>::Input
MpcController<AugmentedLinearizedSystem, p, m>::GetPlantInput(
    const ControlInput& u_control) const {
  Input u = u_offset_;
  for (int i = 0; i < n_control_inputs; i++) {
    u(control_input_index_[i]) += u_control(i);
  }
  return u;
}

/*
 * Calculate a control input based on plant input and offset
 */
template <class AugmentedLinearizedSystem, int p, int m>
const typename MpcController<AugmentedLinearizedSystem, p, m>::ControlInput
MpcController<AugmentedLinearizedSystem, p, m>::GetControlInput(
    const Input& u) const {
  ControlInput u_control;
  for (int i = 0; i < n_control_inputs; i++) {
    u_control(i) =
        u(control_input_index_[i]) - u_offset_(control_input_index_[i]);
  }
  return u_control;
}

/*
 * Solve QP and return optimal input to apply
 */
template <class AugmentedLinearizedSystem, int p, int m>
const typename MpcController<AugmentedLinearizedSystem, p, m>::ControlInput
MpcController<AugmentedLinearizedSystem, p, m>::GetNextInput(const Output& y) {
  ObserveAPosteriori(y);
  auglinsys_.Update(x_, GetPlantInput(u_old_));
  const QP qp = GenerateQP();
  const ControlInput usol = SolveQP(qp);
  ObserveAPriori(usol);

  u_old_ += usol;
  return u_old_;
}

/*
 * Initialize initial states and QP
 */
template <class AugmentedLinearizedSystem, int p, int m>
void MpcController<AugmentedLinearizedSystem, p, m>::SetInitialState(
    const State& x_init, const ControlInput& u_init, const Output& y_init,
    const AugmentedState& dx_init) {
  x_ = x_init;
  u_old_ = u_init;
  // y_old_ = sys_.GetOutput(x_init);
  y_old_ = y_init;
  dx_aug_ = dx_init;

  auglinsys_.Update(x_, GetPlantInput(u_old_));
  const QP qp = GenerateQP();

  // Replicate constraints for number of moves
  const Eigen::Matrix<double, m * n_control_inputs, 1> lb =
      (u_constraints_.lower_bound - u_old_).template replicate<m, 1>();
  const Eigen::Matrix<double, m * n_control_inputs, 1> ub =
      (u_constraints_.upper_bound - u_old_).template replicate<m, 1>();
  const Eigen::Matrix<double, m * n_control_inputs, 1> lbA =
      u_constraints_.lower_rate_bound.template replicate<m, 1>();
  const Eigen::Matrix<double, m * n_control_inputs, 1> ubA =
      u_constraints_.upper_rate_bound.template replicate<m, 1>();
  int n_wsr = n_wsr_max;

  qp_problem_.setPrintLevel(qpOASES::PrintLevel::PL_LOW);  // only print errors
  qp_problem_.init(qp.H.data(), qp.f.data(), Ain_.data(), lb.data(), ub.data(),
                   lbA.data(), ubA.data(), n_wsr, NULL);
}

/*
 * Default values for input constraints
 */
template <class AugmentedLinearizedSystem, int p, int m>
MpcController<AugmentedLinearizedSystem, p,
              m>::InputConstraints::InputConstraints()
    : upper_bound(ControlInput::Constant(std::nan(""))),
      lower_bound(ControlInput::Constant(-std::nan(""))),
      use_rate_constraints(false) {}

#include "mpc_controller_list.h"
