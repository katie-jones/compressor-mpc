#include "distributed_controller.h"

/*
 * Constructor
 */
template <class System, int n_control_inputs, int n_delay_states, int n_disturbance_states, int p, int m, int n_controllers>
DistributedController<System, n_control_inputs, n_delay_states,
                      n_disturbance_states, p, m, n_controllers>::
    DistributedController(const System& sys, const double Ts,
                          const ControlInputIndex& n_delay,
                          const InputConstraints<n_control_inputs>& constraints,
                          const ObserverMatrix& M)
    : observer_(Observer<System, n_delay_states, n_disturbance_states>(
          M, Output::Zero())),
      auglinsys_(
          AugmentedLinearizedSystem<System, n_delay_states,
                                    n_disturbance_states>(sys, Ts, n_delay)),
      qp_solver_(DistributedSolver<n_total_states, n_outputs, n_control_inputs,
                                   p, m, n_controllers>(0, constraints)) {
  static_assert(n_delay_states >= 0,
                "Number of delay states should be positive.");
  static_assert(n_disturbance_states >= 0,
                "Number of disturbance states should be positive.");
  static_assert(n_controllers >= 0,
                "Number of controllers should be positive.");
  static_assert(p >= 0, "Prediction horizon should be positive.");
  static_assert(m >= 0, "Move horizon should be positive.");

  for (int i = 0; i < n_control_inputs; i++) {
    n_delay_[i] = n_delay[i];
  }
}

/*
 * Initialize output, input, state
 */
template <class System, int n_control_inputs, int n_delay_states, int n_disturbance_states, int p, int m, int n_controllers>
void DistributedController<
    System, n_control_inputs, n_delay_states, n_disturbance_states, p, m,
    n_controllers>::Initialize(const Output& y_init, const ControlInput& u_init,
                               const Input& full_u_old, const State& x_init,
                               const AugmentedState& dx_init) {
  // Set arguments
  auglinsys_.Update(x_init, full_u_old);
  u_old_ = u_init;
  x_ = x_init;
  observer_.SetIntialOutput(y_init);
  observer_.SetInitialAugmentedState(dx_init);

  // Initialize QP
  QP qp;

  // Get delta augmented state
  AugmentedState delta_x0;
  delta_x0 << auglinsys_.GetDerivative(),
      observer_.GetStateEstimate().template tail<n_aug_states>();

  int index_delay_states = n_obs_states;
  for (int i = 0; i < n_control_inputs; i++) {
    if (n_delay_[i] != 0) {
      for (int j = 0; j < n_delay_[i]; j++) {
        delta_x0(index_delay_states + j) -= u_old_[i];
      }
      index_delay_states += n_delay_[i];
    }
  }

  Prediction pred;
  pred = auglinsys_.GeneratePrediction(p, m);

  qp_solver_.GenerateDistributedQP(qp, pred.Su, pred.Sx, pred.Sf, delta_x0,
                                   n_aug_states, y_init);
  qp_solver_.InitializeQPProblem(qp, u_old_);
}

/*
 * Get solution based on previous iteration's outputs
 */
template <class System, int n_control_inputs, int n_delay_states, int n_disturbance_states, int p, int m, int n_controllers>
auto DistributedController<
    System, n_control_inputs, n_delay_states, n_disturbance_states, p, m,
    n_controllers>::GenerateInitialQP(const Output& y, const Input& full_u_old)
    -> QP {
  // Observe and update linearization
  x_ += observer_.ObserveAPosteriori(y);
  auglinsys_.Update(x_, full_u_old);

  // Get current derivative (first n_states entries) and augmented state (last
  // n_aug_states entries)
  AugmentedState delta_x0;
  delta_x0 << auglinsys_.GetDerivative(),
      observer_.GetStateEstimate().template tail<n_aug_states>();

  // Adjust delayed inputs to account for linearization before applying
  int index_delay_states = n_obs_states;
  for (int i = 0; i < n_control_inputs; i++) {
    if (n_delay_[i] != 0) {
      for (int j = 0; j < n_delay_[i]; j++) {
        delta_x0(index_delay_states + j) -= u_old_[i];
      }
      index_delay_states += n_delay_[i];
    }
  }

  // Generate QP to solve
  QP qp;
  Prediction pred;
  pred = auglinsys_.GeneratePrediction(p, m);
  qp_solver_.GenerateDistributedQP(qp, pred.Su, pred.Sx, pred.Sf, delta_x0,
                                   n_aug_states, observer_.GetPreviousOutput());
  return qp;
}

#include "distributed_controller_list.h"
