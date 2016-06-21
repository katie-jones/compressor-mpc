#include "distributed_controller.h"

/*
 * Constructor
 */
template <class AugLinSys, typename StateIndices,
          typename ObserverOutputIndices, typename ControlledOutputIndices,
          int p, int m>
DistributedController<AugLinSys, StateIndices, ObserverOutputIndices,
                      ControlledOutputIndices, p, m>::
    DistributedController(const AugLinSys& sys,
                          const InputConstraints<n_control_inputs>& constraints,
                          const ObserverMatrix& M)
    : observer_(Observer<AugLinSys>(M, Output::Zero())),
      auglinsys_(sys),
      qp_solver_(DistributedSolver<n_total_states, n_controlled_outputs,
                                   n_control_inputs, p, m>(0, constraints)) {
  static_assert(n_disturbance_states >= 0,
                "Number of disturbance states should be positive.");
  static_assert(p >= 0, "Prediction horizon should be positive.");
  static_assert(m >= 0, "Move horizon should be positive.");
}

/*
 * Initialize output, input, state
 */
template <class AugLinSys, typename StateIndices,
          typename ObserverOutputIndices, typename ControlledOutputIndices,
          int p, int m>
void DistributedController<AugLinSys, StateIndices, ObserverOutputIndices,
                           ControlledOutputIndices, p,
                           m>::Initialize(const State& x_init,
                                          const FullControlInput& u_init,
                                          const Input& full_u_old,
                                          const Output& y_init,
                                          const AugmentedState& dx_init) {
  // Set arguments
  auglinsys_.Update(x_init, full_u_old);
  u_old_ = u_init;
  x_ = x_init;
  observer_.InitializeSystem(&auglinsys_);
  observer_.SetIntialOutput(y_init);
  observer_.SetInitialAugmentedState(dx_init);

  // Initialize QP
  QP qp;

  // Get delta augmented state
  AugmentedState delta_x0;
  delta_x0 << auglinsys_.GetDerivative(),
      observer_.GetStateEstimate().template tail<n_aug_states>();

  int index_delay_states = n_obs_states;
  for (int i = 0; i < n_full_control_inputs; i++) {
    if (n_delay_[i] != 0) {
      for (int j = 0; j < n_delay_[i]; j++) {
        delta_x0(index_delay_states + j) -= u_old_[i];
      }
      index_delay_states += n_delay_[i];
    }
  }

  Prediction pred;
  pred = auglinsys_.template GeneratePrediction<ControlledOutputIndices>(
      &su_other_, p, m);

  // Take portion of output we are controlling
  ControlOutput y_controlled;
  ControlledOutputIndices::GetSubArray(y_controlled.data(), y_init.data());

  qp_solver_.GenerateDistributedQP(&qp, pred.Su, pred.Sx, pred.Sf, delta_x0,
                                   n_aug_states, y_controlled);
  qp_solver_.InitializeQPProblem(qp, u_old_.template head<n_control_inputs>());
}

/*
 * Get solution based on previous iteration's outputs
 */
template <class AugLinSys, typename StateIndices,
          typename ObserverOutputIndices, typename ControlledOutputIndices,
          int p, int m>
void DistributedController<AugLinSys, StateIndices, ObserverOutputIndices,
                           ControlledOutputIndices, p,
                           m>::GenerateInitialQP(const Output& y,
                                                 const Input& full_u_old) {
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
  for (int i = 0; i < n_full_control_inputs; i++) {
    if (n_delay_[i] != 0) {
      for (int j = 0; j < n_delay_[i]; j++) {
        delta_x0(index_delay_states + j) -= u_old_[i];
      }
      index_delay_states += n_delay_[i];
    }
  }

  // Generate QP to solve
  Prediction pred;
  auglinsys_.template GeneratePrediction<ControlledOutputIndices>(
      &pred.Su, &pred.Sx, &pred.Sf, &su_other_, p, m);

  // Take portion of output we are controlling
  ControlOutput y_controlled;
  ControlledOutputIndices::GetSubArray(y_controlled.data(), y.data());

  qp_solver_.GenerateDistributedQP(&qp_, pred.Su, pred.Sx, pred.Sf, delta_x0,
                                   n_aug_states, y_controlled);
}

#include "distributed_controller_list.h"
