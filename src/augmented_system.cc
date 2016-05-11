#include "augmented_system.h"
#include <iostream>

template <class System, int n_disturbance_states, int n_delay_states>
typename AugmentedSystem<System, n_disturbance_states,
                         n_delay_states>::LinearizedSystem
AugmentedSystem<System, n_disturbance_states, n_delay_states>::DiscretizeRK4(
    LinearizedSystem &sys_continuous) {
  LinearizedSystem sys_discrete;
  Eigen::Matrix<double, System::n_states, System::n_states> A2, A3, Acommon;

  A2 = sys_continuous.A * sys_continuous.A;
  A3 = A2 * sys_continuous.A;
  Acommon = Ts_ * Eigen::Matrix<double, System::n_states,
                                System::n_states>::Identity() +
            Ts_ * Ts_ / 2.0 * sys_continuous.A + Ts_ * Ts_ * Ts_ / 6.0 * A2 +
            Ts_ * Ts_ * Ts_ / 24.0 * A3;

  sys_discrete.A =
      Eigen::Matrix<double, System::n_states, System::n_states>::Identity() +
      Acommon * sys_continuous.A;
  sys_discrete.B = Acommon * sys_continuous.B;
  sys_discrete.C = sys_continuous.C;
  sys_discrete.f = Acommon * sys_continuous.f;
  return sys_discrete;
}

template <class System, int n_disturbance_states, int n_delay_states>
typename AugmentedSystem<System, n_disturbance_states,
                         n_delay_states>::AugmentedLinearizedSystem
AugmentedSystem<System, n_disturbance_states,
                n_delay_states>::LinearizeAndAugment(LinearizedSystem &
                                                         sys_continuous) {
  AugmentedLinearizedSystem sys_out;
  sys_out.A.setZero();

  LinearizedSystem sys_discrete = DiscretizeRK4(sys_continuous);

  sys_out.A.template block<System::n_states, System::n_states>(0, 0) =
      sys_discrete.A;
  sys_out.A.template block<n_disturbance_states, n_disturbance_states>(
      System::n_states, System::n_states) =
      Eigen::Matrix<double, n_disturbance_states,
                    n_disturbance_states>::Identity();

  int index_delay_states = System::n_states + n_disturbance_states;

  for (int i = 0; i < System::n_control_inputs; i++) {
    if (n_delay_[i] == 0) {
      sys_out.B.template block<System::n_states, 1>(0, i) =
          sys_discrete.B.col(i);
    } else {
      sys_out.A.template block<System::n_states, 1>(0, index_delay_states) =
          sys_discrete.B.col(i);
      const int size_block = n_delay_[i] - 1;
      sys_out.A.template block(index_delay_states, index_delay_states + 1,
                               size_block, size_block) =
          Eigen::MatrixXd::Identity(size_block, size_block);
      index_delay_states += n_delay_[i];
    }
  }

  sys_out.C.template leftCols<System::n_states>() = sys_discrete.C;
  sys_out.C.template block<System::n_outputs, n_disturbance_states>(0,
      System::n_states + n_disturbance_states) =
      Eigen::Matrix<double, System::n_outputs,
                    n_disturbance_states>::Identity();

  return sys_out;
}
#include "augmented_system_list.h"
