#ifndef REDUCED_DYNAMIC_SYSTEM_H
#define REDUCED_DYNAMIC_SYSTEM_H

#include "dynamic_system.h"
#include <Eigen/Eigen>

template <int n_states, int n_outputs, typename SubControlInputIndex,
          class OriginalSystem>
class ReducedDynamicSystem
    : public DynamicSystem<n_states, OriginalSystem::n_inputs, n_outputs,
                           typename OriginalSystem::ControlInputIndex::
                               template IndicesSubArray<SubControlInputIndex>> {
 public:
  constexpr static int n_control_inputs = SubControlInputIndex::size;
  using ControlInputIndex =
      typename OriginalSystem::ControlInputIndex::template IndicesSubArray<
          SubControlInputIndex>;

  typedef typename DynamicSystem<n_states, OriginalSystem::n_inputs, n_outputs,
                                 ControlInputIndex>::State State;
  typedef typename DynamicSystem<n_states, OriginalSystem::n_inputs, n_outputs,
                                 ControlInputIndex>::Output Output;
  typedef typename DynamicSystem<n_states, OriginalSystem::n_inputs, n_outputs,
                                 ControlInputIndex>::Input Input;
  typedef typename DynamicSystem<n_states, OriginalSystem::n_inputs, n_outputs,
                                 ControlInputIndex>::ControlInput ControlInput;
  typedef typename DynamicSystem<n_states, OriginalSystem::n_inputs, n_outputs,
                                 ControlInputIndex>::Linearized Linearized;

  // Get state of full system based on x_full_ and input x
  const typename OriginalSystem::State GetFullState(const State& x) const {
    typename OriginalSystem::State x_out = x_full_;
    for (int i = 0; i < n_states; i++) {
      x_out(index_states_(i)) = x(i);
    }
    return x_out;
  }

  // Get state of reduced system based on full state
  State GetReducedState(const typename OriginalSystem::State& x_full) const {
    State x_out;
    for (int i = 0; i < n_states; i++) {
      x_out(i) = x_full(index_states_(i));
    }
    return x_out;
  }

  // Get output of reduced system based on full output
  Output GetReducedOutput(const typename OriginalSystem::Output& y_full) const {
    Output y_out;
    for (int i = 0; i < n_outputs; i++) {
      y_out(i) = y_full(index_outputs_(i));
    }
    return y_out;
  }

  // Get control input of reduced system based on full control input
  ControlInput GetReducedInput(
      const typename OriginalSystem::ControlInput& u_full) const {
    ControlInput u_out;
    for (int i = 0; i < n_control_inputs; i++) {
      u_out(i) = u_full(index_inputs_(i));
    }
    return u_out;
  }

 private:
  const OriginalSystem* p_sys_;
  typename OriginalSystem::State x_full_;
  Input u_offset_;

  const Eigen::Array<int, n_states, 1> index_states_;
  const Eigen::Array<int, n_outputs, 1> index_outputs_;
  const Eigen::Array<int, n_control_inputs, 1> index_inputs_;

 public:
  /// Constructor
  ReducedDynamicSystem(
      const OriginalSystem* p_sys, const typename OriginalSystem::State& x_full,
      const Input& u_offset, const Eigen::Array<int, n_states, 1> index_states,
      const Eigen::Array<int, n_outputs, 1> index_outputs,
      const Eigen::Array<int, n_control_inputs, 1> index_inputs)
      : p_sys_(p_sys),
        x_full_(x_full),
        u_offset_(u_offset),
        index_states_(index_states),
        index_outputs_(index_outputs),
        index_inputs_(index_inputs) {}

  /// Return system linearized about given operating point.
  virtual Linearized GetLinearizedSystem(const State& x, const Input& u) const {
    Linearized sys_out;
    typename OriginalSystem::Linearized sys_orig =
        p_sys_->GetLinearizedSystem(GetFullState(x), u);

    for (int i = 0; i < n_states; i++) {
      sys_out.A.row(i) = GetReducedState(sys_orig.A.row(index_states_(i)));
      sys_out.B.row(i) = GetReducedInput(sys_orig.B.row(index_states_(i)));
      sys_out.C.col(i) =
          GetReducedOutput(sys_orig.C.col(index_states_(i)).transpose());
    }
    sys_out.f = GetReducedState(sys_orig.f);

    return sys_out;
  }

  /// Return derivative of system about given operating point.
  virtual State GetDerivative(const State& x, const Input& u) const {
    return GetReducedState(p_sys_->GetDerivative(GetFullState(x), u));
  }

  /// Return system output at given state.
  virtual Output GetOutput(const State& x) const {
    return GetReducedOutput(p_sys_->GetOutput(GetFullState(x)));
  }
};

#endif
