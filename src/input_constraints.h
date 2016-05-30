#ifndef INPUT_CONSTRAINTS_H
#define INPUT_CONSTRAINTS_H

#include <Eigen/Eigen>

/// Input constraints controller should respect.
template <int n_control_inputs>
struct InputConstraints {
  typedef Eigen::Matrix<double, n_control_inputs, 1> ControlInput;

  ControlInput lower_bound, upper_bound, lower_rate_bound, upper_rate_bound;
  bool use_rate_constraints;
  InputConstraints()
      : upper_bound(ControlInput::Constant(std::nan(""))),
        lower_bound(ControlInput::Constant(-std::nan(""))),
        use_rate_constraints(false) {}
};

#endif
