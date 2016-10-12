#ifndef INPUT_CONSTRAINTS_H
#define INPUT_CONSTRAINTS_H

#include <Eigen/Eigen>

/**
 * Input constraints controller should respect.
 * Template parameters:\n
 * - n_control_inputs: number of control inputs to system
 */
template <int n_control_inputs>
struct InputConstraints {
  /// Control input to system
  typedef Eigen::Matrix<double, n_control_inputs, 1> ControlInput;

  /// Constraints (absolute and rate) on control input
  ControlInput lower_bound, upper_bound, lower_rate_bound, upper_rate_bound;

  /// Activate rate constraints (default: false)
  bool use_rate_constraints;

  /// Constructor
  InputConstraints()
      : upper_bound(ControlInput::Constant(std::nan(""))),
        lower_bound(ControlInput::Constant(-std::nan(""))),
        use_rate_constraints(false) {}
};

#endif
