#ifndef DYNAMIC_SYSTEM_H
#define DYNAMIC_SYSTEM_H

#include <Eigen/Eigen>
#include <iostream>
#include <iomanip>

// class Observer;

/**
 * Abstract class describing a dynamic system.
 * Template with parameters: n_states (number of states), n_inputs (number of
 * inputs), n_outputs (number of outputs), n_control_inputs (number of inputs
 * used for control).
 */
template <int n_states, int n_inputs, int n_outputs, int n_control_inputs>
class DynamicSystem {
  // friend Observer;

 protected:
  /// Type describing system state.
  typedef Eigen::Array<double, n_states, 1> State;

  /// Type describing system input.
  typedef Eigen::Array<double, n_inputs, 1> Input;

  /// Type describing system output.
  typedef Eigen::Matrix<double, n_outputs, 1> Output;

 public:
  /// Linearized form of dynamic system.
  struct Linearized {
    Eigen::Matrix<double, n_states, n_states, Eigen::RowMajor> A;
    Eigen::Matrix<double, n_states, n_control_inputs, Eigen::RowMajor> B;
    Eigen::Matrix<double, n_outputs, n_states, Eigen::RowMajor> C;
    Eigen::Matrix<double, n_outputs, n_control_inputs, Eigen::RowMajor> D;
    Eigen::Matrix<double, n_states, 1> f;
    friend std::ostream& operator<<(std::ostream& os,
                                    const Linearized& linsys) {
      print_matrix(os, linsys.A, "A");
      print_matrix(os, linsys.B, "B");
      print_matrix(os, linsys.C, "C");
      return os;
    }

   private:
    template <typename Derived>
    static std::ostream& print_matrix(std::ostream& os,
                                      const Eigen::MatrixBase<Derived>& x,
                                      const std::string name) {
      os << name << " = " << std::endl;
      for (int i = 0; i < x.rows(); i++) {
        for (int j = 0; j < x.cols(); j++) {
          os << std::setprecision(3) << std::setw(10) << x(i, j) << " ";
        }
        os << std::endl;
      }
      os << std::endl;
      return os;
    }
  };

  virtual ~DynamicSystem() {}

  /// Return system linearized about given operating point.
  virtual Linearized GetLinearizedSystem(const State x,
                                         const Input u) const = 0;

  /// Return derivative of system about given operating point.
  virtual State GetDerivative(const State x, const Input u) const = 0;

  /// Return system output at given state.
  virtual Output GetOutput(const State x) const = 0;
};

#endif