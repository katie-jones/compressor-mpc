#ifndef AUGMENTED_LINEARIZED_SYSTEM_H
#define AUGMENTED_LINEARIZED_SYSTEM_H

#include <Eigen/Eigen>
#include <Eigen/SparseCore>

template <class System, int n_delay_states, int n_disturbance_states>
class AugmentedLinearizedSystem {
  static constexpr int n_control_inputs = System::n_control_inputs;
  static constexpr int n_inputs = System::n_inputs;
  static constexpr int n_outputs = System::n_outputs;
  static constexpr int n_states = System::n_states;

  static constexpr int n_aug_states = n_disturbance_states + n_delay_states;
  static constexpr int n_obs_states = n_states + n_disturbance_states;
  static constexpr int n_total_states = n_aug_states + n_states;

 public:
  /// State of the dynamic system
  typedef Eigen::Matrix<double, n_states, 1> State;

  /// Augmented state including delay and disturbance states
  typedef Eigen::Matrix<double, n_total_states, 1> AugmentedState;

  /// Controlled inputs of dynamic system
  typedef Eigen::Matrix<double, n_control_inputs, 1> ControlInput;

  /// Output of dynamic system
  typedef Eigen::Matrix<double, n_outputs, 1> Output;

  /// Vector of indices for a ControlInput
  typedef std::array<int, n_control_inputs> ControlInputIndex;

 private:
  struct AComposite;
  struct BComposite;

  struct Ctype : public Eigen::Matrix<double, n_outputs, n_total_states> {
    Ctype& operator*=(const AComposite& a);
    Eigen::Matrix<double, System::n_outputs, System::n_control_inputs>
    operator*(const BComposite& b);
  };

  struct AComposite {
    Eigen::Matrix<double, n_states, n_states, Eigen::RowMajor> Aorig;
    Eigen::Matrix<double, n_states, n_control_inputs, Eigen::RowMajor> Adelay;
    Eigen::SparseMatrix<bool, Eigen::RowMajor> Aaug;
    ControlInputIndex n_delay;

    AComposite(const ControlInputIndex& n_delay);
    AugmentedState operator*(const AugmentedState& x) const;
  };

  struct BComposite {
    Eigen::Matrix<double, n_states, n_control_inputs, Eigen::RowMajor> Borig;
    Eigen::SparseMatrix<bool, Eigen::RowMajor> Baug;

    BComposite(const ControlInputIndex& n_delay);
    AugmentedState operator*(const ControlInput& u) const;
  };

  AComposite A;
  BComposite B;
  Eigen::Matrix<double, n_outputs, n_obs_states, Eigen::RowMajor> C;
  State f;

 public:
  AugmentedLinearizedSystem(const ControlInputIndex& n_delay_in);
  void Update(const typename System::Linearized& sys_discrete);
  const AComposite GetA() const { return A; }
  const BComposite GetB() const { return B; }
};

namespace Eigen {
namespace internal {
template <>
struct scalar_product_traits<double, bool> {
  enum { Defined = 1 };
  typedef double ReturnType;
};

template <>
struct scalar_product_traits<bool, double> {
  enum { Defined = 1 };
  typedef double ReturnType;
};
}
}

#endif
