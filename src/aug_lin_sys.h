#ifndef AUGMENTED_LINEARIZED_SYSTEM_H
#define AUGMENTED_LINEARIZED_SYSTEM_H

#include <Eigen/Eigen>
#include <Eigen/SparseCore>

template <class System, int n_delay_states, int n_disturbance_states>
class AugmentedLinearizedSystem {
 public:
  static constexpr int n_control_inputs = System::n_control_inputs;
  static constexpr int n_inputs = System::n_inputs;
  static constexpr int n_outputs = System::n_outputs;
  static constexpr int n_states = System::n_states;

  static constexpr int n_aug_states = n_disturbance_states + n_delay_states;
  static constexpr int n_obs_states = n_states + n_disturbance_states;
  static constexpr int n_total_states = n_aug_states + n_states;

  /// State of the dynamic system
  typedef Eigen::Matrix<double, n_states, 1> State;

  /// Input to dynamic system
  typedef typename System::Input Input;

  /// Augmented state including delay and disturbance states
  typedef Eigen::Matrix<double, n_total_states, 1> AugmentedState;

  /// Controlled inputs of dynamic system
  typedef Eigen::Matrix<double, n_control_inputs, 1> ControlInput;

  /// Output of dynamic system
  typedef Eigen::Matrix<double, n_outputs, 1> Output;

  /// Vector of indices for a ControlInput
  typedef std::array<int, n_control_inputs> ControlInputIndex;

  struct AComposite;
  struct BComposite;

  struct Ctype : public Eigen::Matrix<double, n_outputs, n_total_states> {
    inline Ctype& operator*=(const AComposite& a);
    inline Eigen::Matrix<double, System::n_outputs, System::n_control_inputs>
    operator*(const BComposite& b);
  };

  struct AComposite {
    Eigen::Matrix<double, n_states, n_states> Aorig;
    Eigen::Matrix<double, n_states, n_control_inputs> Adelay;
    Eigen::SparseMatrix<bool> Aaug;
    ControlInputIndex n_delay;

    AComposite(const ControlInputIndex& n_delay);
    inline AugmentedState operator*(const AugmentedState& x) const;
    AugmentedState TimesAugmentedOnly(
        const Eigen::Matrix<double, n_aug_states, 1> x) const;
  };

  struct BComposite {
    Eigen::Matrix<double, n_states, n_control_inputs> Borig;
    Eigen::SparseMatrix<bool> Baug;

    BComposite(const ControlInputIndex& n_delay);
    inline AugmentedState operator*(const ControlInput& u) const;
  };

  /// Structure containing prediction matrices of augmented system
  struct Prediction {
    Eigen::MatrixXd Sx, Sf, Su;
  };

  AComposite A;
  BComposite B;
  Eigen::Matrix<double, n_outputs, n_obs_states> C;
  State f;
  System sys_;
  double sampling_time_;

  AugmentedLinearizedSystem(const System& sys, const double sampling_time,
                            const ControlInputIndex& n_delay_in);
  void Update(const State x, const Input& u);

  // Discretize system using runge-kutta 4 method
  static const typename System::Linearized DiscretizeRK4(
      const typename System::Linearized& sys_continuous, const double Ts);
  const Prediction GeneratePrediction(const int p, const int m) const;

  const AComposite GetA() const { return A; }
  const BComposite GetB() const { return B; }
};

/*
 * Operator x_out = A*x
 */
template <class System, int n_delay_states, int n_disturbance_states>
inline typename AugmentedLinearizedSystem<System, n_delay_states,
                                          n_disturbance_states>::AugmentedState
    AugmentedLinearizedSystem<System, n_delay_states,
                              n_disturbance_states>::AComposite::
    operator*(const AugmentedState& x) const {
  // multiply by augmented part of matrix
  AugmentedState x_out = TimesAugmentedOnly(x.template tail<n_aug_states>());

  // add contribution of Aorig
  x_out.template head<n_states>() += Aorig * x.template head<n_states>();

  return x_out;
}

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
