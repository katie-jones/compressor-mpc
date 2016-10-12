#ifndef TIME_DELAY_H
#define TIME_DELAY_H

#include <Eigen/Eigen>
#include <iostream>

/**
 * Class representing a discrete time delay in a dynamic system. Stores a
 * pre-defined number of delayed values for each input.
 * Template parameters:
 * - Delays: ConstexprArray giving number of delays per input
 */
template <typename Delays>
class TimeDelay {
 public:
  /// Number of inputs
  constexpr static int n_inputs = Delays::size;
  /// Total number of delay states
  constexpr static int n_delay_states = Delays::GetSum();

  /// Input to system
  typedef Eigen::Matrix<double, n_inputs, 1> Input;

  /// Constructor: initializes memory to 0
  TimeDelay() {
    // initialize to zero
    for (int i = 0; i < n_delay_states; i++) {
      u_delay_[i] = 0;
    }

    // initialize current input index
    int sum_delay = 0;
    for (int i = 0; i < n_inputs; i++) {
      current_input_[i] = sum_delay;
      sum_delay += n_delay_[i];
    }
  }

  /// Store the next input values and output the delayed input to be applied to
  /// system
  Input GetDelayedInput(const Input& u_next) {
    Input u_out = Input::Zero();
    int index_delay_states = 0;
    for (int i = 0; i < n_inputs; i++) {
      if (n_delay_[i] == 0) {
        u_out(i) = u_next(i);
      } else {
        index_delay_states += n_delay_[i];
        u_out(i) = u_delay_[current_input_[i]];
        u_delay_[current_input_[i]] = u_next(i);
        current_input_[i]++;
        if (current_input_[i] == index_delay_states) {
          current_input_[i] -= n_delay_[i];
        }
      }
    }
    return u_out;
  }

  /// Debugging function to print the current state stored in memory
  void PrintCurrentState(std::ostream& os) {
    for (int i = 0; i < n_delay_states; i++) {
      os << u_delay_[i] << " ";
    }
    os << std::endl;
    for (int i = 0; i < n_inputs; i++) {
      os << current_input_[i] << " ";
    }
    os << std::endl;
  }

 private:
  double u_delay_[n_delay_states];
  static constexpr Delays n_delay_ = Delays();
  int current_input_[n_inputs];
};

// Definition of constexpr static member
template <typename Delays>
constexpr Delays TimeDelay<Delays>::n_delay_;

#endif
