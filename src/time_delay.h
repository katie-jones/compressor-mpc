#ifndef TIME_DELAY_H
#define TIME_DELAY_H

#include <Eigen/Eigen>
#include <iostream>
#include "mpc_exceptions.h"

template <typename Delays>
class TimeDelay {
 public:
  constexpr static int n_inputs = Delays::size;
  constexpr static int n_delay_states = Delays::GetSum();

  typedef Eigen::Matrix<double, n_inputs, 1> Input;

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
template<typename Delays>
constexpr Delays TimeDelay<Delays>::n_delay_;

#endif
