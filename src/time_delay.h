#ifndef TIME_DELAY_H
#define TIME_DELAY_H

#include <Eigen/Eigen>
#include <iostream>
#include "mpc_exceptions.h"

template <int n_inputs, int n_delay_states>
class TimeDelay {
 public:
  typedef Eigen::Matrix<double, n_inputs, 1> Input;

  TimeDelay(const std::array<int, n_inputs> n_delay) {
    // check if n_delay_states is correct
    int sum_delay = 0;
    for (int i = 0; i < n_inputs; i++) sum_delay += n_delay[i];
    if (sum_delay != n_delay_states) {
      throw delay_states_wrong();
    }

    // initialize to zero
    for (int i=0; i<n_delay_states; i++) {
      u_delay_[i] = 0;
    }

    // insert values into array
    sum_delay = 0;
    for (int i = 0; i < n_inputs; i++) {
      n_delay_[i] = n_delay[i];
      current_input_[i] = sum_delay;
      sum_delay += n_delay_[i];
    }
  }

  Input GetDelayedInput(Input u_next) {
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
    for (int i=0; i<n_delay_states; i++) {
      os << u_delay_[i] << " ";
    }
    os << std::endl;
    for (int i=0; i<n_inputs; i++) {
      os << current_input_[i] << " ";
    }
    os << std::endl;
  }
    

 private:
  double u_delay_[n_delay_states];
  int n_delay_[n_inputs];
  int current_input_[n_inputs];
};

#endif
