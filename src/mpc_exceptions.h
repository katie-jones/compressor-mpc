#ifndef MPC_EXCEPTIONS_H
#define MPC_EXCEPTIONS_H

#include <exception>

class delay_states_wrong : public std::exception {
  virtual const char* what() const throw() {
    return "Total number of delay states given as template parameter is "
           "incorrect.";
  }
};

#endif
