#include <Eigen/Eigen>
#include "time_delay.h"
#include "print_matrix.h"

constexpr int n_inputs = 3;
typedef Eigen::Matrix<double, n_inputs, 1> Input;

int main(void) {
  int n_delay[n_inputs] = {1, 2, 3};
  constexpr int sum_delay = 6;
  Input u_init;
  u_init << 1, 5, 3;

  TimeDelay<n_inputs, sum_delay> delay(n_delay);

  Input u = delay.GetDelayedInput(u_init);
  std::cout << "Delay state: " << std::endl;
  delay.PrintCurrentState(std::cout);
  print_matrix(std::cout, u, "Delay input");

  u = delay.GetDelayedInput(u_init);
  std::cout << "Delay state: " << std::endl;
  delay.PrintCurrentState(std::cout);
  print_matrix(std::cout, u, "Delay input");

  u = delay.GetDelayedInput(u_init);
  std::cout << "Delay state: " << std::endl;
  delay.PrintCurrentState(std::cout);
  print_matrix(std::cout, u, "Delay input");

  u_init = (Input() << 4, 3, 2).finished();

  u = delay.GetDelayedInput(u_init);
  std::cout << "Delay state: " << std::endl;
  delay.PrintCurrentState(std::cout);
  print_matrix(std::cout, u, "Delay input");

  u = delay.GetDelayedInput(u_init);
  std::cout << "Delay state: " << std::endl;
  delay.PrintCurrentState(std::cout);
  print_matrix(std::cout, u, "Delay input");
};
