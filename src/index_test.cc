#include <fstream>
#include <iostream>
#include <tuple>
#include <utility>
#include <boost/timer/timer.hpp>

#include "constexpr_array.h"

int main(void) {
  ConstexprArray<0, 4, 2> x;
  volatile double data[] = {0.2, 0.5, 0.2, -0.7, -3};
  volatile double sub_data[3];

  std::cout << "Vector sub arrays: " << std::endl;

  // Time solution 1
  boost::timer::cpu_timer timer1;
  // for (auto k = 0; k < 1e6; k++) {
  for (auto j = 0; j < 1e2; j++) {
    for (auto i = 0; i < 1e6; i++) {
      x.GetSubArray(sub_data, data);
    }
  }
  // }
  timer1.stop();
  boost::timer::cpu_times int_elapsed = timer1.elapsed();
  boost::timer::nanosecond_type elapsed_ns(int_elapsed.system +
                                           int_elapsed.user);  // system +

  std::cout << sub_data[0] << std::endl;

  std::cout << "Pre-compiled: " << elapsed_ns << std::endl;

  int n_elem = x.size;
  // Time solution 1
  boost::timer::cpu_timer timer2;
  // for (auto k = 0; k < 1e6; k++) {
  for (auto j = 0; j < 1e2; j++) {
    for (auto i = 0; i < 1e6; i++) {
      for (auto n = 0; n < x.size; n++) {
        sub_data[n] = data[x[n]];
      }
    }
  }
  // }
  timer2.stop();
  int_elapsed = timer2.elapsed();
  boost::timer::nanosecond_type elapsed_ns2(int_elapsed.system +
                                            int_elapsed.user);  // system +

  std::cout << "Runtime: " << elapsed_ns2 << std::endl;

  std::cout << "Speed-up: "
            << static_cast<double>(elapsed_ns2) /
                   static_cast<double>(elapsed_ns) << std::endl;

  // ----------- SUB MATRICES --------------
  timer1.stop();
  timer2.stop();

  Eigen::Matrix<volatile double, 3, 3> x_mat;
  for (int i = 0; i < 9; i++) {
    x_mat(i % 3, i / 3) = static_cast<volatile double>(i + 1);
  }
  Eigen::Matrix<volatile double, 2, 2> x_red;

  timer1.start();
  for (auto i = 0; i < 1e7; i++) {
    ConstexprArray<0, 2>::GetSubMatrix(&x_red, x_mat);
  }
  // timer1.stop();
  int_elapsed = timer1.elapsed();
  elapsed_ns = int_elapsed.system + int_elapsed.user;
  // elapsed_ns = int_elapsed.wall;

  timer2.start();
  for (auto i = 0; i < 1e7; i++) {
    for (auto n = 0; n < 2; n++) {
      x_red(n, n) = x_mat(ConstexprArray<0, 2>::GetEntry(n),
                          ConstexprArray<0, 2>::GetEntry(n));
    }
  }
  // timer2.stop();
  int_elapsed = timer2.elapsed();
  elapsed_ns2 = int_elapsed.system + int_elapsed.user;
  // elapsed_ns2 = int_elapsed.wall;

  std::cout << "GetSubMatrix: " << elapsed_ns << std::endl
            << "For loop: " << elapsed_ns2 << std::endl
            << "Speedup: "
            << static_cast<double>(elapsed_ns2) /
                   static_cast<double>(elapsed_ns) << std::endl;
  return 0;
}
