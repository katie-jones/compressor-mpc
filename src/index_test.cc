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

  // Time solution 1
  boost::timer::cpu_timer timer1;
  // for (auto k = 0; k < 1e6; k++) {
  for (auto j = 0; j < 1e2; j++) {
    for (auto i = 0; i < 1e6; i++) {
      x.GetSubArray(sub_data, data);
    }
  }
  // }
  boost::timer::cpu_times int_elapsed = timer1.elapsed();
  boost::timer::nanosecond_type elapsed_ns(int_elapsed.wall);  // system +

  std::cout << sub_data[0] << std::endl;

  std::cout << "Pre-compiled: " << elapsed_ns << std::endl;

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
  int_elapsed = timer1.elapsed();
  boost::timer::nanosecond_type elapsed_ns2(int_elapsed.wall);  // system +

  std::cout << "Runtime: " << elapsed_ns2 << std::endl;

  std::cout << "Speed-up: "
            << static_cast<double>(elapsed_ns2) /
                   static_cast<double>(elapsed_ns) << std::endl;
  x.ExpandArray(data, sub_data);

  return 0;
}
