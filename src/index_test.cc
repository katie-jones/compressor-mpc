#include <fstream>
#include <iostream>
#include <tuple>
#include <utility>
#include <boost/timer/timer.hpp>

#include "constexpr_array.h"

int main(void) {
  ConstexprArray<0, 4, 2> x;
  double data[] = {0.2, 0.5, 0.2, -0.7, -3};
  double sub_data[3];

  // Time solution 1
  boost::timer::cpu_timer timer1;
  for (auto k = 0; k < 1e6; k++) {
    for (auto j = 0; j < 1e6; j++) {
      for (auto i = 0; i < 1e6; i++) {
        // data[0] += 0.01;
        x.GetSubArray(sub_data, data);
      }
    }
  }
  boost::timer::cpu_times int_elapsed = timer1.elapsed();
  boost::timer::nanosecond_type elapsed_ns(int_elapsed.wall);  // system +

  return 0;
}
