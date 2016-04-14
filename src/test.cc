#include <Eigen/Eigen>
#include "gtest/gtest.h"

#include <iostream>

#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>

using namespace std;
using namespace boost::numeric::odeint;

const double sigma = 10.0;
const double R = 28.0;
const double b = 8.0 / 3.0;

typedef boost::array<double, 3> state_type;

void lorenz(const state_type &x, state_type &dxdt, double t) {
  dxdt[0] = sigma * (x[1] - x[0]);
  dxdt[1] = R * x[0] - x[1] - x[0] * x[2];
  dxdt[2] = -b * x[2] + x[0] * x[1];
}

void write_lorenz(const state_type &x, const double t) {
  cout << t << '\t' << x[0] << '\t' << x[1] << '\t' << x[2] << endl;
}

using namespace Eigen;
using namespace std;

TEST(gtest, linking) { EXPECT_EQ(1, 1); }

int main(int argc, char *argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  VectorXd v(10);
  v[0] = 10.0;
  cout << v[0] << endl;
  state_type x = {10.0, 1.0, 1.0};  // initial conditions
  // integrate(lorenz, x, 0.0, 25.0, 0.1, write_lorenz);
  
  return RUN_ALL_TESTS();
}

