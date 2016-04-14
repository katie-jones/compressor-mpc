#include <Eigen/Eigen>

using namespace Eigen;

class comp_state {
  Matrix<double, 5, 1> data;

  template <type_t t>
  comp_state operator+(
