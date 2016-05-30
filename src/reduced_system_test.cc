#include <iostream>
#include <fstream>
#include <Eigen/Eigen>

#include "parallel_compressors.h"
#include "reduced_dynamic_system.h"

constexpr int n_states = 6;
constexpr int n_outputs = 2;
constexpr int n_inputs = 2;
using RedSys =
    ReducedDynamicSystem<n_states, n_outputs, n_inputs, ParallelCompressors>;

int main(void) {
  ParallelCompressors compressor;
  Eigen::Array<int, n_states, 1> index_states;
  Eigen::Array<int, n_outputs, 1> index_outputs;
  Eigen::Array<int, n_inputs, 1> index_inputs;

  index_states << 0, 1, 2, 3, 4, 10;
  index_outputs << 0, 3;
  index_inputs << 0, 1;
  ParallelCompressors::State x_default = ParallelCompressors::GetDefaultState();

  RedSys comp_reduced(&compressor, x_default, compressor.GetDefaultInput(),
                      index_states, index_outputs, index_inputs);
  RedSys::State x;
  x << x_default.head<5>(), x_default(10);

  std::cout << comp_reduced.GetLinearizedSystem(
                   x, ParallelCompressors::GetDefaultInput()) << std::endl;
}
