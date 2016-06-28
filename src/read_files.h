#ifndef READ_FILES_H
#define READ_FILES_H

#include <fstream>

#include "input_constraints.h"

namespace ReadFiles {
template <typename T>
int ReadDataFromFile(T* data_out, const int n_elem, const std::string& filename,
                     const int skip_size = 1) {
  std::ifstream read_file;
  read_file.open(filename);
  for (int i = 0; i < n_elem; i++) {
    if (!(read_file >> data_out[skip_size * i])) {
      std::cerr << "Error reading data from file " << filename << std::endl;
      return 0;
    }
  }
  read_file.close();
  return 1;
}

template <typename T>
int ReadDataFromStream(T* data_out, std::istream& input_stream,
                       const int n_elem, const int skip_size = 1) {
  for (int i = 0; i < n_elem; i++) {
    if (!(input_stream >> data_out[skip_size * i])) {
      std::cerr << "Error reading data from stream." << std::endl;
      return 0;
    }
  }
  return 1;
}

template <int n_control_inputs>
int ReadConstraintsFromFile(InputConstraints<n_control_inputs>* constraints,
                            const std::string& filename) {
  std::ifstream constraints_file;
  constraints_file.open(filename);

  if (!((ReadDataFromStream(constraints->lower_bound.data(), constraints_file,
                           n_control_inputs)) &&
      (ReadDataFromStream(constraints->upper_bound.data(), constraints_file,
                          n_control_inputs)) &&
      (ReadDataFromStream(constraints->lower_rate_bound.data(),
                          constraints_file, n_control_inputs)) &&
      (ReadDataFromStream(constraints->upper_rate_bound.data(),
                          constraints_file, n_control_inputs)))) {
    std::cout << "Error reading from file " << filename << "." << std::endl;
    constraints_file.close();
    return 0;
  }
  return 1;
}
}

#endif
