#ifndef READ_FILES_H
#define READ_FILES_H

#include <exception>
#include <fstream>
#include <sstream>
#include <string>
#include <iostream>

#include "input_constraints.h"

// Read a single string variable from a stream
std::string ReadString(std::ifstream &setup_file, int &line_number) {
  std::string line;
  std::string output;
  std::string remaining_text;

  while (std::getline(setup_file, line)) {
    line_number++;
    // Skip commented and empty lines
    if ((line[0] == '#') || (line.empty()) ||
        (std::all_of(line.begin(), line.end(), isspace)))
      continue;

    std::istringstream line_stream(line);

    if (!(line_stream >> output)) {
      std::string errmsg = "Error reading folder-name at line number" +
                           std::to_string(line_number);
      throw std::runtime_error(errmsg);
    }

    if (line_stream >> remaining_text) {
      std::cerr << "Extra text \"" << remaining_text << "\" in line "
                << line_number << " being ignored." << std::endl;
    }
    return output;
  }

  std::string errmsg =
      "Error reading setup file at line number" + std::to_string(line_number);
  throw std::runtime_error(errmsg);
}

// Read array of numbers from a stream
template <typename T>
int ReadNumbers(T *data_out, const int n_elements, std::ifstream &setup_file,
                int &line_number, bool throw_error = true) {
  std::string line;
  std::istringstream line_stream(line);
  std::string output;
  std::string remaining_text;
  int i = 0;

  while (i < n_elements) {
    if (std::getline(setup_file, line)) {
      line_number++;
      // Skip commented and empty lines
      if ((line[0] == '#') || (line.empty()) ||
          (std::all_of(line.begin(), line.end(), isspace)))
        continue;

      line_stream = std::istringstream(line);
      while (line_stream >> data_out[i]) {
        i++;
        if (i == n_elements) break;
      }
    } else {
      std::string errmsg = "Error reading setup file at line number" +
                           std::to_string(line_number) +
                           " (probably not enough entries given)";
      if (throw_error) throw std::runtime_error(errmsg);
      break;
    }
  }
  if (line_stream >> remaining_text) {
    std::cerr << "Extra text \"" << remaining_text << "\" in line "
              << line_number << " being ignored." << std::endl;
  }
  return i;
}


#endif
