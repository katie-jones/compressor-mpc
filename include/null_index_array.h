#ifndef NULL_INDEX_ARRAY_H
#define NULL_INDEX_ARRAY_H
#include <utility>
#include "constexpr_array.h"

/**
 * Constexpr class representing an array of integers such that x[i] = i.
 * Implements the same functionality as the ConstexprArray class.
 * Template parameters:\n
 * - size_in: size of array
 */
template <int size_in>
class NullIndexArray {
 public:
   /// Constructor
  constexpr NullIndexArray() {
    static_assert(size_in > 0, "Size should be positive");
  }

  /// Get array value at location i
  static constexpr int GetEntry(const int i) { return i; }

  /// Size of array
  static constexpr int size = size_in;

  /// Type of this object
  using type = NullIndexArray<size_in>;

  /// Return std integer sequence with same values as this array
  static constexpr std::make_integer_sequence<int, size> GetIntegerSequence() {
    return std::make_integer_sequence<int, size>();
  }

  /// Return std index sequence with same values as this array
  static constexpr std::make_index_sequence<static_cast<std::size_t>(size)>
  GetIndexSequence() {
    return std::make_index_sequence<static_cast<std::size_t>(size)>();
  }

  /// Return the entries of x_in given in Ints... as smaller array x_out
  template <typename T>
  static void GetSubArray(T* x_out, const T* x_in) {
    for (int i = 0; i < size; i++) {
      x_out[i] = x_in[i];
    }
  }

  /// Get submatrix from an Eigen matrix
  static void GetSubMatrix(Eigen::Matrix<double, size, size>* x_reduced,
                           const Eigen::Matrix<double, size, size>& x) {
    x_reduced = &x;
  }

  /// Get submatrix from an Eigen matrix
  template <int N, typename T>
  static void GetSubMatrix(T* new_data, const T* old_data) {
    int n = size * size;
    for (int i = 0; i < n; i++) {
      *new_data = *old_data;
      new_data++;
      old_data++;
    }
  }

  /// Sub-array of current array with indices SubInts...
  template <int... SubInts>
  using ConstexprSubArray = ConstexprArray<SubInts...>;

  /// Return sub-array of current array based on array given
  template <int... SubInts>
  static constexpr ConstexprSubArray<SubInts...> GetConstexprSubArray(
      ConstexprArray<SubInts...>) {
    return ConstexprSubArray<SubInts...>();
  }

  /// Return sub-array of current array based on array given
  template <int... SubInts>
  static constexpr ConstexprSubArray<SubInts...> GetConstexprSubArray(
      std::integer_sequence<int, SubInts...>) {
    return ConstexprSubArray<SubInts...>();
  }

  /// Sub-array of current array
  template <typename Indices>
  using IndicesSubArray = decltype(GetConstexprSubArray(Indices()));

  /// Add the entries of x_in to entries Ints... in larger array x_out.
  /// Entries of x_out not contained in Ints... are not modified.
  template <typename T>
  static void ExpandArray(T* x_out, const T* x_in) {
    for (int i = 0; i < size; i++) {
      x_out[i] += x_in[i];
    }
  }

  /// Get subvector from an Eigen vector
  static Eigen::Matrix<double, size, 1> GetSubVector(
      const Eigen::Matrix<double, size, 1>& x) {
    return x;
  }

  /// Get subvector from an Eigen vector
  static void GetSubVector(Eigen::Matrix<double, size, 1>* x_out,
                           const Eigen::Matrix<double, size, 1>& x) {
    x_out = &x;
  }

  constexpr int operator[](const int i) const { return i; }
};

#endif
