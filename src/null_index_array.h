#ifndef NULL_INDEX_ARRAY_H
#define NULL_INDEX_ARRAY_H
#include <utility>
#include "constexpr_array.h"

// Type returning original vector
template <int size_in>
class NullIndexArray {
 public:
  constexpr NullIndexArray() {
    static_assert(size_in > 0, "Size should be positive");
  }

  static constexpr int GetEntry(const int i) {
    return i;
  }

  static constexpr int size = size_in;
  using type = NullIndexArray<size_in>;

  static constexpr std::make_integer_sequence<int, size> GetIntegerSequence() {
    return std::make_integer_sequence<int, size>();
  }

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

  // Function to get subvector from an Eigen vector
  static Eigen::Matrix<double, size, 1> GetSubVector(
      const Eigen::Matrix<double, size, 1>& x) {
    return x;
  }

  // Function to get subvector from an Eigen vector
  static void GetSubVector(Eigen::Matrix<double, size, 1>* x_out,
                           const Eigen::Matrix<double, size, 1>& x) {
    *x_out = x;
  }

  constexpr int operator[](const int i) const { return i; }
};

#endif
