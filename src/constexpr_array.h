#ifndef CONSTEXPR_ARRAY_H
#define CONSTEXPR_ARRAY_H
#include <utility>

template <int... Ints>
class ConstexprArray {
 private:
  static constexpr std::array<int, sizeof...(Ints)> data_ = {Ints...};

 public:
  constexpr ConstexprArray() {}

  /// Number of control inputs
  static constexpr int size = static_cast<int>(sizeof...(Ints));

  /// Get sum of entries
  static constexpr int GetSum() {
    int sum = 0;
    int entry = 0;
    for (int i = 0; i < data_.size(); ++i) {
      entry = data_[i];
      sum += entry;
    }
    return sum;
  }

  static constexpr int GetEntry(const int i) { return data_[i]; }

  static constexpr std::integer_sequence<int, Ints...> GetIntegerSequence() {
    return std::integer_sequence<int, Ints...>();
  }

  static constexpr std::index_sequence<static_cast<std::size_t>(Ints)...>
  GetIndexSequence() {
    return std::index_sequence<static_cast<std::size_t>(Ints)...>();
  }

  /// Return the entries of x_in given in Ints... as smaller array x_out
  template <typename T>
  static void GetSubArray(T* x_out, const T* x_in) {
    for (int i = 0; i < size; i++) {
      x_out[i] = x_in[data_[i]];
    }
  }

  /// Return the entries of x_in as entries Ints... in larger array x_out.
  /// Entries of x_out not contained in Ints... are not modified.
  template <typename T>
  static void ExpandArray(T* x_out, const T* x_in) {
    for (int i = 0; i < size; i++) {
      x_out[data_[i]] = x_in[i];
    }
  }

  constexpr int operator[](const int i) const { return data_[i]; }
};

template <int... Ints>
constexpr std::array<int, sizeof...(Ints)> ConstexprArray<Ints...>::data_;

// Return an array from an integer sequence
template <int... Ints>
constexpr ConstexprArray<Ints...> MakeConstexprArray(
    const std::integer_sequence<int, Ints...>&) {
  return ConstexprArray<Ints...>();
}

#endif
