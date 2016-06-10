#ifndef CONSTEXPR_ARRAY_H
#define CONSTEXPR_ARRAY_H
#include <utility>

template <int... Ints>
class ConstexprArray {
 private:
  static constexpr std::array<int, sizeof...(Ints)> data_ = {Ints...};

 public:
  constexpr ConstexprArray() {}
  constexpr ConstexprArray(const std::integer_sequence<int, Ints...>& x) {}

  /// Number of control inputs
  static constexpr int n_inputs = sizeof...(Ints);

  /// Get sum of entries
  static constexpr int GetSum() {
    int sum = 0;
    int entry = 0;
    for (int i = 0; i < data_.size(); ++i) {
      entry = data_[i];
      sum += entry;
    }
    // sum += ConstexprArray<Ints...>::data_[i];
    return sum;
  }

  static constexpr int GetDelay(const int i) {
    return data_[i];
  }

  /// Override square brackets to get entry i
  // constexpr int operator[](const size_t i) const {
  // int x = 0;
  // x = data_[i];
  // return x;
  // }
  constexpr int operator[](const int i) const { return data_[i]; }
};

template <int... Ints>
constexpr std::array<int, sizeof...(Ints)> ConstexprArray<Ints...>::data_;

#endif
