#ifndef CONSTEXPR_ARRAY_H
#define CONSTEXPR_ARRAY_H
#include <utility>
#include <Eigen/Eigen>

template <int... Ints>
class ConstexprArray {
 private:
  static constexpr std::array<int, sizeof...(Ints)> data_ = {Ints...};

 public:
  constexpr ConstexprArray() {}

  using type = ConstexprArray<Ints...>;

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

  /// Get number of nonzero entries
  static constexpr int GetNonzeroEntries() {
    int n = 0;
    for (int i=0; i<data_.size(); ++i) {
      if (data_[i] != 0) n++;
    }
    return n;
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

  /// Sub-array of current array with indices SubInts...
  template <int... SubInts>
  using ConstexprSubArray = ConstexprArray<(data_[SubInts])...>;

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
    using expander = T[];
    (void)expander{(x_out[Ints] += *(x_in++))...};
  }

  // Function to get subvector from an Eigen vector
  template <int N>
  static Eigen::Matrix<double, size, 1> GetSubVector(
      const Eigen::Matrix<double, N, 1>& x) {
    Eigen::Matrix<double, size, 1> x_reduced;
    GetSubArray(x_reduced.data(), x.data());
    return x_reduced;
  }

  // Function to get subvector from an Eigen vector
  template <int N>
  static void GetSubVector(Eigen::Matrix<double, size, 1>* x_reduced,
                           const Eigen::Matrix<double, N, 1>& x) {
    GetSubArray(x_reduced->data(), x.data());
  }

  // Function to get submatrix from an Eigen matrix
  template <int N, typename T>
  static void GetSubMatrix(Eigen::Matrix<T, size, size>* x_reduced,
                           const Eigen::Matrix<T, N, N>& x) {
    T* new_data = x_reduced->data();
    const T* old_data = x.data();

    auto get_sub_row = [](T*& a, const T* b) {
      GetSubArray(a, b);
      a += size;
      return 0;
    };

    int dummy[]{get_sub_row(new_data, old_data + (N * Ints))...};
  }

  // Function to get submatrix from an Eigen matrix
  template <int N, typename T>
  static void GetSubMatrix(T* new_data, const T* old_data) {
    auto get_sub_row = [](T*& a, const T* b) {
      GetSubArray(a, b);
      a += size;
      return 0;
    };

    int dummy[]{get_sub_row(new_data, old_data + (N * Ints))...};
  }

  constexpr int operator[](const int i) const { return data_[i]; }
};

// Declaration of constexpr static member
template <int... Ints>
constexpr std::array<int, sizeof...(Ints)> ConstexprArray<Ints...>::data_;

// Function to return a ConstexprArray from an integer sequence
template <int... Ints>
constexpr ConstexprArray<Ints...> MakeConstexprArray(
    const std::integer_sequence<int, Ints...>&) {
  return ConstexprArray<Ints...>();
}

/// Class containing a list of ConstexprArrays
template <typename... Arrays>
class ConstexprArrayList {
 private:
  static constexpr auto list_ = std::tuple<Arrays...>();

 public:
  constexpr ConstexprArrayList() {}

  /// Number of control inputs
  static constexpr int size = static_cast<int>(sizeof...(Arrays));

  /// Get type of entry i
  template <int i>
  using type = typename std::tuple_element<i, std::tuple<Arrays...>>::type;

  /// Sum of all entries
  static constexpr int GetSum() {
    int sum = 0;
    int dummy[]{(AddEntry<Arrays>())...};
    for (int i = 0; i < size; i++) {
      sum += dummy[i];
    }
    return sum;
  }

  template <typename Array>
  static constexpr int AddEntry() {
    return Array::GetSum();
  }
};

#endif
