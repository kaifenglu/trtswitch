#ifndef __DATAFRAME_LIST__
#define __DATAFRAME_LIST__

#include <algorithm>
#include <cstddef>
#include <cstring>     // std::memcpy
#include <stdexcept>
#include <string>
#include <type_traits>
#include <variant>
#include <vector>
#include <utility>     // std::move
#include <memory>      // std::shared_ptr

#include <Rcpp.h>

#include "ska/flat_hash_map.hpp"

//
// FlatMatrix: contiguous column-major matrix representation (double)
//
struct FlatMatrix {
  std::vector<double> data; // column-major: element (r,c) => data[c * nrow + r]
  std::size_t nrow = 0;
  std::size_t ncol = 0;
  
  FlatMatrix() = default;
  FlatMatrix(std::size_t nr, std::size_t nc) : data(nr * nc),
  nrow(nr), ncol(nc) {}
  
  FlatMatrix(std::vector<double>&& d, std::size_t nr, std::size_t nc)
    : data(std::move(d)), nrow(nr), ncol(nc) {
    if (nr * nc != data.size())
      throw std::runtime_error("FlatMatrix: data size mismatch with dimensions");
  }
  
  inline void resize(std::size_t nr, std::size_t nc) {
    nrow = nr; ncol = nc;
    data.resize(nr * nc);
  }
  
  inline void fill(double v) { std::fill(data.begin(), data.end(), v); }
  
  inline bool empty() const noexcept { return data.empty() || nrow == 0 || ncol == 0; }
  
  inline std::size_t size() const noexcept { return data.size(); }
  
  // column-major index helper (use size_t everywhere)
  inline static std::size_t idx_col(std::size_t row, std::size_t col, std::size_t nrows) noexcept {
    return col * nrows + row;
  }
  
  inline double& operator()(std::size_t r, std::size_t c) {
    return data[idx_col(r, c, nrow)];
  }
  inline double operator()(std::size_t r, std::size_t c) const {
    return data[idx_col(r, c, nrow)];
  }
  
  // raw pointer accessors for parallel-friendly use
  inline const double* data_ptr() const noexcept { return data.empty() ? nullptr : data.data(); }
  inline double* data_ptr() noexcept { return data.empty() ? nullptr : data.data(); }
};

//
// IntMatrix: contiguous column-major integer matrix representation (int)
//
struct IntMatrix {
  std::vector<int> data; // column-major: element (r,c) => data[c * nrow + r]
  std::size_t nrow = 0;
  std::size_t ncol = 0;
  
  IntMatrix() = default;
  IntMatrix(std::size_t nr, std::size_t nc) : data(nr * nc),
  nrow(nr), ncol(nc) {}
  
  IntMatrix(std::vector<int>&& d, std::size_t nr, std::size_t nc)
    : data(std::move(d)), nrow(nr), ncol(nc) {
    if (nr * nc != data.size())
      throw std::runtime_error("IntMatrix: data size mismatch with dimensions");
  }
  
  inline void resize(std::size_t nr, std::size_t nc) {
    nrow = nr; ncol = nc;
    data.resize(nr * nc);
  }
  
  inline void fill(int v) { std::fill(data.begin(), data.end(), v); }
  
  inline bool empty() const noexcept { return data.empty() || nrow == 0 || ncol == 0; }
  
  inline std::size_t size() const noexcept { return data.size(); }
  
  // column-major index helper
  inline static std::size_t idx_col(std::size_t row, std::size_t col, std::size_t nrows) noexcept {
    return col * nrows + row;
  }
  
  inline int& operator()(std::size_t r, std::size_t c) {
    return data[idx_col(r, c, nrow)];
  }
  inline int operator()(std::size_t r, std::size_t c) const {
    return data[idx_col(r, c, nrow)];
  }
  
  // raw pointer accessors for parallel-friendly use
  inline const int* data_ptr() const noexcept { return data.empty() ? nullptr : data.data(); }
  inline int* data_ptr() noexcept { return data.empty() ? nullptr : data.data(); }
};

//
// DataFrameCpp using ska::flat_hash_map for columns
//
struct DataFrameCpp {
  std::vector<std::string> names_; // insertion order
  
  // use ska::flat_hash_map for better performance than unordered_map
  ska::flat_hash_map<std::string, std::vector<double>> numeric_cols;
  ska::flat_hash_map<std::string, std::vector<int>> int_cols;
  ska::flat_hash_map<std::string, std::vector<bool>> bool_cols;
  ska::flat_hash_map<std::string, std::vector<std::string>> string_cols;
  
  DataFrameCpp() = default;
  
  inline static std::size_t idx_col(std::size_t row, std::size_t col, std::size_t nrows) noexcept {
    return FlatMatrix::idx_col(row, col, nrows);
  }
  
  std::size_t nrows() const {
    if (!names_.empty()) {
      const std::string& nm = names_.front();
      if (numeric_cols.count(nm)) return numeric_cols.at(nm).size();
      if (int_cols.count(nm)) return int_cols.at(nm).size();
      if (bool_cols.count(nm)) return bool_cols.at(nm).size();
      if (string_cols.count(nm)) return string_cols.at(nm).size();
    }
    return 0;
  }
  
  std::size_t size() const { return names_.size(); }
  const std::vector<std::string>& names() const { return names_; }
  
  bool containElementNamed(const std::string& name) const {
    return numeric_cols.count(name) || int_cols.count(name) ||
      bool_cols.count(name) || string_cols.count(name);
  }
  
  void check_row_size(std::size_t size, const std::string& name) const {
    std::size_t cur = nrows();
    if (cur > 0 && size != cur) throw std::runtime_error("Column '" + name + "' has inconsistent number of rows");
  }
  
  // Push single column overloads (double)
  void push_back(const std::vector<double>& col, const std::string& name);
  void push_back(std::vector<double>&& col, const std::string& name);
  
  // int column overloads (const& + &&)
  void push_back(const std::vector<int>& col, const std::string& name);
  void push_back(std::vector<int>&& col, const std::string& name);
  
  // bool column overloads (const& + &&)
  void push_back(const std::vector<bool>& col, const std::string& name);
  void push_back(std::vector<bool>&& col, const std::string& name);
  
  // string column overloads (const& + &&)
  void push_back(const std::vector<std::string>& col, const std::string& name);
  void push_back(std::vector<std::string>&& col, const std::string& name);
  
  // Scalar expansions
  void push_back(double value, const std::string& name);
  void push_back(int value, const std::string& name);
  void push_back(bool value, const std::string& name);
  void push_back(const std::string& value, const std::string& name);
  
  // Efficient: push_back_flat accepts a column-major flattened buffer containing nrows * p values
  // and will create p new columns named base_name, base_name.1, ..., base_name.p (if p>1)
  void push_back_flat(const std::vector<double>& flat_col_major, std::size_t nrows, const std::string& base_name);
  
  // Push front variants (const& + && for all types)
  void push_front(const std::vector<double>& col, const std::string& name);
  void push_front(std::vector<double>&& col, const std::string& name);
  
  void push_front(const std::vector<int>& col, const std::string& name);
  void push_front(std::vector<int>&& col, const std::string& name);
  
  void push_front(const std::vector<bool>& col, const std::string& name);
  void push_front(std::vector<bool>&& col, const std::string& name);
  
  void push_front(const std::vector<std::string>& col, const std::string& name);
  void push_front(std::vector<std::string>&& col, const std::string& name);
  
  // Scalar expansions for push_front (requested)
  void push_front(double value, const std::string& name);
  void push_front(int value, const std::string& name);
  void push_front(bool value, const std::string& name);
  void push_front(const std::string& value, const std::string& name);
  
  // Erase column
  void erase(const std::string& name);
  
  // Accessors with type checking
  template <typename T>
  std::vector<T>& get(const std::string& name) {
    if constexpr (std::is_same_v<T, double>) {
      if (numeric_cols.count(name)) return numeric_cols.at(name);
    } else if constexpr (std::is_same_v<T, int>) {
      if (int_cols.count(name)) return int_cols.at(name);
    } else if constexpr (std::is_same_v<T, bool>) {
      if (bool_cols.count(name)) return bool_cols.at(name);
    } else if constexpr (std::is_same_v<T, std::string>) {
      if (string_cols.count(name)) return string_cols.at(name);
    }
    throw std::runtime_error("Column '" + name + "' not found or type mismatch.");
  }
  
  template <typename T>
  const std::vector<T>& get(const std::string& name) const {
    if constexpr (std::is_same_v<T, double>) {
      if (numeric_cols.count(name)) return numeric_cols.at(name);
    } else if constexpr (std::is_same_v<T, int>) {
      if (int_cols.count(name)) return int_cols.at(name);
    } else if constexpr (std::is_same_v<T, bool>) {
      if (bool_cols.count(name)) return bool_cols.at(name);
    } else if constexpr (std::is_same_v<T, std::string>) {
      if (string_cols.count(name)) return string_cols.at(name);
    }
    throw std::runtime_error("Column '" + name + "' not found or type mismatch.");
  }
  
  // Convenience parallel-friendly accessors
  const double* numeric_col_ptr(const std::string& name) const noexcept {
    auto it = numeric_cols.find(name);
    return it == numeric_cols.end() ? nullptr : it->second.data();
  }
  std::size_t numeric_col_nrows(const std::string& name) const noexcept {
    auto it = numeric_cols.find(name);
    return it == numeric_cols.end() ? 0 : it->second.size();
  }
  
  const int* int_col_ptr(const std::string& name) const noexcept {
    auto it = int_cols.find(name);
    return it == int_cols.end() ? nullptr : it->second.data();
  }
  
  std::size_t int_col_nrows(const std::string& name) const noexcept {
    auto it = int_cols.find(name);
    return it == int_cols.end() ? 0 : it->second.size();
  }
  
  // reserve map buckets (useful to avoid rehashing)
  void reserve_columns(std::size_t expected);
};

//
// Forward declaration for recursive ListCpp usage
//
struct ListCpp;
using ListPtr = std::shared_ptr<ListCpp>;

//
// ListCpp: heterogeneous container using ska::flat_hash_map for storage
//
struct ListCpp {
  std::vector<std::string> names_; // insertion order
  ska::flat_hash_map<std::string, std::variant<
    bool, int, double, std::string,
    std::vector<bool>,
    std::vector<int>,
    std::vector<double>,
    std::vector<std::string>,
    FlatMatrix,
    IntMatrix,
    DataFrameCpp,
    ListPtr,                                 // pointer to nested ListCpp
    std::vector<DataFrameCpp>,
    std::vector<ListPtr>                     // vector of pointers to nested ListCpp
    >> data;
    
    ListCpp() = default;
    
    std::size_t size() const { return data.size(); }
    
    bool containsElementNamed(const std::string& name) const { return data.count(name) > 0; }
    
    // Generic push_back for value-like types that match the variant alternatives
    template<typename T>
    void push_back(T value, const std::string& name) {
      if (containsElementNamed(name)) throw std::runtime_error("Element '" + name + "' already exists.");
      data.emplace(name, std::move(value));
      names_.push_back(name);
    }
    
    // Convenience overloads for adding ListCpp by value or by shared_ptr
    void push_back(const ListCpp& l, const std::string& name);
    void push_back(ListCpp&& l, const std::string& name);
    void push_back(const ListPtr& p, const std::string& name);
    void push_back(ListPtr&& p, const std::string& name);
    
    // Convenience overloads for FlatMatrix / IntMatrix / DataFrameCpp
    void push_back(const FlatMatrix& fm, const std::string& name);
    void push_back(FlatMatrix&& fm, const std::string& name);
    void push_front(const FlatMatrix& fm, const std::string& name);
    void push_front(FlatMatrix&& fm, const std::string& name);
    
    void push_back(const IntMatrix& im, const std::string& name);
    void push_back(IntMatrix&& im, const std::string& name);
    void push_front(const IntMatrix& im, const std::string& name);
    void push_front(IntMatrix&& im, const std::string& name);
    
    void push_back(const DataFrameCpp& df, const std::string& name);
    void push_back(DataFrameCpp&& df, const std::string& name);
    void push_front(const DataFrameCpp& df, const std::string& name);
    void push_front(DataFrameCpp&& df, const std::string& name);
    
    // push_front variants
    template<typename T>
    void push_front(T value, const std::string& name) {
      if (containsElementNamed(name)) throw std::runtime_error("Element '" + name + "' already exists.");
      data.emplace(name, std::move(value));
      names_.insert(names_.begin(), name);
    }
    
    void push_front(const ListCpp& l, const std::string& name);
    void push_front(ListCpp&& l, const std::string& name);
    void push_front(const ListPtr& p, const std::string& name);
    void push_front(ListPtr&& p, const std::string& name);
    
    std::vector<std::string> names() const { return names_; }
    
    template <typename T>
    T& get(const std::string& name) {
      if (!containsElementNamed(name)) throw std::runtime_error("Element with name '" + name + "' not found.");
      return std::get<T>(data.at(name));
    }
    
    template <typename T>
    const T& get(const std::string& name) const {
      if (!containsElementNamed(name)) throw std::runtime_error("Element with name '" + name + "' not found.");
      return std::get<T>(data.at(name));
    }
    
    // Convenience: get nested ListCpp by reference (throws if not present)
    ListCpp& get_list(const std::string& name);
    const ListCpp& get_list(const std::string& name) const;
    
    void erase(const std::string& name);
};

// --------------------------- Converters between R and C++ types (declarations) ----

DataFrameCpp convertRDataFrameToCpp(const Rcpp::DataFrame& r_df);
Rcpp::DataFrame convertDataFrameCppToR(const DataFrameCpp& df);
Rcpp::List convertListCppToR(const ListCpp& L);

// Rcpp wrap specializations
namespace Rcpp {
template <> inline SEXP wrap(const DataFrameCpp& df) {
  return Rcpp::wrap(convertDataFrameCppToR(df));
}
template <> inline SEXP wrap(const ListCpp& l) {
  return Rcpp::wrap(convertListCppToR(l));
}
template <> inline SEXP wrap(const FlatMatrix& fm) {
  if (fm.nrow == 0 || fm.ncol == 0) return R_NilValue;
  Rcpp::NumericMatrix M(static_cast<int>(fm.nrow), static_cast<int>(fm.ncol));
  std::memcpy(REAL(M), fm.data.data(), fm.data.size() * sizeof(double));
  return Rcpp::wrap(M);
}
template <> inline SEXP wrap(const IntMatrix& im) {
  if (im.nrow == 0 || im.ncol == 0) return R_NilValue;
  Rcpp::IntegerMatrix M(static_cast<int>(im.nrow), static_cast<int>(im.ncol));
  std::memcpy(INTEGER(M), im.data.data(), im.data.size() * sizeof(int));
  return Rcpp::wrap(M);
}
} // namespace Rcpp

// Declarations for helpers implemented in dataframe_list.cpp
FlatMatrix flatmatrix_from_Rmatrix(const Rcpp::NumericMatrix& M);
IntMatrix intmatrix_from_Rmatrix(const Rcpp::IntegerMatrix& M);
DataFrameCpp dataframe_from_flatmatrix(const FlatMatrix& fm, const std::vector<std::string>& names = {});
FlatMatrix flatten_numeric_columns(const DataFrameCpp& df, const std::vector<std::string>& cols = {});
const double* numeric_column_ptr(const DataFrameCpp& df, const std::string& name) noexcept;
void move_numeric_column(DataFrameCpp& df, std::vector<double>&& col, const std::string& name);
const int* int_column_ptr(const DataFrameCpp& df, const std::string& name) noexcept;
void move_int_column(DataFrameCpp& df, std::vector<int>&& col, const std::string& name);
DataFrameCpp subset_rows(const DataFrameCpp& df, const std::vector<std::size_t>& row_idx);
std::shared_ptr<ListCpp> listcpp_from_rlist(const Rcpp::List& rlist);

#endif // __DATAFRAME_LIST__