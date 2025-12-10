#ifndef __DATAFRAME_LIST__
#define __DATAFRAME_LIST__

// DataFrameCpp / ListCpp that uses ska::flat_hash_map (flat hash map)
// for the internal maps to improve lookup performance and cache locality.
// Uses FlatMatrix (column-major flattened std::vector<double>)
// for efficient R interop (single memcpy).

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
// FlatMatrix: contiguous column-major matrix representation
//
struct FlatMatrix {
  std::vector<double> data; // column-major: element (r,c) => data[c * nrow + r]
  int nrow = 0;
  int ncol = 0;
  
  FlatMatrix() = default;
  FlatMatrix(int nr, int nc) : data(static_cast<size_t>(nr) * static_cast<size_t>(nc)),
  nrow(nr), ncol(nc) {}
  
  FlatMatrix(std::vector<double>&& d, int nr, int nc)
    : data(std::move(d)), nrow(nr), ncol(nc) {
    if (static_cast<size_t>(nr) * static_cast<size_t>(nc) != data.size())
      throw std::runtime_error("FlatMatrix: data size mismatch with dimensions");
  }
  
  inline void resize(int nr, int nc) {
    nrow = nr; ncol = nc;
    data.resize(static_cast<size_t>(nr) * static_cast<size_t>(nc));
  }
  
  inline void fill(double v) { std::fill(data.begin(), data.end(), v); }
  
  inline bool empty() const noexcept { return data.empty() || nrow == 0 || ncol == 0; }
  
  inline size_t size() const noexcept { return data.size(); }
  
  // column-major index helper
  inline static size_t idx_col(int row, int col, int nrows) noexcept {
    return static_cast<size_t>(col) * static_cast<size_t>(nrows) + static_cast<size_t>(row);
  }
  
  inline double& operator()(int r, int c) {
    return data[idx_col(r, c, nrow)];
  }
  inline double operator()(int r, int c) const {
    return data[idx_col(r, c, nrow)];
  }
  
  // raw pointer accessors for parallel-friendly use
  inline const double* data_ptr() const noexcept { return data.empty() ? nullptr : data.data(); }
  inline double* data_ptr() noexcept { return data.empty() ? nullptr : data.data(); }
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
  
  inline static size_t idx_col(int row, int col, int nrows) noexcept {
    return FlatMatrix::idx_col(row, col, nrows);
  }
  
  size_t nrows() const {
    if (!names_.empty()) {
      const std::string& nm = names_.front();
      if (numeric_cols.count(nm)) return numeric_cols.at(nm).size();
      if (int_cols.count(nm)) return int_cols.at(nm).size();
      if (bool_cols.count(nm)) return bool_cols.at(nm).size();
      if (string_cols.count(nm)) return string_cols.at(nm).size();
    }
    return 0;
  }
  
  size_t size() const { return names_.size(); }
  const std::vector<std::string>& names() const { return names_; }
  
  bool containElementNamed(const std::string& name) const {
    return numeric_cols.count(name) || int_cols.count(name) ||
      bool_cols.count(name) || string_cols.count(name);
  }
  
  void check_row_size(size_t size, const std::string& name) const {
    size_t cur = nrows();
    if (cur > 0 && size != cur) throw std::runtime_error("Column '" + name + "' has inconsistent number of rows");
  }
  
  // Push single column overloads (double)
  void push_back(const std::vector<double>& col, const std::string& name) {
    if (containElementNamed(name)) throw std::runtime_error("Column '" + name + "' already exists.");
    check_row_size(col.size(), name);
    numeric_cols[name] = col;
    names_.push_back(name);
  }
  
  // rvalue overload: move into the map to avoid copy when caller can move
  void push_back(std::vector<double>&& col, const std::string& name) {
    if (containElementNamed(name)) throw std::runtime_error("Column '" + name + "' already exists.");
    check_row_size(col.size(), name);
    numeric_cols.emplace(name, std::move(col));
    names_.push_back(name);
  }
  
  // int column overloads (const& + &&)
  void push_back(const std::vector<int>& col, const std::string& name) {
    if (containElementNamed(name)) throw std::runtime_error("Column '" + name + "' already exists.");
    check_row_size(col.size(), name);
    int_cols[name] = col;
    names_.push_back(name);
  }
  void push_back(std::vector<int>&& col, const std::string& name) {
    if (containElementNamed(name)) throw std::runtime_error("Column '" + name + "' already exists.");
    check_row_size(col.size(), name);
    int_cols.emplace(name, std::move(col));
    names_.push_back(name);
  }
  
  // bool column overloads (const& + &&)
  void push_back(const std::vector<bool>& col, const std::string& name) {
    if (containElementNamed(name)) throw std::runtime_error("Column '" + name + "' already exists.");
    check_row_size(col.size(), name);
    bool_cols[name] = col;
    names_.push_back(name);
  }
  void push_back(std::vector<bool>&& col, const std::string& name) {
    if (containElementNamed(name)) throw std::runtime_error("Column '" + name + "' already exists.");
    check_row_size(col.size(), name);
    bool_cols.emplace(name, std::move(col));
    names_.push_back(name);
  }
  
  // string column overloads (const& + &&)
  void push_back(const std::vector<std::string>& col, const std::string& name) {
    if (containElementNamed(name)) throw std::runtime_error("Column '" + name + "' already exists.");
    check_row_size(col.size(), name);
    string_cols[name] = col;
    names_.push_back(name);
  }
  void push_back(std::vector<std::string>&& col, const std::string& name) {
    if (containElementNamed(name)) throw std::runtime_error("Column '" + name + "' already exists.");
    check_row_size(col.size(), name);
    string_cols.emplace(name, std::move(col));
    names_.push_back(name);
  }
  
  // Scalar expansions
  void push_back(double value, const std::string& name) {
    size_t cur = nrows();
    if (cur == 0 && !names_.empty()) throw std::runtime_error("Cannot push scalar when DataFrame has 0 rows");
    if (cur == 0) cur = 1;
    std::vector<double> col(cur, value);
    push_back(std::move(col), name);
  }
  void push_back(int value, const std::string& name) {
    size_t cur = nrows();
    if (cur == 0 && !names_.empty()) throw std::runtime_error("Cannot push scalar when DataFrame has 0 rows");
    if (cur == 0) cur = 1;
    std::vector<int> col(cur, value);
    push_back(std::move(col), name);
  }
  void push_back(bool value, const std::string& name) {
    size_t cur = nrows();
    if (cur == 0 && !names_.empty()) throw std::runtime_error("Cannot push scalar when DataFrame has 0 rows");
    if (cur == 0) cur = 1;
    std::vector<bool> col(cur, value);
    push_back(std::move(col), name);
  }
  void push_back(const std::string& value, const std::string& name) {
    size_t cur = nrows();
    if (cur == 0 && !names_.empty()) throw std::runtime_error("Cannot push scalar when DataFrame has 0 rows");
    if (cur == 0) cur = 1;
    std::vector<std::string> col(cur, value);
    push_back(std::move(col), name);
  }
  
  // Efficient: push_back_flat accepts a column-major flattened buffer containing nrows * p values
  // and will create p new columns named base_name, base_name.1, ..., base_name.p (if p>1)
  void push_back_flat(const std::vector<double>& flat_col_major, int nrows, const std::string& base_name) {
    if (flat_col_major.empty()) return;
    if (containElementNamed(base_name)) throw std::runtime_error("Column '" + base_name + "' already exists.");
    if (nrows <= 0) throw std::runtime_error("nrows must be > 0");
    if (flat_col_major.size() % static_cast<size_t>(nrows) != 0)
      throw std::runtime_error("flattened data size is not divisible by nrows");
    int p = static_cast<int>(flat_col_major.size() / static_cast<size_t>(nrows));
    check_row_size(static_cast<size_t>(nrows), base_name);
    
    if (p == 1) {
      std::vector<double> col(nrows);
      std::copy_n(flat_col_major.begin(), static_cast<size_t>(nrows), col.begin());
      numeric_cols[base_name] = std::move(col);
      names_.push_back(base_name);
    } else {
      for (int c = 0; c < p; ++c) {
        std::string col_name = base_name + "." + std::to_string(c + 1);
        if (containElementNamed(col_name)) throw std::runtime_error("Column '" + col_name + "' already exists.");
        std::vector<double> col(nrows);
        size_t offset = FlatMatrix::idx_col(0, c, nrows);
        for (int r = 0; r < nrows; ++r) {
          col[r] = flat_col_major[offset + static_cast<size_t>(r)];
        }
        numeric_cols.emplace(col_name, std::move(col));
        names_.push_back(col_name);
      }
    }
  }
  
  // Push front variants (const& + && for all types)
  void push_front(const std::vector<double>& col, const std::string& name) {
    if (containElementNamed(name)) throw std::runtime_error("Column '" + name + "' already exists.");
    check_row_size(col.size(), name);
    numeric_cols[name] = col;
    names_.insert(names_.begin(), name);
  }
  void push_front(std::vector<double>&& col, const std::string& name) {
    if (containElementNamed(name)) throw std::runtime_error("Column '" + name + "' already exists.");
    check_row_size(col.size(), name);
    numeric_cols.emplace(name, std::move(col));
    names_.insert(names_.begin(), name);
  }
  
  void push_front(const std::vector<int>& col, const std::string& name) {
    if (containElementNamed(name)) throw std::runtime_error("Column '" + name + "' already exists.");
    check_row_size(col.size(), name);
    int_cols[name] = col;
    names_.insert(names_.begin(), name);
  }
  void push_front(std::vector<int>&& col, const std::string& name) {
    if (containElementNamed(name)) throw std::runtime_error("Column '" + name + "' already exists.");
    check_row_size(col.size(), name);
    int_cols.emplace(name, std::move(col));
    names_.insert(names_.begin(), name);
  }
  
  void push_front(const std::vector<bool>& col, const std::string& name) {
    if (containElementNamed(name)) throw std::runtime_error("Column '" + name + "' already exists.");
    check_row_size(col.size(), name);
    bool_cols[name] = col;
    names_.insert(names_.begin(), name);
  }
  void push_front(std::vector<bool>&& col, const std::string& name) {
    if (containElementNamed(name)) throw std::runtime_error("Column '" + name + "' already exists.");
    check_row_size(col.size(), name);
    bool_cols.emplace(name, std::move(col));
    names_.insert(names_.begin(), name);
  }
  
  void push_front(const std::vector<std::string>& col, const std::string& name) {
    if (containElementNamed(name)) throw std::runtime_error("Column '" + name + "' already exists.");
    check_row_size(col.size(), name);
    string_cols[name] = col;
    names_.insert(names_.begin(), name);
  }
  void push_front(std::vector<std::string>&& col, const std::string& name) {
    if (containElementNamed(name)) throw std::runtime_error("Column '" + name + "' already exists.");
    check_row_size(col.size(), name);
    string_cols.emplace(name, std::move(col));
    names_.insert(names_.begin(), name);
  }
  
  // Erase column
  void erase(const std::string& name) {
    if (numeric_cols.erase(name)) { }
    else if (int_cols.erase(name)) { }
    else if (bool_cols.erase(name)) { }
    else string_cols.erase(name);
    names_.erase(std::remove(names_.begin(), names_.end(), name), names_.end());
  }
  
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
  int numeric_col_nrows(const std::string& name) const noexcept {
    auto it = numeric_cols.find(name);
    return it == numeric_cols.end() ? 0 : static_cast<int>(it->second.size());
  }
  
  // reserve map buckets (useful to avoid rehashing)
  void reserve_columns(size_t expected) {
    numeric_cols.reserve(expected);
    int_cols.reserve(expected);
    bool_cols.reserve(expected);
    string_cols.reserve(expected);
  }
};

//
// Forward declaration for recursive ListCpp usage
//
struct ListCpp;
using ListPtr = std::shared_ptr<ListCpp>;

//
// ListCpp: heterogeneous container using ska::flat_hash_map for storage
// Recursive occurrences of ListCpp are represented with std::shared_ptr<ListCpp>
// to avoid illegal by-value recursion.
//
struct ListCpp {
  std::vector<std::string> names_; // insertion order
  ska::flat_hash_map<std::string, std::variant<
    bool, int, double, std::string,
    std::vector<bool>,
    std::vector<int>,
    std::vector<double>,
    std::vector<std::string>,
    std::vector<std::vector<bool>>,
    std::vector<std::vector<int>>,
    std::vector<std::vector<double>>,
    std::vector<std::vector<std::string>>,
    FlatMatrix,
    DataFrameCpp,
    ListPtr,                                 // pointer to nested ListCpp
    std::vector<DataFrameCpp>,
    std::vector<ListPtr>                     // vector of pointers to nested ListCpp
    >> data;
    
    ListCpp() = default;
    
    size_t size() const { return data.size(); }
    
    bool containsElementNamed(const std::string& name) const { return data.count(name) > 0; }
    
    // Generic push_back for value-like types that match the variant alternatives
    template<typename T>
    void push_back(T value, const std::string& name) {
      if (containsElementNamed(name)) throw std::runtime_error("Element '" + name + "' already exists.");
      data.emplace(name, std::move(value));
      names_.push_back(name);
    }
    
    // Convenience overloads for adding ListCpp by value or by shared_ptr
    void push_back(const ListCpp& l, const std::string& name) {
      push_back(std::make_shared<ListCpp>(l), name);
    }
    void push_back(ListCpp&& l, const std::string& name) {
      push_back(std::make_shared<ListCpp>(std::move(l)), name);
    }
    void push_back(const ListPtr& p, const std::string& name) {
      if (containsElementNamed(name)) throw std::runtime_error("Element '" + name + "' already exists.");
      data.emplace(name, p);
      names_.push_back(name);
    }
    void push_back(ListPtr&& p, const std::string& name) {
      if (containsElementNamed(name)) throw std::runtime_error("Element '" + name + "' already exists.");
      data.emplace(name, std::move(p));
      names_.push_back(name);
    }
    
    // push_front variants
    template<typename T>
    void push_front(T value, const std::string& name) {
      if (containsElementNamed(name)) throw std::runtime_error("Element '" + name + "' already exists.");
      data.emplace(name, std::move(value));
      names_.insert(names_.begin(), name);
    }
    
    void push_front(const ListCpp& l, const std::string& name) {
      push_front(std::make_shared<ListCpp>(l), name);
    }
    void push_front(ListCpp&& l, const std::string& name) {
      push_front(std::make_shared<ListCpp>(std::move(l)), name);
    }
    void push_front(const ListPtr& p, const std::string& name) {
      if (containsElementNamed(name)) throw std::runtime_error("Element '" + name + "' already exists.");
      data.emplace(name, p);
      names_.insert(names_.begin(), name);
    }
    void push_front(ListPtr&& p, const std::string& name) {
      if (containsElementNamed(name)) throw std::runtime_error("Element '" + name + "' already exists.");
      data.emplace(name, std::move(p));
      names_.insert(names_.begin(), name);
    }
    
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
    ListCpp& get_list(const std::string& name) {
      if (!containsElementNamed(name)) throw std::runtime_error("Element with name '" + name + "' not found.");
      auto& var = data.at(name);
      if (auto p = std::get_if<ListPtr>(&var)) {
        if (!*p) throw std::runtime_error("List pointer is null for element '" + name + "'");
        return **p;
      }
      throw std::runtime_error("Element '" + name + "' is not a ListCpp");
    }
    const ListCpp& get_list(const std::string& name) const {
      if (!containsElementNamed(name)) throw std::runtime_error("Element with name '" + name + "' not found.");
      const auto& var = data.at(name);
      if (auto p = std::get_if<ListPtr>(&var)) {
        if (!*p) throw std::runtime_error("List pointer is null for element '" + name + "'");
        return **p;
      }
      throw std::runtime_error("Element '" + name + "' is not a ListCpp");
    }
    
    void erase(const std::string& name) {
      data.erase(name);
      names_.erase(std::remove(names_.begin(), names_.end(), name), names_.end());
    }
};

//
// Converters between R and C++ types
//

inline DataFrameCpp convertRDataFrameToCpp(const Rcpp::DataFrame& r_df) {
  DataFrameCpp df;
  Rcpp::CharacterVector cn = r_df.names();
  size_t nc = r_df.size();
  if (nc == 0) return df;
  for (size_t j = 0; j < nc; ++j) {
    std::string name = Rcpp::as<std::string>(cn[j]);
    Rcpp::RObject col = r_df[j];
    if (Rcpp::is<Rcpp::NumericVector>(col)) {
      Rcpp::NumericVector nv = col;
      df.push_back(Rcpp::as<std::vector<double>>(nv), name);
    } else if (Rcpp::is<Rcpp::IntegerVector>(col)) {
      Rcpp::IntegerVector iv = col;
      df.push_back(Rcpp::as<std::vector<int>>(iv), name);
    } else if (Rcpp::is<Rcpp::LogicalVector>(col)) {
      Rcpp::LogicalVector lv = col;
      df.push_back(Rcpp::as<std::vector<bool>>(lv), name);
    } else if (Rcpp::is<Rcpp::CharacterVector>(col)) {
      Rcpp::CharacterVector cv = col;
      df.push_back(Rcpp::as<std::vector<std::string>>(cv), name);
    } else {
      Rcpp::warning("Unsupported column type in DataFrame conversion: " + name);
    }
  }
  return df;
}

inline Rcpp::DataFrame convertDataFrameCppToR(const DataFrameCpp& df) {
  Rcpp::List cols;
  Rcpp::CharacterVector names;
  for (const auto& nm : df.names_) {
    names.push_back(nm);
    if (df.numeric_cols.count(nm)) {
      cols.push_back(Rcpp::wrap(df.numeric_cols.at(nm)));
    } else if (df.int_cols.count(nm)) {
      cols.push_back(Rcpp::wrap(df.int_cols.at(nm)));
    } else if (df.bool_cols.count(nm)) {
      cols.push_back(Rcpp::wrap(df.bool_cols.at(nm)));
    } else if (df.string_cols.count(nm)) {
      cols.push_back(Rcpp::wrap(df.string_cols.at(nm)));
    } else {
      cols.push_back(R_NilValue);
    }
  }
  cols.names() = names;
  return Rcpp::as<Rcpp::DataFrame>(cols);
}

//
// Visitor for ListCpp -> R conversion, with special handling for FlatMatrix and nested ListCpp
//
struct RcppVisitor {
  template <typename T>
  Rcpp::RObject operator()(const T& x) const {
    return Rcpp::wrap(x);
  }
  
  Rcpp::RObject operator()(const FlatMatrix& fm) const {
    if (fm.nrow <= 0 || fm.ncol <= 0) return R_NilValue;
    Rcpp::NumericMatrix M(fm.nrow, fm.ncol);
    std::memcpy(REAL(M), fm.data.data(), fm.data.size() * sizeof(double));
    return M;
  }
  
  Rcpp::RObject operator()(const DataFrameCpp& df) const {
    return convertDataFrameCppToR(df);
  }
  
  Rcpp::RObject operator()(const ListPtr& p) const {
    if (!p) return R_NilValue;
    return convertListCppToR(*p);
  }
  
  Rcpp::RObject operator()(const std::vector<DataFrameCpp>& dfs) const {
    Rcpp::List out;
    for (const auto& df : dfs) out.push_back(convertDataFrameCppToR(df));
    return out;
  }
  
  Rcpp::RObject operator()(const std::vector<ListPtr>& lists) const {
    Rcpp::List out;
    for (const auto& p : lists) {
      if (!p) out.push_back(R_NilValue);
      else out.push_back(convertListCppToR(*p));
    }
    return out;
  }
};

inline Rcpp::List convertListCppToR(const ListCpp& L) {
  Rcpp::List out;
  Rcpp::CharacterVector names;
  for (const auto& nm : L.names_) names.push_back(nm);
  out.names() = names;
  for (const auto& nm : L.names_) {
    const auto& var = L.data.at(nm);
    Rcpp::RObject obj = std::visit(RcppVisitor{}, var);
    out.push_back(obj);
  }
  return out;
}

namespace Rcpp {
template <> inline SEXP wrap(const DataFrameCpp& df) {
  return Rcpp::wrap(convertDataFrameCppToR(df));
}
template <> inline SEXP wrap(const ListCpp& l) {
  return Rcpp::wrap(convertListCppToR(l));
}
template <> inline SEXP wrap(const FlatMatrix& fm) {
  if (fm.nrow <= 0 || fm.ncol <= 0) return R_NilValue;
  Rcpp::NumericMatrix M(fm.nrow, fm.ncol);
  std::memcpy(REAL(M), fm.data.data(), fm.data.size() * sizeof(double));
  return Rcpp::wrap(M);
}
} // namespace Rcpp

#endif // __DATAFRAME_LIST__