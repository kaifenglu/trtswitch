#include "dataframe_list.h"

#include <cstring>   // std::memcpy
#include <algorithm>
#include <stdexcept>
#include <memory>

#include <Rcpp.h> // ensure Rcpp types available in this TU

using std::string;
using std::vector;

// --------------------------- DataFrameCpp members (small) -------------------

void DataFrameCpp::push_back(const std::vector<double>& col, const std::string& name) {
  if (containElementNamed(name)) throw std::runtime_error("Column '" + name + "' already exists.");
  check_row_size(col.size(), name);
  numeric_cols[name] = col;
  names_.push_back(name);
}
void DataFrameCpp::push_back(std::vector<double>&& col, const std::string& name) {
  if (containElementNamed(name)) throw std::runtime_error("Column '" + name + "' already exists.");
  check_row_size(col.size(), name);
  numeric_cols.emplace(name, std::move(col));
  names_.push_back(name);
}

void DataFrameCpp::push_back(const std::vector<int>& col, const std::string& name) {
  if (containElementNamed(name)) throw std::runtime_error("Column '" + name + "' already exists.");
  check_row_size(col.size(), name);
  int_cols[name] = col;
  names_.push_back(name);
}
void DataFrameCpp::push_back(std::vector<int>&& col, const std::string& name) {
  if (containElementNamed(name)) throw std::runtime_error("Column '" + name + "' already exists.");
  check_row_size(col.size(), name);
  int_cols.emplace(name, std::move(col));
  names_.push_back(name);
}

void DataFrameCpp::push_back(const std::vector<bool>& col, const std::string& name) {
  if (containElementNamed(name)) throw std::runtime_error("Column '" + name + "' already exists.");
  check_row_size(col.size(), name);
  bool_cols[name] = col;
  names_.push_back(name);
}
void DataFrameCpp::push_back(std::vector<bool>&& col, const std::string& name) {
  if (containElementNamed(name)) throw std::runtime_error("Column '" + name + "' already exists.");
  check_row_size(col.size(), name);
  bool_cols.emplace(name, std::move(col));
  names_.push_back(name);
}

void DataFrameCpp::push_back(const std::vector<std::string>& col, const std::string& name) {
  if (containElementNamed(name)) throw std::runtime_error("Column '" + name + "' already exists.");
  check_row_size(col.size(), name);
  string_cols[name] = col;
  names_.push_back(name);
}
void DataFrameCpp::push_back(std::vector<std::string>&& col, const std::string& name) {
  if (containElementNamed(name)) throw std::runtime_error("Column '" + name + "' already exists.");
  check_row_size(col.size(), name);
  string_cols.emplace(name, std::move(col));
  names_.push_back(name);
}

// Scalar expansions (push_back)
void DataFrameCpp::push_back(double value, const std::string& name) {
  size_t cur = nrows();
  if (cur == 0 && !names_.empty()) throw std::runtime_error("Cannot push scalar when DataFrame has 0 rows");
  if (cur == 0) cur = 1;
  std::vector<double> col(cur, value);
  push_back(std::move(col), name);
}
void DataFrameCpp::push_back(int value, const std::string& name) {
  size_t cur = nrows();
  if (cur == 0 && !names_.empty()) throw std::runtime_error("Cannot push scalar when DataFrame has 0 rows");
  if (cur == 0) cur = 1;
  std::vector<int> col(cur, value);
  push_back(std::move(col), name);
}
void DataFrameCpp::push_back(bool value, const std::string& name) {
  size_t cur = nrows();
  if (cur == 0 && !names_.empty()) throw std::runtime_error("Cannot push scalar when DataFrame has 0 rows");
  if (cur == 0) cur = 1;
  std::vector<bool> col(cur, value);
  push_back(std::move(col), name);
}
void DataFrameCpp::push_back(const std::string& value, const std::string& name) {
  size_t cur = nrows();
  if (cur == 0 && !names_.empty()) throw std::runtime_error("Cannot push scalar when DataFrame has 0 rows");
  if (cur == 0) cur = 1;
  std::vector<std::string> col(cur, value);
  push_back(std::move(col), name);
}

// Efficient: push_back_flat accepts a column-major flattened buffer containing nrows * p values
// and will create p new columns named base_name, base_name.1, ..., base_name.p (if p>1)
void DataFrameCpp::push_back_flat(const std::vector<double>& flat_col_major, int nrows, const std::string& base_name) {
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

// Push front variants (vectors)
void DataFrameCpp::push_front(const std::vector<double>& col, const std::string& name) {
  if (containElementNamed(name)) throw std::runtime_error("Column '" + name + "' already exists.");
  check_row_size(col.size(), name);
  numeric_cols[name] = col;
  names_.insert(names_.begin(), name);
}
void DataFrameCpp::push_front(std::vector<double>&& col, const std::string& name) {
  if (containElementNamed(name)) throw std::runtime_error("Column '" + name + "' already exists.");
  check_row_size(col.size(), name);
  numeric_cols.emplace(name, std::move(col));
  names_.insert(names_.begin(), name);
}
void DataFrameCpp::push_front(const std::vector<int>& col, const std::string& name) {
  if (containElementNamed(name)) throw std::runtime_error("Column '" + name + "' already exists.");
  check_row_size(col.size(), name);
  int_cols[name] = col;
  names_.insert(names_.begin(), name);
}
void DataFrameCpp::push_front(std::vector<int>&& col, const std::string& name) {
  if (containElementNamed(name)) throw std::runtime_error("Column '" + name + "' already exists.");
  check_row_size(col.size(), name);
  int_cols.emplace(name, std::move(col));
  names_.insert(names_.begin(), name);
}
void DataFrameCpp::push_front(const std::vector<bool>& col, const std::string& name) {
  if (containElementNamed(name)) throw std::runtime_error("Column '" + name + "' already exists.");
  check_row_size(col.size(), name);
  bool_cols[name] = col;
  names_.insert(names_.begin(), name);
}
void DataFrameCpp::push_front(std::vector<bool>&& col, const std::string& name) {
  if (containElementNamed(name)) throw std::runtime_error("Column '" + name + "' already exists.");
  check_row_size(col.size(), name);
  bool_cols.emplace(name, std::move(col));
  names_.insert(names_.begin(), name);
}
void DataFrameCpp::push_front(const std::vector<std::string>& col, const std::string& name) {
  if (containElementNamed(name)) throw std::runtime_error("Column '" + name + "' already exists.");
  check_row_size(col.size(), name);
  string_cols[name] = col;
  names_.insert(names_.begin(), name);
}
void DataFrameCpp::push_front(std::vector<std::string>&& col, const std::string& name) {
  if (containElementNamed(name)) throw std::runtime_error("Column '" + name + "' already exists.");
  check_row_size(col.size(), name);
  string_cols.emplace(name, std::move(col));
  names_.insert(names_.begin(), name);
}

// Scalar expansions for push_front (new implementations)
void DataFrameCpp::push_front(double value, const std::string& name) {
  size_t cur = nrows();
  if (cur == 0 && !names_.empty()) throw std::runtime_error("Cannot push scalar when DataFrame has 0 rows");
  if (cur == 0) cur = 1;
  std::vector<double> col(cur, value);
  push_front(std::move(col), name);
}
void DataFrameCpp::push_front(int value, const std::string& name) {
  size_t cur = nrows();
  if (cur == 0 && !names_.empty()) throw std::runtime_error("Cannot push scalar when DataFrame has 0 rows");
  if (cur == 0) cur = 1;
  std::vector<int> col(cur, value);
  push_front(std::move(col), name);
}
void DataFrameCpp::push_front(bool value, const std::string& name) {
  size_t cur = nrows();
  if (cur == 0 && !names_.empty()) throw std::runtime_error("Cannot push scalar when DataFrame has 0 rows");
  if (cur == 0) cur = 1;
  std::vector<bool> col(cur, value);
  push_front(std::move(col), name);
}
void DataFrameCpp::push_front(const std::string& value, const std::string& name) {
  size_t cur = nrows();
  if (cur == 0 && !names_.empty()) throw std::runtime_error("Cannot push scalar when DataFrame has 0 rows");
  if (cur == 0) cur = 1;
  std::vector<std::string> col(cur, value);
  push_front(std::move(col), name);
}

void DataFrameCpp::erase(const std::string& name) {
  if (numeric_cols.erase(name)) { }
  else if (int_cols.erase(name)) { }
  else if (bool_cols.erase(name)) { }
  else string_cols.erase(name);
  names_.erase(std::remove(names_.begin(), names_.end(), name), names_.end());
}

void DataFrameCpp::reserve_columns(size_t expected) {
  numeric_cols.reserve(expected);
  int_cols.reserve(expected);
  bool_cols.reserve(expected);
  string_cols.reserve(expected);
}



// --------------------------- ListCpp small members --------------------------

void ListCpp::push_back(const ListCpp& l, const std::string& name) {
  push_back(std::make_shared<ListCpp>(l), name);
}
void ListCpp::push_back(ListCpp&& l, const std::string& name) {
  push_back(std::make_shared<ListCpp>(std::move(l)), name);
}
void ListCpp::push_back(const ListPtr& p, const std::string& name) {
  if (containsElementNamed(name)) throw std::runtime_error("Element '" + name + "' already exists.");
  data.emplace(name, p);
  names_.push_back(name);
}
void ListCpp::push_back(ListPtr&& p, const std::string& name) {
  if (containsElementNamed(name)) throw std::runtime_error("Element '" + name + "' already exists.");
  data.emplace(name, std::move(p));
  names_.push_back(name);
}

// FlatMatrix overloads
void ListCpp::push_back(const FlatMatrix& fm, const std::string& name) {
  if (containsElementNamed(name)) throw std::runtime_error("Element '" + name + "' already exists.");
  data.emplace(name, fm);
  names_.push_back(name);
}
void ListCpp::push_back(FlatMatrix&& fm, const std::string& name) {
  if (containsElementNamed(name)) throw std::runtime_error("Element '" + name + "' already exists.");
  data.emplace(name, std::move(fm));
  names_.push_back(name);
}
void ListCpp::push_front(const FlatMatrix& fm, const std::string& name) {
  if (containsElementNamed(name)) throw std::runtime_error("Element '" + name + "' already exists.");
  data.emplace(name, fm);
  names_.insert(names_.begin(), name);
}
void ListCpp::push_front(FlatMatrix&& fm, const std::string& name) {
  if (containsElementNamed(name)) throw std::runtime_error("Element '" + name + "' already exists.");
  data.emplace(name, std::move(fm));
  names_.insert(names_.begin(), name);
}

// IntMatrix overloads
void ListCpp::push_back(const IntMatrix& im, const std::string& name) {
  if (containsElementNamed(name)) throw std::runtime_error("Element '" + name + "' already exists.");
  data.emplace(name, im);
  names_.push_back(name);
}
void ListCpp::push_back(IntMatrix&& im, const std::string& name) {
  if (containsElementNamed(name)) throw std::runtime_error("Element '" + name + "' already exists.");
  data.emplace(name, std::move(im));
  names_.push_back(name);
}
void ListCpp::push_front(const IntMatrix& im, const std::string& name) {
  if (containsElementNamed(name)) throw std::runtime_error("Element '" + name + "' already exists.");
  data.emplace(name, im);
  names_.insert(names_.begin(), name);
}
void ListCpp::push_front(IntMatrix&& im, const std::string& name) {
  if (containsElementNamed(name)) throw std::runtime_error("Element '" + name + "' already exists.");
  data.emplace(name, std::move(im));
  names_.insert(names_.begin(), name);
}

// DataFrameCpp overloads
void ListCpp::push_back(const DataFrameCpp& df, const std::string& name) {
  if (containsElementNamed(name)) throw std::runtime_error("Element '" + name + "' already exists.");
  data.emplace(name, df);
  names_.push_back(name);
}
void ListCpp::push_back(DataFrameCpp&& df, const std::string& name) {
  if (containsElementNamed(name)) throw std::runtime_error("Element '" + name + "' already exists.");
  data.emplace(name, std::move(df));
  names_.push_back(name);
}
void ListCpp::push_front(const DataFrameCpp& df, const std::string& name) {
  if (containsElementNamed(name)) throw std::runtime_error("Element '" + name + "' already exists.");
  data.emplace(name, df);
  names_.insert(names_.begin(), name);
}
void ListCpp::push_front(DataFrameCpp&& df, const std::string& name) {
  if (containsElementNamed(name)) throw std::runtime_error("Element '" + name + "' already exists.");
  data.emplace(name, std::move(df));
  names_.insert(names_.begin(), name);
}

void ListCpp::push_front(const ListCpp& l, const std::string& name) {
  push_front(std::make_shared<ListCpp>(l), name);
}
void ListCpp::push_front(ListCpp&& l, const std::string& name) {
  push_front(std::make_shared<ListCpp>(std::move(l)), name);
}
void ListCpp::push_front(const ListPtr& p, const std::string& name) {
  if (containsElementNamed(name)) throw std::runtime_error("Element '" + name + "' already exists.");
  data.emplace(name, p);
  names_.insert(names_.begin(), name);
}
void ListCpp::push_front(ListPtr&& p, const std::string& name) {
  if (containsElementNamed(name)) throw std::runtime_error("Element '" + name + "' already exists.");
  data.emplace(name, std::move(p));
  names_.insert(names_.begin(), name);
}

ListCpp& ListCpp::get_list(const std::string& name) {
  if (!containsElementNamed(name)) throw std::runtime_error("Element with name '" + name + "' not found.");
  auto& var = data.at(name);
  if (auto p = std::get_if<ListPtr>(&var)) {
    if (!*p) throw std::runtime_error("List pointer is null for element '" + name + "'");
    return **p;
  }
  throw std::runtime_error("Element '" + name + "' is not a ListCpp");
}
const ListCpp& ListCpp::get_list(const std::string& name) const {
  if (!containsElementNamed(name)) throw std::runtime_error("Element with name '" + name + "' not found.");
  const auto& var = data.at(name);
  if (auto p = std::get_if<ListPtr>(&var)) {
    if (!*p) throw std::runtime_error("List pointer is null for element '" + name + "'");
    return **p;
  }
  throw std::runtime_error("Element '" + name + "' is not a ListCpp");
}

void ListCpp::erase(const std::string& name) {
  data.erase(name);
  names_.erase(std::remove(names_.begin(), names_.end(), name), names_.end());
}

// --------------------------- Converters implementations ---------------------

// convert R data.frame -> DataFrameCpp
DataFrameCpp convertRDataFrameToCpp(const Rcpp::DataFrame& r_df) {
  DataFrameCpp df;
  Rcpp::CharacterVector cn = r_df.names();
  R_xlen_t nc = r_df.size();
  if (nc == 0) return df;
  
  for (R_xlen_t j = 0; j < nc; ++j) {
    std::string name;
    if (cn.size() > 0 && static_cast<R_xlen_t>(cn.size()) > j && !Rcpp::StringVector::is_na(cn[j]) && std::string(cn[j]) != "") {
      name = Rcpp::as<std::string>(cn[j]);
    } else {
      name = "V" + std::to_string(static_cast<long long>(j + 1));
    }
    
    Rcpp::RObject col = r_df[j];
    if (Rcpp::is<Rcpp::NumericVector>(col)) {
      df.push_back(Rcpp::as<std::vector<double>>(col), name);
    } else if (Rcpp::is<Rcpp::IntegerVector>(col)) {
      df.push_back(Rcpp::as<std::vector<int>>(col), name);
    } else if (Rcpp::is<Rcpp::LogicalVector>(col)) {
      df.push_back(Rcpp::as<std::vector<bool>>(col), name);
    } else if (Rcpp::is<Rcpp::CharacterVector>(col)) {
      df.push_back(Rcpp::as<std::vector<std::string>>(col), name);
    } else {
      Rcpp::warning("Unsupported column type in DataFrame conversion: " + name);
    }
  }
  return df;
}

// convert DataFrameCpp -> R data.frame
Rcpp::DataFrame convertDataFrameCppToR(const DataFrameCpp& df) {
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

// --------------------------- Visitor (complete) -----------------------------
// Define RcppVisitor before convertListCppToR so the type is complete where used.
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
  
  Rcpp::RObject operator()(const IntMatrix& im) const {
    if (im.nrow <= 0 || im.ncol <= 0) return R_NilValue;
    Rcpp::IntegerMatrix M(im.nrow, im.ncol);
    std::memcpy(INTEGER(M), im.data.data(), im.data.size() * sizeof(int));
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

// convert ListCpp -> R list
Rcpp::List convertListCppToR(const ListCpp& L) {
  size_t p = L.names_.size();
  Rcpp::List out(static_cast<int>(p));
  Rcpp::CharacterVector names(static_cast<int>(p));
  RcppVisitor visitor; // now the type is complete
  for (size_t i = 0; i < p; ++i) {
    const auto& nm = L.names_[i];
    const auto& var = L.data.at(nm);
    Rcpp::RObject obj = std::visit(visitor, var);
    out[static_cast<int>(i)] = obj;
    names[static_cast<int>(i)] = nm;
  }
  out.names() = names;
  return out;
}

// --------------------------- R <-> FlatMatrix / IntMatrix / DataFrame helpers ------------

// Convert Rcpp::NumericMatrix to FlatMatrix (efficient memcpy when possible)
FlatMatrix flatmatrix_from_Rmatrix(const Rcpp::NumericMatrix& M) {
  int nr = M.nrow(), nc = M.ncol();
  if (nr == 0 || nc == 0) return FlatMatrix();
  FlatMatrix fm(nr, nc);
  std::memcpy(fm.data.data(), REAL(M), static_cast<size_t>(nr) * static_cast<size_t>(nc) * sizeof(double));
  return fm;
}

// Convert Rcpp::IntegerMatrix to IntMatrix (efficient memcpy)
IntMatrix intmatrix_from_Rmatrix(const Rcpp::IntegerMatrix& M) {
  int nr = M.nrow(), nc = M.ncol();
  if (nr == 0 || nc == 0) return IntMatrix();
  IntMatrix im(nr, nc);
  std::memcpy(im.data.data(), INTEGER(M), static_cast<size_t>(nr) * static_cast<size_t>(nc) * sizeof(int));
  return im;
}

// Build a DataFrameCpp from a FlatMatrix
DataFrameCpp dataframe_from_flatmatrix(const FlatMatrix& fm, const std::vector<std::string>& names) {
  DataFrameCpp out;
  if (fm.ncol <= 0 || fm.nrow < 0) return out;
  if (!names.empty() && static_cast<int>(names.size()) != fm.ncol) throw std::invalid_argument("dataframe_from_flatmatrix: names size must equal fm.ncol");
  for (int c = 0; c < fm.ncol; ++c) {
    std::vector<double> col(static_cast<size_t>(fm.nrow));
    size_t offset = FlatMatrix::idx_col(0, c, fm.nrow);
    std::memcpy(col.data(), fm.data.data() + offset, static_cast<size_t>(fm.nrow) * sizeof(double));
    std::string nm = names.empty() ? ("V" + std::to_string(c + 1)) : names[static_cast<size_t>(c)];
    out.push_back(std::move(col), nm);
  }
  return out;
}

// Flatten numeric columns from DataFrameCpp into a contiguous FlatMatrix
FlatMatrix flatten_numeric_columns(const DataFrameCpp& df, const std::vector<std::string>& cols_in) {
  std::vector<std::string> cols;
  if (cols_in.empty()) {
    // include both numeric (double) and integer columns, and booleans
    for (const auto& nm : df.names()) {
      if (df.numeric_cols.count(nm) || df.int_cols.count(nm) || df.bool_cols.count(nm)) cols.push_back(nm);
    }
  } else cols = cols_in;
  if (cols.empty()) return FlatMatrix();
  int ncol = static_cast<int>(cols.size()), nrow = static_cast<int>(df.nrows());
  if (nrow == 0) return FlatMatrix();
  FlatMatrix out(nrow, ncol);
  for (int c = 0; c < ncol; ++c) {
    const std::string& name = cols[static_cast<size_t>(c)];
    size_t offset = FlatMatrix::idx_col(0, c, nrow);
    
    if (df.numeric_cols.count(name)) {
      const std::vector<double>& src = df.numeric_cols.at(name);
      if (static_cast<int>(src.size()) != nrow) throw std::runtime_error("flatten_numeric_columns: inconsistent row count for column '" + name + "'");
      std::memcpy(out.data.data() + offset, src.data(), static_cast<size_t>(nrow) * sizeof(double));
    } else if (df.int_cols.count(name)) {
      const std::vector<int>& src = df.int_cols.at(name);
      if (static_cast<int>(src.size()) != nrow) throw std::runtime_error("flatten_numeric_columns: inconsistent row count for column '" + name + "'");
      for (int r = 0; r < nrow; ++r) {
        out.data[offset + static_cast<size_t>(r)] = static_cast<double>(src[static_cast<size_t>(r)]);
      }
    } else if (df.bool_cols.count(name)) {
      const std::vector<bool>& src = df.bool_cols.at(name);
      if (static_cast<int>(src.size()) != nrow) throw std::runtime_error("flatten_numeric_columns: inconsistent row count for column '" + name + "'");
      for (int r = 0; r < nrow; ++r) {
        out.data[offset + static_cast<size_t>(r)] = src[static_cast<size_t>(r)] ? 1.0 : 0.0;
      }
    } else {
      throw std::runtime_error("flatten_numeric_columns: column '" + name + "' not found or not numeric/integer/bool");
    }
  }
  return out;
}

// numeric_column_ptr
const double* numeric_column_ptr(const DataFrameCpp& df, const std::string& name) noexcept {
  return df.numeric_col_ptr(name);
}

// int_column_ptr
const int* int_column_ptr(const DataFrameCpp& df, const std::string& name) noexcept {
  return df.int_col_ptr(name);
}

// move_numeric_column
void move_numeric_column(DataFrameCpp& df, std::vector<double>&& col, const std::string& name) {
  df.push_back(std::move(col), name);
}

// move_int_column
void move_int_column(DataFrameCpp& df, std::vector<int>&& col, const std::string& name) {
  df.push_back(std::move(col), name);
}

// subset_rows
DataFrameCpp subset_rows(const DataFrameCpp& df, const std::vector<int>& row_idx) {
  DataFrameCpp out;
  if (row_idx.empty()) {
    for (const auto& nm : df.names()) {
      if (df.numeric_cols.count(nm)) out.push_back(std::vector<double>{}, nm);
      else if (df.int_cols.count(nm)) out.push_back(std::vector<int>{}, nm);
      else if (df.bool_cols.count(nm)) out.push_back(std::vector<bool>{}, nm);
      else if (df.string_cols.count(nm)) out.push_back(std::vector<std::string>{}, nm);
    }
    return out;
  }
  int max_index = *std::max_element(row_idx.begin(), row_idx.end());
  if (max_index < 0 || static_cast<size_t>(max_index) >= df.nrows()) throw std::runtime_error("subset_rows: row index out of range");
  for (const auto& nm : df.names()) {
    if (df.numeric_cols.count(nm)) {
      const auto& src = df.numeric_cols.at(nm);
      std::vector<double> dest; dest.reserve(row_idx.size());
      for (int idx : row_idx) dest.push_back(src[static_cast<size_t>(idx)]);
      out.push_back(std::move(dest), nm);
    } else if (df.int_cols.count(nm)) {
      const auto& src = df.int_cols.at(nm);
      std::vector<int> dest; dest.reserve(row_idx.size());
      for (int idx : row_idx) dest.push_back(src[static_cast<size_t>(idx)]);
      out.push_back(std::move(dest), nm);
    } else if (df.bool_cols.count(nm)) {
      const auto& src = df.bool_cols.at(nm);
      std::vector<bool> dest; dest.reserve(row_idx.size());
      for (int idx : row_idx) dest.push_back(src[static_cast<size_t>(idx)]);
      out.push_back(std::move(dest), nm);
    } else if (df.string_cols.count(nm)) {
      const auto& src = df.string_cols.at(nm);
      std::vector<std::string> dest; dest.reserve(row_idx.size());
      for (int idx : row_idx) dest.push_back(src[static_cast<size_t>(idx)]);
      out.push_back(std::move(dest), nm);
    }
  }
  return out;
}

// --------------------------- ListCpp construction from R list -----------------
// Helper: use Rcpp::as<T>(obj) to perform conversions from RObject safely.

static bool convert_r_object_to_listcpp_variant(ListCpp& target, const std::string& name, const Rcpp::RObject& obj) {
  // NumericMatrix -> FlatMatrix
  if (Rcpp::is<Rcpp::NumericMatrix>(obj)) {
    Rcpp::NumericMatrix M = Rcpp::as<Rcpp::NumericMatrix>(obj);
    FlatMatrix fm = flatmatrix_from_Rmatrix(M);
    target.push_back(std::move(fm), name);
    return true;
  }
  
  // IntegerMatrix -> IntMatrix
  if (Rcpp::is<Rcpp::IntegerMatrix>(obj)) {
    Rcpp::IntegerMatrix M = Rcpp::as<Rcpp::IntegerMatrix>(obj);
    IntMatrix im = intmatrix_from_Rmatrix(M);
    target.push_back(std::move(im), name);
    return true;
  }
  
  // DataFrame -> DataFrameCpp
  if (Rcpp::is<Rcpp::DataFrame>(obj)) {
    Rcpp::DataFrame rdf = Rcpp::as<Rcpp::DataFrame>(obj);
    DataFrameCpp cppdf = convertRDataFrameToCpp(rdf);
    target.push_back(std::move(cppdf), name);
    return true;
  }
  
  // List -> nested ListCpp (as shared_ptr)
  if (Rcpp::is<Rcpp::List>(obj)) {
    Rcpp::List rl = Rcpp::as<Rcpp::List>(obj);
    std::shared_ptr<ListCpp> nested = listcpp_from_rlist(rl);
    target.push_back(std::move(nested), name);
    return true;
  }
  
  // NumericVector -> std::vector<double>
  if (Rcpp::is<Rcpp::NumericVector>(obj)) {
    Rcpp::NumericVector nv = Rcpp::as<Rcpp::NumericVector>(obj);
    std::vector<double> v = Rcpp::as<std::vector<double>>(nv);
    target.push_back(std::move(v), name);
    return true;
  }
  
  // IntegerVector -> std::vector<int>
  if (Rcpp::is<Rcpp::IntegerVector>(obj)) {
    Rcpp::IntegerVector iv = Rcpp::as<Rcpp::IntegerVector>(obj);
    std::vector<int> v = Rcpp::as<std::vector<int>>(iv);
    target.push_back(std::move(v), name);
    return true;
  }
  
  // LogicalVector -> std::vector<bool>
  if (Rcpp::is<Rcpp::LogicalVector>(obj)) {
    Rcpp::LogicalVector lv = Rcpp::as<Rcpp::LogicalVector>(obj);
    std::vector<bool> v = Rcpp::as<std::vector<bool>>(lv);
    target.push_back(std::move(v), name);
    return true;
  }
  
  // CharacterVector -> std::vector<string>
  if (Rcpp::is<Rcpp::CharacterVector>(obj)) {
    Rcpp::CharacterVector cv = Rcpp::as<Rcpp::CharacterVector>(obj);
    std::vector<std::string> v = Rcpp::as<std::vector<std::string>>(cv);
    target.push_back(std::move(v), name);
    return true;
  }
  
  // Fallback: unsupported type
  return false;
}

std::shared_ptr<ListCpp> listcpp_from_rlist(const Rcpp::List& rlist) {
  auto out = std::make_shared<ListCpp>();
  Rcpp::CharacterVector rn = rlist.names();
  R_xlen_t n = rlist.size();
  for (R_xlen_t i = 0; i < n; ++i) {
    std::string name;
    if (rn.size() > 0 && static_cast<R_xlen_t>(rn.size()) > i && !Rcpp::StringVector::is_na(rn[i]) && std::string(rn[i]) != "")
      name = Rcpp::as<std::string>(rn[i]);
    else
      name = "V" + std::to_string(static_cast<long long>(i + 1));
    Rcpp::RObject obj = rlist[i];
    bool ok = convert_r_object_to_listcpp_variant(*out, name, obj);
    if (!ok) Rcpp::warning("listcpp_from_rlist: skipping unsupported element type for name '" + name + "'");
  }
  return out;
}