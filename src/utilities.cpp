#include "utilities.h"
#include "dataframe_list.h"

#include <Rcpp.h>

#include <algorithm>  // fill, lower_bound, max_element, memmove, 
// min_element, sort, swap, upper_bound
#include <cmath>      // copysign, exp, fabs, isinf, isnan, log, pow, sqrt
#include <cstddef>    // size_t
#include <cstring>    // memcpy
#include <functional> // function
#include <limits>     // numeric_limits
#include <memory>     // make_shared, shared_ptr
#include <numeric>    // accumulate, inner_product, iota 
#include <sstream>    // ostringstream
#include <stdexcept>  // invalid_argument, runtime_error
#include <string>     // string
#include <vector>     // vector

#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/logistic.hpp>
#include <boost/math/distributions/extreme_value.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/students_t.hpp>

using std::size_t;

// --------------------------- Distribution helpers --------------------------

double boost_pnorm(double q, double mean, double sd, bool lower_tail) {
  if (std::isnan(q)) return std::numeric_limits<double>::quiet_NaN();
  if (sd <= 0) throw std::invalid_argument("Standard deviation must be positive.");
  
  double z = (q - mean) / sd;
  if (lower_tail) {
    if (z <= -EXTREME_Z) return 0.0;
    if (z >= EXTREME_Z) return 1.0;
  } else {
    if (z >= EXTREME_Z) return 0.0;
    if (z <= -EXTREME_Z) return 1.0;
  }
  
  boost::math::normal_distribution<> dist(mean, sd);
  if (lower_tail) return boost::math::cdf(dist, q);
  else return boost::math::cdf(boost::math::complement(dist, q));
}

double boost_qnorm(double p, double mean, double sd, bool lower_tail) {
  if (std::isnan(p)) return std::numeric_limits<double>::quiet_NaN();
  if (sd <= 0) throw std::invalid_argument("Standard deviation must be positive.");
  if (p < 0.0 || p > 1.0) throw std::invalid_argument(
      "Probability must be between 0 and 1.");
  
  // Clamp extreme probabilities to avoid overflow in boost::math::quantile
  bool at_extreme = false;
  double extreme_quantile = 0.0;
  
  if (lower_tail) {
    if (p <= MIN_PROB) {
      at_extreme = true;
      extreme_quantile = MIN_NORMAL_QUANTILE;
    } else if (p >= MAX_PROB) {
      at_extreme = true;
      extreme_quantile = MAX_NORMAL_QUANTILE;
    }
  } else {
    // When lower_tail = false, we compute quantile(1 - p)
    if (p <= MIN_PROB) {
      at_extreme = true;
      extreme_quantile = MAX_NORMAL_QUANTILE;
    } else if (p >= MAX_PROB) {
      at_extreme = true;
      extreme_quantile = MIN_NORMAL_QUANTILE;
    }
  }
  
  // If at extreme, return scaled value directly
  if (at_extreme) {
    return mean + sd * extreme_quantile;
  }
  
  // Safe to call boost quantile
  boost::math::normal_distribution<> dist(mean, sd);
  return lower_tail ? boost::math::quantile(dist, p) :
    boost::math::quantile(dist, 1.0 - p);
}

double boost_dnorm(double x, double mean, double sd) {
  if (std::isnan(x)) return std::numeric_limits<double>::quiet_NaN();
  if (sd <= 0) throw std::invalid_argument("Standard deviation must be positive.");
  boost::math::normal_distribution<> dist(mean, sd);
  return boost::math::pdf(dist, x);
}

double boost_plogis(double q, double location, double scale, bool lower_tail) {
  if (std::isnan(q)) return std::numeric_limits<double>::quiet_NaN();
  if (scale <= 0) throw std::invalid_argument("Scale must be positive.");
  boost::math::logistic_distribution<> dist(location, scale);
  if (lower_tail) return boost::math::cdf(dist, q);
  else return boost::math::cdf(boost::math::complement(dist, q));
}

double boost_qlogis(double p, double location, double scale, bool lower_tail) {
  if (std::isnan(p)) return std::numeric_limits<double>::quiet_NaN();
  if (scale <= 0) throw std::invalid_argument("Scale must be positive.");
  if (p < 0.0 || p > 1.0) throw std::invalid_argument(
      "Probability must be between 0 and 1.");
  boost::math::logistic_distribution<> dist(location, scale);
  return lower_tail ? boost::math::quantile(dist, p) : 
    boost::math::quantile(dist, 1.0 - p);
}

double boost_dlogis(double x, double location, double scale) {
  if (std::isnan(x)) return std::numeric_limits<double>::quiet_NaN();
  if (scale <= 0) throw std::invalid_argument("Scale must be positive.");
  boost::math::logistic_distribution<> dist(location, scale);
  return boost::math::pdf(dist, x);
}

double boost_pextreme(double q, double location, double scale, bool lower_tail) {
  if (std::isnan(q)) return std::numeric_limits<double>::quiet_NaN();
  if (scale <= 0) throw std::invalid_argument("Scale must be positive.");
  boost::math::extreme_value_distribution<> dist(location, scale);
  // keep semantics consistent with complementary log-log link
  if (lower_tail) return boost::math::cdf(complement(dist, 2.0 * location - q));
  else return boost::math::cdf(dist, 2.0 * location - q);
}

double boost_qextreme(double p, double location, double scale, bool lower_tail) {
  if (std::isnan(p)) return std::numeric_limits<double>::quiet_NaN();
  if (scale <= 0) throw std::invalid_argument("Scale must be positive.");
  if (p < 0.0 || p > 1.0) throw std::invalid_argument(
      "Probability must be between 0 and 1.");
  boost::math::extreme_value_distribution<> dist(location, scale);
  return lower_tail? -boost::math::quantile(complement(dist, p)) : 
    -boost::math::quantile(complement(dist, 1.0 - p));
}

double boost_dextreme(double x, double location, double scale) {
  if (std::isnan(x)) return std::numeric_limits<double>::quiet_NaN();
  if (scale <= 0) throw std::invalid_argument("Scale must be positive.");
  boost::math::extreme_value_distribution<> dist(location, scale);
  return boost::math::pdf(dist, -x);
}

double boost_pchisq(double q, double df, bool lower_tail) {
  if (std::isnan(q)) return std::numeric_limits<double>::quiet_NaN();
  if (df <= 0) throw std::invalid_argument("Degrees of freedom must be positive.");
  if (std::isinf(q)) {
    if (q > 0.0) return lower_tail ? 1.0 : 0.0;
    else return lower_tail ? 0.0 : 1.0;
  }
  boost::math::chi_squared_distribution<> dist(df);
  if (lower_tail) return boost::math::cdf(dist, q);
  else return boost::math::cdf(boost::math::complement(dist, q));
}

double boost_qchisq(double p, double df, bool lower_tail) {
  if (std::isnan(p)) return std::numeric_limits<double>::quiet_NaN();
  if (df <= 0) throw std::invalid_argument("Degrees of freedom must be positive.");
  if (p < 0.0 || p > 1.0) throw std::invalid_argument(
      "Probability must be between 0 and 1.");
  boost::math::chi_squared_distribution<> dist(df);
  return lower_tail ? boost::math::quantile(dist, p) : 
    boost::math::quantile(dist, 1.0 - p);
}

double boost_pt(double q, double df, bool lower_tail) {
  if (std::isnan(q)) return std::numeric_limits<double>::quiet_NaN();
  if (df <= 0) throw std::invalid_argument("Degrees of freedom must be positive.");
  boost::math::students_t_distribution<> dist(df);
  if (lower_tail) return boost::math::cdf(dist, q);
  else return boost::math::cdf(boost::math::complement(dist, q));
}

double boost_qt(double p, double df, bool lower_tail) {
  if (std::isnan(p)) return std::numeric_limits<double>::quiet_NaN();
  if (df <= 0) throw std::invalid_argument("Degrees of freedom must be positive.");
  if (p < 0.0 || p > 1.0) throw std::invalid_argument(
      "Probability must be between 0 and 1.");
  boost::math::students_t_distribution<> dist(df);
  return lower_tail ? boost::math::quantile(dist, p) : 
    boost::math::quantile(dist, 1.0 - p);
}

// --------------------------- Small utilities --------------------------------

std::vector<size_t> seqcpp(size_t start, size_t end) {
  if (start > end) throw std::invalid_argument(
      "start must be less than or equal to end for the sequence function.");
  size_t size = end - start + 1;
  std::vector<size_t> result(size);
  std::iota(result.begin(), result.end(), start);
  return result;
}

std::vector<size_t> which(const std::vector<unsigned char>& vec) {
  std::vector<size_t> indices;
  indices.reserve(vec.size());
  for (size_t i = 0; i < vec.size(); ++i) {
    if (vec[i] != 0 && vec[i] != 255) indices.push_back(i);
  }
  return indices;
}

size_t findInterval1(const double x,
                     const std::vector<double>& v,
                     bool rightmost_closed,
                     bool all_inside,
                     bool left_open) {
  
  const double* v_begin = v.data();
  const double* v_end   = v_begin + v.size();
  const size_t nv = v.size();
  
  const double* pos = left_open ? std::lower_bound(v_begin, v_end, x) :
    std::upper_bound(v_begin, v_end, x);
  size_t idx = static_cast<size_t>(pos - v_begin);
  if (rightmost_closed) {
    if (left_open) {
      if (x == v[0]) idx = 1;
    } else {
      if (x == v[nv - 1]) idx = nv - 1;
    }
  }
  if (all_inside) {
    if (idx == 0) idx = 1;
    else if (idx == nv) idx = nv - 1;
  }
  
  return idx;
}


std::vector<size_t> findInterval3(const std::vector<double>& x,
                                  const std::vector<double>& v,
                                  bool rightmost_closed,
                                  bool all_inside,
                                  bool left_open) {
  std::vector<size_t> out(x.size());
  const double* v_begin = v.data();
  const double* v_end   = v_begin + v.size();
  const size_t n = x.size();
  const size_t nv = v.size();
  
  for (size_t i = 0; i < n; ++i) {
    double xi = x[i];
    const double* pos = left_open ? std::lower_bound(v_begin, v_end, xi) :
      std::upper_bound(v_begin, v_end, xi);
    size_t idx = static_cast<size_t>(pos - v_begin);
    if (rightmost_closed) {
      if (left_open) {
        if (nv > 0 && xi == v[0]) idx = 1;
      } else {
        if (nv > 0 && xi == v[nv - 1]) idx = nv - 1;
      }
    }
    if (all_inside) {
      if (idx == 0) idx = 1;
      else if (idx == nv) idx = nv - 1;
    }
    out[i] = idx;
  }
  return out;
}

// --------------------------- Quantiles -------------------------------------

double quantilecpp(const std::vector<double>& x, double p) {
  size_t n = x.size();
  if (n == 0) throw std::invalid_argument("Empty vector");
  if (p < 0.0 || p > 1.0) throw std::invalid_argument("p must be in [0,1]");
  std::vector<double> y(x);
  std::sort(y.begin(), y.end());
  double h = (n - 1) * p + 1;
  size_t j = static_cast<size_t>(std::floor(h));
  double g = h - j;
  if (j <= 0) return y.front();
  if (j >= n) return y.back();
  return (1 - g) * y[j - 1] + g * y[j];
}



// in-place truncation of vector
void truncate_in_place(std::vector<double>& v, bool trunc_upper_only, double trunc) {
  if (v.empty()) return;
  if (trunc < 0.0 || trunc >= 0.5) {
    throw std::invalid_argument("trunc must lie in [0, 0.5)");
  }
  if (trunc_upper_only) {
    double upper = quantilecpp(v, 1 - trunc);
    for (double &x : v) if (x > upper) x = upper;
  } else {
    double lower = quantilecpp(v, trunc);
    double upper = quantilecpp(v, 1 - trunc);
    for (double &x : v) x = std::clamp(x, lower, upper);
  }
};


// bygroup: group-by helper that builds lookup tables and combined indices
ListCpp bygroup(const DataFrameCpp& data,
                const std::vector<std::string>& variables) {
  
  size_t n = data.nrows();
  size_t p = variables.size();
  ListCpp result;
  std::vector<size_t> nlevels(p);
  
  // IntMatrix for indices (n rows, p cols), column-major storage
  IntMatrix indices(n, p);
  
  // Flattened lookup buffers and per-variable metadata
  struct VarLookupInfo {
    int type; // 0=int, 1=double, 2=bool, 3=string, 4=size_t
    size_t offset;
  };
  std::vector<VarLookupInfo> var_info(p);
  
  std::vector<int> int_flat;
  std::vector<double> dbl_flat;
  std::vector<unsigned char> bool_flat;
  std::vector<std::string> str_flat;
  std::vector<size_t> size_t_flat;
  
  ListCpp lookups_per_variable; // will contain a std::vector for each variable
  
  for (size_t i = 0; i < p; ++i) {
    const std::string& var = variables[i];
    if (!data.containElementNamed(var))
      throw std::invalid_argument("Data must contain variable: " + var);
    
    if (data.int_cols.count(var)) {
      const auto& col = data.int_cols.at(var);
      auto w = unique_sorted(col);
      nlevels[i] = w.size();
      auto idx = matchcpp(col, w); // indices 0..(levels-1)
      
      // append w to flat buffer and record metadata
      size_t off = int_flat.size();
      int_flat.insert(int_flat.end(), w.begin(), w.end());
      var_info[i] = VarLookupInfo{0, off};
      
      // fill indices column i (column-major layout)
      intmatrix_set_column(indices, i, idx);
      lookups_per_variable.push_back(std::move(w), var);
    } else if (data.numeric_cols.count(var)) {
      const auto& col = data.numeric_cols.at(var);
      auto w = unique_sorted(col);
      nlevels[i] = w.size();
      auto idx = matchcpp(col, w);
      
      size_t off = dbl_flat.size();
      dbl_flat.insert(dbl_flat.end(), w.begin(), w.end());
      var_info[i] = VarLookupInfo{1, off};
      
      intmatrix_set_column(indices, i, idx);
      lookups_per_variable.push_back(std::move(w), var);
    } else if (data.bool_cols.count(var)) {
      const auto& col = data.bool_cols.at(var);
      auto w = unique_sorted(col);
      nlevels[i] = w.size();
      auto idx = matchcpp(col, w);
      
      size_t off = bool_flat.size();
      bool_flat.insert(bool_flat.end(), w.begin(), w.end());
      var_info[i] = VarLookupInfo{2, off};
      
      intmatrix_set_column(indices, i, idx);
      lookups_per_variable.push_back(std::move(w), var);
    } else if (data.string_cols.count(var)) {
      const auto& col = data.string_cols.at(var);
      auto w = unique_sorted(col);
      nlevels[i] = w.size();
      auto idx = matchcpp(col, w);
      
      size_t off = str_flat.size();
      str_flat.insert(str_flat.end(), w.begin(), w.end());
      var_info[i] = VarLookupInfo{3, off};
      
      intmatrix_set_column(indices, i, idx);
      lookups_per_variable.push_back(std::move(w), var);
    } else if (data.size_t_cols.count(var)) {
      const auto& col = data.size_t_cols.at(var);
      auto w = unique_sorted(col);
      nlevels[i] = w.size();
      auto idx = matchcpp(col, w);
      
      size_t off = size_t_flat.size();
      size_t_flat.insert(size_t_flat.end(), w.begin(), w.end());
      var_info[i] = VarLookupInfo{4, off};
      
      intmatrix_set_column(indices, i, idx);
      lookups_per_variable.push_back(std::move(w), var);
    } else {
      throw std::invalid_argument("Unsupported variable type in bygroup: " + var);
    }
  } // end for variables
  
  // compute combined index
  std::vector<int> combined_index(n, 0);
  size_t orep = 1;
  for (size_t i = 0; i < p; ++i) orep *= nlevels[i];
  size_t lookup_nrows = orep;
  
  for (size_t i = 0; i < p; ++i) {
    orep /= nlevels[i];
    const int* col_ptr = indices.data_ptr() + i * n;
    for (size_t j = 0; j < n; ++j) {
      combined_index[j] += col_ptr[j] * orep;
    }
  }
  
  // Build lookup_df with columns repeated in the same pattern as before.
  DataFrameCpp lookup_df;
  size_t repeat_each = lookup_nrows;
  for (size_t i = 0; i < p; ++i) {
    const std::string& var = variables[i];
    size_t nlevels_i = nlevels[i];
    repeat_each /= nlevels_i;
    size_t times = lookup_nrows / ( nlevels_i * repeat_each );
    
    VarLookupInfo info = var_info[i];
    if (info.type == 0) {
      const int* base = int_flat.data() + info.offset;
      std::vector<int> col(lookup_nrows);
      size_t idxw = 0;
      for (size_t t = 0; t < times; ++t) {
        for (size_t level = 0; level < nlevels_i; ++level)
          for (size_t r = 0; r < repeat_each; ++r) col[idxw++] = base[level];
      }
      lookup_df.push_back(std::move(col), var);
    } else if (info.type == 1) {
      const double* base = dbl_flat.data() + info.offset;
      std::vector<double> col(lookup_nrows);
      size_t idxw = 0;
      for (size_t t = 0; t < times; ++t) {
        for (size_t level = 0; level < nlevels_i; ++level)
          for (size_t r = 0; r < repeat_each; ++r) col[idxw++] = base[level];
      }
      lookup_df.push_back(std::move(col), var);
    } else if (info.type == 2) {
      const unsigned char* base = bool_flat.data() + info.offset;
      std::vector<unsigned char> col(lookup_nrows);
      size_t idxw = 0;
      for (size_t t = 0; t < times; ++t) {
        for (size_t level = 0; level < nlevels_i; ++level) {
          for (size_t r = 0; r < repeat_each; ++r) col[idxw++] = base[level];
        }
      }
      lookup_df.push_back(std::move(col), var);
    } else if (info.type == 3) { // string
      const std::string* base = str_flat.data() + info.offset;
      std::vector<std::string> col(lookup_nrows);
      size_t idxw = 0;
      for (size_t t = 0; t < times; ++t) {
        for (size_t level = 0; level < nlevels_i; ++level)
          for (size_t r = 0; r < repeat_each; ++r) col[idxw++] = base[level];
      }
      lookup_df.push_back(std::move(col), var);
    } else { // size_t
      const size_t* base = size_t_flat.data() + info.offset;
      std::vector<size_t> col(lookup_nrows);
      size_t idxw = 0;
      for (size_t t = 0; t < times; ++t) {
        for (size_t level = 0; level < nlevels_i; ++level)
          for (size_t r = 0; r < repeat_each; ++r) col[idxw++] = base[level];
      }
      lookup_df.push_back(std::move(col), var);
    }
  }
  
  result.push_back(std::move(nlevels), "nlevels");
  result.push_back(std::move(indices), "indices");
  result.push_back(std::move(lookups_per_variable), "lookups_per_variable");
  result.push_back(std::move(combined_index), "index");
  result.push_back(std::move(lookup_df), "lookup");
  return result;
}

// --------------------------- Matrix utilities (FlatMatrix) ------------------

std::vector<double> mat_vec_mult(const FlatMatrix& A, const std::vector<double>& x) {
  size_t m = A.nrow;
  size_t p = A.ncol;
  if (x.size() != p)
    throw std::invalid_argument("Vector size mismatch");
  std::vector<double> result(m, 0.0);
  for (size_t c = 0; c < p; ++c) {
    double xc = x[c];
    size_t offset = c * m;
    const double* colptr = A.data_ptr() + offset;
    for (size_t r = 0; r < m; ++r) result[r] += colptr[r] * xc;
  }
  return result;
}

FlatMatrix mat_mat_mult(const FlatMatrix& A, const FlatMatrix& B) {
  size_t m = A.nrow;
  size_t k = A.ncol;
  size_t k2 = B.nrow;
  size_t n = B.ncol;
  if (k != k2) throw std::invalid_argument("Matrix dimensions mismatch");
  if (m == 0 || k == 0 || n == 0) return FlatMatrix();
  FlatMatrix C(m, n);
  // Column-major: For each column j in B/C, compute
  // C[:,j] = sum_{t=0..k-1} A[:,t] * B[t,j]
  for (size_t j = 0; j < n; ++j) {
    const double* bcol = B.data_ptr() + j * k;
    double* ccol = C.data_ptr() + j * m;
    for (size_t t = 0; t < k; ++t) {
      const double* acol = A.data_ptr() + t * m;
      double scale = bcol[t];
      if (scale == 0.0) continue;
      for (size_t i = 0; i < m; ++i) {
        ccol[i] += acol[i] * scale;
      }
    }
  }
  return C;
}

// Transpose a FlatMatrix (double)
FlatMatrix transpose(const FlatMatrix& M) {
  if (M.nrow == 0 || M.ncol == 0) return FlatMatrix();
  
  const size_t src_nrow = M.nrow;
  const size_t src_ncol = M.ncol;
  FlatMatrix out(src_ncol, src_nrow); // swapped dims
  
  const double* src = M.data_ptr();
  double* dst = out.data_ptr();
  
  for (size_t c = 0; c < src_ncol; ++c) {
    const double* src_col = src + c * src_nrow;
    for (size_t r = 0; r < src_nrow; ++r) {
      dst[r * src_ncol + c] = src_col[r];
    }
  }
  
  return out;
}

// Transpose an IntMatrix (int)
IntMatrix transpose(const IntMatrix& M) {
  if (M.nrow == 0 || M.ncol == 0) return IntMatrix();
  
  const size_t src_nrow = M.nrow;
  const size_t src_ncol = M.ncol;
  IntMatrix out(src_ncol, src_nrow); // swapped dims
  
  const int* src = M.data_ptr();
  int* dst = out.data_ptr();
  
  for (size_t c = 0; c < src_ncol; ++c) {
    const int* src_col = src + c * src_nrow;
    for (size_t r = 0; r < src_nrow; ++r) {
      dst[r * src_ncol + c] = src_col[r];
    }
  }
  
  return out;
}

// Transpose an SztMatrix (std::size_t)
SztMatrix transpose(const SztMatrix& M) {
  if (M.nrow == 0 || M.ncol == 0) return SztMatrix();
  
  const size_t src_nrow = M.nrow;
  const size_t src_ncol = M.ncol;
  SztMatrix out(src_ncol, src_nrow); // swapped dims
  
  const size_t* src = M.data_ptr();
  size_t* dst = out.data_ptr();
  
  for (size_t c = 0; c < src_ncol; ++c) {
    const size_t* src_col = src + c * src_nrow;
    for (size_t r = 0; r < src_nrow; ++r) {
      dst[r * src_ncol + c] = src_col[r];
    }
  }
  
  return out;
}

// Transpose a BoolMatrix (unsigned char)
BoolMatrix transpose(const BoolMatrix& M) {
  if (M.nrow == 0 || M.ncol == 0) return BoolMatrix();
  
  const size_t src_nrow = M.nrow;
  const size_t src_ncol = M.ncol;
  BoolMatrix out(src_ncol, src_nrow); // swapped dims
  
  const unsigned char* src = M.data_ptr();
  unsigned char* dst = out.data_ptr();
  
  for (size_t c = 0; c < src_ncol; ++c) {
    const unsigned char* src_col = src + c * src_nrow;
    for (size_t r = 0; r < src_nrow; ++r) {
      dst[r * src_ncol + c] = src_col[r];
    }
  }
  
  return out;
}

double quadsym(const std::vector<double>& u, const FlatMatrix& v) {
  size_t p = u.size();
  const double* vptr = v.data_ptr();
  const double* uptr = u.data();
  double sum = 0.0;
  
  for (size_t j = 0; j < p; ++j) {
    const double* col = vptr + j * p;
    // diagonal term
    sum += uptr[j] * uptr[j] * col[j]; // col[j] == v(j,j)
    // off-diagonals i < j. Access column j contiguous for i = 0..j-1
    double s = 0.0;
    for (size_t i = 0; i < j; ++i) s += col[i] * uptr[i];
    sum += 2.0 * uptr[j] * s; // account for symmetric pair (i,j) and (j,i)
  }
  return sum;
}

// --------------------------- Linear algebra helpers (FlatMatrix-backed) ----
// cholesky2: in-place working on FlatMatrix (n x n), returns rank * nonneg
int cholesky2(FlatMatrix& matrix, size_t n, double toler) {
  double* base = matrix.data_ptr();
  double eps = 0.0;
  for (size_t i = 0; i < n; ++i) {
    double val = matrix(i, i);
    if (val > eps) eps = val;
  }
  if (eps == 0.0) eps = toler; else eps *= toler;
  int nonneg = 1;
  int rank = 0;
  
  for (size_t i = 0; i < n; ++i) {
    double* col_i = base + i * n;
    double pivot = col_i[i];
    if (std::isinf(pivot) || pivot < eps) {
      col_i[i] = 0.0;
      if (pivot < -8.0 * eps) nonneg = -1;
    } else {
      ++rank;
      for (size_t j = i + 1; j < n; ++j) {
        double* col_j = base + j * n;
        double temp = col_i[j] / pivot;
        col_i[j] = temp;
        col_j[j] -= temp * temp * pivot;
        for (size_t k = j + 1; k < n; ++k) {
          col_j[k] -= temp * col_i[k];
        }
      }
    }
  }
  return rank * nonneg;
}


// chsolve2 assumes matrix holds the representation produced by cholesky2
void chsolve2(FlatMatrix& matrix, size_t n, double* y) {
  // Forward substitution L * z = y
  double* base = matrix.data_ptr();
  for (size_t j = 0; j < n-1; ++j) {
    double yj = y[j];
    if (yj == 0.0) continue;
    double* col_j = base + j * n;
    for (size_t i = j + 1; i < n; ++i) {
      y[i] -= yj * col_j[i];
    }
  }
  // Now y holds z; solve L^T * x = z
  if (n == 0) return;
  for (size_t i = n; i-- > 0; ) {
    double* col_i = base + i * n;
    double diag = col_i[i];
    if (diag == 0.0) {
      y[i] = 0.0;
    } else {
      double temp = y[i] / diag;
      for (size_t j = i + 1; j < n; ++j) temp -= y[j] * col_i[j];
      y[i] = temp;
    }
  }
}


// invsympd: returns the inverse of a symmetric positive definite matrix
FlatMatrix invsympd(const FlatMatrix& matrix, size_t n, double toler) {
  FlatMatrix v = matrix; // copy
  cholesky2(v, n, toler);
  FlatMatrix iv(n, n);
  for (size_t i = 0; i < n; ++i) {
    iv(i,i) = 1.0;
    double* ycol = iv.data_ptr() + i * n;
    chsolve2(v, n, ycol);
  }
  return iv;
}


// -------------------------- Survival helpers --------------------------------

DataFrameCpp survsplitcpp(const std::vector<double>& tstart,
                          const std::vector<double>& tstop,
                          const std::vector<double>& cut) {
  size_t n = tstart.size();
  size_t ncut = cut.size();
  size_t extra = 0;
  for (size_t i = 0; i < n; ++i) {
    if (std::isnan(tstart[i]) || std::isnan(tstop[i])) continue;
    for (size_t j = 0; j < ncut; ++j) {
      if (cut[j] > tstart[i] && cut[j] < tstop[i]) ++extra;
    }
  }
  size_t n2 = n + extra;
  std::vector<size_t> row(n2);
  std::vector<size_t> interval(n2);
  std::vector<double> start(n2);
  std::vector<double> end(n2);
  std::vector<int> censor(n2);
  size_t k = 0;
  for (size_t i = 0; i < n; ++i) {
    if (std::isnan(tstart[i]) || std::isnan(tstop[i])) {
      start[k] = tstart[i];
      end[k] = tstop[i];
      row[k] = i;
      interval[k] = 1;
      ++k;
    } else {
      size_t j = 0;
      while (j < ncut && cut[j] <= tstart[i]) ++j;
      start[k] = tstart[i];
      row[k] = i;
      interval[k] = j;
      for (; j < ncut && cut[j] < tstop[i]; ++j) {
        if (cut[j] > tstart[i]) {
          end[k] = cut[j];
          censor[k] = 1;
          ++k;
          start[k] = cut[j];
          row[k] = i;
          interval[k] = j + 1;
        }
      }
      end[k] = tstop[i];
      censor[k] = 0;
      ++k;
    }
  }
  DataFrameCpp df;
  df.push_back(std::move(row), "row");
  df.push_back(std::move(start), "start");
  df.push_back(std::move(end), "end");
  df.push_back(std::move(censor), "censor");
  df.push_back(std::move(interval), "interval");
  return df;
}


//' @title Split a survival data set at specified cut points
//' @description For a given survival dataset and specified cut times, 
//' each record is split into multiple subrecords at each cut time. 
//' The resulting dataset is in counting process format, with each 
//' subrecord containing a start time, stop time, and event status.
//' This is adapted from the survsplit.c function from the survival package.
//'
//' @param tstart The starting time of the time interval for 
//'   counting-process data.
//' @param tstop The stopping time of the time interval for 
//'   counting-process data.
//' @param cut The vector of cut points.
//'
//' @return A data frame with the following variables:
//'
//' * \code{row}: The row number of the observation in the input data 
//'   (starting from 0).
//'
//' * \code{start}: The starting time of the resulting subrecord.
//'
//' * \code{end}: The ending time of the resulting subrecord.
//'
//' * \code{censor}: Whether the subrecord lies strictly within a record
//'   in the input data (1 for all but the last interval and 0 for the 
//'   last interval).
//'
//' * \code{interval}: The interval number derived from cut (starting 
//'   from 0 if the interval lies to the left of the first cutpoint).
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @keywords internal
//'
//' @examples
//'
//' survsplit(15, 60, c(10, 30, 40))
//'
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame survsplit(const std::vector<double>& tstart,
                          const std::vector<double>& tstop,
                          const std::vector<double>& cut) {
  auto dfcpp = survsplitcpp(tstart, tstop, cut);
  return Rcpp::wrap(dfcpp);
}


// ------------------------- QR and other helpers -----------------------------
double sumsq(const std::vector<double>& x) {
  double s = 0.0;
  for (double xi : x) s += xi * xi;
  return s;
}

// Householder vector
// Given an n-vector x, this function computes an n-vector v with v(1) = 1
// such that (I - 2*v*t(v)/t(v)*v)*x is zero in all but the first component.
std::vector<double> house(const std::vector<double>& x) {
  size_t n = x.size();
  double sumxx = std::inner_product(x.begin(), x.end(), x.begin(), 0.0);
  double mu = std::sqrt(sumxx);
  std::vector<double> v = x;
  if (mu > 0.0) {
    double beta = x[0] + std::copysign(1.0, x[0])*mu;
    for (size_t i=1; i<n; ++i) {
      v[i] /= beta;
    }
  }
  v[0] = 1.0;
  return v;
}

// Householder pre-multiplication
// Given an m-by-n matrix A and a nonzero m-vector v with v(1) = 1,
// the following algorithm overwrites A with P*A where
// P = I - 2*v*t(v)/t(v)*v.
void row_house(FlatMatrix& A, const size_t i1, const size_t i2, const size_t j1,
               const size_t j2, const std::vector<double>& v) {
  if (i1 < 0 || i1 > i2 || i2 >= A.nrow) {
    throw std::invalid_argument("Invalid row indices i1 and i2");
  }
  if (j1 < 0 || j1 > j2 || j2 >= A.ncol) {
    throw std::invalid_argument("Invalid column indices j1 and j2");
  }
  
  size_t m = i2-i1+1, n = j2-j1+1;
  double sumvv = std::inner_product(v.begin(), v.end(), v.begin(), 0.0);
  double beta = -2.0 / sumvv;
  std::vector<double> w(n);
  for (size_t j=0; j<n; ++j) {
    for (size_t i=0; i<m; ++i) {
      w[j] += A(i+i1,j+j1)*v[i];
    }
    w[j] *= beta;
  }
  
  for (size_t j=0; j<n; ++j) {
    for (size_t i=0; i<m; ++i) {
      A(i+i1,j+j1) += v[i]*w[j];
    }
  }
}

ListCpp qrcpp(const FlatMatrix& X, double tol) {
  size_t m = X.nrow;
  size_t n = X.ncol;
  FlatMatrix A = X;
  std::vector<double> c(n, 0.0);
  for (size_t j = 0; j < n; ++j) {
    double* acol = A.data_ptr() + j * m;
    double s = 0.0;
    for (size_t i = 0; i < m; ++i) {
      double v = acol[i];
      s += v * v;
    }
    c[j] = s;
  }
  
  size_t r = 0;
  std::vector<int> piv(n);
  std::iota(piv.begin(), piv.end(), 0);
  double tau = *std::max_element(c.begin(), c.end());
  
  while (tau > tol) {
    size_t k = r;
    for (; k < n; ++k) {
      if (c[k] > tol) break;
    }
    if (k != r) {
      size_t off_r = r * m;
      size_t off_k = k * m;
      for (size_t i = 0; i < m; ++i) {
        std::swap(A.data[off_r + i], A.data[off_k + i]);
      }
      std::swap(piv[r], piv[k]);
      std::swap(c[r], c[k]);
    }
    
    size_t msub = m - r;
    std::vector<double> x(msub);
    for (size_t i = 0; i < msub; ++i) x[i] = A(r + i, r);
    std::vector<double> v = house(x);
    
    if (msub > 0 && r < n) row_house(A, r, m - 1, r, n - 1, v);
    
    for (size_t i = 1; i < msub; ++i) A(r + i, r) = v[i];
    for (size_t i = r + 1; i < n; ++i) {
      double val = A(r, i); c[i] -= val * val;
    }
    if (r < n - 1) {
      tau = *std::max_element(c.begin() + (r + 1), c.end());
    } else tau = 0.0;
    
    r++;
  }
  
  FlatMatrix Qf(m, m);
  for (size_t c = 0; c < m; ++c) Qf(c, c) = 1.0;
  
  if (r > 0) {
    for (size_t k = r; k-- > 0; ) {
      size_t msub_k = m - k;
      std::vector<double> vks(msub_k);
      vks[0] = 1.0;
      for (size_t i = 1; i < msub_k; ++i) {
        vks[i] = A(k + i, k);
      }
      row_house(Qf, k, m - 1, k, m - 1, vks);
    }
  }
  
  FlatMatrix Rf(m, n);
  std::fill(Rf.data.begin(), Rf.data.end(), 0.0);
  for (size_t j = 0; j < n; ++j) {
    for (size_t i = 0; i <= j && i < m; ++i) {
      Rf(i, j) = A(i, j);
    }
  }
  
  ListCpp result;
  result.push_back(std::move(A), "qr");
  result.push_back(r, "rank");
  result.push_back(std::move(piv), "pivot");
  result.push_back(std::move(Qf), "Q");
  result.push_back(std::move(Rf), "R");
  return result;
}


// -------------------------- Matching and other helpers ----------------------

std::vector<int> match3(const std::vector<int>& id1,
                        const std::vector<double>& v1,
                        const std::vector<int>& id2,
                        const std::vector<double>& v2) {
  std::vector<int> result;
  result.reserve(id1.size());
  size_t i = 0, j = 0;
  size_t n1 = id1.size(), n2 = id2.size();
  while (i < n1 && j < n2) {
    if (id1[i] < id2[j] || (id1[i] == id2[j] && v1[i] < v2[j])) {
      result.push_back(-1); ++i;
    } else if (id1[i] > id2[j] || (id1[i] == id2[j] && v1[i] > v2[j])) {
      ++j;
    } else {
      result.push_back(j); ++i; ++j;
    }
  }
  while (i < n1) { result.push_back(-1); ++i; }
  return result;
}


// --------------------------- Misc helpers -----------------------------------

std::string sanitize(const std::string& s) {
  std::string out = s;
  for (char &c : out) {
    if (!std::isalnum(static_cast<unsigned char>(c)) && 
        static_cast<unsigned char>(c) != '_') 
      c = '.';
  }
  return out;
}


double getAccrualDurationFromN1(
    const double nsubjects,
    const std::vector<double>& accrualTime,
    const std::vector<double>& accrualIntensity) {
  
  size_t J = accrualTime.size();
  std::vector<double> p(J);
  p[0] = 0;
  for (size_t j = 0; j < J - 1; ++j) {
    p[j + 1] = p[j] + accrualIntensity[j] * (accrualTime[j + 1] - accrualTime[j]);
  }
  
  size_t m = findInterval1(nsubjects, p) - 1;
  double t = accrualTime[m] + (nsubjects - p[m]) / accrualIntensity[m];
  return t;
}


// -------------------------- Counterfactual helpers --------------------------

DataFrameCpp untreated(double psi,
                       const std::vector<int>& id,
                       const std::vector<double>& time,
                       const std::vector<int>& event,
                       const std::vector<int>& treat,
                       const std::vector<double>& rx,
                       const std::vector<double>& censor_time,
                       bool recensor,
                       bool autoswitch) {
  size_t n = id.size();
  double a = std::exp(psi);
  std::vector<double> u_star(n), t_star(n);
  std::vector<int> d_star = event;
  for (size_t i = 0; i < n; ++i) { 
    u_star[i] = time[i] * ((1.0 - rx[i]) + rx[i] * a); t_star[i] = u_star[i]; }
  if (recensor) {
    std::vector<double> c_star = censor_time;
    for (size_t i = 0; i < n; ++i) c_star[i] *= std::min(1.0, a);
    if (autoswitch) {
      bool all_rx1 = true, all_rx0 = true;
      for (size_t i = 0; i < n; ++i) {
        if (treat[i] == 1 && rx[i] != 1.0) all_rx1 = false;
        if (treat[i] == 0 && rx[i] != 0.0) all_rx0 = false;
      }
      for (size_t i = 0; i < n; ++i) {
        if (treat[i] == 1 && all_rx1) 
          c_star[i] = std::numeric_limits<double>::infinity();
        if (treat[i] == 0 && all_rx0) 
          c_star[i] = std::numeric_limits<double>::infinity();
      }
    }
    for (size_t i = 0; i < n; ++i) {
      if (c_star[i] < u_star[i]) { t_star[i] = c_star[i]; d_star[i] = 0; }
    }
  }
  DataFrameCpp df;
  df.push_back(id, "uid");
  df.push_back(std::move(t_star), "t_star");
  df.push_back(std::move(d_star), "d_star");
  df.push_back(treat, "treated");
  return df;
}

DataFrameCpp unswitched(double psi,
                        const std::vector<int>& id,
                        const std::vector<double>& time,
                        const std::vector<int>& event,
                        const std::vector<int>& treat,
                        const std::vector<double>& rx,
                        const std::vector<double>& censor_time,
                        bool recensor,
                        bool autoswitch) {
  size_t n = id.size();
  double a0 = std::exp(psi);
  double a1 = std::exp(-psi);
  std::vector<double> u_star(n), t_star(n);
  std::vector<int> d_star = event;
  for (size_t i = 0; i < n; ++i) {
    if (treat[i] == 0) u_star[i] = time[i] * ((1.0 - rx[i]) + rx[i] * a0);
    else u_star[i] = time[i] * (rx[i] + (1.0 - rx[i]) * a1);
    t_star[i] = u_star[i];
  }
  if (recensor) {
    std::vector<double> c_star(n);
    for (size_t i = 0; i < n; ++i) 
      c_star[i] = treat[i] == 0 ? censor_time[i] * std::min(1.0, a0) : 
      censor_time[i] * std::min(1.0, a1);
    if (autoswitch) {
      bool all_rx1 = true, all_rx0 = true;
      for (size_t i = 0; i < n; ++i) {
        if (treat[i] == 1 && rx[i] != 1.0) all_rx1 = false;
        if (treat[i] == 0 && rx[i] != 0.0) all_rx0 = false;
      }
      for (size_t i = 0; i < n; ++i) {
        if (treat[i] == 1 && all_rx1) 
          c_star[i] = std::numeric_limits<double>::infinity();
        if (treat[i] == 0 && all_rx0) 
          c_star[i] = std::numeric_limits<double>::infinity();
      }
    }
    for (size_t i = 0; i < n; ++i) {
      if (c_star[i] < u_star[i]) { t_star[i] = c_star[i]; d_star[i] = 0; }
    }
  }
  DataFrameCpp df;
  df.push_back(id, "uid");
  df.push_back(std::move(t_star), "t_star");
  df.push_back(std::move(d_star), "d_star");
  df.push_back(treat, "treated");
  return df;
}


// --------------------------- Root selection and helpers ---------------------

ListCpp getpsiest(double target,
                  const std::vector<double>& psi,
                  const std::vector<double>& Z,
                  int direction) {
  size_t n = psi.size();
  if (n != Z.size()) 
    throw std::invalid_argument("psi and Z must have the same length");
  if (n < 2) 
    throw std::invalid_argument("Need at least two points to find roots");
  std::vector<double> Zt(n);
  for (size_t i = 0; i < n; ++i) Zt[i] = Z[i] - target;
  std::vector<double> roots;
  for (size_t i = 1; i < n; ++i) {
    double z1 = Zt[i - 1];
    double z2 = Zt[i];
    if (std::isnan(z1) || std::isnan(z2) || z1 == z2) continue;
    if (z1 == 0.0) roots.push_back(psi[i - 1]);
    else if (z1 * z2 < 0.0) {
      double psi_root = psi[i - 1] - z1 * (psi[i] - psi[i - 1]) / (z2 - z1);
      roots.push_back(psi_root);
    }
  }
  double root = NaN;
  if (!roots.empty()) {
    if (direction == -1) root = roots.front();
    else if (direction == 1) root = roots.back();
    else {
      root = roots[0];
      double minabs = std::abs(roots[0]);
      size_t m = roots.size();
      for (size_t j = 1; j < m; ++j) {
        double a = std::abs(roots[j]);
        if (a < minabs) { minabs = a; root = roots[j]; }
      }
    }
  }
  ListCpp result;
  result.push_back(std::move(roots), "all_roots");
  result.push_back(root, "selected_root");
  return result;
}

double getpsiend(const std::function<double(double)>& f,
                 bool lowerend, double initialend) {
  double psiend = initialend;
  double zend = f(initialend);
  const double LIMIT = 10.0;
  if (lowerend) {
    if ((std::isinf(zend) && zend > 0) || std::isnan(zend)) {
      while (((std::isinf(zend) && zend > 0) || std::isnan(zend)) && 
             psiend <= LIMIT) {
        psiend += 1; zend = f(psiend);
      }
      if (psiend > LIMIT) return NaN;
    }
    if (zend < 0) {
      while (!std::isinf(zend) && zend < 0 && psiend >= -LIMIT) {
        psiend -= 1; zend = f(psiend);
      }
      if (std::isinf(zend) || std::isnan(zend) || psiend < -LIMIT) return NaN;
    }
  } else {
    if ((std::isinf(zend) && zend < 0) || std::isnan(zend)) {
      while (((std::isinf(zend) && zend < 0) || std::isnan(zend)) && 
             psiend >= -LIMIT) {
        psiend -= 1; zend = f(psiend);
      }
      if (psiend < -LIMIT) return NaN;
    }
    if (zend > 0) {
      while (!std::isinf(zend) && zend > 0 && psiend <= LIMIT) {
        psiend += 1; zend = f(psiend);
      }
      if (std::isinf(zend) || std::isnan(zend) || psiend > LIMIT) return NaN;
    }
  }
  return psiend;
}
