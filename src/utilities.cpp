#include "utilities.h"
#include "dataframe_list.h"

#include <Rcpp.h>

#include <cstddef>    // size_t
#include <vector>     // vector
#include <string>     // string
#include <algorithm>  // fill, lower_bound, max_element, memmove, 
                      // min_element, sort, swap, upper_bound
#include <numeric>    // accumulate, inner_product, iota 
#include <functional> // function
#include <cmath>      // copysign, exp, fabs, isinf, isnan, log, pow, sqrt
#include <limits>     // numeric_limits
#include <sstream>    // ostringstream
#include <stdexcept>  // invalid_argument, runtime_error
#include <cstring>    // memcpy
#include <memory>     // make_shared, shared_ptr

#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/logistic.hpp>
#include <boost/math/distributions/extreme_value.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/students_t.hpp>

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

std::vector<int> seqcpp(int start, int end) {
  if (start > end) throw std::invalid_argument(
      "start must be less than or equal to end for the sequence function.");
  int size = end - start + 1;
  std::vector<int> result(size);
  std::iota(result.begin(), result.end(), start);
  return result;
}

std::vector<int> which(const std::vector<unsigned char>& vec) {
  std::vector<int> indices;
  indices.reserve(vec.size());
  int n = vec.size();
  for (int i = 0; i < n; ++i) 
    if (vec[i] != 0 && vec[i] != 255) indices.push_back(i);
  return indices;
}

std::vector<int> findInterval3(const std::vector<double>& x,
                               const std::vector<double>& v,
                               bool rightmost_closed,
                               bool all_inside,
                               bool left_open) {
  std::vector<int> out(x.size());
  const double* v_begin = v.data();
  const double* v_end   = v_begin + v.size();
  const int n = x.size();
  const int nv = v.size();
  
  for (int i = 0; i < n; ++i) {
    double xi = x[i];
    if (std::isnan(xi)) { out[i] = -1; continue; }
    const double* pos = left_open ? std::lower_bound(v_begin, v_end, xi) : 
      std::upper_bound(v_begin, v_end, xi);
    int idx = static_cast<int>(pos - v_begin);
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

// --------------------------- Root finders -----------------------------------

static const double EPS = 3.0e-8;
inline double SIGN(double a, double b) { 
  return (b >= 0.0 ? std::fabs(a) : -std::fabs(a)); }

double brent(const std::function<double(double)>& f, 
             double x1, double x2, double tol, int maxiter) {
  double a = x1, b = x2, c = x2;
  double fa = f(a), fb = f(b), fc = fb;
  double d = 0.0, d1 = 0.0;
  
  if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0)) 
    throw std::invalid_argument("Root must be bracketed in brent");
  
  for (int iter = 1; iter <= maxiter; ++iter) {
    if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
      c = a; fc = fa; d = b - a; d1 = d;
    }
    
    if (std::fabs(fc) < std::fabs(fb)) {
      a = b; b = c; c = a;
      fa = fb; fb = fc; fc = fa;
    }
    
    double tol1 = 2.0 * EPS * std::fabs(b) + 0.5 * tol;
    double xm = 0.5 * (c - b);
    if (std::fabs(xm) <= tol1 || fb == 0.0) return b;
    
    double p, q, r, s;
    if (std::fabs(d1) >= tol1 && std::fabs(fa) > std::fabs(fb)) {
      s = fb / fa;
      if (a == c) {
        p = 2.0 * xm * s;
        q = 1.0 - s;
      } else {
        q = fa / fc;
        r = fb / fc;
        p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0));
        q = (q - 1.0) * (r - 1.0) * (s - 1.0);
      }
      if (p > 0.0) q = -q;
      p = std::fabs(p);
      double min1 = 3.0 * xm * q - std::fabs(tol1 * q);
      double min2 = std::fabs(d1 * q);
      if (2.0 * p < (min1 < min2 ? min1 : min2)) {
        d1 = d; d = p / q;
      } else {
        d = xm; d1 = d;
      }
    } else {
      d = xm; d1 = d;
    }
    
    a = b; fa = fb;
    if (std::fabs(d) > tol1) b += d; else b += SIGN(tol1, xm);
    fb = f(b);
  }
  throw std::runtime_error("Maximum iterations exceeded in brent");
}

double bisect(const std::function<double(double)>& f,
              double x1, double x2, double tol, int maxiter) {
  double a = x1, b = x2;
  double fa = f(a), fb = f(b);
  if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0)) 
    throw std::invalid_argument("Root must be bracketed in bisect");
  if (std::fabs(fa) < tol) return a;
  if (std::fabs(fb) < tol) return b;
  double xmid, fmid;
  for (int j = 1; j <= maxiter; ++j) {
    xmid = a + 0.5 * (b - a);
    fmid = f(xmid);
    if (std::fabs(fmid) < tol || (b - a) < tol) return xmid;
    if ((fa > 0.0 && fmid < 0.0) || (fa < 0.0 && fmid > 0.0)) { 
      b = xmid; fb = fmid; }
    else { a = xmid; fa = fmid; }
  }
  throw std::runtime_error("Maximum number of iterations exceeded in bisect");
}

// --------------------------- Quantiles -------------------------------------

double quantilecpp(const std::vector<double>& x, double p) {
  int n = static_cast<int>(x.size());
  if (n == 0) throw std::invalid_argument("Empty vector");
  if (p < 0.0 || p > 1.0) throw std::invalid_argument("p must be in [0,1]");
  std::vector<double> y(x);
  std::sort(y.begin(), y.end());
  double h = (n - 1) * p + 1;
  int j = static_cast<int>(std::floor(h));
  double g = h - j;
  if (j <= 0) return y.front();
  if (j >= n) return y.back();
  return (1 - g) * y[j - 1] + g * y[j];
}

double squantilecpp(const std::function<double(double)>& S, double p, double tol) {
  if (p < 0.0 || p > 1.0) throw std::invalid_argument("p must be in [0,1]");
  double lower = 0.0, upper = 1.0;
  double Su = S(upper);
  while (Su > p) {
    lower = upper; upper *= 2.0; Su = S(upper);
    if (upper > 1e12) throw std::runtime_error(
        "Cannot find suitable upper bound for quantile search");
  }
  auto f = [&](double t) -> double { return S(t) - p; };
  return brent(f, lower, upper, tol);
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
  int n = data.nrows();
  int p = variables.size();
  ListCpp result;
  std::vector<int> nlevels(p);
  
  // IntMatrix for indices (n rows, p cols), column-major storage
  IntMatrix indices(n, p);
  
  // Flattened lookup buffers and per-variable metadata
  struct VarLookupInfo {
    int type; // 0=int, 1=double, 2=bool, 3=string
    int offset;
  };
  std::vector<VarLookupInfo> var_info(p);
  
  std::vector<int> int_flat;
  std::vector<double> dbl_flat;
  std::vector<unsigned char> bool_flat;
  std::vector<std::string> str_flat;
  
  ListCpp lookups_per_variable; // will contain a std::vector for each variable
  
  for (int i = 0; i < p; ++i) {
    const std::string& var = variables[i];
    if (!data.containElementNamed(var)) 
      throw std::invalid_argument("Data must contain variable: " + var);
    
    if (data.int_cols.count(var)) {
      const auto& col = data.int_cols.at(var);
      auto w = unique_sorted(col);
      nlevels[i] = w.size();
      auto idx = matchcpp(col, w); // indices 0..(levels-1)
      
      // append w to flat buffer and record metadata
      int off = int_flat.size();
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
      
      int off = dbl_flat.size();
      dbl_flat.insert(dbl_flat.end(), w.begin(), w.end());
      var_info[i] = VarLookupInfo{1, off};
      
      intmatrix_set_column(indices, i, idx);
      lookups_per_variable.push_back(std::move(w), var);
    } else if (data.bool_cols.count(var)) {
      const auto& col = data.bool_cols.at(var);
      auto w = unique_sorted(col);
      nlevels[i] = w.size();
      auto idx = matchcpp(col, w);
      
      int off = bool_flat.size();
      bool_flat.insert(bool_flat.end(), w.begin(), w.end());
      var_info[i] = VarLookupInfo{2, off};
      
      intmatrix_set_column(indices, i, idx);
      lookups_per_variable.push_back(std::move(w), var);
    } else if (data.string_cols.count(var)) {
      const auto& col = data.string_cols.at(var);
      auto w = unique_sorted(col);
      nlevels[i] = w.size();
      auto idx = matchcpp(col, w);
      
      int off = str_flat.size();
      str_flat.insert(str_flat.end(), w.begin(), w.end());
      var_info[i] = VarLookupInfo{3, off};
      
      intmatrix_set_column(indices, i, idx);
      lookups_per_variable.push_back(std::move(w), var);
    } else {
      throw std::invalid_argument("Unsupported variable type in bygroup: " + var);
    }
  } // end for variables
  
  // compute combined index
  std::vector<int> combined_index(n, 0);
  int orep = 1;
  for (int i = 0; i < p; ++i) orep *= nlevels[i];
  int lookup_nrows = orep;
  
  for (int i = 0; i < p; ++i) {
    orep /= nlevels[i];
    const int* col_ptr = indices.data_ptr() + i * n;
    for (int j = 0; j < n; ++j) {
      combined_index[j] += col_ptr[j] * orep;
    }
  }
  
  // Build lookup_df with columns repeated in the same pattern as before.
  DataFrameCpp lookup_df;
  int repeat_each = lookup_nrows;
  for (int i = 0; i < p; ++i) {
    const std::string& var = variables[i];
    int nlevels_i = nlevels[i];
    repeat_each /= nlevels_i;
    int times = lookup_nrows / ( nlevels_i * repeat_each );
    
    VarLookupInfo info = var_info[i];
    if (info.type == 0) {
      const int* base = int_flat.data() + info.offset;
      std::vector<int> col(lookup_nrows);
      int idxw = 0;
      for (int t = 0; t < times; ++t) {
        for (int level = 0; level < nlevels_i; ++level) 
          for (int r = 0; r < repeat_each; ++r) col[idxw++] = base[level];
      }
      lookup_df.push_back(std::move(col), var);
    } else if (info.type == 1) {
      const double* base = dbl_flat.data() + info.offset;
      std::vector<double> col(lookup_nrows);
      int idxw = 0;
      for (int t = 0; t < times; ++t) {
        for (int level = 0; level < nlevels_i; ++level) 
          for (int r = 0; r < repeat_each; ++r) col[idxw++] = base[level];
      }
      lookup_df.push_back(std::move(col), var);
    } else if (info.type == 2) {
      const unsigned char* base = bool_flat.data() + info.offset;
      std::vector<unsigned char> col(lookup_nrows);
      int idxw = 0;
      for (int t = 0; t < times; ++t) {
        for (int level = 0; level < nlevels_i; ++level) {
          for (int r = 0; r < repeat_each; ++r) col[idxw++] = base[level];
        }
      }
      lookup_df.push_back(std::move(col), var);
    } else { // string
      const std::string* base = str_flat.data() + info.offset;
      std::vector<std::string> col(lookup_nrows);
      int idxw = 0;
      for (int t = 0; t < times; ++t) {
        for (int level = 0; level < nlevels_i; ++level) 
          for (int r = 0; r < repeat_each; ++r) col[idxw++] = base[level];
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
  int m = A.nrow;
  int p = A.ncol;
  if (static_cast<int>(x.size()) != p) 
    throw std::invalid_argument("Vector size mismatch");
  std::vector<double> result(m, 0.0);
  for (int c = 0; c < p; ++c) {
    double xc = x[c];
    int offset = c * m;
    const double* colptr = A.data_ptr() + offset;
    for (int r = 0; r < m; ++r) result[r] += colptr[r] * xc;
  }
  return result;
}

FlatMatrix mat_mat_mult(const FlatMatrix& A, const FlatMatrix& B) {
  int m = A.nrow;
  int k = A.ncol;
  int k2 = B.nrow;
  int n = B.ncol;
  if (k != k2) throw std::invalid_argument("Matrix dimensions mismatch");
  if (m == 0 || k == 0 || n == 0) return FlatMatrix();
  FlatMatrix C(m, n);
  // Column-major: For each column j in B/C, compute 
  // C[:,j] = sum_{t=0..k-1} A[:,t] * B[t,j]
  for (int j = 0; j < n; ++j) {
    const double* bcol = B.data_ptr() + j * k;
    double* ccol = C.data_ptr() + j * m;
    for (int t = 0; t < k; ++t) {
      const double* acol = A.data_ptr() + t * m;
      double scale = bcol[t];
      if (scale == 0.0) continue;
      for (int i = 0; i < m; ++i) {
        ccol[i] += acol[i] * scale;
      }
    }
  }
  return C;
}

// Transpose a FlatMatrix (double)
FlatMatrix transpose(const FlatMatrix& M) {
  if (M.nrow == 0 || M.ncol == 0) return FlatMatrix();
  
  const int src_nrow = M.nrow;
  const int src_ncol = M.ncol;
  FlatMatrix out(src_ncol, src_nrow); // swapped dims
  
  const double* src = M.data_ptr();
  double* dst = out.data_ptr();
  
  for (int c = 0; c < src_ncol; ++c) {
    const double* src_col = src + c * src_nrow;
    for (int r = 0; r < src_nrow; ++r) {
      dst[r * src_ncol + c] = src_col[r];
    }
  }
  
  return out;
}

// Transpose an IntMatrix (int)
IntMatrix transpose(const IntMatrix& M) {
  if (M.nrow == 0 || M.ncol == 0) return IntMatrix();
  
  const int src_nrow = M.nrow;
  const int src_ncol = M.ncol;
  IntMatrix out(src_ncol, src_nrow); // swapped dims
  
  const int* src = M.data_ptr();
  int* dst = out.data_ptr();
  
  for (int c = 0; c < src_ncol; ++c) {
    const int* src_col = src + c * src_nrow;
    for (int r = 0; r < src_nrow; ++r) {
      dst[r * src_ncol + c] = src_col[r];
    }
  }
  
  return out;
}

// Transpose a BoolMatrix (unsigned char)
BoolMatrix transpose(const BoolMatrix& M) {
  if (M.nrow == 0 || M.ncol == 0) return BoolMatrix();
  
  const int src_nrow = M.nrow;
  const int src_ncol = M.ncol;
  BoolMatrix out(src_ncol, src_nrow); // swapped dims
  
  const unsigned char* src = M.data_ptr();
  unsigned char* dst = out.data_ptr();
  
  for (int c = 0; c < src_ncol; ++c) {
    const unsigned char* src_col = src + c * src_nrow;
    for (int r = 0; r < src_nrow; ++r) {
      dst[r * src_ncol + c] = src_col[r];
    }
  }
  
  return out;
}

double quadsym(const std::vector<double>& u, const FlatMatrix& v) {
  int p = u.size();
  const double* vptr = v.data_ptr();
  const double* uptr = u.data();
  double sum = 0.0;
  
  for (int j = 0; j < p; ++j) {
    const double* col = vptr + j * p;
    // diagonal term
    sum += uptr[j] * uptr[j] * col[j]; // col[j] == v(j,j)
    // off-diagonals i < j. Access column j contiguous for i = 0..j-1
    double s = 0.0;
    for (int i = 0; i < j; ++i) s += col[i] * uptr[i];
    sum += 2.0 * uptr[j] * s; // account for symmetric pair (i,j) and (j,i)
  }
  return sum;
}

// --------------------------- Linear algebra helpers (FlatMatrix-backed) ----
// cholesky2: in-place working on FlatMatrix (n x n), returns rank * nonneg
int cholesky2(FlatMatrix& matrix, int n, double toler) {
  double* base = matrix.data_ptr();
  double eps = 0.0;
  for (int i = 0; i < n; ++i) {
    double val = matrix(i, i);
    if (val > eps) eps = val;
  }
  if (eps == 0.0) eps = toler; else eps *= toler;
  int nonneg = 1;
  int rank = 0;
  
  for (int i = 0; i < n; ++i) {
    double* col_i = base + i * n;
    double pivot = col_i[i];
    if (std::isinf(pivot) || pivot < eps) {
      col_i[i] = 0.0;
      if (pivot < -8.0 * eps) nonneg = -1;
    } else {
      ++rank;
      for (int j = i + 1; j < n; ++j) {
        double* col_j = base + j * n;
        double temp = col_i[j] / pivot;
        col_i[j] = temp;
        col_j[j] -= temp * temp * pivot;
        for (int k = j + 1; k < n; ++k) {
          col_j[k] -= temp * col_i[k];
        }
      }
    }
  }
  return rank * nonneg;
}

// chsolve2 assumes matrix holds the representation produced by cholesky2
void chsolve2(FlatMatrix& matrix, int n, double* y) {
  // Forward substitution L * z = y
  double* base = matrix.data_ptr();
  for (int j = 0; j < n-1; ++j) {
    double yj = y[j];
    if (yj == 0.0) continue;
    double* col_j = base + j * n;
    for (int i = j + 1; i < n; ++i) {
      y[i] -= yj * col_j[i];
    }
  }
  // Now y holds z; solve L^T * x = z
  if (n == 0) return;
  for (int i = n - 1; i >= 0; --i) {
    double* col_i = base + i * n;
    double diag = col_i[i];
    if (diag == 0.0) {
      y[i] = 0.0;
    } else {
      double temp = y[i] / diag;
      for (int j = i + 1; j < n; ++j) temp -= y[j] * col_i[j];
      y[i] = temp;
    }
  }
}

FlatMatrix invsympd(const FlatMatrix& matrix, int n, double toler) {
  FlatMatrix v = matrix; // copy
  cholesky2(v, n, toler);
  FlatMatrix iv(n, n);
  for (int i = 0; i < n; ++i) {
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
  int n = tstart.size();
  int ncut = cut.size();
  int extra = 0;
  for (int i = 0; i < n; ++i) {
    if (std::isnan(tstart[i]) || std::isnan(tstop[i])) continue;
    for (int j = 0; j < ncut; ++j) {
      if (cut[j] > tstart[i] && cut[j] < tstop[i]) ++extra;
    }
  }
  int n2 = n + extra;
  std::vector<int> row(n2);
  std::vector<int> interval(n2);
  std::vector<double> start(n2);
  std::vector<double> end(n2);
  std::vector<int> censor(n2);
  int k = 0;
  for (int i = 0; i < n; ++i) {
    if (std::isnan(tstart[i]) || std::isnan(tstop[i])) {
      start[k] = tstart[i];
      end[k] = tstop[i];
      row[k] = i;
      interval[k] = 1;
      ++k;
    } else {
      int j = 0;
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

std::vector<double> house(const std::vector<double>& x) {
  double mu = std::sqrt(sumsq(x));
  std::vector<double> v = x; // copy
  if (mu > 0.0) {
    double beta = x[0] + std::copysign(mu, x[0]);
    int n = x.size();
    for (int i = 1; i < n; ++i) v[i] /= beta;
  }
  v[0] = 1.0;
  return v;
}

void row_house(FlatMatrix& A, int i1, int i2, int j1, int j2, 
               const std::vector<double>& v) {
  int m_total = A.nrow;
  if (m_total == 0) return;
  int n_total = A.ncol;
  if (i1 > i2 || i2 >= m_total) 
    throw std::invalid_argument("Invalid row indices i1 and i2");
  if (j1 > j2 || j2 >= n_total) 
    throw std::invalid_argument("Invalid column indices j1 and j2");
  int m = i2 - i1 + 1;
  int n = j2 - j1 + 1;
  double beta = -2.0 / sumsq(v);
  std::vector<double> w(n, 0.0);
  for (int j = 0; j < n; ++j) {
    double* acol = A.data_ptr() + (j1 + j) * m_total;
    double acc = 0.0;
    for (int i = 0; i < m; ++i) acc += acol[i1 + i] * v[i];
    w[j] = acc * beta;
    for (int i = 0; i < m; ++i) acol[i1 + i] += v[i] * w[j];
  }
}

ListCpp qrcpp(const FlatMatrix& X, double tol) {
  int m = X.nrow;
  int n = X.ncol;
  FlatMatrix A = X;
  std::vector<double> c(n, 0.0);
  for (int j = 0; j < n; ++j) {
    double* acol = A.data_ptr() + j * m;
    double s = 0.0;
    for (int i = 0; i < m; ++i) {
      double v = acol[i];
      s += v * v;
    }
    c[j] = s;
  }
  
  int r = -1;
  std::vector<int> piv(n);
  std::iota(piv.begin(), piv.end(), 0);
  double tau = *std::max_element(c.begin(), c.end());
  
  while (tau > tol) {
    ++r;
    int k = r;
    for (; k < n; ++k) {
      if (c[k] > tol) break;
    }
    if (k != r) {
      int off_r = r * m;
      int off_k = k * m;
      for (int i = 0; i < m; ++i) {
        std::swap(A.data[off_r + i], A.data[off_k + i]);
      }
      std::swap(piv[r], piv[k]);
      std::swap(c[r], c[k]);
    }
    
    int msub = m - r;
    std::vector<double> x(msub);
    for (int i = 0; i < msub; ++i) x[i] = A(r + i, r);
    std::vector<double> v = house(x);
    
    if (msub > 0 && r < n) row_house(A, r, m - 1, r, n - 1, v);
    
    for (int i = 1; i < msub; ++i) A(r + i, r) = v[i];
    for (int i = r + 1; i < n; ++i) {
      double val = A(r, i); c[i] -= val * val;
    }
    if (r < n - 1) {
      tau = *std::max_element(c.begin() + (r + 1), c.end());
    } else tau = 0.0;
  }
  
  FlatMatrix Qf(m, m);
  for (int c = 0; c < m; ++c) Qf(c, c) = 1.0;
  
  if (r >= 0) {
    for (int k = r; k >= 0; --k) {
      int msub_k = m - k;
      std::vector<double> vks(msub_k);
      vks[0] = 1.0;
      for (int i = 1; i < msub_k; ++i) {
        vks[i] = A(k + i, k);
      }
      row_house(Qf, k, m - 1, k, m - 1, vks);
    }
  }
  
  FlatMatrix Rf(m, n);
  std::fill(Rf.data.begin(), Rf.data.end(), 0.0);
  for (int j = 0; j < n; ++j) {
    for (int i = 0; i <= j && i < m; ++i) {
      Rf(i, j) = A(i, j);
    }
  }
  
  ListCpp result;
  result.push_back(std::move(A), "qr");
  result.push_back(r + 1, "rank");
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
  int i = 0, j = 0;
  int n1 = id1.size(), n2 = id2.size();
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
  int n = id.size();
  double a = std::exp(psi);
  std::vector<double> u_star(n), t_star(n);
  std::vector<int> d_star = event;
  for (int i = 0; i < n; ++i) { 
    u_star[i] = time[i] * ((1.0 - rx[i]) + rx[i] * a); t_star[i] = u_star[i]; }
  if (recensor) {
    std::vector<double> c_star = censor_time;
    for (int i = 0; i < n; ++i) c_star[i] *= std::min(1.0, a);
    if (autoswitch) {
      bool all_rx1 = true, all_rx0 = true;
      for (int i = 0; i < n; ++i) {
        if (treat[i] == 1 && rx[i] != 1.0) all_rx1 = false;
        if (treat[i] == 0 && rx[i] != 0.0) all_rx0 = false;
      }
      for (int i = 0; i < n; ++i) {
        if (treat[i] == 1 && all_rx1) 
          c_star[i] = std::numeric_limits<double>::infinity();
        if (treat[i] == 0 && all_rx0) 
          c_star[i] = std::numeric_limits<double>::infinity();
      }
    }
    for (int i = 0; i < n; ++i) {
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
  int n = id.size();
  double a0 = std::exp(psi);
  double a1 = std::exp(-psi);
  std::vector<double> u_star(n), t_star(n);
  std::vector<int> d_star = event;
  for (int i = 0; i < n; ++i) {
    if (treat[i] == 0) u_star[i] = time[i] * ((1.0 - rx[i]) + rx[i] * a0);
    else u_star[i] = time[i] * (rx[i] + (1.0 - rx[i]) * a1);
    t_star[i] = u_star[i];
  }
  if (recensor) {
    std::vector<double> c_star(n);
    for (int i = 0; i < n; ++i) 
      c_star[i] = treat[i] == 0 ? censor_time[i] * std::min(1.0, a0) : 
      censor_time[i] * std::min(1.0, a1);
    if (autoswitch) {
      bool all_rx1 = true, all_rx0 = true;
      for (int i = 0; i < n; ++i) {
        if (treat[i] == 1 && rx[i] != 1.0) all_rx1 = false;
        if (treat[i] == 0 && rx[i] != 0.0) all_rx0 = false;
      }
      for (int i = 0; i < n; ++i) {
        if (treat[i] == 1 && all_rx1) 
          c_star[i] = std::numeric_limits<double>::infinity();
        if (treat[i] == 0 && all_rx0) 
          c_star[i] = std::numeric_limits<double>::infinity();
      }
    }
    for (int i = 0; i < n; ++i) {
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

double qtpwexpcpp1(const double p,
                   const std::vector<double>& piecewiseSurvivalTime,
                   const std::vector<double>& lambda,
                   const double lowerBound,
                   const bool lowertail,
                   const bool logp) {
  int m = static_cast<int>(piecewiseSurvivalTime.size());
  double u = logp ? std::exp(p) : p;
  if (!lowertail) u = 1.0 - u;
  if (u <= 0.0) return lowerBound;
  if (u >= 1.0) return std::numeric_limits<double>::infinity();
  double v1 = -log1p(-u);
  int j = 0;
  while (j < m && piecewiseSurvivalTime[j] <= lowerBound) ++j;
  int j1 = (j == 0) ? 0 : (j - 1);
  double v = 0.0;
  if (j1 == m - 1) {
    double lj = lambda[j1];
    if (lj <= 0.0) return std::numeric_limits<double>::infinity();
    return lowerBound + v1 / lj;
  }
  for (j = j1; j < m - 1; ++j) {
    double dt = (j == j1) ? piecewiseSurvivalTime[j + 1] - lowerBound : 
    piecewiseSurvivalTime[j + 1] - piecewiseSurvivalTime[j];
    double lj = lambda[j];
    if (lj > 0.0) v += lj * dt;
    if (v >= v1) break;
  }
  double lj = lambda[j];
  if (lj <= 0.0) return std::numeric_limits<double>::infinity();
  if (j == m - 1) {
    double dt = (v1 - v) / lj;
    return piecewiseSurvivalTime[j] + dt;
  }
  double dt = (v - v1) / lj;
  return piecewiseSurvivalTime[j + 1] - dt;
}

// --------------------------- Root selection and helpers ---------------------

ListCpp getpsiest(double target,
                  const std::vector<double>& psi,
                  const std::vector<double>& Z,
                  int direction) {
  int n = psi.size();
  if (n != static_cast<int>(Z.size())) 
    throw std::invalid_argument("psi and Z must have the same length");
  if (n < 2) 
    throw std::invalid_argument("Need at least two points to find roots");
  std::vector<double> Zt(n);
  for (int i = 0; i < n; ++i) Zt[i] = Z[i] - target;
  std::vector<double> roots;
  for (int i = 1; i < n; ++i) {
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
      int m = roots.size();
      for (int j = 1; j < m; ++j) {
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
