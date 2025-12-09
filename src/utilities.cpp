#include "utilities.h"


// Normal distribution functions using Boost.Math
double boost_pnorm(double q, double mean, double sd, bool lower_tail) {
  if (sd <= 0) {
    throw std::invalid_argument("Standard deviation must be positive.");
  }
  
  boost::math::normal_distribution<> dist(mean, sd);
  double p = boost::math::cdf(dist, q);
  
  if (!lower_tail) {
    p = 1.0 - p;
  }
  
  return p;
}

double boost_qnorm(double p, double mean, double sd, bool lower_tail) {
  if (sd <= 0) {
    throw std::invalid_argument("Standard deviation must be positive.");
  }
  if (p < 0.0 || p > 1.0) {
    throw std::invalid_argument("Probability must be between 0 and 1.");
  }
  
  boost::math::normal_distribution<> dist(mean, sd);
  
  double q;
  if (lower_tail) {
    q = boost::math::quantile(dist, p);
  } else {
    q = boost::math::quantile(dist, 1.0 - p);
  }
  
  return q;
}

double boost_dnorm(double x, double mean, double sd) {
  if (sd <= 0) {
    throw std::invalid_argument("Standard deviation must be positive.");
  }
  
  boost::math::normal_distribution<> dist(mean, sd);
  return boost::math::pdf(dist, x);
}


// Logistic distribution functions using Boost.Math
double boost_plogis(double q, double location, double scale, bool lower_tail) {
  if (scale <= 0) {
    throw std::invalid_argument("Scale must be positive.");
  }
  
  boost::math::logistic_distribution<> dist(location, scale);
  double p = boost::math::cdf(dist, q);
  
  if (!lower_tail) {
    p = 1.0 - p;
  }
  
  return p;
}

double boost_qlogis(double p, double location, double scale, bool lower_tail) {
  if (scale <= 0) {
    throw std::invalid_argument("Scale must be positive.");
  }
  if (p < 0.0 || p > 1.0) {
    throw std::invalid_argument("Probability must be between 0 and 1.");
  }
  
  boost::math::logistic_distribution<> dist(location, scale);
  
  double q;
  if (lower_tail) {
    q = boost::math::quantile(dist, p);
  } else {
    q = boost::math::quantile(dist, 1.0 - p);
  }
  
  return q;
}

double boost_dlogis(double x, double location, double scale) {
  if (scale <= 0) {
    throw std::invalid_argument("Scale must be positive.");
  }
  
  boost::math::logistic_distribution<> dist(location, scale);
  return boost::math::pdf(dist, x);
}


// Extreme value distribution functions using Boost.Math
double boost_pextreme(double q, double location, double scale, bool lower_tail) {
  if (scale <= 0) {
    throw std::invalid_argument("Scale must be positive.");
  }
  
  boost::math::extreme_value_distribution<> dist(location, scale);
  double p = boost::math::cdf(complement(dist, -q));
  
  if (!lower_tail) {
    p = 1.0 - p;
  }
  
  return p;
}

double boost_qextreme(double p, double location, double scale, bool lower_tail) {
  if (scale <= 0) {
    throw std::invalid_argument("Scale must be positive.");
  }
  if (p < 0.0 || p > 1.0) {
    throw std::invalid_argument("Probability must be between 0 and 1.");
  }
  
  boost::math::extreme_value_distribution<> dist(location, scale);
  
  double q;
  if (lower_tail) {
    q = -boost::math::quantile(complement(dist, p));
  } else {
    q = -boost::math::quantile(complement(dist, 1.0 - p));
  }
  
  return q;
}

double boost_dextreme(double x, double location, double scale) {
  if (scale <= 0) {
    throw std::invalid_argument("Scale must be positive.");
  }
  
  boost::math::extreme_value_distribution<> dist(location, scale);
  return boost::math::pdf(dist, -x);
}


// Chi-squared distribution functions using Boost.Math
double boost_pchisq(double q, double df, bool lower_tail) {
  if (df <= 0) {
    throw std::invalid_argument("Degrees of freedom must be positive.");
  }
  
  boost::math::chi_squared_distribution<> dist(df);
  
  double p = boost::math::cdf(dist, q);
  
  if (!lower_tail) {
    p = 1.0 - p;
  }
  
  return p;
}


double boost_qchisq(double p, double df, bool lower_tail) {
  if (df <= 0) {
    throw std::invalid_argument("Degrees of freedom must be positive.");
  }
  if (p < 0.0 || p > 1.0) {
    throw std::invalid_argument("Probability must be between 0 and 1.");
  }
  
  boost::math::chi_squared_distribution<> dist(df);
  
  double q;
  if (lower_tail) {
    q = boost::math::quantile(dist, p);
  } else {
    q = boost::math::quantile(dist, 1.0 - p);
  }
  
  return q;
}


// Function to generate a sequence of integers from start to end (inclusive)
std::vector<int> seqcpp(int start, int end) {
  if (start > end) {
    throw std::invalid_argument("start must be less than or equal to end for the sequence function.");
  }
  
  // Calculate the size required for the vector (inclusive range: end - start + 1)
  size_t size = static_cast<size_t>(end - start + 1);
  
  // Create the vector and resize it
  std::vector<int> result(size);
  
  // Use std::iota to fill the vector sequentially starting from i0
  std::iota(result.begin(), result.end(), start);
  
  return result;
}


// Function to find the indices of all true elements in a boolean vector
std::vector<int> which(const std::vector<bool>& vec) {
  std::vector<int> true_indices;
  true_indices.reserve(vec.size()); // optional, improves performance
  
  for (size_t i = 0; i < vec.size(); ++i) {
    if (vec[i]) true_indices.push_back(static_cast<int>(i));
  }
  
  true_indices.shrink_to_fit(); // optional, reduce memory
  return true_indices;
}



// Return interval indices for each element in x relative to breakpoints v.
// Both x and v must be sorted in non-decreasing order.
std::vector<int> findInterval3(const std::vector<double>& x,
                               const std::vector<double>& v,
                               bool rightmost_closed,
                               bool all_inside,
                               bool left_open) {
  std::vector<int> out;
  out.resize(x.size());
  
  const double* v_begin = v.data();
  const double* v_end   = v_begin + v.size();
  const int nv = static_cast<int>(v.size());
  
  for (size_t i = 0; i < x.size(); ++i) {
    
    double xi = x[i];
    
    // Handle NaN in x (R equivalent: NA → NA_INTEGER)
    // You may instead choose: out[i] = -1 or continue
    if (std::isnan(xi)) {
      out[i] = -1; // placeholder for NA handling
      continue;
    }
    
    const double* pos;
    if (left_open) {
      pos = std::lower_bound(v_begin, v_end, xi);
    } else {
      pos = std::upper_bound(v_begin, v_end, xi);
    }
    
    int idx = static_cast<int>(pos - v_begin);
    
    if (rightmost_closed) {
      if (left_open) {
        if (nv > 0 && xi == v[0]) idx = 1;
      } else {
        if (nv > 0 && xi == v[nv - 1]) idx = nv - 1;
      }
    }
    
    if (all_inside) {
      if (idx == 0) {
        idx = 1;
      } else if (idx == nv) {
        idx = nv - 1;
      }
    }
    
    out[i] = idx;
  }
  
  return out;
}



static const double EPS = 3.0e-8;

inline double SIGN(double a, double b) {
  return (b >= 0.0 ? std::fabs(a) : -std::fabs(a));
}

// Brent's method to find root of f in [x1, x2] with tolerance tol
double brent(const std::function<double(double)>& f,
             double x1, double x2, double tol, int maxiter) {
  double a = x1, b = x2, c = x2;
  double fa = f(a), fb = f(b), fc = fb;
  double d = 0.0, d1 = 0.0;
  
  if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0)) {
    throw std::invalid_argument("Root must be bracketed in brent");
  }
  
  for (int iter = 1; iter <= maxiter; ++iter) {
    // Adjust bounds
    if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
      c = a;
      fc = fa;
      d = b - a;
      d1 = d;
    }
    
    // Ensure best guess is b
    if (std::fabs(fc) < std::fabs(fb)) {
      a = b; b = c; c = a;
      fa = fb; fb = fc; fc = fa;
    }
    
    // Convergence check
    double tol1 = 2.0 * EPS * std::fabs(b) + 0.5 * tol;
    double xm = 0.5 * (c - b);
    if (std::fabs(xm) <= tol1 || fb == 0.0) {
      return b;
    }
    
    double p, q, r, s;
    if (std::fabs(d1) >= tol1 && std::fabs(fa) > std::fabs(fb)) {
      // Inverse quadratic interpolation
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
        d1 = d;
        d = p / q;
      } else {
        d = xm;
        d1 = d;
      }
    } else {
      // Bisection fallback
      d = xm;
      d1 = d;
    }
    
    a = b;
    fa = fb;
    if (std::fabs(d) > tol1) {
      b += d;
    } else {
      b += SIGN(tol1, xm);
    }
    
    fb = f(b);
  }
  
  throw std::runtime_error("Maximum iterations exceeded in brent");
}


// Bisection method to find root of f in [x1, x2] with tolerance tol
double bisect(const std::function<double(double)>& f,
              double x1, double x2, double tol, int maxiter) {
  double a = x1, b = x2;
  double fa = f(a), fb = f(b);
  
  if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0)) {
    throw std::invalid_argument("Root must be bracketed in bisect");
  }
  
  double xmid, fmid;
  
  // Check if endpoints are already the root within tolerance
  if (std::fabs(fa) < tol) return a;
  if (std::fabs(fb) < tol) return b;
  
  for (int j = 1; j <= maxiter; ++j) {
    xmid = a + (b - a) * 0.5;
    fmid = f(xmid);
    
    // Convergence check: root found within tolerance
    if (std::fabs(fmid) < tol || (b - a) < tol) {
      return xmid;
    }
    
    if ((fa > 0.0 && fmid < 0.0) || (fa < 0.0 && fmid > 0.0)) {
      b = xmid;
      fb = fmid;
    } else {
      a = xmid;
      fa = fmid;
    }
  }
  
  throw std::runtime_error("Maximum number of iterations exceeded in bisect");
}


// Compute the p-th quantile of vector x using R type 7 method
double quantilecpp(const std::vector<double>& x, double p) {
  int n = static_cast<int>(x.size());
  if (n == 0) throw std::invalid_argument("Empty vector");
  if (p < 0.0 || p > 1.0) throw std::invalid_argument("p must be in [0,1]");
  
  // Copy input vector and sort it to find elements reliably
  std::vector<double> y(x);
  std::sort(y.begin(), y.end());
  
  // Compute fractional index (R type 7 formula: h = (n-1)*p + 1)
  double h = (n - 1) * p + 1;
  int j = static_cast<int>(std::floor(h)); // Integer part (1-based index)
  double g = h - j;                      // Fractional part
  
  // Interpolate between the j-th and (j+1)-th element (1-based indices).
  // In C++, these are 0-based indices j-1 and j.
  
  // Handle boundaries
  if (j <= 0) return y.front(); // If h <= 1, return min element
  if (j >= n) return y.back();  // If h >= n, return max element
  
  double lower_val = y[j - 1]; // Value at 0-based index j-1 (j-th element)
  double upper_val = y[j];     // Value at 0-based index j (j+1)-th element
  
  // Interpolate: (1 - g) * lower_val + g * upper_val
  return (1 - g) * lower_val + g * upper_val;
}


// Compute the quantile corresponding to survival function S at probability p
double squantilecpp(const std::function<double(double)>& S, double p, double tol) {
  if (p < 0.0 || p > 1.0) throw std::invalid_argument("p must be in [0,1]");
  
  double lower = 0.0;
  double upper = 1.0;
  double Su = S(upper);
  
  // Efficiently find an upper bound where S(upper) <= p
  while (Su > p) {
    lower = upper;
    upper *= 2.0;
    Su = S(upper);
    if (upper > 1e12)  // avoid infinite loop
      throw std::runtime_error("Cannot find suitable upper bound for quantile search");
  }
  
  // Define the function to find the root for
  auto f = [&S, p](double t) -> double {
    return S(t) - p;
  };
  
  return brent(f, lower, upper, tol);
}


// Main function to perform by-group indexing and return ListCpp
ListCpp bygroup(const DataFrameCpp& data, const std::vector<std::string>& variables) {
  int n = static_cast<int>(data.nrows());
  int p = static_cast<int>(variables.size());
  
  ListCpp result;
  std::vector<int> nlevels(p);
  std::vector<std::vector<int>> indices(n, std::vector<int>(p, 0));
  ListCpp lookups_per_variable;
  
  std::vector<std::vector<int>> int_lookups;
  std::vector<std::vector<double>> numeric_lookups;
  std::vector<std::vector<bool>> bool_lookups;
  std::vector<std::vector<std::string>> string_lookups;
  
  for (int i = 0; i < p; ++i) {
    const std::string& var = variables[i];
    if (!data.containElementNamed(var)) 
      throw std::invalid_argument("Data must contain variable: " + var);
    
    if (data.int_cols.count(var)) {
      auto col = data.int_cols.at(var);
      auto w = unique_sorted(col);
      nlevels[i] = static_cast<int>(w.size());
      auto idx = matchcpp(col, w);
      for (int j = 0; j < n; ++j) indices[j][i] = idx[j];
      
      DataFrameCpp df_uv;
      df_uv.push_back(w, var);
      lookups_per_variable.push_back(df_uv, var);
      
      int_lookups.push_back(w);
    } else if (data.numeric_cols.count(var)) {
      auto col = data.numeric_cols.at(var);
      auto w = unique_sorted(col);
      nlevels[i] = static_cast<int>(w.size());
      auto idx = matchcpp(col, w);
      for (int j = 0; j < n; ++j) indices[j][i] = idx[j];
      
      DataFrameCpp df_uv;
      df_uv.push_back(w, var);
      lookups_per_variable.push_back(df_uv, var);
      
      numeric_lookups.push_back(w);
    } else if (data.bool_cols.count(var)) {
      auto col = data.bool_cols.at(var);
      auto w = unique_sorted(col);
      nlevels[i] = static_cast<int>(w.size());
      auto idx = matchcpp(col, w);
      for (int j = 0; j < n; ++j) indices[j][i] = idx[j];
      
      DataFrameCpp df_uv;
      df_uv.push_back(w, var);
      lookups_per_variable.push_back(df_uv, var);
      
      bool_lookups.push_back(w);
    } else if (data.string_cols.count(var)) {
      auto col = data.string_cols.at(var);
      auto w = unique_sorted(col);
      nlevels[i] = static_cast<int>(w.size());
      auto idx = matchcpp(col, w);
      for (int j = 0; j < n; ++j) indices[j][i] = idx[j];
      
      DataFrameCpp df_uv;
      df_uv.push_back(w, var);
      lookups_per_variable.push_back(df_uv, var);
      
      string_lookups.push_back(w);
    } else {
      throw std::invalid_argument("Unsupported variable type in bygroup: " + var);
    }
  }
  
  std::vector<int> combined_index(n, 0);
  int orep = 1;
  for (int i = 0; i < p; ++i) orep *= nlevels[i];
  
  for (int i = 0; i < p; ++i) {
    orep /= nlevels[i];
    for (int j = 0; j < n; ++j) {
      combined_index[j] += indices[j][i] * orep;
    }
  }
  
  DataFrameCpp lookup_df;
  int lookup_nrows = 1;
  for (int i = 0; i < p; ++i) lookup_nrows *= nlevels[i];
  
  for (int i = 0; i < p; ++i) {
    const std::string& var = variables[i];
    int nlevels_i = nlevels[i];
    int repeat_each = 1;
    for (int j = i + 1; j < p; ++j) repeat_each *= nlevels[j];
    int times = lookup_nrows / (nlevels_i * repeat_each);
    
    if (data.int_cols.count(var)) {
      auto& w = int_lookups.front();
      std::vector<int> col(lookup_nrows);
      int idx = 0;
      for (int t = 0; t < times; ++t) {
        for (int level = 0; level < nlevels_i; ++level) {
          for (int r = 0; r < repeat_each; ++r) {
            col[idx++] = w[level];
          }
        }
      }
      lookup_df.push_back(col, var);
      int_lookups.erase(int_lookups.begin());
    } else if (data.numeric_cols.count(var)) {
      auto& w = numeric_lookups.front();
      std::vector<double> col(lookup_nrows);
      int idx = 0;
      for (int t = 0; t < times; ++t) {
        for (int level = 0; level < nlevels_i; ++level) {
          for (int r = 0; r < repeat_each; ++r) {
            col[idx++] = w[level];
          }
        }
      }
      lookup_df.push_back(col, var);
      numeric_lookups.erase(numeric_lookups.begin());
    } else if (data.bool_cols.count(var)) {
      auto& w = bool_lookups.front();
      std::vector<bool> col(lookup_nrows);
      int idx = 0;
      for (int t = 0; t < times; ++t) {
        for (int level = 0; level < nlevels_i; ++level) {
          for (int r = 0; r < repeat_each; ++r) {
            col[idx++] = w[level];
          }
        }
      }
      lookup_df.push_back(col, var);
      bool_lookups.erase(bool_lookups.begin());
    } else if (data.string_cols.count(var)) {
      auto& w = string_lookups.front();
      std::vector<std::string> col(lookup_nrows);
      int idx = 0;
      for (int t = 0; t < times; ++t) {
        for (int level = 0; level < nlevels_i; ++level) {
          for (int r = 0; r < repeat_each; ++r) {
            col[idx++] = w[level];
          }
        }
      }
      lookup_df.push_back(col, var);
      string_lookups.erase(string_lookups.begin());
    }
  }
  
  result.push_back(nlevels, "nlevels");
  result.push_back(indices, "indices");
  result.push_back(lookups_per_variable, "lookups_per_variable");
  result.push_back(combined_index, "index");
  result.push_back(lookup_df, "lookup");
  
  return result;
}


// Helper to perform matrix-vector multiplication.
std::vector<double> mat_vec_mult(const std::vector<std::vector<double>>& A,
                                 const std::vector<double>& x) {
  size_t n = A.size();
  if (n == 0) return {};
  size_t p = (n > 0) ? A[0].size() : 0;
  if (x.size() != p) throw std::invalid_argument("Vector size mismatch");
  
  std::vector<double> result(n, 0.0);
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < p; ++j) {
      result[i] += A[i][j] * x[j];
    }
  }
  return result;
}

// Helper to perform matrix-matrix multiplication.
std::vector<std::vector<double>> mat_mat_mult(const std::vector<std::vector<double>>& A,
                                              const std::vector<std::vector<double>>& B) {
  size_t rows_A = A.size();
  if (rows_A == 0) return {};
  size_t cols_A = A[0].size();
  size_t rows_B = B.size();
  if (rows_B == 0) return {};
  size_t cols_B = B[0].size();
  if (cols_A != rows_B) 
    throw std::invalid_argument("Matrix dimensions mismatch");
  
  std::vector<std::vector<double>> C(rows_A, std::vector<double>(cols_B, 0.0));
  for (size_t i = 0; i < rows_A; ++i) {
    for (size_t j = 0; j < cols_B; ++j) {
      for (size_t k = 0; k < cols_A; ++k) {
        C[i][j] += A[i][k] * B[k][j];
      }
    }
  }
  return C;
}

// Helper to transpose a matrix.
std::vector<std::vector<double>> transpose(const std::vector<std::vector<double>>& A) {
  size_t rows = A.size();
  if (rows == 0) return {};
  size_t cols = A[0].size();
  std::vector<std::vector<double>> At(cols, std::vector<double>(rows));
  for (size_t i = 0; i < rows; ++i) {
    for (size_t j = 0; j < cols; ++j) {
      At[j][i] = A[i][j];
    }
  }
  return At;
}


// The following three utilities functions are from the survival package

// Cholesky decomposition with tolerance, returns rank*nonneg
int cholesky2(std::vector<std::vector<double>>& matrix, int n, double toler) {
  double eps = 0.0;
  
  // Find largest diagonal element
  for (int i = 0; i < n; ++i) {
    if (matrix[i][i] > eps) eps = matrix[i][i];
  }
  
  if (eps == 0.0) eps = toler;
  else eps *= toler;
  
  int nonneg = 1;
  int rank = 0;
  
  for (int i = 0; i < n; ++i) {
    double pivot = matrix[i][i];
    
    if (std::isinf(pivot) || pivot < eps) {
      matrix[i][i] = 0.0;
      if (pivot < -8.0 * eps) nonneg = -1;
    } else {
      ++rank;
      for (int j = i + 1; j < n; ++j) {
        double temp = matrix[i][j] / pivot;
        matrix[i][j] = temp;
        matrix[j][j] -= temp * temp * pivot;
        for (int k = j + 1; k < n; ++k) {
          matrix[j][k] -= temp * matrix[i][k];
        }
      }
    }
  }
  
  return rank * nonneg;
}


// Solve system after cholesky2 decomposition
void chsolve2(std::vector<std::vector<double>>& matrix, int n, std::vector<double>& y) {
  // Forward substitution: solve L * z = y
  for (int i = 0; i < n; ++i) {
    double temp = y[i];
    for (int j = 0; j < i; ++j) {
      temp -= y[j] * matrix[j][i];
    }
    y[i] = temp;
  }
  
  // Backward substitution: solve L^T * x = z
  for (int i = n - 1; i >= 0; --i) {
    if (matrix[i][i] == 0.0) {
      y[i] = 0.0;
    } else {
      double temp = y[i] / matrix[i][i];
      for (int j = i + 1; j < n; ++j) {
        temp -= y[j] * matrix[i][j];
      }
      y[i] = temp;
    }
  }
}


// Invert a matrix after cholesky2 decomposition
void chinv2(std::vector<std::vector<double>>& matrix, int n) {
  // Step 1: invert diagonal and apply sweep operator
  for (int i = 0; i < n; ++i) {
    if (matrix[i][i] > 0.0) {
      matrix[i][i] = 1.0 / matrix[i][i]; // invert D
      for (int j = i + 1; j < n; ++j) {
        matrix[i][j] = -matrix[i][j];
        for (int k = 0; k < i; ++k) {
          matrix[k][j] += matrix[i][j] * matrix[k][i];
        }
      }
    }
  }
  
  // Step 2: finalize inverse
  for (int i = 0; i < n; ++i) {
    if (matrix[i][i] == 0.0) { // singular row
      for (int j = 0; j < i; ++j) matrix[i][j] = 0.0;
      for (int j = i; j < n; ++j) matrix[j][i] = 0.0;
    } else {
      for (int j = i + 1; j < n; ++j) {
        double temp = matrix[i][j] * matrix[j][j];
        matrix[j][i] = temp;
        for (int k = i; k < j; ++k) {
          matrix[k][i] += temp * matrix[k][j];
        }
      }
    }
  }
}


// Invert a symmetric positive definite matrix
std::vector<std::vector<double>> invsympd(const std::vector<std::vector<double>>& matrix, int n, double toler) {
  // Clone the input matrix
  std::vector<std::vector<double>> v = matrix;
  
  // Step 1: Cholesky decomposition (in-place)
  cholesky2(v, n, toler);
  
  // Step 2: Invert using sweep operator
  chinv2(v, n);
  
  // Step 3: Fill lower triangle from upper triangle to make symmetric
  for (int i = 1; i < n; ++i) {
    for (int j = 0; j < i; ++j) {
      v[j][i] = v[i][j];
    }
  }
  
  return v;
}


// Split survival data into intervals based on cut points
DataFrameCpp survsplit(const std::vector<double>& tstart,
                       const std::vector<double>& tstop,
                       const std::vector<double>& cut) {
  int n = static_cast<int>(tstart.size());
  int ncut = static_cast<int>(cut.size());
  
  // Count extra rows
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
  std::vector<int> censor(n2, 0); // use int instead of bool for consistency
  
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
  
  // Construct the DataFrameCpp output
  DataFrameCpp df;
  df.push_back(row, "row");
  df.push_back(start, "start");
  df.push_back(end, "end");
  df.push_back(censor, "censor");
  df.push_back(interval, "interval");
  
  return df;
}




// Compute the Euclidean norm squared
double sumsq(const std::vector<double>& x) {
  double s = 0.0;
  for (double xi : x) s += xi*xi;
  return s;
}

// Householder vector
std::vector<double> house(const std::vector<double>& x) {
  int n = static_cast<int>(x.size());
  double mu = std::sqrt(sumsq(x));
  std::vector<double> v = x; // copy
  if (mu > 0.0) {
    double beta = x[0] + std::copysign(mu, x[0]);
    for (int i = 1; i < n; ++i) {
      v[i] /= beta;
    }
  }
  v[0] = 1.0;
  return v;
}

// Row Householder pre-multiplication
// A: m-by-n matrix represented as vector of vectors
void row_house(std::vector<std::vector<double>>& A,
               int i1, int i2,
               int j1, int j2,
               const std::vector<double>& v) {
  int m_total = static_cast<int>(A.size());
  if (m_total == 0) return;
  int n_total = static_cast<int>(A[0].size());
  
  if (i1 < 0 || i1 > i2 || i2 >= m_total) 
    throw std::invalid_argument("Invalid row indices i1 and i2");
  if (j1 < 0 || j1 > j2 || j2 >= n_total)
    throw std::invalid_argument("Invalid column indices j1 and j2");
  
  int m = i2 - i1 + 1;
  int n = j2 - j1 + 1;
  
  double beta = -2.0 / sumsq(v);
  
  std::vector<double> w(n, 0.0);
  
  for (int j = 0; j < n; ++j) {
    for (int i = 0; i < m; ++i) {
      w[j] += A[i + i1][j + j1] * v[i];
    }
    w[j] *= beta;
  }
  
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < n; ++j) {
      A[i + i1][j + j1] += v[i] * w[j];
    }
  }
}


// Maximum element in a vector
double max_elem(const std::vector<double>& x, int start, int end) {
  if (start > end || start < 0 || end >= static_cast<int>(x.size())) {
    throw std::invalid_argument("Invalid start or end indices in max_elem");
  }
  return *std::max_element(x.begin() + start, x.begin() + end + 1);
}


// QR decomposition with Householder and column pivoting
ListCpp qrcpp(const std::vector<std::vector<double>>& X, double tol) {
  int m = static_cast<int>(X.size());
  int n = static_cast<int>(X[0].size());
  std::vector<std::vector<double>> A = X;
  
  // Compute squared column norms
  std::vector<double> c(n);
  for (int j = 0; j < n; ++j) {
    double s = 0.0;
    for (int i = 0; i < m; ++i) s += A[i][j]*A[i][j];
    c[j] = s;
  }
  
  // Initial pivoting
  int k = 0;
  for (; k < n; ++k) if (c[k] > tol) break;
  
  int r = -1;
  std::vector<int> piv(n);
  std::iota(piv.begin(), piv.end(), 0);
  
  double tau = max_elem(c, 0, n-1);
  
  while (tau > tol) {
    ++r;
    
    // Exchange column r with column k
    std::swap(piv[r], piv[k]);
    for (int i = 0; i < m; ++i) std::swap(A[i][r], A[i][k]);
    std::swap(c[r], c[k]);
    
    // Householder vector
    std::vector<double> v(m-r);
    for (int i = 0; i < m-r; ++i) v[i] = A[i+r][r];
    v = house(v);
    
    // Apply Householder
    row_house(A, r, m-1, r, n-1, v);
    
    // Store Householder vector in sub-diagonal
    for (int i = 1; i < m-r; ++i) A[i+r][r] = v[i];
    
    // Update squared norms
    for (int i = r+1; i < n; ++i) c[i] -= A[r][i]*A[r][i];
    
    // Next pivot column
    if (r < n-1) {
      tau = max_elem(c, r+1, n-1);
      for (k = r+1; k < n; ++k) if (c[k] > tol) break;
    } else {
      tau = 0.0;
    }
  }
  
  // Recover Q
  std::vector<std::vector<double>> Q(m, std::vector<double>(m, 0.0));
  for (int i = 0; i < m; ++i) Q[i][i] = 1.0;
  for (int k = r; k >= 0; --k) {
    std::vector<double> v(m-k);
    v[0] = 1.0;
    for (int i = 1; i < m-k; ++i) v[i] = A[i+k][k];
    row_house(Q, k, m-1, k, m-1, v);
  }
  
  // Recover R
  std::vector<std::vector<double>> R(m, std::vector<double>(n, 0.0));
  for (int j = 0; j < n; ++j) {
    for (int i = 0; i <= j && i < m; ++i) R[i][j] = A[i][j];
  }
  
  // Create and populate the ListCpp object
  ListCpp result;
  result.push_back(A, "qr");
  result.push_back(r + 1, "rank");
  result.push_back(piv, "pivot");
  result.push_back(Q, "Q");
  result.push_back(R, "R");
  
  return result;
}



// Return -1 if no match is found
std::vector<int> match3(const std::vector<int>& id1,
                        const std::vector<double>& v1,
                        const std::vector<int>& id2,
                        const std::vector<double>& v2) {
  std::vector<int> result;
  result.reserve(id1.size());
  
  size_t i = 0, j = 0;
  size_t n1 = id1.size(), n2 = id2.size();
  
  while (i < n1 && j < n2) {
    // Condition 1: item1 < item2 (move i forward, mark as -1)
    if (id1[i] < id2[j] || (id1[i] == id2[j] && v1[i] < v2[j])) {
      result.push_back(-1);
      ++i;
    } 
    // Condition 2: item1 > item2 (move j forward)
    else if (id1[i] > id2[j] || (id1[i] == id2[j] && v1[i] > v2[j])) {
      ++j;
    } 
    // Condition 3: item1 == item2 (match found, move both i and j forward)
    else {
      result.push_back(static_cast<int>(j));
      ++i;
      ++j;
    }
  }
  
  // Handle remaining i elements which have no possible match in the rest of j
  while (i < n1) {
    result.push_back(-1);
    ++i;
  }
  
  return result;
}



// Untreated counterfactual time calculation
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
    u_star[i] = time[i] * ((1.0 - rx[i]) + rx[i] * a);
    t_star[i] = u_star[i];
  }
  
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
        // If the whole group had the same RX value, 
        // ignore censoring for that group
        if (treat[i] == 1 && all_rx1) 
          c_star[i] = std::numeric_limits<double>::infinity();
        if (treat[i] == 0 && all_rx0) 
          c_star[i] = std::numeric_limits<double>::infinity();
      }
    }
    
    for (size_t i = 0; i < n; ++i) {
      // Apply the new censoring if c* < u*
      if (c_star[i] < u_star[i]) {
        t_star[i] = c_star[i];
        d_star[i] = 0; // Event indicator becomes 0 (censored)
      }
    }
  }
  
  DataFrameCpp df;
  df.push_back(id, "uid");
  df.push_back(t_star, "t_star");
  df.push_back(d_star, "d_star");
  df.push_back(treat, "treated");
  
  return df;
}


// Unswitched counterfactual time calculation
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
  
  // Calculate the underlying counterfactual time u* based on 
  // assigned treatment group
  for (size_t i = 0; i < n; ++i) {
    if (treat[i] == 0) // Control group counterfactual time
      u_star[i] = time[i] * ((1.0 - rx[i]) + rx[i] * a0);
    else // Treated group counterfactual time
      u_star[i] = time[i] * (rx[i] + (1.0 - rx[i]) * a1);
    
    t_star[i] = u_star[i];
  }
  
  if (recensor) {
    std::vector<double> c_star(n);
    for (size_t i = 0; i < n; ++i) {
      // Scale censor times based on group-specific 'a' (a0 or a1)
      c_star[i] = treat[i] == 0 ? censor_time[i] * std::min(1.0, a0)
        : censor_time[i] * std::min(1.0, a1);
    }
    
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
      // Apply the new censoring if c* < u*
      if (c_star[i] < u_star[i]) {
        t_star[i] = c_star[i];
        d_star[i] = 0; // Event indicator becomes 0 (censored)
      }
    }
  }
  
  DataFrameCpp df;
  df.push_back(id, "uid");
  df.push_back(t_star, "t_star");
  df.push_back(d_star, "d_star");
  df.push_back(treat, "treated");
  
  return df;
}


// Sanitize a string: replace non-alphanumeric and non-underscore with dot
std::string sanitize(const std::string& s) {
  std::string out = s;
  for (char &c : out) {
    if (!std::isalnum(static_cast<unsigned char>(c)) && 
        static_cast<unsigned char>(c) != '_') {
      c = '.';
    }
  }
  return out;
}


// Piecewise exponential quantile function
double qtpwexpcpp1(const double p,
                   const std::vector<double>& piecewiseSurvivalTime,
                   const std::vector<double>& lambda,
                   const double lowerBound,
                   const bool lowertail,
                   const bool logp) {
  int m = static_cast<int>(piecewiseSurvivalTime.size());
  
  if (m == 0 || static_cast<int>(lambda.size()) != m)
    throw std::invalid_argument("Invalid piecewise model inputs.");
  
  // Convert p -> cumulative probability above lowerBound
  double u;
  if (logp)
    u = std::exp(p);
  else
    u = p;
  
  if (!lowertail)
    u = 1.0 - u;
  
  // Bound u away from 0 and 1 for numerical stability
  if (u <= 0.0) return lowerBound;
  if (u >= 1.0) return std::numeric_limits<double>::infinity();
  
  double v1 = -log1p(-u);  // better numerics than -log(1-u)
  
  // identify interval for lowerBound
  int j = 0;
  while (j < m && piecewiseSurvivalTime[j] <= lowerBound) j++;
  int j1 = std::max(0, j - 1);
  
  double v = 0.0;
  
  // If starting in the last interval
  if (j1 == m - 1) {
    double lj = lambda[j1];
    if (lj <= 0.0) return std::numeric_limits<double>::infinity();
    return lowerBound + v1 / lj;
  }
  
  // accumulate hazards until reaching target v1
  for (j = j1; j < m - 1; ++j) {
    double dt = (j == j1)
    ? piecewiseSurvivalTime[j + 1] - lowerBound
    : piecewiseSurvivalTime[j + 1] - piecewiseSurvivalTime[j];
    
    double lj = lambda[j];
    if (lj > 0.0)
      v += lj * dt;
    
    if (v >= v1) break;
  }
  
  // solve within the interval j
  double lj = lambda[j];
  if (lj <= 0.0)
    return std::numeric_limits<double>::infinity();
  
  if (j == m - 1) {
    // last interval
    double dt = (v1 - v) / lj;
    return piecewiseSurvivalTime[j] + dt;
  }
  
  // inside interval j before end of grid
  double dt = (v - v1) / lj;
  return piecewiseSurvivalTime[j + 1] - dt;
}


// find the roots of an implicit function given a set of discrete sample points
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
  for (int i = 0; i < n; ++i)
    Zt[i] = Z[i] - target;
  
  std::vector<double> roots;
  
  for (int i = 1; i < n; ++i) {
    double z1 = Zt[i - 1];
    double z2 = Zt[i];
    
    if (std::isnan(z1) || std::isnan(z2) || z1 == z2)
      continue;
    
    if (z1 == 0.0)
      // If the exact point is a zero
      roots.push_back(psi[i - 1]);
    else if (z1 * z2 < 0.0) {
      // If a sign change (zero crossing) occurs
      double psi_root = psi[i - 1] - z1 * (psi[i] - psi[i - 1]) / (z2 - z1);
      roots.push_back(psi_root);
    }
  }
  
  double root = NAN;
  
  if (!roots.empty()) {
    if (direction == -1) {
      // Leftmost root
      root = roots.front();
    }
    else if (direction == 1) {
      // Rightmost
      root = roots.back();
    }
    else {
      // Closest to zero
      root = roots[0];
      double minabs = std::abs(roots[0]);
      for (size_t j = 1; j < roots.size(); ++j) {
        double a = std::abs(roots[j]);
        if (a < minabs) {
          minabs = a;
          root = roots[j];
        }
      }
    }
  }
  
  ListCpp result;
  result.push_back(roots, "all_roots");
  result.push_back(root, "selected_root");
  return result;
}


// Find an endpoint psi value where f(psiend) has the desired sign
double getpsiend(const std::function<double(double)>& f,
                 bool lowerend,
                 double initialend) {
  double psiend = initialend;
  double zend = f(initialend);
  const double LIMIT = 10.0;
  
  if (lowerend) {
    if ((std::isinf(zend) && zend > 0) || std::isnan(zend)) {
      while (((std::isinf(zend) && zend > 0) || std::isnan(zend)) &&
             psiend <= LIMIT) {
        psiend += 1;
        zend = f(psiend);
      }
      if (psiend > LIMIT)
        return NAN;
    }
    
    if (zend < 0) {
      while (!std::isinf(zend) && zend < 0 && psiend >= -LIMIT) {
        psiend -= 1;
        zend = f(psiend);
      }
      if (std::isinf(zend) || std::isnan(zend) || psiend < -LIMIT)
        return NAN;
    }
  }
  else { // upper end
    if ((std::isinf(zend) && zend < 0) || std::isnan(zend)) {
      while (((std::isinf(zend) && zend < 0) || std::isnan(zend)) &&
             psiend >= -LIMIT) {
        psiend -= 1;
        zend = f(psiend);
      }
      if (psiend < -LIMIT)
        return NAN;
    }
    
    if (zend > 0) {
      while (!std::isinf(zend) && zend > 0 && psiend <= LIMIT) {
        psiend += 1;
        zend = f(psiend);
      }
      if (std::isinf(zend) || std::isnan(zend) || psiend > LIMIT)
        return NAN;
    }
  }
  
  return psiend;
}