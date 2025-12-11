#include "utilities.h"
#include "dataframe_list.h" // DataFrameCpp, ListCpp

#include <vector>
#include <string>
#include <algorithm>   // sort, upper_bound, max_element, etc.
#include <numeric>     // iota, accumulate, inner_product
#include <functional>  // std::function
#include <cmath>       // isnan, isinf, fabs, pow, log, exp, sqrt
#include <limits>      // numeric_limits
#include <stdexcept>   // exceptions
#include <cstring>     // std::memcpy
#include <memory>

#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/logistic.hpp>
#include <boost/math/distributions/extreme_value.hpp>
#include <boost/math/distributions/chi_squared.hpp>

// --------------------------- Distribution helpers --------------------------

double boost_pnorm(double q, double mean, double sd, bool lower_tail) {
  if (sd <= 0) throw std::invalid_argument("Standard deviation must be positive.");
  boost::math::normal_distribution<> dist(mean, sd);
  double p = boost::math::cdf(dist, q);
  return lower_tail ? p : (1.0 - p);
}

double boost_qnorm(double p, double mean, double sd, bool lower_tail) {
  if (sd <= 0) throw std::invalid_argument("Standard deviation must be positive.");
  if (p < 0.0 || p > 1.0) throw std::invalid_argument("Probability must be between 0 and 1.");
  boost::math::normal_distribution<> dist(mean, sd);
  return lower_tail ? boost::math::quantile(dist, p) : boost::math::quantile(dist, 1.0 - p);
}

double boost_dnorm(double x, double mean, double sd) {
  if (sd <= 0) throw std::invalid_argument("Standard deviation must be positive.");
  boost::math::normal_distribution<> dist(mean, sd);
  return boost::math::pdf(dist, x);
}

double boost_plogis(double q, double location, double scale, bool lower_tail) {
  if (scale <= 0) throw std::invalid_argument("Scale must be positive.");
  boost::math::logistic_distribution<> dist(location, scale);
  double p = boost::math::cdf(dist, q);
  return lower_tail ? p : (1.0 - p);
}

double boost_qlogis(double p, double location, double scale, bool lower_tail) {
  if (scale <= 0) throw std::invalid_argument("Scale must be positive.");
  if (p < 0.0 || p > 1.0) throw std::invalid_argument("Probability must be between 0 and 1.");
  boost::math::logistic_distribution<> dist(location, scale);
  return lower_tail ? boost::math::quantile(dist, p) : boost::math::quantile(dist, 1.0 - p);
}

double boost_dlogis(double x, double location, double scale) {
  if (scale <= 0) throw std::invalid_argument("Scale must be positive.");
  boost::math::logistic_distribution<> dist(location, scale);
  return boost::math::pdf(dist, x);
}

double boost_pextreme(double q, double location, double scale, bool lower_tail) {
  if (scale <= 0) throw std::invalid_argument("Scale must be positive.");
  boost::math::extreme_value_distribution<> dist(location, scale);
  // note: original code used complement and -q; keep semantics consistent with previous implementation
  double p = boost::math::cdf(complement(dist, -q));
  return lower_tail ? p : (1.0 - p);
}

double boost_qextreme(double p, double location, double scale, bool lower_tail) {
  if (scale <= 0) throw std::invalid_argument("Scale must be positive.");
  if (p < 0.0 || p > 1.0) throw std::invalid_argument("Probability must be between 0 and 1.");
  boost::math::extreme_value_distribution<> dist(location, scale);
  double q;
  if (lower_tail) q = -boost::math::quantile(complement(dist, p));
  else q = -boost::math::quantile(complement(dist, 1.0 - p));
  return q;
}

double boost_dextreme(double x, double location, double scale) {
  if (scale <= 0) throw std::invalid_argument("Scale must be positive.");
  boost::math::extreme_value_distribution<> dist(location, scale);
  return boost::math::pdf(dist, -x);
}

double boost_pchisq(double q, double df, bool lower_tail) {
  if (df <= 0) throw std::invalid_argument("Degrees of freedom must be positive.");
  boost::math::chi_squared_distribution<> dist(df);
  double p = boost::math::cdf(dist, q);
  return lower_tail ? p : (1.0 - p);
}

double boost_qchisq(double p, double df, bool lower_tail) {
  if (df <= 0) throw std::invalid_argument("Degrees of freedom must be positive.");
  if (p < 0.0 || p > 1.0) throw std::invalid_argument("Probability must be between 0 and 1.");
  boost::math::chi_squared_distribution<> dist(df);
  return lower_tail ? boost::math::quantile(dist, p) : boost::math::quantile(dist, 1.0 - p);
}

// --------------------------- Small utilities --------------------------------

std::vector<int> seqcpp(int start, int end) {
  if (start > end) throw std::invalid_argument("start must be less than or equal to end for the sequence function.");
  size_t size = static_cast<size_t>(end - start + 1);
  std::vector<int> result(size);
  std::iota(result.begin(), result.end(), start);
  return result;
}

std::vector<int> which(const std::vector<bool>& vec) {
  std::vector<int> true_indices;
  true_indices.reserve(vec.size());
  for (size_t i = 0; i < vec.size(); ++i) if (vec[i]) true_indices.push_back(static_cast<int>(i));
  true_indices.shrink_to_fit();
  return true_indices;
}

// findInterval3: adapted from previous implementation
std::vector<int> findInterval3(const std::vector<double>& x,
                               const std::vector<double>& v,
                               bool rightmost_closed,
                               bool all_inside,
                               bool left_open) {
  std::vector<int> out(x.size());
  const double* v_begin = v.data();
  const double* v_end   = v_begin + v.size();
  const int nv = static_cast<int>(v.size());
  
  for (size_t i = 0; i < x.size(); ++i) {
    double xi = x[i];
    if (std::isnan(xi)) { out[i] = -1; continue; }
    const double* pos = left_open ? std::lower_bound(v_begin, v_end, xi) : std::upper_bound(v_begin, v_end, xi);
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
inline double SIGN(double a, double b) { return (b >= 0.0 ? std::fabs(a) : -std::fabs(a)); }

double brent(const std::function<double(double)>& f,
             double x1, double x2, double tol, int maxiter) {
  double a = x1, b = x2, c = x2;
  double fa = f(a), fb = f(b), fc = fb;
  double d = 0.0, d1 = 0.0;
  
  if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0)) throw std::invalid_argument("Root must be bracketed in brent");
  
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
  if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0)) throw std::invalid_argument("Root must be bracketed in bisect");
  if (std::fabs(fa) < tol) return a;
  if (std::fabs(fb) < tol) return b;
  double xmid, fmid;
  for (int j = 1; j <= maxiter; ++j) {
    xmid = a + 0.5 * (b - a);
    fmid = f(xmid);
    if (std::fabs(fmid) < tol || (b - a) < tol) return xmid;
    if ((fa > 0.0 && fmid < 0.0) || (fa < 0.0 && fmid > 0.0)) { b = xmid; fb = fmid; }
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
  double lower_val = y[j - 1];
  double upper_val = y[j];
  return (1 - g) * lower_val + g * upper_val;
}

double squantilecpp(const std::function<double(double)>& S, double p, double tol) {
  if (p < 0.0 || p > 1.0) throw std::invalid_argument("p must be in [0,1]");
  double lower = 0.0, upper = 1.0;
  double Su = S(upper);
  while (Su > p) {
    lower = upper; upper *= 2.0; Su = S(upper);
    if (upper > 1e12) throw std::runtime_error("Cannot find suitable upper bound for quantile search");
  }
  auto f = [&S, p](double t) -> double { return S(t) - p; };
  return brent(f, lower, upper, tol);
}

// --------------------------- Matrix utilities (FlatMatrix) ------------------

// Multiply matrix A (m x p) by vector x (length p), A is column-major FlatMatrix.
std::vector<double> mat_vec_mult(const FlatMatrix& A, const std::vector<double>& x) {
  int m = A.nrow;
  int p = A.ncol;
  if (static_cast<int>(x.size()) != p) throw std::invalid_argument("Vector size mismatch");
  std::vector<double> result(static_cast<size_t>(m), 0.0);
  for (int c = 0; c < p; ++c) {
    double xc = x[static_cast<size_t>(c)];
    size_t offset = FlatMatrix::idx_col(0, c, m);
    const double* colptr = A.data_ptr() ? A.data_ptr() + offset : nullptr;
    if (!colptr) continue;
    for (int r = 0; r < m; ++r) result[static_cast<size_t>(r)] += colptr[static_cast<size_t>(r)] * xc;
  }
  return result;
}

// Multiply A (m x k) * B (k x n) -> C (m x n) using column-major arithmetic and contiguous blocks.
FlatMatrix mat_mat_mult(const FlatMatrix& A, const FlatMatrix& B) {
  int m = A.nrow;
  int k = A.ncol;
  int k2 = B.nrow;
  int n = B.ncol;
  if (k != k2) throw std::invalid_argument("Matrix dimensions mismatch");
  if (m == 0 || k == 0 || n == 0) return FlatMatrix();
  FlatMatrix C(m, n);
  // Column-major: For each column j in B/C, compute C[:,j] = sum_{t=0..k-1} A[:,t] * B[t,j]
  for (int j = 0; j < n; ++j) {
    size_t coff = FlatMatrix::idx_col(0, j, m);
    for (int t = 0; t < k; ++t) {
      size_t aoff = FlatMatrix::idx_col(0, t, m);
      double scale = B.data[static_cast<size_t>(FlatMatrix::idx_col(t, j, k))];
      if (scale == 0.0) continue;
      for (int i = 0; i < m; ++i) {
        C.data[coff + static_cast<size_t>(i)] += A.data[aoff + static_cast<size_t>(i)] * scale;
      }
    }
  }
  return C;
}

// Transpose a FlatMatrix
FlatMatrix transpose(const FlatMatrix& A) {
  if (A.nrow == 0 || A.ncol == 0) return FlatMatrix();
  FlatMatrix At(A.ncol, A.nrow);
  for (int c = 0; c < A.ncol; ++c) {
    for (int r = 0; r < A.nrow; ++r) {
      At.data[FlatMatrix::idx_col(c, r, A.ncol)] = A.data[FlatMatrix::idx_col(r, c, A.nrow)];
    }
  }
  return At;
}

// --------------------------- Linear algebra helpers (FlatMatrix-backed) ----

// sumsq
double sumsq(const std::vector<double>& x) {
  double s = 0.0;
  for (double xi : x) s += xi * xi;
  return s;
}

// Householder vector (same as original, operates on std::vector<double>)
std::vector<double> house(const std::vector<double>& x) {
  int n = static_cast<int>(x.size());
  double mu = std::sqrt(sumsq(x));
  std::vector<double> v = x; // copy
  if (mu > 0.0) {
    double beta = x[0] + std::copysign(mu, x[0]);
    for (int i = 1; i < n; ++i) v[i] /= beta;
  }
  v[0] = 1.0;
  return v;
}

// Apply a row Householder (v) to submatrix A[i1..i2, j1..j2]
// A is represented as FlatMatrix (column-major)
void row_house(FlatMatrix& A, int i1, int i2, int j1, int j2, const std::vector<double>& v) {
  int m_total = A.nrow;
  if (m_total == 0) return;
  int n_total = A.ncol;
  if (i1 < 0 || i1 > i2 || i2 >= m_total) throw std::invalid_argument("Invalid row indices i1 and i2");
  if (j1 < 0 || j1 > j2 || j2 >= n_total) throw std::invalid_argument("Invalid column indices j1 and j2");
  int m = i2 - i1 + 1;
  int n = j2 - j1 + 1;
  double beta = -2.0 / sumsq(v);
  std::vector<double> w(static_cast<size_t>(n), 0.0);
  for (int jj = 0; jj < n; ++jj) {
    size_t coloff = FlatMatrix::idx_col(0, j1 + jj, m_total);
    double acc = 0.0;
    for (int ii = 0; ii < m; ++ii) acc += A.data[coloff + static_cast<size_t>(i1 + ii)] * v[static_cast<size_t>(ii)];
    w[static_cast<size_t>(jj)] = acc * beta;
  }
  for (int ii = 0; ii < m; ++ii) {
    for (int jj = 0; jj < n; ++jj) {
      size_t idx = FlatMatrix::idx_col(i1 + ii, j1 + jj, m_total);
      A.data[idx] += v[static_cast<size_t>(ii)] * w[static_cast<size_t>(jj)];
    }
  }
}

// cholesky2: in-place working on FlatMatrix (n x n), returns rank * nonneg
int cholesky2(FlatMatrix& matrix, int n, double toler) {
  double eps = 0.0;
  for (int i = 0; i < n; ++i) {
    double val = matrix.data[FlatMatrix::idx_col(i, i, n)];
    if (val > eps) eps = val;
  }
  if (eps == 0.0) eps = toler; else eps *= toler;
  int nonneg = 1;
  int rank = 0;
  
  for (int i = 0; i < n; ++i) {
    double pivot = matrix.data[FlatMatrix::idx_col(i, i, n)];
    if (std::isinf(pivot) || pivot < eps) {
      matrix.data[FlatMatrix::idx_col(i, i, n)] = 0.0;
      if (pivot < -8.0 * eps) nonneg = -1;
    } else {
      ++rank;
      for (int j = i + 1; j < n; ++j) {
        double temp = matrix.data[FlatMatrix::idx_col(i, j, n)] / pivot;
        matrix.data[FlatMatrix::idx_col(i, j, n)] = temp;
        matrix.data[FlatMatrix::idx_col(j, j, n)] -= temp * temp * pivot;
        for (int k = j + 1; k < n; ++k) {
          matrix.data[FlatMatrix::idx_col(j, k, n)] -= temp * matrix.data[FlatMatrix::idx_col(i, k, n)];
        }
      }
    }
  }
  return rank * nonneg;
}

// chsolve2 assumes matrix holds the representation produced by cholesky2
void chsolve2(FlatMatrix& matrix, int n, std::vector<double>& y) {
  // Forward substitution L * z = y
  for (int i = 0; i < n; ++i) {
    double temp = y[static_cast<size_t>(i)];
    for (int j = 0; j < i; ++j) temp -= y[static_cast<size_t>(j)] * matrix.data[FlatMatrix::idx_col(j, i, n)];
    y[static_cast<size_t>(i)] = temp;
  }
  // Backward substitution L^T * x = z
  for (int ii = n - 1; ii >= 0; --ii) {
    int i = ii;
    double diag = matrix.data[FlatMatrix::idx_col(i, i, n)];
    if (diag == 0.0) {
      y[static_cast<size_t>(i)] = 0.0;
    } else {
      double temp = y[static_cast<size_t>(i)] / diag;
      for (int j = i + 1; j < n; ++j) temp -= y[static_cast<size_t>(j)] * matrix.data[FlatMatrix::idx_col(i, j, n)];
      y[static_cast<size_t>(i)] = temp;
    }
  }
}

// chinv2: invert after decomposition in-place (FlatMatrix)
void chinv2(FlatMatrix& matrix, int n) {
  // Step 1: invert diagonal and apply sweep operator
  for (int i = 0; i < n; ++i) {
    double mii = matrix.data[FlatMatrix::idx_col(i, i, n)];
    if (mii > 0.0) {
      matrix.data[FlatMatrix::idx_col(i, i, n)] = 1.0 / mii;
      for (int j = i + 1; j < n; ++j) {
        size_t idx_ij = FlatMatrix::idx_col(i, j, n);
        matrix.data[idx_ij] = -matrix.data[idx_ij];
        for (int k = 0; k < i; ++k) {
          matrix.data[FlatMatrix::idx_col(k, j, n)] += matrix.data[idx_ij] * matrix.data[FlatMatrix::idx_col(k, i, n)];
        }
      }
    }
  }
  // Step 2: finalize inverse and symmetrize
  for (int i = 0; i < n; ++i) {
    double mii = matrix.data[FlatMatrix::idx_col(i, i, n)];
    if (mii == 0.0) {
      for (int j = 0; j < i; ++j) matrix.data[FlatMatrix::idx_col(i, j, n)] = 0.0;
      for (int j = i; j < n; ++j) matrix.data[FlatMatrix::idx_col(j, i, n)] = 0.0;
    } else {
      for (int j = i + 1; j < n; ++j) {
        double temp = matrix.data[FlatMatrix::idx_col(i, j, n)] * matrix.data[FlatMatrix::idx_col(j, j, n)];
        matrix.data[FlatMatrix::idx_col(j, i, n)] = temp;
        for (int k = i; k < j; ++k) matrix.data[FlatMatrix::idx_col(k, i, n)] += temp * matrix.data[FlatMatrix::idx_col(k, j, n)];
      }
    }
  }
}

// invsympd: invert symmetric positive definite matrix (returns FlatMatrix)
FlatMatrix invsympd(const FlatMatrix& matrix, int n, double toler) {
  FlatMatrix v = matrix; // copy
  cholesky2(v, n, toler);
  chinv2(v, n);
  // fill symmetric entries (upper -> lower)
  for (int i = 1; i < n; ++i) {
    for (int j = 0; j < i; ++j) {
      v.data[FlatMatrix::idx_col(j, i, n)] = v.data[FlatMatrix::idx_col(i, j, n)];
    }
  }
  return v;
}

// --------------------------- Survival helpers --------------------------------

DataFrameCpp survsplit(const std::vector<double>& tstart,
                       const std::vector<double>& tstop,
                       const std::vector<double>& cut) {
  int n = static_cast<int>(tstart.size());
  int ncut = static_cast<int>(cut.size());
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
  std::vector<int> censor(n2, 0);
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
  df.push_back(row, "row");
  df.push_back(start, "start");
  df.push_back(end, "end");
  df.push_back(censor, "censor");
  df.push_back(interval, "interval");
  return df;
}

// --------------------------- Misc math helpers --------------------------------

double max_elem(const std::vector<double>& x, int start, int end) {
  if (start > end || start < 0 || end >= static_cast<int>(x.size())) throw std::invalid_argument("Invalid start or end indices in max_elem");
  return *std::max_element(x.begin() + start, x.begin() + end + 1);
}

// Householder-based QR with column pivoting.
// Accepts X as FlatMatrix (m x n) and returns ListCpp containing results.
// Implementation reuses original algorithm but converts columns to a working
// row-major container internally for clarity; outputs Q/R as FlatMatrix.
ListCpp qrcpp(const FlatMatrix& X, double tol) {
  int m = X.nrow;
  int n = X.ncol;
  if (m == 0 || n == 0) {
    ListCpp empty_res;
    empty_res.push_back(0, "rank");
    return empty_res;
  }
  
  // Work on a local copy of X (we will store Householder vectors into it)
  FlatMatrix A = X;
  
  // squared column norms
  std::vector<double> cvec(static_cast<size_t>(n), 0.0);
  for (int j = 0; j < n; ++j) {
    size_t off = FlatMatrix::idx_col(0, j, m);
    double s = 0.0;
    for (int i = 0; i < m; ++i) {
      double v = A.data[off + static_cast<size_t>(i)];
      s += v * v;
    }
    cvec[static_cast<size_t>(j)] = s;
  }
  
  int r = -1;
  std::vector<int> piv(static_cast<size_t>(n));
  std::iota(piv.begin(), piv.end(), 0);
  double tau = 0.0;
  if (n > 0) tau = max_elem(cvec, 0, n - 1);
  
  while (tau > tol) {
    ++r;
    // find next pivot index kidx (first column with cvec > tol starting at r)
    int kidx = r;
    for (int kk = r; kk < n; ++kk) {
      if (cvec[static_cast<size_t>(kk)] > tol) { kidx = kk; break; }
    }
    
    // swap columns r and kidx in A, swap pivot indices and cvec
    if (kidx != r) {
      // swap column data
      size_t off_r = FlatMatrix::idx_col(0, r, m);
      size_t off_k = FlatMatrix::idx_col(0, kidx, m);
      for (int ii = 0; ii < m; ++ii) {
        std::swap(A.data[off_r + static_cast<size_t>(ii)], A.data[off_k + static_cast<size_t>(ii)]);
      }
      std::swap(piv[static_cast<size_t>(r)], piv[static_cast<size_t>(kidx)]);
      std::swap(cvec[static_cast<size_t>(r)], cvec[static_cast<size_t>(kidx)]);
    }
    
    // Householder on column r (subcolumn rows r..m-1)
    int msub = m - r;
    std::vector<double> x(static_cast<size_t>(msub));
    for (int i = 0; i < msub; ++i) {
      x[static_cast<size_t>(i)] = A.data[FlatMatrix::idx_col(r + i, r, m)];
    }
    std::vector<double> v = house(x); // v[0] == 1
    
    // apply Householder to A[r..m-1, r..n-1]
    if (msub > 0 && r < n) row_house(A, r, m - 1, r, n - 1, v);
    
    // store Householder vector in sub-diagonal of A (positions r+1..m-1, column r)
    for (int i = 1; i < msub; ++i) {
      A.data[FlatMatrix::idx_col(r + i, r, m)] = v[static_cast<size_t>(i)];
    }
    
    // update squared norms for remaining columns
    for (int j = r + 1; j < n; ++j) {
      double val = A.data[FlatMatrix::idx_col(r, j, m)];
      cvec[static_cast<size_t>(j)] -= val * val;
      if (cvec[static_cast<size_t>(j)] < 0.0) cvec[static_cast<size_t>(j)] = 0.0;
    }
    
    if (r < n - 1) tau = max_elem(cvec, r + 1, n - 1);
    else tau = 0.0;
  } // end while
  
  // rank r (last pivot index found). Following previous convention return r (and push r+1 as rank)
  // Recover Q (m x m) from stored Householder vectors in A
  FlatMatrix Qf(m, m);
  // initialize identity
  for (int c = 0; c < m; ++c) {
    size_t off = FlatMatrix::idx_col(0, c, m);
    for (int i = 0; i < m; ++i) Qf.data[off + static_cast<size_t>(i)] = 0.0;
    Qf.data[FlatMatrix::idx_col(c, c, m)] = 1.0;
  }
  
  if (r >= 0) {
    for (int kk = r; kk >= 0; --kk) {
      int msub_k = m - kk;
      std::vector<double> vks(static_cast<size_t>(msub_k));
      vks[0] = 1.0;
      for (int ii = 1; ii < msub_k; ++ii) {
        vks[static_cast<size_t>(ii)] = A.data[FlatMatrix::idx_col(kk + ii, kk, m)];
      }
      
      // apply to Qf[kk..m-1, kk..m-1]
      std::vector<double> w(static_cast<size_t>(msub_k), 0.0);
      for (int jj = kk; jj < m; ++jj) {
        double acc = 0.0;
        for (int ii = 0; ii < msub_k; ++ii) {
          acc += Qf.data[FlatMatrix::idx_col(kk + ii, jj, m)] * vks[static_cast<size_t>(ii)];
        }
        w[static_cast<size_t>(jj - kk)] = acc * (-2.0 / sumsq(vks));
      }
      for (int ii = 0; ii < msub_k; ++ii) {
        for (int jj = kk; jj < m; ++jj) {
          size_t idx = FlatMatrix::idx_col(kk + ii, jj, m);
          Qf.data[idx] += vks[static_cast<size_t>(ii)] * w[static_cast<size_t>(jj - kk)];
        }
      }
    }
  }
  
  // Recover R (m x n) using upper triangle of A
  FlatMatrix Rf(m, n);
  // initialize zeros
  std::fill(Rf.data.begin(), Rf.data.end(), 0.0);
  for (int j = 0; j < n; ++j) {
    for (int i = 0; i <= j && i < m; ++i) {
      Rf.data[FlatMatrix::idx_col(i, j, m)] = A.data[FlatMatrix::idx_col(i, j, m)];
    }
  }
  
  ListCpp result;
  result.push_back(A, "qr");         // internal transformed A (contains Householder data)
  result.push_back(r + 1, "rank");   // number of pivots found
  result.push_back(piv, "pivot");
  result.push_back(Qf, "Q");
  result.push_back(Rf, "R");
  return result;
}

// --------------------------- Matching and other helpers ----------------------

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
      result.push_back(static_cast<int>(j)); ++i; ++j;
    }
  }
  while (i < n1) { result.push_back(-1); ++i; }
  return result;
}

// --------------------------- Counterfactual helpers --------------------------

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
  for (size_t i = 0; i < n; ++i) { u_star[i] = time[i] * ((1.0 - rx[i]) + rx[i] * a); t_star[i] = u_star[i]; }
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
        if (treat[i] == 1 && all_rx1) c_star[i] = std::numeric_limits<double>::infinity();
        if (treat[i] == 0 && all_rx0) c_star[i] = std::numeric_limits<double>::infinity();
      }
    }
    for (size_t i = 0; i < n; ++i) {
      if (c_star[i] < u_star[i]) { t_star[i] = c_star[i]; d_star[i] = 0; }
    }
  }
  DataFrameCpp df;
  df.push_back(id, "uid");
  df.push_back(t_star, "t_star");
  df.push_back(d_star, "d_star");
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
    for (size_t i = 0; i < n; ++i) c_star[i] = treat[i] == 0 ? censor_time[i] * std::min(1.0, a0) : censor_time[i] * std::min(1.0, a1);
    if (autoswitch) {
      bool all_rx1 = true, all_rx0 = true;
      for (size_t i = 0; i < n; ++i) {
        if (treat[i] == 1 && rx[i] != 1.0) all_rx1 = false;
        if (treat[i] == 0 && rx[i] != 0.0) all_rx0 = false;
      }
      for (size_t i = 0; i < n; ++i) {
        if (treat[i] == 1 && all_rx1) c_star[i] = std::numeric_limits<double>::infinity();
        if (treat[i] == 0 && all_rx0) c_star[i] = std::numeric_limits<double>::infinity();
      }
    }
    for (size_t i = 0; i < n; ++i) {
      if (c_star[i] < u_star[i]) { t_star[i] = c_star[i]; d_star[i] = 0; }
    }
  }
  DataFrameCpp df;
  df.push_back(id, "uid");
  df.push_back(t_star, "t_star");
  df.push_back(d_star, "d_star");
  df.push_back(treat, "treated");
  return df;
}

// --------------------------- Misc helpers -----------------------------------

std::string sanitize(const std::string& s) {
  std::string out = s;
  for (char &c : out) {
    if (!std::isalnum(static_cast<unsigned char>(c)) && static_cast<unsigned char>(c) != '_') c = '.';
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
  if (m == 0 || static_cast<int>(lambda.size()) != m) throw std::invalid_argument("Invalid piecewise model inputs.");
  double u = logp ? std::exp(p) : p;
  if (!lowertail) u = 1.0 - u;
  if (u <= 0.0) return lowerBound;
  if (u >= 1.0) return std::numeric_limits<double>::infinity();
  double v1 = -log1p(-u);
  int j = 0;
  while (j < m && piecewiseSurvivalTime[j] <= lowerBound) ++j;
  int j1 = std::max(0, j - 1);
  double v = 0.0;
  if (j1 == m - 1) {
    double lj = lambda[j1];
    if (lj <= 0.0) return std::numeric_limits<double>::infinity();
    return lowerBound + v1 / lj;
  }
  for (j = j1; j < m - 1; ++j) {
    double dt = (j == j1) ? piecewiseSurvivalTime[j + 1] - lowerBound : piecewiseSurvivalTime[j + 1] - piecewiseSurvivalTime[j];
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
  int n = static_cast<int>(psi.size());
  if (n != static_cast<int>(Z.size())) throw std::invalid_argument("psi and Z must have the same length");
  if (n < 2) throw std::invalid_argument("Need at least two points to find roots");
  std::vector<double> Zt(n);
  for (int i = 0; i < n; ++i) Zt[static_cast<size_t>(i)] = Z[static_cast<size_t>(i)] - target;
  std::vector<double> roots;
  for (int i = 1; i < n; ++i) {
    double z1 = Zt[static_cast<size_t>(i - 1)];
    double z2 = Zt[static_cast<size_t>(i)];
    if (std::isnan(z1) || std::isnan(z2) || z1 == z2) continue;
    if (z1 == 0.0) roots.push_back(psi[static_cast<size_t>(i - 1)]);
    else if (z1 * z2 < 0.0) {
      double psi_root = psi[static_cast<size_t>(i - 1)] - z1 * (psi[static_cast<size_t>(i)] - psi[static_cast<size_t>(i - 1)]) / (z2 - z1);
      roots.push_back(psi_root);
    }
  }
  double root = NAN;
  if (!roots.empty()) {
    if (direction == -1) root = roots.front();
    else if (direction == 1) root = roots.back();
    else {
      root = roots[0];
      double minabs = std::abs(roots[0]);
      for (size_t j = 1; j < roots.size(); ++j) {
        double a = std::abs(roots[j]);
        if (a < minabs) { minabs = a; root = roots[j]; }
      }
    }
  }
  ListCpp result;
  result.push_back(roots, "all_roots");
  result.push_back(root, "selected_root");
  return result;
}

double getpsiend(const std::function<double(double)>& f,
                 bool lowerend,
                 double initialend) {
  double psiend = initialend;
  double zend = f(initialend);
  const double LIMIT = 10.0;
  if (lowerend) {
    if ((std::isinf(zend) && zend > 0) || std::isnan(zend)) {
      while (((std::isinf(zend) && zend > 0) || std::isnan(zend)) && psiend <= LIMIT) {
        psiend += 1; zend = f(psiend);
      }
      if (psiend > LIMIT) return NAN;
    }
    if (zend < 0) {
      while (!std::isinf(zend) && zend < 0 && psiend >= -LIMIT) {
        psiend -= 1; zend = f(psiend);
      }
      if (std::isinf(zend) || std::isnan(zend) || psiend < -LIMIT) return NAN;
    }
  } else {
    if ((std::isinf(zend) && zend < 0) || std::isnan(zend)) {
      while (((std::isinf(zend) && zend < 0) || std::isnan(zend)) && psiend >= -LIMIT) {
        psiend -= 1; zend = f(psiend);
      }
      if (psiend < -LIMIT) return NAN;
    }
    if (zend > 0) {
      while (!std::isinf(zend) && zend > 0 && psiend <= LIMIT) {
        psiend += 1; zend = f(psiend);
      }
      if (std::isinf(zend) || std::isnan(zend) || psiend > LIMIT) return NAN;
    }
  }
  return psiend;
}

// --------------------------- bygroup (rewritten to use DataFrameCpp/ListCpp) ----
ListCpp bygroup(const DataFrameCpp& data, const std::vector<std::string>& variables) {
  int n = static_cast<int>(data.nrows());
  int p = static_cast<int>(variables.size());
  ListCpp result;
  std::vector<int> nlevels(p);
  std::vector<std::vector<int>> indices(static_cast<size_t>(n), std::vector<int>(static_cast<size_t>(p), 0));
  ListCpp lookups_per_variable;
  
  // store unique lookups by type
  std::vector<std::vector<int>> int_lookups;
  std::vector<std::vector<double>> numeric_lookups;
  std::vector<std::vector<bool>> bool_lookups;
  std::vector<std::vector<std::string>> string_lookups;
  
  for (int i = 0; i < p; ++i) {
    const std::string& var = variables[static_cast<size_t>(i)];
    if (!data.containElementNamed(var)) throw std::invalid_argument("Data must contain variable: " + var);
    
    if (data.int_cols.count(var)) {
      auto col = data.int_cols.at(var);
      auto w = unique_sorted(col);
      nlevels[static_cast<size_t>(i)] = static_cast<int>(w.size());
      auto idx = matchcpp(col, w);
      for (int j = 0; j < n; ++j) indices[static_cast<size_t>(j)][static_cast<size_t>(i)] = idx[static_cast<size_t>(j)];
      DataFrameCpp df_uv;
      df_uv.push_back(w, var);
      lookups_per_variable.push_back(df_uv, var);
      int_lookups.push_back(w);
    } else if (data.numeric_cols.count(var)) {
      auto col = data.numeric_cols.at(var);
      auto w = unique_sorted(col);
      nlevels[static_cast<size_t>(i)] = static_cast<int>(w.size());
      auto idx = matchcpp(col, w);
      for (int j = 0; j < n; ++j) indices[static_cast<size_t>(j)][static_cast<size_t>(i)] = idx[static_cast<size_t>(j)];
      DataFrameCpp df_uv;
      df_uv.push_back(w, var);
      lookups_per_variable.push_back(df_uv, var);
      numeric_lookups.push_back(w);
    } else if (data.bool_cols.count(var)) {
      auto col = data.bool_cols.at(var);
      auto w = unique_sorted(col);
      nlevels[static_cast<size_t>(i)] = static_cast<int>(w.size());
      auto idx = matchcpp(col, w);
      for (int j = 0; j < n; ++j) indices[static_cast<size_t>(j)][static_cast<size_t>(i)] = idx[static_cast<size_t>(j)];
      DataFrameCpp df_uv;
      df_uv.push_back(w, var);
      lookups_per_variable.push_back(df_uv, var);
      bool_lookups.push_back(w);
    } else if (data.string_cols.count(var)) {
      auto col = data.string_cols.at(var);
      auto w = unique_sorted(col);
      nlevels[static_cast<size_t>(i)] = static_cast<int>(w.size());
      auto idx = matchcpp(col, w);
      for (int j = 0; j < n; ++j) indices[static_cast<size_t>(j)][static_cast<size_t>(i)] = idx[static_cast<size_t>(j)];
      DataFrameCpp df_uv;
      df_uv.push_back(w, var);
      lookups_per_variable.push_back(df_uv, var);
      string_lookups.push_back(w);
    } else {
      throw std::invalid_argument("Unsupported variable type in bygroup: " + var);
    }
  }
  
  // compute combined index
  std::vector<int> combined_index(static_cast<size_t>(n), 0);
  int orep = 1;
  for (int i = 0; i < p; ++i) orep *= nlevels[static_cast<size_t>(i)];
  for (int i = 0; i < p; ++i) {
    orep /= nlevels[static_cast<size_t>(i)];
    for (int j = 0; j < n; ++j) combined_index[static_cast<size_t>(j)] += indices[static_cast<size_t>(j)][static_cast<size_t>(i)] * orep;
  }
  
  int lookup_nrows = 1;
  for (int i = 0; i < p; ++i) lookup_nrows *= nlevels[static_cast<size_t>(i)];
  
  DataFrameCpp lookup_df;
  for (int i = 0; i < p; ++i) {
    const std::string& var = variables[static_cast<size_t>(i)];
    int nlevels_i = nlevels[static_cast<size_t>(i)];
    int repeat_each = 1;
    for (int j = i + 1; j < p; ++j) repeat_each *= nlevels[static_cast<size_t>(j)];
    int times = lookup_nrows / (nlevels_i * repeat_each);
    
    if (data.int_cols.count(var)) {
      auto& w = int_lookups.front();
      std::vector<int> col(static_cast<size_t>(lookup_nrows));
      int idx = 0;
      for (int t = 0; t < times; ++t) {
        for (int level = 0; level < nlevels_i; ++level) for (int r = 0; r < repeat_each; ++r) col[static_cast<size_t>(idx++)] = w[static_cast<size_t>(level)];
      }
      lookup_df.push_back(std::move(col), var);
      int_lookups.erase(int_lookups.begin());
    } else if (data.numeric_cols.count(var)) {
      auto& w = numeric_lookups.front();
      std::vector<double> col(static_cast<size_t>(lookup_nrows));
      int idx = 0;
      for (int t = 0; t < times; ++t) {
        for (int level = 0; level < nlevels_i; ++level) for (int r = 0; r < repeat_each; ++r) col[static_cast<size_t>(idx++)] = w[static_cast<size_t>(level)];
      }
      lookup_df.push_back(std::move(col), var);
      numeric_lookups.erase(numeric_lookups.begin());
    } else if (data.bool_cols.count(var)) {
      auto& w = bool_lookups.front();
      std::vector<bool> col(static_cast<size_t>(lookup_nrows));
      int idx = 0;
      for (int t = 0; t < times; ++t) {
        for (int level = 0; level < nlevels_i; ++level) for (int r = 0; r < repeat_each; ++r) col[static_cast<size_t>(idx++)] = w[static_cast<size_t>(level)];
      }
      lookup_df.push_back(std::move(col), var);
      bool_lookups.erase(bool_lookups.begin());
    } else if (data.string_cols.count(var)) {
      auto& w = string_lookups.front();
      std::vector<std::string> col(static_cast<size_t>(lookup_nrows));
      int idx = 0;
      for (int t = 0; t < times; ++t) {
        for (int level = 0; level < nlevels_i; ++level) for (int r = 0; r < repeat_each; ++r) col[static_cast<size_t>(idx++)] = w[static_cast<size_t>(level)];
      }
      lookup_df.push_back(std::move(col), var);
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
