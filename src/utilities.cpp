#include <Rcpp.h>
#include <R_ext/Applic.h>
#include "utilities.h"

using namespace Rcpp;


void set_seed(int seed) {
  Environment base_env("package:base");
  Function set_seed_r = base_env["set.seed"];
  set_seed_r(seed);
}


// Function to find the indices of all TRUE elements in a logical vector
IntegerVector which(const LogicalVector& vector) {
  IntegerVector true_indices;
  for (int i = 0; i < vector.size(); ++i) {
    if (vector[i]) {
      true_indices.push_back(i);
    }
  }
  return true_indices;
}


//' @title Find Interval Numbers of Indices
//' @description The implementation of \code{findInterval()} in R from
//' Advanced R by Hadley Wickham. Given a vector of non-decreasing
//' breakpoints in v, find the interval containing each element of x; i.e.,
//' if \code{i <- findInterval3(x,v)}, for each index \code{j} in \code{x},
//' \code{v[i[j]] <= x[j] < v[i[j] + 1]}, where \code{v[0] := -Inf},
//' \code{v[N+1] := +Inf}, and \code{N = length(v)}.
//'
//' @param x The numeric vector of interest.
//' @param v The vector of break points.
//' @param rightmost_closed Logical; if true, the rightmost interval, 
//'   `vec[N-1] .. vec[N]` is treated as closed if `left_open` is false, 
//'   and the leftmost interval, `vec[1] .. vec[2]` is treated as 
//'   closed if left_open is true.
//' @param all_inside Logical; if true, the returned indices are coerced
//'   into `1, ..., N-1`, i.e., `0` is mapped to `1` and 
//'   `N` is mapped to `N-1`.
//' @param left_open Logical; if true, all intervals are open at left and 
//'   closedat right. This may be useful, .e.g., in survival analysis.   
//' @return A vector of \code{length(x)} with values in \code{0:N} where
//'   \code{N = length(v)}.
//'
//' @keywords internal
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//' x <- 2:18
//' v <- c(5, 10, 15) # create two bins [5,10) and [10,15)
//' cbind(x, findInterval3(x, v))
//'
//' @export
// [[Rcpp::export]]
IntegerVector findInterval3(NumericVector x, NumericVector v, 
                            bool rightmost_closed = false, 
                            bool all_inside = false, 
                            bool left_open = false) {
  IntegerVector out(x.size());
  
  NumericVector::iterator it, pos;
  IntegerVector::iterator out_it;
  
  NumericVector::iterator x_begin=x.begin(), x_end=x.end();
  NumericVector::iterator v_begin=v.begin(), v_end=v.end();
  
  int nv = v.size();
  
  for(it = x_begin, out_it = out.begin(); it != x_end; ++it, ++out_it) {
    // Handle NA in x
    if (NumericVector::is_na(*it)) {
      *out_it = NA_INTEGER;
      continue;
    }
    
    // Choose bound depending on left_open
    if (left_open) {
      pos = std::lower_bound(v_begin, v_end, *it);
    } else {
      pos = std::upper_bound(v_begin, v_end, *it);
    }
    
    int idx = static_cast<int>(std::distance(v_begin, pos));
    
    if (rightmost_closed) {
      if (left_open) {
        if (*it == v[0]) idx = 1;
      } else {
        if (*it == v[nv-1]) idx = nv-1;
      }
    }
    
    if (all_inside) {
      if (idx == 0) {
        idx = 1;
      } else if (idx == nv) {
        idx = nv-1;
      }
    }
    
    *out_it = idx;
  }
  
  return out;
}


#include <algorithm>
#define ITMAX 100
#define EPS 3.0e-8
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
double brent(const std::function<double(double)>& f,
             double x1, double x2, double tol) {

  double a=x1, b=x2, c=x2, d, d1 = 0.0, min1, min2;
  double fa=f(a), fb=f(b), fc, p, q, r, s, tol1, xm;

  if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0)) {
    stop("Root must be bracketed in brent");
  }

  fc = fb;
  for (int iter=1; iter<=ITMAX; ++iter) {
    if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
      c = a;     // Rename a, b, c and adjust bounding interval d
      fc = fa;
      d = b - a;
      d1 = d;
    }
    if (fabs(fc) < fabs(fb)) {
      a = b;
      b = c;
      c = a;
      fa = fb;
      fb = fc;
      fc = fa;
    }
    // Convergence check
    tol1 = 2.0*EPS*fabs(b) + 0.5*tol;
    xm = 0.5*(c-b);
    if (fabs(xm) <= tol1 || fb == 0.0) {
      return b;
    }

    if (fabs(d1) >= tol1 && fabs(fa) > fabs(fb)) {
      s = fb/fa; // Attempt inverse quadratic interpolation
      if (a == c) {
        p = 2.0*xm*s;
        q = 1.0-s;
      } else {
        q = fa/fc;
        r = fb/fc;
        p = s*(2.0*xm*q*(q-r) - (b-a)*(r-1.0));
        q = (q-1.0)*(r-1.0)*(s-1.0);
      }
      if (p > 0.0) {
        q = -q;  // Check whether in bounds
      }
      p = fabs(p);
      min1 = 3.0*xm*q - fabs(tol1*q);
      min2 = fabs(d1)*fabs(q);
      if (2.0*p < (min1 < min2 ? min1 : min2)) {
        d1 = d;  // Accept interpolation
        d = p/q;
      } else {  // Interpolation failed, use bisection
        d = xm;
        d1 = d;
      }
    } else {  // Bounds decreasing too slowly, use bisection
      d = xm;
      d1 = d;
    }
    a = b;  // Move last best guess to a
    fa = fb;
    if (fabs(d) > tol1) { // Evaluate new trial root
      b += d;
    } else {
      b += SIGN(tol1, xm);
    }
    fb = f(b);
  }
  stop("Maximum number of iterations exceeded in brent");
  return 0.0; // Never get here
}


double bisect(const std::function<double(double)>& f,
              double x1, double x2, double tol) {
  double f1 = f(x1);
  double f2 = f(x2);
  
  if ((f1 > 0.0 && f2 > 0.0) || (f1 < 0.0 && f2 < 0.0)) {
    stop("Root must be bracketed in bisect");
  }
  
  // rtb will hold the endpoint where f is negative
  // dx is the signed interval width
  double rtb, dx;
  if (f1 < 0.0) {
    dx  = x2 - x1;
    rtb = x1;
  } else {
    dx  = x1 - x2;
    rtb = x2;
  }
  
  for (int j = 1; j <= ITMAX; ++j) {
    dx *= 0.5;
    double xmid = rtb + dx;
    double fmid = f(xmid);
    
    if (fmid <= 0.0) {
      rtb = xmid;
    }
    
    if (std::fabs(dx) < tol || fmid == 0.0) {
      return rtb;
    }
  }
  
  stop("Maximum number of iterations exceeded in bisect");
  return 0.0; // never reached
}


// [[Rcpp::export]]
bool hasVariable(DataFrame df, std::string varName) {
  StringVector names = df.names();
  for (int i = 0; i < names.size(); ++i) {
    if (names[i] == varName) {
      return true;
    }
  }
  return false;
}


double quantilecpp(const NumericVector& x, const double p) {
  int n = static_cast<int>(x.size());
  NumericVector y = clone(x);
  y.sort();
  double u = n*p + 1 - p;
  int j = static_cast<int>(std::floor(u));
  double g = u - j;
  double result = (1-g)*y[j-1] + g*y[j];
  return result;
}



double squantilecpp(const std::function<double(double)>& S, double p) {
  double lower = 0;
  double upper = 1;
  while (S(upper) > p) {
    lower = upper;
    upper = 2*upper;
  }
  
  auto f = [S, p](double t)->double{
    return S(t) - p;
  };
  
  return brent(f, lower, upper, 1e-6);
}


IntegerVector c_vectors_i(IntegerVector vec1, IntegerVector vec2) {
  IntegerVector result(vec1.size() + vec2.size());
  std::copy(vec1.begin(), vec1.end(), result.begin());
  std::copy(vec2.begin(), vec2.end(), result.begin() + vec1.size());
  return result;
}


NumericVector c_vectors(NumericVector vec1, NumericVector vec2) {
  NumericVector result(vec1.size() + vec2.size());
  std::copy(vec1.begin(), vec1.end(), result.begin());
  std::copy(vec2.begin(), vec2.end(), result.begin() + vec1.size());
  return result;
}


NumericMatrix subset_matrix_by_row(NumericMatrix a, IntegerVector q) {
  int n = static_cast<int>(q.size()), p = a.ncol();
  NumericMatrix b(n,p);
  for (int j=0; j<p; ++j) {
    for (int i=0; i<n; ++i) {
      b(i,j) = a(q[i],j);
    }
  }
  return b;
}


NumericMatrix c_matrices(NumericMatrix a1, NumericMatrix a2) {
  int n1 = a1.nrow(), n2 = a2.nrow(), p = a1.ncol();
  NumericMatrix b(n1+n2, p);
  for (int i=0; i<n1; ++i) {
    for (int j=0; j<p; ++j) {
      b(i,j) = a1(i,j);
    }
  }

  for (int i=0; i<n2; ++i) {
    int h = i+n1;
    for (int j=0; j<p; ++j) {
      b(h,j) = a2(i,j);
    }
  }

  return b;
}


List bygroup(DataFrame data, const StringVector& variables) {
  int n = data.nrows();
  int p = static_cast<int>(variables.size());
  
  IntegerVector d(p);   // the number of unique values
  List u(p);            // the vector of unique values
  IntegerMatrix x(n,p); // indices of original values in unique values
  for (int i=0; i<p; ++i) {
    std::string s = as<std::string>(variables[i]);
    if (!hasVariable(data, s)) {
      stop("data must contain the variables");
    }
    SEXP col = data[s];
    SEXPTYPE col_type = TYPEOF(col);
    if (col_type == LGLSXP || col_type == INTSXP) {
      IntegerVector v = col;
      IntegerVector w = unique(v);
      w.sort();
      d[i] = static_cast<int>(w.size());
      u[i] = w;
      x(_,i) = match(v,w) - 1;
    } else if (col_type == REALSXP) {
      NumericVector v = col;
      NumericVector w = unique(v);
      w.sort();
      d[i] = static_cast<int>(w.size());
      u[i] = w;
      x(_,i) = match(v,w) - 1;
    } else if (col_type == STRSXP) {
      StringVector v = col;
      StringVector w = unique(v);
      w.sort();
      d[i] = static_cast<int>(w.size());
      u[i] = w;
      x(_,i) = match(v,w) - 1;
    } else {
      stop("Unsupported variable type in bygroup" + s);
    }
  }
  
  int frac = 1;
  int orep = 1;
  for (int i=0; i<p; ++i) {
    orep *= d[i];
  }
  
  IntegerVector index(n);
  List lookup;
  for (int i=0; i<p; ++i) {
    orep /= d[i];
    index += x(_,i)*orep;
    
    IntegerVector j = rep(rep_each(seq(0, d[i]-1), orep), frac);
    std::string s = as<std::string>(variables[i]);
    SEXP data_col;
    SEXP col = u[i];
    SEXPTYPE col_type = TYPEOF(col);
    if (col_type == LGLSXP || col_type == INTSXP) {
      IntegerVector w = col;
      data_col = w[j];
    } else if (col_type == REALSXP) {
      NumericVector w = col;
      data_col = w[j];
    } else if (col_type == STRSXP) {
      StringVector w = col;
      data_col = w[j];
    } else {
      stop("Unsupported variable type in bygroup" + s);
    }
    lookup.push_back(data_col, s);
    
    frac = frac*d[i];
  }
  
  return List::create(
    Named("nlevels") = d,
    Named("indices") = x,
    Named("lookups") = u,
    Named("index") = index,
    Named("lookup") = as<DataFrame>(lookup));
}


// The following three utilities functions are from the survival package
int cholesky2(NumericMatrix matrix, int n, double toler) {
  double eps = 0;
  for (int i=0; i<n; ++i) {
    if (matrix(i,i) > eps) eps = matrix(i,i);
  }
  if (eps==0) eps = toler; // no positive diagonals!
  else eps *= toler;

  int nonneg = 1;
  int rank = 0;
  for (int i=0; i<n; ++i) {
    double pivot = matrix(i,i);
    if (std::isinf(pivot) == 1 || pivot < eps) {
      matrix(i,i) = 0;
      if (pivot < -8*eps) nonneg = -1;
    }
    else  {
      ++rank;
      for (int j=i+1; j<n; ++j) {
        double temp = matrix(i,j)/pivot;
        matrix(i,j) = temp;
        matrix(j,j) -= temp*temp*pivot;
        for (int k=j+1; k<n; ++k) matrix(j,k) -= temp*matrix(i,k);
      }
    }
  }

  return(rank*nonneg);
}


void chsolve2(NumericMatrix matrix, int n, NumericVector y) {
  for (int i=0; i<n; ++i) {
    double temp = y[i];
    for (int j=0; j<i; ++j)
      temp -= y[j]*matrix(j,i);
    y[i] = temp;
  }

  for (int i=n-1; i>=0; i--) {
    if (matrix(i,i) == 0) y[i] = 0;
    else {
      double temp = y[i]/matrix(i,i);
      for (int j=i+1; j<n; ++j)
        temp -= y[j]*matrix(i,j);
      y[i] = temp;
    }
  }
}


void chinv2(NumericMatrix matrix, int n) {
  for (int i=0; i<n; ++i){
    if (matrix(i,i) > 0) {
      matrix(i,i) = 1/matrix(i,i);   // this line inverts D
      for (int j=i+1; j<n; ++j) {
        matrix(i,j) = -matrix(i,j);
        for (int k=0; k<i; ++k)     // sweep operator
          matrix(k,j) += matrix(i,j)*matrix(k,i);
      }
    }
  }

  for (int i=0; i<n; ++i) {
    if (matrix(i,i) == 0) {  // singular row
      for (int j=0; j<i; ++j) matrix(i,j) = 0;
      for (int j=i; j<n; ++j) matrix(j,i) = 0;
    }
    else {
      for (int j=i+1; j<n; ++j) {
        double temp = matrix(i,j)*matrix(j,j);
        matrix(j,i) = temp;
        for (int k=i; k<j; ++k)
          matrix(k,i) += temp*matrix(k,j);
      }
    }
  }
}


NumericMatrix invsympd(NumericMatrix matrix, int n, double toler) {
  NumericMatrix v = clone(matrix);
  cholesky2(v, n, toler);
  chinv2(v, n);
  for (int i=1; i<n; ++i) {
    for (int j=0; j<i; ++j) {
      v(j,i) = v(i,j);
    }
  }

  return v;
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
DataFrame survsplit(NumericVector tstart,
                    NumericVector tstop,
                    NumericVector cut) {
  int extra;
  int n = static_cast<int>(tstart.size());
  int ncut = static_cast<int>(cut.size());

  // Each cut point strictly within an interval generates an extra line.
  // NA inputs are left alone.
  extra = 0;
  for (int i=0; i<n; ++i) {
    for (int j=0; j<ncut; ++j) {
      if (!std::isnan(tstart[i]) && !std::isnan(tstop[i]) &&
          cut[j] > tstart[i] && cut[j] < tstop[i]) ++extra;
    }
  }

  int n2 = n + extra;
  IntegerVector row(n2), interval(n2);
  NumericVector start(n2), end(n2);
  LogicalVector censor(n2);

  int k = 0;
  for (int i=0; i<n; ++i) {
    if (std::isnan(tstart[i]) || std::isnan(tstop[i])) {
      start[k] = tstart[i];
      end[k] = tstop[i];
      row[k] = i;           // row in the original data
      interval[k] = 1;
      ++k;
    } else {
      // find the first cut point after tstart
      int j = 0;
      for (; j < ncut && cut[j] <= tstart[i]; ++j);
      start[k] = tstart[i];
      row[k] = i;
      interval[k] = j;
      for (; j < ncut && cut[j] < tstop[i]; ++j) {
        if (cut[j] > tstart[i]) {
          end[k] = cut[j];
          censor[k] = 1;
          ++k; // create the next sub-interval
          start[k] = cut[j];
          row[k] = i;
          interval[k] = j+1;
        }
      }
      end[k] = tstop[i]; // finish the last sub-interval
      censor[k] = 0;
      ++k;
    }
  }

  DataFrame result = DataFrame::create(
    Named("row") = row,
    Named("start") = start,
    Named("end") = end,
    Named("censor") = censor,
    Named("interval") = interval);

  return result;
}


bool is_sorted(NumericVector x) {
  int n = x.size();
  
  // Loop through the vector and check if it is sorted
  for (int i = 1; i < n; ++i) {
    if (x[i] < x[i - 1]) {
      return 0;  // Return false if any element is smaller than the previous
    }
  }
  
  return 1;  // If no violations, the vector is sorted
}


// Householder vector
// Given an n-vector x, this function computes an n-vector v with v(1) = 1
// such that (I - 2*v*t(v)/t(v)*v)*x is zero in all but the first component.
NumericVector house(const NumericVector& x) {
  int n = static_cast<int>(x.size());
  double mu = sqrt(sum(x*x));
  NumericVector v = clone(x);
  if (mu > 0.0) {
    double beta = x[0] + std::copysign(1.0, x[0])*mu;
    for (int i=1; i<n; ++i) {
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
void row_house(NumericMatrix& A, const int i1, const int i2,
               const int j1, const int j2, const NumericVector& v) {
  if (i1 < 0 || i1 > i2 || i2 >= A.nrow()) {
    stop("Invalid row indices i1 and i2");
  }
  if (j1 < 0 || j1 > j2 || j2 >= A.ncol()) {
    stop("Invalid column indices j1 and j2");
  }
  
  int m = i2-i1+1, n = j2-j1+1;
  double beta = -2.0/sum(v*v);
  NumericVector w(n);
  for (int j=0; j<n; ++j) {
    for (int i=0; i<m; ++i) {
      w[j] += A(i+i1,j+j1)*v[i];
    }
    w[j] *= beta;
  }
  
  for (int i=0; i<m; ++i) {
    for (int j=0; j<n; ++j) {
      A(i+i1,j+j1) += v[i]*w[j];
    }
  }
}


//' @title QR Decomposition of a Matrix
//' @description Computes the QR decomposition of a matrix.
//'
//' @param X A numeric matrix whose QR decomposition is to be computed.
//' @param tol The tolerance for detecting linear dependencies in the
//'   columns of \code{X}.
//'
//' @details
//' This function performs Householder QR with column pivoting:
//' Given an \eqn{m}-by-\eqn{n} matrix \eqn{A} with \eqn{m \geq n},
//' the following algorithm computes \eqn{r = \textrm{rank}(A)} and
//' the factorization \eqn{Q^T A P} equal to
//' \tabular{ccccc}{
//' | \tab \eqn{R_{11}} \tab \eqn{R_{12}} \tab | \tab \eqn{r} \cr
//' | \tab 0 \tab 0 \tab | \tab \eqn{m-r} \cr
//'   \tab \eqn{r} \tab \eqn{n-r} \tab \tab
//' }
//' with \eqn{Q = H_1 \cdots H_r} and \eqn{P = P_1 \cdots P_r}.
//' The upper triangular part of \eqn{A}
//' is overwritten by the upper triangular part of \eqn{R} and
//' components \eqn{(j+1):m} of
//' the \eqn{j}th Householder vector are stored in \eqn{A((j+1):m, j)}.
//' The permutation \eqn{P} is encoded in an integer vector \code{pivot}.
//'
//' @return A list with the following components:
//'
//' * \code{qr}: A matrix with the same dimensions as \code{X}. The upper
//'   triangle contains the \code{R} of the decomposition and the lower
//'   triangle contains Householder vectors (stored in compact form).
//'
//' * \code{rank}: The rank of \code{X} as computed by the decomposition.
//'
//' * \code{pivot}: The column permutation for the pivoting strategy used
//'   during the decomposition.
//'
//' * \code{Q}: The complete \eqn{m}-by-\eqn{m} orthogonal matrix \eqn{Q}.
//'
//' * \code{R}: The complete \eqn{m}-by-\eqn{n} upper triangular
//'   matrix \eqn{R}.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @references
//' Gene N. Golub and Charles F. Van Loan.
//' Matrix Computations, second edition. Baltimore, Maryland:
//' The John Hopkins University Press, 1989, p.235.
//'
//' @examples
//'
//' hilbert <- function(n) { i <- 1:n; 1 / outer(i - 1, i, `+`) }
//' h9 <- hilbert(9)
//' qrcpp(h9)
//'
//' @export
// [[Rcpp::export]]
List qrcpp(const NumericMatrix& X, double tol = 1e-12) {
  int m = X.nrow(), n = X.ncol();
  NumericMatrix A = clone(X);
  NumericVector c(n);
  for (int j=0; j<n; ++j) {
    c[j] = sum(A(_,j)*A(_,j));
  }
  
  double tau = max(c);
  int k = 0;
  for (; k<n; ++k) {
    if (c[k] > tol) break;
  }
  
  int r = -1;
  IntegerVector piv = seq(0,n-1);
  double u;
  while (tau > tol) {
    ++r;
    
    // exchange column r with column k
    int l = piv[r];
    piv[r] = piv[k];
    piv[k] = l;
    
    for (int i=0; i<m; ++i) {
      u = A(i,r);
      A(i,r) = A(i,k);
      A(i,k) = u;
    }
    
    u = c[r];
    c[r] = c[k];
    c[k] = u;
    
    // find the Householder vector
    NumericVector v(m-r);
    for (int i=0; i<m-r; ++i) {
      v[i] = A(i+r,r);
    }
    v = house(v);
    
    // pre-multiply by the Householder matrix
    row_house(A, r, m-1, r, n-1, v);
    
    // update the sub-diagonal elements of column r
    for (int i=1; i<m-r; ++i) {
      A(i+r,r) = v[i];
    }
    
    // go to the next column and update the squared norm
    for (int i=r+1; i<n; ++i) {
      c[i] -= A(r,i)*A(r,i);
    }
    
    // identify the pivot column
    if (r < n-1) {
      tau = max(c[Range(r+1,n-1)]);
      for (k=r+1; k<n; ++k) {
        if (c[k] > tol) break;
      }
    } else {
      tau = 0.0;
    }
  }
  
  // recover the Q matrix
  NumericMatrix Q = NumericMatrix::diag(m, 1.0);
  for (int k=r; k>=0; k--) {
    NumericVector v(m-k);
    v[0] = 1.0;
    for (int i=1; i<m-k; ++i) {
      v[i] = A(i+k,k);
    }
    
    row_house(Q, k, m-1, k, m-1, v);
  }
  
  // recover the R matrix
  NumericMatrix R(m,n);
  for (int j=0; j<n; ++j) {
    for (int i=0; i<=j; ++i) {
      R(i,j) = A(i,j);
    }
  }
  
  List result = List::create(
    Named("qr") = A,
    Named("rank") = r+1,
    Named("pivot") = piv+1,
    Named("Q") = Q,
    Named("R") = R
  );
  
  return result;
}


// [[Rcpp::export]]
IntegerVector match3(
  const IntegerVector id1, const NumericVector v1,
  const IntegerVector id2, const NumericVector v2) {
  
  IntegerVector result;
  int i=0, j=0;
  while (i < id1.size() && j < id2.size()) {
    if (id1[i] < id2[j] || (id1[i] == id2[j] && v1[i] < v2[j])) {
      result.push_back(-1);
      ++i;
    } else if (id1[i] > id2[j] || (id1[i] == id2[j] && v1[i] > v2[j])) {
      ++j;
    } else {
      result.push_back(j);
      ++i;
      ++j;
    }
  }
  
  while (i < id1.size()) {
    result.push_back(-1);
    ++i;
  }
  
  return result;
}

// counterfactual untreated survival times and event indicators
DataFrame untreated(
    const double psi,
    const IntegerVector& id,
    const NumericVector& time,
    const IntegerVector& event,
    const IntegerVector& treat,
    const NumericVector& rx,
    const NumericVector& censor_time,
    const bool recensor,
    const bool autoswitch) {
  
  double a = exp(psi);
  NumericVector u_star = time*((1 - rx) + rx*a);
  NumericVector t_star = clone(u_star);
  IntegerVector d_star = clone(event);
  
  if (recensor) {
    NumericVector c_star = censor_time*std::min(1.0, a);
    
    if (autoswitch) {
      NumericVector rx1 = rx[treat == 1];
      NumericVector rx0 = rx[treat == 0];
      if (is_true(all(rx1 == 1.0))) c_star[treat == 1] = R_PosInf;
      if (is_true(all(rx0 == 0.0))) c_star[treat == 0] = R_PosInf;
    }
    
    t_star = pmin(u_star, c_star);
    d_star[c_star < u_star] = 0;
  }
  
  DataFrame result = DataFrame::create(
    Named("uid") = id,
    Named("t_star") = t_star,
    Named("d_star") = d_star,
    Named("treated") = treat
  );
  
  return result;
}


// counterfactual unswitched survival times and event indicators
DataFrame unswitched(
    const double psi,
    const int n,
    const IntegerVector& id,
    const NumericVector& time,
    const IntegerVector& event,
    const IntegerVector& treat,
    const NumericVector& rx,
    const NumericVector& censor_time,
    const bool recensor,
    const bool autoswitch) {
  
  double a = exp(psi), a1 = exp(-psi);
  NumericVector u_star(n), t_star(n);
  IntegerVector d_star(n);
  for (int i=0; i<n; ++i) {
    if (treat[i] == 0) {
      u_star[i] = time[i]*((1 - rx[i]) + rx[i]*a);
    } else {
      u_star[i] = time[i]*(rx[i] + (1 - rx[i])*a1);
    }
    t_star[i] = u_star[i];
    d_star[i] = event[i];
  }
  
  if (recensor) {
    double c0 = std::min(1.0, a), c1 = std::min(1.0, a1);
    NumericVector c_star(n);
    for (int i=0; i<n; ++i) {
      if (treat[i] == 0) {
        c_star[i] = censor_time[i]*c0;
      } else {
        c_star[i] = censor_time[i]*c1;
      }
    }
    
    if (autoswitch) {
      NumericVector rx1 = rx[treat == 1];
      NumericVector rx0 = rx[treat == 0];
      if (is_true(all(rx1 == 1.0))) c_star[treat == 1] = R_PosInf;
      if (is_true(all(rx0 == 0.0))) c_star[treat == 0] = R_PosInf;
    }
    
    t_star = pmin(u_star, c_star);
    d_star[c_star < u_star] = 0;
  }
  
  DataFrame result = DataFrame::create(
    Named("uid") = id,
    Named("t_star") = t_star,
    Named("d_star") = d_star,
    Named("treated") = treat
  );
  
  return result;
}


std::string sanitize(const std::string& s) {
  std::string out = s;
  for (char &c : out) {
    if (!std::isalnum(static_cast<unsigned char>(c)) && c != '_') {
      c = '.';
    }
  }
  return out;
}

// [[Rcpp::export]]
double qtpwexpcpp1(const double p,
                   const NumericVector& piecewiseSurvivalTime,
                   const NumericVector& lambda,
                   const double lowerBound,
                   const bool lowertail,
                   const bool logp) {
  int m = static_cast<int>(piecewiseSurvivalTime.size());
  
  // cumulative hazard from lowerBound until the quantile
  double u = p;
  if (logp) u = exp(p);
  if (!lowertail) u = 1.0 - u;
  
  // identify the time interval containing the lowerBound
  int j = 0;
  for (; j<m; ++j) {
    if (piecewiseSurvivalTime[j] > lowerBound) break;
  }
  int j1 = (j==0 ? 0 : j-1); // to handle floating point precision

  double q;
  double v1 = -log(1.0 - u);
  if (j1 == m-1) { // in the last interval
    q = (lambda[j1]==0.0 ? 1.0e+8 : v1/lambda[j1] + lowerBound);
  } else {
    // accumulate the pieces on the cumulative hazard scale
    double v = 0;
    for (j=j1; j<m-1; ++j) {
      if (j==j1) {
        v += lambda[j]*(piecewiseSurvivalTime[j+1] - lowerBound);
      } else {
        v += lambda[j]*(piecewiseSurvivalTime[j+1] -
          piecewiseSurvivalTime[j]);
      }
      if (v >= v1) break;
    }
    
    if (j == m-1) { // in the last interval
      q = (lambda[j]==0.0 ? 1.0e+8 :
             (v1 - v)/lambda[j] + piecewiseSurvivalTime[j]);
    } else {
      q = (lambda[j]==0.0 ? 1.0e+8 :
             piecewiseSurvivalTime[j+1] - (v - v1)/lambda[j]);
    }
  }
  
  return q;
}


List getpsiest(const double target, const NumericVector& psi, 
                 const NumericVector& Z, const int direction = 0) {
  
  int n = psi.size();
  if (n != Z.size()) stop("psi and Z must have the same length");
  if (n < 2) stop("Need at least two points to find roots");
  NumericVector Zt = Z - target; // shifted Z values
  
  std::vector<double> roots;
  for (int i = 1; i < n; ++i) {
    double z1 = Zt[i-1];
    double z2 = Zt[i];
    
    // skip missing values and identical values
    if (std::isnan(z1) || std::isnan(z2) || z1 == z2)
      continue;
    
    // exact zero
    if (z1 == 0.0)
      roots.push_back(psi[i - 1]);
    else if (z1 * z2 < 0.0) {
      // linear interpolation for zero crossing
      double psi_root = psi[i - 1] - z1 * (psi[i] - psi[i - 1]) / (z2 - z1);
      roots.push_back(psi_root);
    }
  }

  NumericVector out_roots(roots.begin(), roots.end());
  
  double root = NA_REAL;
  if (!roots.empty()) {
    if (direction == -1) { // leftmost
      root = out_roots[0]; 
    } else if (direction == 1) { // rightmost
      root = out_roots[out_roots.size() - 1];
    } else { // closest to zero
      root = NA_REAL;
      double minabs = std::abs(roots[0]);
      root = roots[0];
      for (size_t j = 1; j < roots.size(); ++j) {
        double a = std::abs(roots[j]);
        if (a < minabs) {
          minabs = a;
          root = roots[j];
        }
      }
    }
  }
  
  return List::create(
    _["roots"] = out_roots,
    _["root"] = root
  );
}


double getpsiend(const std::function<double(double)>& f,
                 const bool lowerend, const double initialend) {
  double psiend = initialend, zend = f(initialend);
  if (lowerend) {
    if ((std::isinf(zend) && zend > 0) || std::isnan(zend)) {
      while (((std::isinf(zend) && zend > 0) || std::isnan(zend)) && 
             psiend <= 10) {
        psiend = psiend + 1; zend = f(psiend);
      }
      if (psiend > 10) {
        psiend = NA_REAL; zend = NA_REAL;
      }
    }
    
    if (zend < 0) {
      while (!std::isinf(zend) && zend < 0 && psiend >= -10) {
        psiend = psiend - 1; zend = f(psiend);
      }
      if (std::isinf(zend) || std::isnan(zend) || psiend < -10) {
        psiend = NA_REAL;
      }
    }
  } else { // upper end
    if ((std::isinf(zend) && zend < 0) || std::isnan(zend)) {
      while (((std::isinf(zend) && zend < 0) || std::isnan(zend)) && 
             psiend >= -10) {
        psiend = psiend - 1; zend = f(psiend);
      }
      if (psiend < -10) {
        psiend = NA_REAL; zend = NA_REAL;
      }
    }
    
    if (zend > 0) {
      while (!std::isinf(zend) && zend > 0 && psiend <= 10) {
        psiend = psiend + 1; zend = f(psiend);
      }
      if (std::isinf(zend) || std::isnan(zend) || psiend > 10) {
        psiend = NA_REAL;
      }
    }
  }
  
  return psiend;
}
