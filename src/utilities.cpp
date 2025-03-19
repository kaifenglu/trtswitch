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
  for (int i = 0; i < vector.size(); i++) {
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
IntegerVector findInterval3(NumericVector x, NumericVector v) {
  IntegerVector out(x.size());

  NumericVector::iterator it, pos;
  IntegerVector::iterator out_it;

  NumericVector::iterator x_begin=x.begin(), x_end=x.end();
  NumericVector::iterator v_begin=v.begin(), v_end=v.end();

  for(it = x_begin, out_it = out.begin(); it != x_end; ++it, ++out_it) {
    pos = std::upper_bound(v_begin, v_end, *it);
    *out_it = static_cast<int>(std::distance(v_begin, pos));
  }

  return out;
}


#include <algorithm>
#define ITMAX 100
#define EPS 3.0e-8
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

//' @title Brent's Method for Root-Finding
//' @description Using Brent's method, find the root of a function known to
//' lie between x1 and x2. Program based on the book - Numerical Recipes in C
//' The Art of Scientific Computing - Second Edition, by William H. Press,
//' Saul A. Teukolsky, William T. Vetterling, and Brian P. Flannery.
//' It mimics the uniroot() function in R.
//'
//' @param f Name of the univariate objective function.
//' @param x1 One end of the interval bracket.
//' @param x2 The other end of the interval bracket.
//' @param tol The tolerance limit for stopping the iteration.
//'
//' @return The root x between x1 and x2 such that f(x) = 0.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//' brent(sin, -1, 1, 0.0001)
//' @export
//'
// [[Rcpp::plugins(cpp11)]]
double brent(const std::function<double(double)>& f,
             double x1, double x2, double tol) {
  int iter;
  double a=x1, b=x2, c=x2, d, d1 = 0.0, min1, min2;
  double fa=f(a), fb=f(b), fc, p, q, r, s, tol1, xm;

  if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0)) {
    stop("Root must be bracketed in brent");
  }

  fc = fb;
  for (iter=1; iter<=ITMAX; iter++) {
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


// [[Rcpp::export]]
bool hasVariable(DataFrame df, std::string varName) {
  StringVector names = df.names();
  for (int i = 0; i < names.size(); i++) {
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


// [[Rcpp::plugins(cpp11)]]
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
  int i, j, n = static_cast<int>(q.size()), p = a.ncol();
  NumericMatrix b(n,p);
  for (j=0; j<p; j++) {
    for (i=0; i<n; i++) {
      b(i,j) = a(q[i],j);
    }
  }
  return b;
}


NumericMatrix c_matrices(NumericMatrix a1, NumericMatrix a2) {
  int h, i, j, n1 = a1.nrow(), n2 = a2.nrow(), p = a1.ncol();
  NumericMatrix b(n1+n2, p);
  for (i=0; i<n1; i++) {
    for (j=0; j<p; j++) {
      b(i,j) = a1(i,j);
    }
  }

  for (i=0; i<n2; i++) {
    h = i+n1;
    for (j=0; j<p; j++) {
      b(h,j) = a2(i,j);
    }
  }

  return b;
}


List bygroup(DataFrame data, const StringVector& variables) {
  int i;
  int n = data.nrows();
  int p = static_cast<int>(variables.size());

  IntegerVector d(p);   // the number of unique values
  List u(p);            // the vector of unique values
  IntegerMatrix x(n,p); // indices of original values in unique values
  for (i=0; i<p; i++) {
    String s = variables[i];
    if (!hasVariable(data, s)) {
      stop("data must contain the variables");
    }

    if (TYPEOF(data[s]) == LGLSXP || TYPEOF(data[s]) == INTSXP) {
      IntegerVector v = data[s];
      IntegerVector w = unique(v);
      w.sort();
      d[i] = static_cast<int>(w.size());
      u[i] = w;
      x(_,i) = match(v,w) - 1;
    } if (TYPEOF(data[s]) == REALSXP) {
      NumericVector v = data[s];
      NumericVector w = unique(v);
      w.sort();
      d[i] = static_cast<int>(w.size());
      u[i] = w;
      x(_,i) = match(v,w) - 1;
    } if (TYPEOF(data[s]) == STRSXP) {
      StringVector v = data[s];
      StringVector w = unique(v);
      w.sort();
      d[i] = static_cast<int>(w.size());
      u[i] = w;
      x(_,i) = match(v,w) - 1;
    }
  }

  int frac = 1;
  int orep = 1;
  for (i=0; i<p; i++) {
    orep = orep*d[i];
  }

  IntegerVector index(n);
  DataFrame lookup;
  for (i=0; i<p; i++) {
    orep = orep/d[i];
    index = index + x(_,i)*orep;

    IntegerVector j = rep(rep_each(seq(0, d[i]-1), orep), frac);
    String s = variables[i];
    if (TYPEOF(data[s]) == LGLSXP || TYPEOF(data[s]) == INTSXP) {
      IntegerVector w = u[i];
      lookup.push_back(w[j],s);
    } else if (TYPEOF(data[s]) == REALSXP) {
      NumericVector w = u[i];
      lookup.push_back(w[j],s);
    } else if (TYPEOF(data[s]) == STRSXP) {
      StringVector w = u[i];
      lookup.push_back(w[j],s);
    }

    frac = frac*d[i];
  }

  return List::create(
    Named("nlevels") = d,
    Named("indices") = x+1,
    Named("lookups") = u,
    Named("index") = index+1,
    Named("lookup") = lookup);
}


// The following three utilities functions are from the survival package
int cholesky2(NumericMatrix matrix, int n, double toler) {
  double temp;
  int i, j, k;
  double eps, pivot;
  int rank;
  int nonneg;

  nonneg = 1;
  eps = 0;
  for (i=0; i<n; i++) {
    if (matrix(i,i) > eps) eps = matrix(i,i);
  }
  if (eps==0) eps = toler; // no positive diagonals!
  else eps *= toler;

  rank = 0;
  for (i=0; i<n; i++) {
    pivot = matrix(i,i);
    if (std::isinf(pivot) == 1 || pivot < eps) {
      matrix(i,i) = 0;
      if (pivot < -8*eps) nonneg = -1;
    }
    else  {
      rank++;
      for (j=i+1; j<n; j++) {
        temp = matrix(i,j)/pivot;
        matrix(i,j) = temp;
        matrix(j,j) -= temp*temp*pivot;
        for (k=j+1; k<n; k++) matrix(j,k) -= temp*matrix(i,k);
      }
    }
  }

  return(rank*nonneg);
}


void chsolve2(NumericMatrix matrix, int n, NumericVector y) {
  int i, j;
  double temp;

  for (i=0; i<n; i++) {
    temp = y[i];
    for (j=0; j<i; j++)
      temp -= y[j]*matrix(j,i);
    y[i] = temp;
  }

  for (i=n-1; i>=0; i--) {
    if (matrix(i,i) == 0) y[i] = 0;
    else {
      temp = y[i]/matrix(i,i);
      for (j=i+1; j<n; j++)
        temp -= y[j]*matrix(i,j);
      y[i] = temp;
    }
  }
}


void chinv2(NumericMatrix matrix, int n) {
  double temp;
  int i, j, k;

  for (i=0; i<n; i++){
    if (matrix(i,i) > 0) {
      matrix(i,i) = 1/matrix(i,i);   // this line inverts D
      for (j=i+1; j<n; j++) {
        matrix(i,j) = -matrix(i,j);
        for (k=0; k<i; k++)     // sweep operator
          matrix(k,j) += matrix(i,j)*matrix(k,i);
      }
    }
  }

  for (i=0; i<n; i++) {
    if (matrix(i,i) == 0) {  // singular row
      for (j=0; j<i; j++) matrix(i,j) = 0;
      for (j=i; j<n; j++) matrix(j,i) = 0;
    }
    else {
      for (j=i+1; j<n; j++) {
        temp = matrix(i,j)*matrix(j,j);
        matrix(j,i) = temp;
        for (k=i; k<j; k++)
          matrix(k,i) += temp*matrix(k,j);
      }
    }
  }
}


NumericMatrix invsympd(NumericMatrix matrix, int n, double toler) {
  int i, j;
  NumericMatrix v = clone(matrix);
  i = cholesky2(v, n, toler);
  chinv2(v, n);
  for (i=1; i<n; i++) {
    for (j=0; j<i; j++) {
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
//' * \code{stop}: The stopping time of the resulting subrecord.
//'
//' * \code{censor}: Whether the subrecord lies strictly within a record
//'   in the input data.
//'
//' * \code{interval}: The interval number.
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
  int i, j, k, extra;
  int n = static_cast<int>(tstart.size());
  int ncut = static_cast<int>(cut.size());

  // Each cut point strictly within an interval generates an extra line.
  // NA inputs are left alone.
  extra = 0;
  for (i=0; i<n; i++) {
    for (j=0; j<ncut; j++) {
      if (!std::isnan(tstart[i]) && !std::isnan(tstop[i]) &&
          cut[j] > tstart[i] && cut[j] < tstop[i]) extra++;
    }
  }

  int n2 = n + extra;
  IntegerVector row(n2), interval(n2);
  NumericVector start(n2), end(n2);
  LogicalVector censor(n2);

  k = 0;
  for (i=0; i<n; i++) {
    if (std::isnan(tstart[i]) || std::isnan(tstop[i])) {
      start[k] = tstart[i];
      end[k] = tstop[i];
      row[k] = i;           // row in the original data
      interval[k] = 1;
      k++;
    } else {
      // find the first cut point after tstart
      for (j=0; j < ncut && cut[j] <= tstart[i]; j++);
      start[k] = tstart[i];
      row[k] = i;
      interval[k] = j;
      for (; j < ncut && cut[j] < tstop[i]; j++) {
        if (cut[j] > tstart[i]) {
          end[k] = cut[j];
          censor[k] = 1;
          k++; // create the next sub-interval
          start[k] = cut[j];
          row[k] = i;
          interval[k] = j+1;
        }
      }
      end[k] = tstop[i]; // finish the last sub-interval
      censor[k] = 0;
      k++;
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
    for (int i=1; i<n; i++) {
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
  
  int i, j, m = i2-i1+1, n = j2-j1+1;
  double beta = -2.0/sum(v*v);
  NumericVector w(n);
  for (j=0; j<n; j++) {
    for (i=0; i<m; i++) {
      w[j] += A(i+i1,j+j1)*v[i];
    }
    w[j] *= beta;
  }
  
  for (i=0; i<m; i++) {
    for (j=0; j<n; j++) {
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
  int i, j, k, l, m = X.nrow(), n = X.ncol();
  NumericMatrix A = clone(X);
  NumericVector c(n);
  for (j=0; j<n; j++) {
    c[j] = sum(A(_,j)*A(_,j));
  }
  
  double tau = max(c);
  for (k=0; k<n; k++) {
    if (c[k] > tol) break;
  }
  
  int r = -1;
  IntegerVector piv = seq(0,n-1);
  double u;
  while (tau > tol) {
    r++;
    
    // exchange column r with column k
    l = piv[r];
    piv[r] = piv[k];
    piv[k] = l;
    
    for (i=0; i<m; i++) {
      u = A(i,r);
      A(i,r) = A(i,k);
      A(i,k) = u;
    }
    
    u = c[r];
    c[r] = c[k];
    c[k] = u;
    
    // find the Householder vector
    NumericVector v(m-r);
    for (i=0; i<m-r; i++) {
      v[i] = A(i+r,r);
    }
    v = house(v);
    
    // pre-multiply by the Householder matrix
    row_house(A, r, m-1, r, n-1, v);
    
    // update the sub-diagonal elements of column r
    for (i=1; i<m-r; i++) {
      A(i+r,r) = v[i];
    }
    
    // go to the next column and update the squared norm
    for (i=r+1; i<n; i++) {
      c[i] -= A(r,i)*A(r,i);
    }
    
    // identify the pivot column
    if (r < n-1) {
      tau = max(c[Range(r+1,n-1)]);
      for (k=r+1; k<n; k++) {
        if (c[k] > tol) break;
      }
    } else {
      tau = 0.0;
    }
  }
  
  // recover the Q matrix
  NumericMatrix Q = NumericMatrix::diag(m, 1.0);
  for (k=r; k>=0; k--) {
    NumericVector v(m-k);
    v[0] = 1.0;
    for (i=1; i<m-k; i++) {
      v[i] = A(i+k,k);
    }
    
    row_house(Q, k, m-1, k, m-1, v);
  }
  
  // recover the R matrix
  NumericMatrix R(m,n);
  for (j=0; j<n; j++) {
    for (i=0; i<=j; i++) {
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
  
  int i;
  double a = exp(psi), a1 = exp(-psi);
  NumericVector u_star(n), t_star(n);
  IntegerVector d_star(n);
  for (i=0; i<n; i++) {
    if (treat[i] == 0) {
      u_star[i] = time[i]*((1 - rx[i]) + rx[i]*a);
    } else {
      u_star[i] = time[i]*(rx[i] + (1 - rx[i])*a1);
    }
    t_star[i] = u_star[i];
    d_star[i] = event[i];
  }
  
  if (recensor) {
    NumericVector c_star(n);
    for (i=0; i<n; i++) {
      if (treat[i] == 0) {
        c_star[i] = censor_time[i]*std::min(1.0, a);
      } else {
        c_star[i] = censor_time[i]*std::min(1.0, a1);
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
