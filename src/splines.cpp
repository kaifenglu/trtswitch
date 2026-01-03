#include <RcppThread.h>

#include "splines.h"
#include "utilities.h"
#include "dataframe_list.h"
#include "thread_utils.h"

#include <algorithm> // any_of, is_sorted, max_element, min_element, none_of, sort
#include <cmath>     // isfinite, isnan, pow
#include <stdexcept> // invalid_argument
#include <limits>    // numeric_limits
#include <vector>    // vector
#include <string>    // string, to_string
#include <utility>   // move, swap
#include <cstddef>   // size_t

// Helper: indices of non-NaN entries
static std::vector<int> indices_not_nan(const std::vector<double>& v) {
  std::vector<int> idx;
  idx.reserve(v.size());
  int n = v.size();
  for (int i = 0; i < n; ++i) {
    if (!std::isnan(v[i])) idx.push_back(i);
  }
  return idx;
}

// Create FlatMatrix filled with NaN (column-major)
static FlatMatrix make_nan_matrix(int rows, int cols) {
  if (rows <= 0 || cols <= 0) return FlatMatrix();
  FlatMatrix fm(rows, cols);
  fm.fill(std::numeric_limits<double>::quiet_NaN());
  return fm;
}

// Helper: findSpan (The NURBS Book)
int findSpan(int n, int p, double u, const std::vector<double>& U) {
  if (u == U[n+1]) return n;
  int low = p;
  int high = n + 1;
  int mid = (low + high) / 2;
  while (u < U[mid] || u >= U[mid+1]) {
    if (u < U[mid]) high = mid;
    else low = mid;
    mid = (low + high) / 2;
  }
  return mid;
}

// Helper: dersBasisFuns (algorithm A2.3) computes derivatives up to order n
// ders(k,j) where k=0..n, j=0..p (p = degree)
FlatMatrix dersBasisFuns(int i, double u, int p, int n, const std::vector<double>& U) {
  FlatMatrix ders(n+1, p+1);

  // Allocate small temporary arrays
  FlatMatrix ndu(p+1, p+1);
  std::vector<double> left(p+1), right(p+1);
  
  ndu(0,0) = 1.0;
  for (int j = 1; j <= p; ++j) {
    left[j] = u - U[i + 1 - j];
    right[j] = U[i + j] - u;
    double saved = 0.0;
    for (int r = 0; r < j; ++r) {
      ndu(j,r) = right[r+1] + left[j-r];
      double temp = ndu(r,j-1) / ndu(j,r);
      ndu(r,j) = saved + right[r+1] * temp;
      saved = left[j-r] * temp;
    }
    ndu(j,j) = saved;
  }
  // load basis functions
  for (int j = 0; j <= p; ++j) ders(0,j) = ndu(j,p);
  
  // compute derivatives
  FlatMatrix a(2, p+1);
  for (int r = 0; r <= p; ++r) {
    int s1 = 0, s2 = 1;
    a(0,0) = 1.0;
    // compute the kth derivative
    for (int k = 1; k <= n; ++k) {
      double d = 0.0;
      int rk = r - k;
      int pk = p - k;
      if (r >= k) {
        a(s2,0) = a(s1,0) / ndu(pk+1,rk);
        d = a(s2,0) * ndu(rk,pk);
      }
      
      int j1 = std::max(1, -rk);
      int j2 = std::min(k-1, p-r);
      for (int j = j1; j <= j2; ++j) {
        a(s2,j) = (a(s1,j) - a(s1,j-1)) / ndu(pk+1,rk + j);
        d += a(s2,j) * ndu(rk + j, pk);
      }
      
      if (r <= pk) {
        a(s2,k) = -a(s1,k-1) / ndu(pk+1,r);
        d += a(s2,k) * ndu(r,pk);
      }
      ders(k,r) = d;
      std::swap(s1, s2);
    }
  }
  
  // multiply through by the correct factors
  int r = p;
  for (int k = 1; k <= n; ++k) {
    for (int j = 0; j <= p; ++j) ders(k,j) *= r;
    r *= (p - k);
  }
  
  return ders;
};


// splineDesigncpp using findSpan + basisFuns + dersBasisFuns (The NURBS Book)
// - ord = p+1 where p is degree
// - computes k-th derivatives up to requested order efficiently
FlatMatrix splineDesigncpp(
    const std::vector<double>& knots_in,
    const std::vector<double>& x,
    int ord,
    const std::vector<int>& derivs)
{
  std::vector<double> knots = knots_in;
  int nk = knots.size();
  if (nk <= 0) throw std::invalid_argument("must have at least 'ord' knots");
  if (!std::is_sorted(knots.begin(), knots.end())) 
    std::sort(knots.begin(), knots.end());
  
  int nx = x.size();
  int nd = derivs.size();
  if (nd > nx) throw std::invalid_argument(
      "length of 'derivs' is larger than length of 'x'");
  if (nd < 1) throw std::invalid_argument("empty 'derivs'");
  
  if (ord < 1 || ord > nk) 
    throw std::invalid_argument(
        "'ord' must be a positive integer, at most the number of knots");
  if (nk < 2*ord - 1) throw std::invalid_argument("need at least 2*ord - 1 knots");
  
  int p = ord - 1;              // degree
  int ncol = nk - ord;          // number of basis functions
  if (ncol <= 0) ncol = 0;
  
  // Check x within allowable outer range (same as before)
  double xmin_allowed = knots[p];
  double xmax_allowed = knots[nk - p - 1];
  for (int i = 0; i < nx; ++i) {
    if (std::isnan(x[i])) continue;
    if (x[i] < xmin_allowed || x[i] > xmax_allowed) {
      thread_utils::push_thread_warning("the 'x' values must be in the range "
                                          "of knots[ord], knots[nk - degree]");
      throw std::out_of_range("the 'x' values must be in the range "
                                "of knots[ord], knots[nk - degree]");
    }
  }
  
  // Prepare derivs recycled
  std::vector<int> derivs2(nx);
  for (int i = 0; i < nx; ++i) {
    derivs2[i] = derivs[i % nd];
  }

  // Precompute n parameter for findSpan: n = ncol - 1
  int n = ncol - 1;
  
  // Prepare output design matrix as column-major FlatMatrix (nx rows x ncol columns)
  FlatMatrix design(nx, ncol);

  // Temporary container for ders
  FlatMatrix ders; // will be (deriv+1) x (p+1) per evaluation
  
  // Main loop over x values
  for (int i = 0; i < nx; ++i) {
    double xv = x[i];
    if (std::isnan(xv)) continue;
    
    int deriv_order = derivs2[i];
    if (deriv_order > p) deriv_order = p;
    
    // find span
    int span = findSpan(n, p, xv, knots);
    // compute ders up to deriv_order
    ders = dersBasisFuns(span, xv, p, deriv_order, knots);
    
    // ders(k,j) corresponds to derivative order k of basis function N_{span-p+j,p}(x)
    int start = span - p; // start index of nonzero basis functions
    for (int j = 0; j <= p; ++j) {
      int global_i = start + j; // global basis function index (column)
      if (global_i >= 0 && global_i < ncol) {
        design(i, global_i) = ders(deriv_order, j);
      }
    }
  }
  
  return design;
}


// ---------------------------------------------------------------------------
// wrappers that assemble a ListCpp result with the requested metadata
// ---------------------------------------------------------------------------

// bscpp now returns a ListCpp with the requested fields
ListCpp bscpp(
    const std::vector<double>& x,
    int df, // number of columns of the basis matrix
    const std::vector<double>& knots_in,
    int degree,
    bool intercept,
    const std::vector<double>& boundary_knots_in,
    bool warn_outside)
{
  // compute basis matrix 
  if (x.empty() || (x.size() == 1 && std::isnan(x[0]))) {
    throw std::invalid_argument("input 'x' is empty or NA");
  }
  
  int m = x.size();
  int ord = 1 + degree;
  if (ord <= 1) throw std::invalid_argument("'degree' must be a positive integer");
  
  // NA handling
  std::vector<int> non_na_idx = indices_not_nan(x);
  std::vector<double> z = subset(x, non_na_idx);
  int n = z.size();
  
  // boundary knots
  std::vector<double> boundary_knots = boundary_knots_in;
  
  if (boundary_knots.empty() || 
      (boundary_knots.size() == 1 && std::isnan(boundary_knots[0]))) {
    if (n == 1) {
      boundary_knots = { z[0]*7.0/8.0, z[0]*9.0/8.0 };
    } else {
      double mn = *std::min_element(z.begin(), z.end());
      double mx = *std::max_element(z.begin(), z.end());
      boundary_knots = {mn, mx};
    }
  } else {
    if (boundary_knots.size() != 2) 
      throw std::invalid_argument("boundary_knots must have length 2");
  }
  
  
  // compute K and knots
  std::vector<double> knots = knots_in;
  bool mk_knots = (df > 0) && (knots.empty() || 
                   (knots.size() == 1 && std::isnan(knots[0])));
  int K = 0; // number of inner knots
  if (mk_knots) {
    K = df - ord + (1 - (intercept ? 1 : 0));
    if (K < 0) {
      K = 0;
      thread_utils::push_thread_warning("'df' is too small");
    }
    std::vector<double> knots1(K);
    std::vector<double> z_no_outside;
    for (double zi : z) 
      if (!(zi < boundary_knots[0]) && !(zi > boundary_knots[1])) 
        z_no_outside.push_back(zi);
    for (int k = 0; k < K; ++k) 
      knots1[k] = quantilecpp(z_no_outside, (k+1.0)/(K+1.0));
    knots = std::move(knots1);
  } else {
    for (double kv : knots) if (!std::isfinite(kv)) 
      throw std::invalid_argument("non-finite knots");
    K = knots.size();
  }
  
  // If mk_knots && K>0, handle interior knots equalling boundary knots- shove inside
  if (mk_knots && K > 0) {
    bool lrEq0 = std::any_of(knots.begin(), knots.end(), 
                             [&](double v){ return v == boundary_knots[0]; });
    bool lrEq1 = std::any_of(knots.begin(), knots.end(), 
                             [&](double v){ return v == boundary_knots[1]; });
    bool anyWarning = false;
    if (lrEq0) {
      double piv = boundary_knots[0];
      std::vector<char> sub(K);
      bool aE0 = true;
      for (int i = 0; i < K; ++i) { 
        sub[i] = (knots[i] == piv); 
        if (!sub[i]) aE0 = false; 
      }
      if (aE0) {
        thread_utils::push_thread_warning(
          "all interior knots match left boundary knot");
      } else {
        std::vector<double> knots2;
        for (double kv : knots) if (kv > piv) knots2.push_back(kv);
        double shift = (*std::min_element(knots2.begin(), knots2.end()) - piv) / 8.0;
        for (int i = 0; i < K; ++i) if (sub[i]) knots[i] = knots[i] + shift;
      }
      anyWarning = true;
    }
    if (lrEq1) {
      double piv = boundary_knots[1];
      std::vector<char> sub(K);
      bool aE1 = true;
      for (int i = 0; i < K; ++i) { 
        sub[i] = (knots[i] == piv); 
        if (!sub[i]) aE1 = false; 
      }
      if (aE1) {
        thread_utils::push_thread_warning(
          "all interior knots match right boundary knot");
      } else {
        std::vector<double> knots2;
        for (double kv : knots) if (kv < piv) knots2.push_back(kv);
        double shift = (piv - *std::max_element(knots2.begin(), knots2.end())) / 8.0;
        for (int i = 0; i < K; ++i) if (sub[i]) knots[i] = knots[i] - shift;
      }
      anyWarning = true;
    }
    if (anyWarning && warn_outside) 
      thread_utils::push_thread_warning(
        "shoving 'interior' knots matching boundary knots to inside");
  }
  
  // Build full knot vector u with length K + 2*ord
  int nU = K + 2*ord;
  std::vector<double> u(nU);
  for (int k = 0; k < ord; ++k) u[k] = boundary_knots[0];
  for (int k = ord; k < K + ord; ++k) u[k] = knots[k-ord];
  for (int k = K + ord; k < K + 2 * ord; ++k) u[k] = boundary_knots[1];
  
  // Prepare deriv0 vector
  std::vector<int> deriv0(1, 0);
  
  // Design matrix (for non-NA points) as FlatMatrix: n rows x (K+ord) cols
  int Lfull = K + ord;
  FlatMatrix design(n, Lfull);

  // Identify outside points
  std::vector<unsigned char> outleft(n, 0), outright(n, 0), outside(n, 0);
  for (int i = 0; i < n; ++i) {
    outleft[i] = (z[i] < boundary_knots[0]);
    outright[i] = (z[i] > boundary_knots[1]);
    outside[i] = (outleft[i] || outright[i]);
  }
  
  if (std::any_of(outside.begin(), outside.end(), [](char v){ return v; })) {
    if (warn_outside) {
      thread_utils::push_thread_warning("some 'x' values beyond boundary knots "
                                          "may cause ill-conditioned bases");
    }
    
    // construct factorial scaling vector scalef
    std::vector<double> scalef(ord);
    scalef[0] = 1.0;
    for (int i = 1; i < ord; ++i) scalef[i] = scalef[i-1] * i;
    
    double e = 0.25;
    
    // Left extrapolation
    if (std::any_of(outleft.begin(), outleft.end(), [](char v){ return v; })) {
      double k_pivot = (1 - e) * boundary_knots[0] + e * u[ord];
      std::vector<int> idx = which(outleft);
      std::vector<double> zol = subset(z, idx);
      int r = zol.size();
      
      FlatMatrix zl(r, ord);
      double* zptr = zl.data_ptr();
      for (int i = 0; i < r; ++i) zptr[i] = 1.0;
      for (int j = 1; j <= degree; ++j) {
        double* zcol = zptr + j * r;
        for (int i = 0; i < r; ++i) zcol[i] = std::pow(zol[i] - k_pivot, j);
      }
      
      std::vector<double> kp(ord, k_pivot);
      std::vector<int> derivs(ord);
      for (int i = 0; i < ord; ++i) derivs[i] = i;
      FlatMatrix tt = splineDesigncpp(u, kp, ord, derivs);

      for (int j = 0; j < Lfull; ++j) {
        double* design_col = design.data_ptr() + j * n;
        for (int l = 0; l < ord; ++l) {
          const double* zl_col = zl.data_ptr() + l * r;
          double tt_val = tt(l, j) / scalef[l];
          for (int i = 0; i < r; ++i) {
            design_col[idx[i]] += zl_col[i] * tt_val;
          }
        }
      }
    }
    
    // Right extrapolation
    if (std::any_of(outright.begin(), outright.end(), [](char v){ return v; })) {
      double k_pivot = (1 - e) * boundary_knots[1] + e * u[K+ord-1];
      std::vector<int> idx = which(outright);
      std::vector<double> zor = subset(z, idx);
      int r = zor.size();
      
      FlatMatrix zr(r, ord);
      double* zptr = zr.data_ptr();
      for (int i = 0; i < r; ++i) zptr[i] = 1.0;
      for (int j = 1; j <= degree; ++j) {
        double* zcol = zptr + j * r;
        for (int i = 0; i < r; ++i) zcol[i] = std::pow(zor[i] - k_pivot, j);
      }
      
      std::vector<double> kp(ord, k_pivot);
      std::vector<int> derivs(ord);
      for (int i = 0; i < ord; ++i) derivs[i] = i;
      FlatMatrix tt = splineDesigncpp(u, kp, ord, derivs);
      
      for (int j = 0; j < Lfull; ++j) {
        double* design_col = design.data_ptr() + j * n;
        for (int l = 0; l < ord; ++l) {
          const double* zr_col = zr.data_ptr() + l * r;
          double tt_val = tt(l, j) / scalef[l];
          for (int i = 0; i < r; ++i) {
            design_col[idx[i]] += zr_col[i] * tt_val;
          }
        }
      }
    }
    
    // Inside values
    std::vector<unsigned char> inside(n);
    for (int i = 0; i < n; ++i) inside[i] = !outside[i];
    if (std::any_of(inside.begin(), inside.end(), [](char v){ return v; })) {
      std::vector<int> idx = which(inside);
      std::vector<double> z_inside = subset(z, idx);
      FlatMatrix v = splineDesigncpp(u, z_inside, ord, deriv0);
      // column-major copy v -> design
      for (int j = 0; j < Lfull; ++j) {
        const double* v_col = v.data_ptr() + j * v.nrow;
        double* design_col = design.data_ptr() + j * design.nrow;
        for (int i = 0; i < v.nrow; ++i) {
          design_col[idx[i]] = v_col[i];
        }
      }
    }
  } else {
    // All inside: compute design directly
    design = splineDesigncpp(u, z, ord, deriv0);
  }
  
  // If not intercept, drop first column
  int ncol_full = Lfull;
  int col_from = (intercept ? 0 : 1);
  int ncol = ncol_full - col_from;
  if (ncol <= 0) ncol = 0;
  
  // Prepare basis FlatMatrix with NA rows for original input length m
  FlatMatrix basis = make_nan_matrix(m, ncol);
  
  // Fill basis rows for non-NA positions
  // Column-major copy: for each output column, copy from the corresponding column
  for (int j = 0; j < ncol; ++j) {
    double* basis_col_ptr = basis.data_ptr() + j * m;
    const double* design_col_ptr = design.data_ptr() + (j + col_from) * n;
    for (int i = 0, k = 0; i < m; ++i) {
      if (std::isnan(x[i])) continue;
      basis_col_ptr[i] = design_col_ptr[k++];
    }
  }

  // Assemble result ListCpp using push_back(value, name)
  ListCpp out;
  out.push_back(std::move(basis), "basis");
  out.push_back(degree, "degree");
  out.push_back(knots, "knots");
  out.push_back(boundary_knots, "boundary_knots");
  out.push_back(intercept, "intercept");
  
  return out;
}


ListCpp nscpp(
    const std::vector<double>& x,
    int df,
    const std::vector<double>& knots_in,
    bool intercept,
    const std::vector<double>& boundary_knots_in)
{
  // compute basis matrix 
  if (x.empty() || (x.size() == 1 && std::isnan(x[0]))) {
    throw std::invalid_argument("input 'x' is empty or NA");
  }
  
  int m = x.size();
  
  // NA handling
  std::vector<int> non_na_idx = indices_not_nan(x);
  std::vector<double> z = subset(x, non_na_idx);
  int n = z.size();
  
  // boundary knots
  std::vector<double> boundary_knots = boundary_knots_in;
  if (boundary_knots.empty() || 
      (boundary_knots.size() == 1 && std::isnan(boundary_knots[0]))) {
    if (n == 1) boundary_knots = { z[0]*7.0/8.0, z[0]*9.0/8.0 };
    else {
      double mn = *std::min_element(z.begin(), z.end());
      double mx = *std::max_element(z.begin(), z.end());
      boundary_knots = {mn, mx};
    }
  } else {
    if (boundary_knots.size() != 2) 
      throw std::invalid_argument("boundary_knots must have length 2");
  }
  
  // compute K and knots
  std::vector<double> knots = knots_in;
  bool mk_knots = (df > 0) && (knots.empty() || 
                   (knots.size() == 1 && std::isnan(knots[0])));
  int K = 0; // number of interior knots
  if (mk_knots) {
    K = df - 1 - (intercept ? 1 : 0);
    if (K < 0) {
      K = 0;
      thread_utils::push_thread_warning("'df' is too small");
    }
    std::vector<double> knots1(K);
    std::vector<double> z_no_outside;
    for (double zi : z) 
      if (!(zi < boundary_knots[0]) && !(zi > boundary_knots[1])) 
        z_no_outside.push_back(zi);
    for (int k = 0; k < K; ++k) 
      knots1[k] = quantilecpp(z_no_outside, (k + 1.0)/(K + 1.0));
    knots = std::move(knots1);
  } else {
    for (double kv : knots) if (!std::isfinite(kv)) 
      throw std::invalid_argument("non-finite knots");
    K = knots.size();
  }
  
  // adjust knots matching boundaries (like original)
  if (mk_knots && K > 0) {
    bool lrEq0 = std::any_of(knots.begin(), knots.end(), 
                             [&](double v){ return v == boundary_knots[0]; });
    bool lrEq1 = std::any_of(knots.begin(), knots.end(), 
                             [&](double v){ return v == boundary_knots[1]; });
    bool anyWarning = false;
    if (lrEq0) {
      double piv = boundary_knots[0];
      std::vector<char> sub(K);
      bool aE0 = true;
      for (int i = 0; i < K; ++i) { 
        sub[i] = (knots[i] == piv); 
        if (!sub[i]) aE0 = false; 
      }
      if (aE0) {
        thread_utils::push_thread_warning(
          "all interior knots match left boundary knot");
      } else {
        std::vector<double> knots2; 
        for (double kv : knots) if (kv > piv) knots2.push_back(kv);
        double shift = ( *std::min_element(knots2.begin(), knots2.end()) - piv ) / 8.0;
        for (int i = 0; i < K; ++i) if (sub[i]) knots[i] = knots[i] + shift;
      }
    }
    if (lrEq1) {
      double piv = boundary_knots[1];
      std::vector<char> sub(K);
      bool aE1 = true;
      for (int i = 0; i < K; ++i) { 
        sub[i] = (knots[i] == piv); 
        if (!sub[i]) aE1 = false; 
      }
      if (aE1) {
        thread_utils::push_thread_warning(
          "all interior knots match right boundary knot");
      } else {
        std::vector<double> knots2; 
        for (double kv : knots) if (kv < piv) knots2.push_back(kv);
        double shift = ( piv - *std::max_element(knots2.begin(), knots2.end()) ) / 8.0;
        for (int i = 0; i < K; ++i) if (sub[i]) knots[i] = knots[i] - shift;
      }
    }
    if (anyWarning) 
      thread_utils::push_thread_warning(
        "shoving 'interior' knots matching boundary knots to inside");
  }
  
  // Build u for cubic case (ord=4)
  std::vector<double> u(K + 8);
  for (int k = 0; k < 4; ++k) u[k] = boundary_knots[0];
  for (int k = 4; k < K + 4; ++k) u[k] = knots[k - 4];
  for (int k = K + 4; k < K + 8; ++k) u[k] = boundary_knots[1];
  
  // deriv0
  std::vector<int> deriv0(1,0);
  
  // compute design for interior points and extrapolants using splineDesigncpp
  int L = K + 4;
  FlatMatrix design(n, L);
  
  // Identify outside points (on z)
  std::vector<unsigned char> outleft(n, 0), outright(n, 0), outside(n, 0);
  for (int i = 0; i < n; ++i) {
    outleft[i] = (z[i] < boundary_knots[0]);
    outright[i] = (z[i] > boundary_knots[1]);
    outside[i] = (outleft[i] || outright[i]);
  }
  
  if (std::any_of(outside.begin(), outside.end(), [](char v){ return v; })) {

    // Left extrapolation
    if (std::any_of(outleft.begin(), outleft.end(), [](char v){ return v; })) {
      double k_pivot = boundary_knots[0];
      std::vector<int> idx = which(outleft);
      std::vector<double> zol = subset(z, idx);
      int r = zol.size();
      
      FlatMatrix zl(r,2);
      for (int i = 0; i < r; ++i) zl(i,0) = 1.0;
      for (int i = 0; i < r; ++i) zl(i,1) = zol[i] - k_pivot;
      
      std::vector<double> kp(2, k_pivot);
      std::vector<int> derivs = {0, 1};
      FlatMatrix tt = splineDesigncpp(u, kp, 4, derivs);
      
      for (int j = 0; j < L; ++j) {
        double* design_col = design.data_ptr() + j * n;
        for (int l = 0; l < 2; ++l) {
          const double* zl_col = zl.data_ptr() + l * r;
          double tt_val = tt(l, j);
          for (int i = 0; i < r; ++i) {
            design_col[idx[i]] += zl_col[i] * tt_val;
          }
        }
      }
    }
    
    // Right extrapolation
    if (std::any_of(outright.begin(), outright.end(), [](char v){ return v; })) {
      double k_pivot = boundary_knots[1];
      std::vector<int> idx = which(outright);
      std::vector<double> zor = subset(z, idx);
      int r = zor.size();
      
      FlatMatrix zr(r,2);
      for (int i = 0; i < r; ++i) zr(i,0) = 1.0;
      for (int i = 0; i < r; ++i) zr(i,1) = zor[i] - k_pivot;
      
      std::vector<double> kp(2, k_pivot);
      std::vector<int> derivs = {0,1};
      FlatMatrix tt = splineDesigncpp(u, kp, 4, derivs);
      
      for (int j = 0; j < L; ++j) {
        double* design_col = design.data_ptr() + j * n;
        for (int l = 0; l < 2; ++l) {
          const double* zr_col = zr.data_ptr() + l * r;
          double tt_val = tt(l, j);
          for (int i = 0; i < r; ++i) {
            design_col[idx[i]] += zr_col[i] * tt_val;
          }
        }
      }
    }
    
    // Inside values
    std::vector<unsigned char> inside(n);
    for (int i = 0; i < n; ++i) inside[i] = !outside[i];
    if (std::any_of(inside.begin(), inside.end(), [](char v){ return v; })) {
      std::vector<int> idx = which(inside);
      std::vector<double> z_inside = subset(z, idx);
      FlatMatrix v = splineDesigncpp(u, z_inside, 4, deriv0);
      int r = z_inside.size();
      
      for (int j = 0; j < L; ++j) {
        double* design_col = design.data_ptr() + j * n;
        const double* v_col = v.data_ptr() + j * r;
        for (int i = 0; i < r; ++i) {
          design_col[idx[i]] = v_col[i];
        }
      }
    }
  } else {
    // All inside: compute design directly
    design = splineDesigncpp(u, z, 4, deriv0);
  }
  
  // Build constraint matrix con = splineDesigncpp(u, boundary_knots, 4, deriv2)
  std::vector<int> deriv2(1,2);
  FlatMatrix con = splineDesigncpp(u, boundary_knots, 4, deriv2);
  
  if (!intercept) {
    // drop first column of con and design
    int nr = con.nrow, nc = con.ncol;
    double* base = con.data_ptr();
    std::memmove(base, base + nr, sizeof(double) * nr * (nc - 1));
    con.resize(nr, nc-1);
    
    nr = design.nrow; nc = design.ncol;
    base = design.data_ptr();
    std::memmove(base, base + nr, sizeof(double) * nr * (nc - 1));
    design.resize(nr, nc-1);
  }
  
  
  FlatMatrix tcon = transpose(con);
  ListCpp qr = qrcpp(tcon, 1e-12);
  FlatMatrix Q = qr.get<FlatMatrix>("Q");
  
  L = design.ncol;
  FlatMatrix basis2(n, L - 2);
  for (int j = 0; j < L - 2; ++j) {
    double* basis2_col = basis2.data_ptr() + j * n;
    for (int k = 0; k < L; ++k) {
      const double* design_col = design.data_ptr() + k * n;
      double q_val = Q(k, j + 2);
      for (int i = 0; i < n; ++i) {
        basis2_col[i] += design_col[i] * q_val;
      }
    }
  }
  
  // Prepare basis FlatMatrix with NA rows for original input length m
  FlatMatrix basis = make_nan_matrix(m, L-2);
  
  // Fill basis rows for non-NA positions
  // Column-major copy: for each output column, copy from the corresponding column
  for (int j = 0; j < L-2; ++j) {
    double* basis_col_ptr = basis.data_ptr() + j * m;
    const double* basis2_col_ptr = basis2.data_ptr() + j * n;
    for (int i = 0, k = 0; i < m; ++i) {
      if (std::isnan(x[i])) continue;
      basis_col_ptr[i] = basis2_col_ptr[k++];
    }
  }
  
  // Assemble ListCpp result using push_back(value, name)
  ListCpp out;
  out.push_back(std::move(basis), "basis");
  out.push_back(3, "degree"); // cubic
  out.push_back(knots, "knots");
  out.push_back(boundary_knots, "boundary_knots");
  out.push_back(intercept, "intercept");
  
  return out;
}



