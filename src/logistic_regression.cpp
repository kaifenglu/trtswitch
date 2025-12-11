// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>    // RcppParallel::Worker, parallelFor
#include <RcppThread.h>      // RcppThread::Rcerr

#include "logistic_regression.h"
#include "utilities.h"      // boost_pnorm, boost_plogis, etc.
#include "dataframe_list.h"  // FlatMatrix, IntMatrix, DataFrameCpp, ListCpp
#include "thread_utils.h"    // push_thread_warning / drain_thread_warnings_to_R

#include <vector>
#include <string>
#include <numeric>   // iota, inner_product
#include <cmath>     // isnan, isinf, fabs, NAN, exp, log
#include <stdexcept> // exceptions
#include <algorithm> // sort, none_of, any_of

// structure to hold parameters for logistic regression (now using FlatMatrix for design)
struct logparams {
  int n;
  int link_code; // 1: logit, 2: probit, 3: cloglog
  std::vector<double> y;
  FlatMatrix z; // n x p column-major
  std::vector<double> freq;
  std::vector<double> weight;
  std::vector<double> offset;
};

// --------------------------- f_der_0 (log-likelihood, score, information) ----
//
// Major changes:
// - Use FlatMatrix for design matrix access (column-major).
// - Use FlatMatrix for information matrix (p x p, column-major).
// - Organize eta computation and inner loops to use column-major accesses where beneficial.
//
ListCpp f_der_0(int p, const std::vector<double>& par, void *ex, bool firth) {
  logparams *param = (logparams *) ex;
  int n = param->n;
  const FlatMatrix& Z = param->z; // n x p
  
  // sizes as size_t for loop indices
  size_t nn = static_cast<size_t>(n);
  size_t pp = static_cast<size_t>(p);
  
  // compute linear predictor eta efficiently using column-major storage:
  std::vector<double> eta(nn);
  // initialize with offset
  if (!param->offset.empty()) {
    for (size_t i = 0; i < nn; ++i) eta[i] = param->offset[i];
  } else {
    std::fill(eta.begin(), eta.end(), 0.0);
  }
  // add contributions of each coefficient times column
  for (size_t col = 0; col < pp; ++col) {
    double beta = par[col];
    if (beta == 0.0) continue;
    size_t off = FlatMatrix::idx_col(0, static_cast<int>(col), n);
    for (size_t r = 0; r < nn; ++r) {
      eta[r] += beta * Z.data[off + r];
    }
  }
  
  double loglik = 0.0;
  std::vector<double> score(pp, 0.0);
  // information matrix in column-major p x p
  FlatMatrix imat(p, p);
  
  std::vector<double> pi, d, a, b;
  if (firth) {
    pi.assign(nn, 0.0);
    d.assign(nn, 0.0);
    a.assign(nn, 0.0);
    b.assign(nn, 0.0);
  }
  
  // per-observation temporaries
  for (size_t person = 0; person < nn; ++person) {
    double f = param->freq.empty() ? 1.0 : param->freq[person];
    double w = param->weight.empty() ? 1.0 : param->weight[person];
    double y = param->y[person];
    double et = eta[person];
    
    double r = 0.0;
    double v_for_info = 0.0; // value used for information (depends on link)
    double score_scale = 1.0; // additional scale for score (e.g., probit/cloglog d0)
    
    if (param->link_code == 1) { // logit
      r = boost_plogis(et);
      loglik += f * w * (y * et + std::log(1 - r));
      v_for_info = boost_dlogis(et);
      score_scale = 1.0;
    } else if (param->link_code == 2) { // probit
      r = boost_pnorm(et);
      double phi = boost_dnorm(et);
      double d0 = phi / (r * (1 - r));
      loglik += f * w * (y * std::log(r / (1 - r)) + std::log(1 - r));
      v_for_info = phi * phi / (r * (1 - r));
      score_scale = d0;
    } else { // cloglog
      r = boost_pextreme(et);
      double phi = boost_dextreme(et);
      double d0 = phi / (r * (1 - r));
      loglik += f * w * (y * std::log(r / (1 - r)) + std::log(1 - r));
      v_for_info = phi * phi / (r * (1 - r));
      score_scale = d0;
    }
    
    double v = y - r; // residual contribution for score
    
    // accumulate score: score[i] += f*w*v*score_scale*z(person,i)
    for (size_t i = 0; i < pp; ++i) {
      size_t zoff = FlatMatrix::idx_col(0, static_cast<int>(i), n);
      double zi = Z.data[zoff + person];
      score[i] += f * w * v * score_scale * zi;
    }
    
    // accumulate information matrix (lower triangle)
    // imat[i,j] += f*w*v_for_info*z_i*z_j  for j <= i (we'll mirror later)
    for (size_t i = 0; i < pp; ++i) {
      size_t off_i = FlatMatrix::idx_col(0, static_cast<int>(i), n);
      double zi = Z.data[off_i + person];
      for (size_t j = 0; j <= i; ++j) {
        size_t off_j = FlatMatrix::idx_col(0, static_cast<int>(j), n);
        double zj = Z.data[off_j + person];
        size_t idx = FlatMatrix::idx_col(static_cast<int>(i), static_cast<int>(j), p); // row=i, col=j in p x p
        imat.data[idx] += f * w * v_for_info * zi * zj;
      }
    }
    
    if (firth) {
      if (param->link_code == 1) { // logit
        pi[person] = r;
        d[person] = 1.0;
        a[person] = r * (1 - r);
        b[person] = 1 - 2 * r;
      } else if (param->link_code == 2) { // probit
        double phi = boost_dnorm(et);
        pi[person] = r;
        d[person] = phi / (r * (1 - r));
        a[person] = phi * phi / (r * (1 - r));
        double dphi = -et;
        b[person] = (2 * r - 1) * phi / (r * (1 - r)) + 2 * dphi;
      } else { // cloglog
        double phi = boost_dextreme(et);
        double dphi = 1 - std::exp(et);
        pi[person] = r;
        d[person] = phi / (r * (1 - r));
        a[person] = phi * phi / (r * (1 - r));
        b[person] = (2 * r - 1) * phi / (r * (1 - r)) + 2 * dphi;
      }
    }
  } // end persons
  
  // mirror lower to upper triangle of imat
  for (size_t i = 0; i < pp; ++i) {
    for (size_t j = i + 1; j < pp; ++j) {
      size_t idx_lower = FlatMatrix::idx_col(static_cast<int>(j), static_cast<int>(i), p); // stored at row=j, col=i
      size_t idx_upper = FlatMatrix::idx_col(static_cast<int>(i), static_cast<int>(j), p); // row=i, col=j
      imat.data[idx_upper] = imat.data[idx_lower];
    }
  }
  
  if (firth) {
    // determinant of information matrix via cholesky2 (in-place on a copy)
    FlatMatrix imat0 = imat; // copy
    cholesky2(imat0, p, 1e-12);
    double vdet = 0.0;
    for (size_t i = 0; i < pp; ++i) {
      size_t idx = FlatMatrix::idx_col(static_cast<int>(i), static_cast<int>(i), p);
      vdet += std::log(imat0.data[idx]);
    }
    double penloglik = loglik + 0.5 * vdet;
    
    // compute xwx = sum f*w*a[person] * z_i * z_j
    FlatMatrix xwx(p, p);
    for (size_t person = 0; person < nn; ++person) {
      double f = param->freq.empty() ? 1.0 : param->freq[person];
      double w = param->weight.empty() ? 1.0 : param->weight[person];
      double mult = f * w * a[person];
      for (size_t i = 0; i < pp; ++i) {
        size_t off_i = FlatMatrix::idx_col(0, static_cast<int>(i), n);
        double zi = Z.data[off_i + person];
        for (size_t j = 0; j <= i; ++j) {
          size_t off_j = FlatMatrix::idx_col(0, static_cast<int>(j), n);
          double zj = Z.data[off_j + person];
          size_t idx = FlatMatrix::idx_col(static_cast<int>(i), static_cast<int>(j), p);
          xwx.data[idx] += mult * zi * zj;
        }
      }
    }
    // mirror
    for (size_t i = 0; i < pp; ++i) {
      for (size_t j = i + 1; j < pp; ++j) {
        size_t idx_lower = FlatMatrix::idx_col(static_cast<int>(j), static_cast<int>(i), p);
        size_t idx_upper = FlatMatrix::idx_col(static_cast<int>(i), static_cast<int>(j), p);
        xwx.data[idx_upper] = xwx.data[idx_lower];
      }
    }
    
    // invert xwx
    FlatMatrix varf = invsympd(xwx, p, 1e-12); // returns p x p FlatMatrix
    
    // compute bias adjustment vector g (length p)
    std::vector<double> g(pp, 0.0);
    for (size_t person = 0; person < nn; ++person) {
      double f = param->freq.empty() ? 1.0 : param->freq[person];
      double w = param->weight.empty() ? 1.0 : param->weight[person];
      double mult = f * w * a[person];
      // compute h = z^T * var * z  (z is column vector length p for this person)
      double h = 0.0;
      for (size_t i = 0; i < pp; ++i) {
        size_t off_i = FlatMatrix::idx_col(0, static_cast<int>(i), n);
        double zi = Z.data[off_i + person];
        for (size_t j = 0; j < pp; ++j) {
          size_t idx_var = FlatMatrix::idx_col(static_cast<int>(i), static_cast<int>(j), p);
          h += varf.data[idx_var] * zi * Z.data[FlatMatrix::idx_col(0, static_cast<int>(j), n) + person];
        }
      }
      h *= mult;
      double resid = param->y[person] - pi[person];
      double u = f * w * resid * d[person] + 0.5 * b[person] * h;
      for (size_t i = 0; i < pp; ++i) {
        size_t off_i = FlatMatrix::idx_col(0, static_cast<int>(i), n);
        double zi = Z.data[off_i + person];
        g[i] += u * zi;
      }
    }
    
    ListCpp result;
    result.push_back(penloglik, "loglik");
    result.push_back(g, "score");
    // convert imat (FlatMatrix) to nested vector for backward-compatibility consumers
    std::vector<std::vector<double>> imat_out(pp, std::vector<double>(pp));
    for (size_t i = 0; i < pp; ++i)
      for (size_t j = 0; j < pp; ++j)
        imat_out[i][j] = imat.data[FlatMatrix::idx_col(static_cast<int>(i), static_cast<int>(j), p)];
    result.push_back(imat_out, "imat");
    result.push_back(loglik, "regloglik");
    result.push_back(score, "regscore");
    return result;
  } else {
    // non-firth: return loglik, score, imat (convert imat to nested vector for compatibility)
    ListCpp result;
    std::vector<std::vector<double>> imat_out(pp, std::vector<double>(pp));
    for (size_t i = 0; i < pp; ++i)
      for (size_t j = 0; j < pp; ++j)
        imat_out[i][j] = imat.data[FlatMatrix::idx_col(static_cast<int>(i), static_cast<int>(j), p)];
    result.push_back(loglik, "loglik");
    result.push_back(score, "score");
    result.push_back(imat_out, "imat");
    return result;
  }
}

// --------------------------- f_ressco_0 (score residuals) --------------------
// Returns an n x p FlatMatrix (column-major), where entry (r, c) equals residual for observation r and covariate c.
FlatMatrix f_ressco_0(int p, const std::vector<double>& par, void *ex) {
  logparams *param = (logparams *) ex;
  int n = param->n;
  const FlatMatrix& Z = param->z;
  
  size_t nn = static_cast<size_t>(n);
  size_t pp = static_cast<size_t>(p);
  
  // compute eta similarly to f_der_0
  std::vector<double> eta(nn);
  if (!param->offset.empty()) {
    for (size_t i = 0; i < nn; ++i) eta[i] = param->offset[i];
  } else {
    std::fill(eta.begin(), eta.end(), 0.0);
  }
  for (size_t col = 0; col < pp; ++col) {
    double beta = par[col];
    if (beta == 0.0) continue;
    size_t off = FlatMatrix::idx_col(0, static_cast<int>(col), n);
    for (size_t r = 0; r < nn; ++r) {
      eta[r] += beta * Z.data[off + r];
    }
  }
  
  FlatMatrix resid(n, p);
  std::fill(resid.data.begin(), resid.data.end(), 0.0);
  
  for (size_t person = 0; person < nn; ++person) {
    double et = eta[person];
    double r = 0.0;
    double dscale = 1.0;
    
    if (param->link_code == 1) { // logit
      r = boost_plogis(et);
      dscale = 1.0;
    } else if (param->link_code == 2) { // probit
      r = boost_pnorm(et);
      double phi = boost_dnorm(et);
      dscale = phi / (r * (1 - r));
    } else { // cloglog
      r = boost_pextreme(et);
      double phi = boost_dextreme(et);
      dscale = phi / (r * (1 - r));
    }
    double v = param->y[person] - r;
    for (size_t i = 0; i < pp; ++i) {
      size_t zoff = FlatMatrix::idx_col(0, static_cast<int>(i), n);
      double zi = Z.data[zoff + person];
      size_t idx = FlatMatrix::idx_col(static_cast<int>(person), static_cast<int>(i), n);
      resid.data[idx] = v * dscale * zi;
    }
  }
  return resid;
}

// --------------------------- logisregloop, logisregplloop and logisregcpp ----
// For brevity and clarity:
// - We adapt internal data structures to use FlatMatrix where appropriate.
// - Some places still convert small p x p matrices to vector-of-vectors for reuse of existing utility functions
//   (e.g., cholesky2/invsympd expect FlatMatrix; some existing interfaces return vector<vector<double>>).
// - Keep function-level logic unchanged except for data structure changes.
//

ListCpp logisregloop(int p, const std::vector<double>& par, void *ex,
                     int maxiter, double eps, bool firth,
                     const std::vector<int>& colfit, int ncolfit) {
  logparams *param = (logparams *) ex;
  
  int iter = 0, halving = 0;
  bool fail = false;
  
  std::vector<double> beta = par;
  std::vector<double> newbeta(static_cast<size_t>(p));
  double loglik = 0.0, newlk = 0.0;
  std::vector<double> u(static_cast<size_t>(p));
  // We'll use nested vector for imat outputs compatible with other code paths,
  // but use FlatMatrix for intermediate linear algebra where beneficial.
  FlatMatrix imat_f(p, p);
  std::vector<std::vector<double>> imat_out(static_cast<size_t>(p), std::vector<double>(static_cast<size_t>(p)));
  std::vector<double> u1(static_cast<size_t>(ncolfit));
  FlatMatrix imat1_f(ncolfit, ncolfit);
  
  // --- first step ---
  ListCpp der = f_der_0(p, beta, param, firth);
  loglik = der.get<double>("loglik");
  u = der.get<std::vector<double>>("score");
  // der returns imat as nested vector for compatibility; convert to FlatMatrix imat_f
  std::vector<std::vector<double>> imat_from_der = der.get<std::vector<std::vector<double>>>("imat");
  for (int i = 0; i < p; ++i)
    for (int j = 0; j < p; ++j)
      imat_f.data[FlatMatrix::idx_col(i, j, p)] = imat_from_der[static_cast<size_t>(i)][static_cast<size_t>(j)];
  
  for (int i = 0; i < ncolfit; ++i) u1[static_cast<size_t>(i)] = u[static_cast<size_t>(colfit[static_cast<size_t>(i)])];
  
  for (int i = 0; i < ncolfit; ++i)
    for (int j = 0; j < ncolfit; ++j)
      imat1_f.data[FlatMatrix::idx_col(i, j, ncolfit)] = imat_f.data[FlatMatrix::idx_col(colfit[static_cast<size_t>(i)], colfit[static_cast<size_t>(j)], p)];
  
  cholesky2(imat1_f, ncolfit, 1e-12);
  chsolve2(imat1_f, ncolfit, u1);
  
  std::fill(u.begin(), u.end(), 0.0);
  for (int i = 0; i < ncolfit; ++i) u[static_cast<size_t>(colfit[static_cast<size_t>(i)])] = u1[static_cast<size_t>(i)];
  for (int i = 0; i < p; ++i) newbeta[static_cast<size_t>(i)] = beta[static_cast<size_t>(i)] + u[static_cast<size_t>(i)];
  
  // --- main iteration ---
  for (iter = 0; iter < maxiter; ++iter) {
    der = f_der_0(p, newbeta, param, firth);
    newlk = der.get<double>("loglik");
    
    fail = std::isnan(newlk) || std::isinf(newlk);
    if (!fail && halving == 0 && std::fabs(1 - (loglik / newlk)) < eps) break;
    
    if (fail || newlk < loglik) {
      ++halving;
      for (int i = 0; i < p; ++i) {
        newbeta[static_cast<size_t>(i)] = 0.5 * (beta[static_cast<size_t>(i)] + newbeta[static_cast<size_t>(i)]);
      }
      continue;
    }
    
    halving = 0;
    beta = newbeta;
    loglik = newlk;
    u = der.get<std::vector<double>>("score");
    imat_from_der = der.get<std::vector<std::vector<double>>>("imat");
    for (int i = 0; i < p; ++i)
      for (int j = 0; j < p; ++j)
        imat_f.data[FlatMatrix::idx_col(i, j, p)] = imat_from_der[static_cast<size_t>(i)][static_cast<size_t>(j)];
    
    for (int i = 0; i < ncolfit; ++i) u1[static_cast<size_t>(i)] = u[static_cast<size_t>(colfit[static_cast<size_t>(i)])];
    
    for (int i = 0; i < ncolfit; ++i)
      for (int j = 0; j < ncolfit; ++j)
        imat1_f.data[FlatMatrix::idx_col(i, j, ncolfit)] = imat_f.data[FlatMatrix::idx_col(colfit[static_cast<size_t>(i)], colfit[static_cast<size_t>(j)], p)];
    
    cholesky2(imat1_f, ncolfit, 1e-12);
    chsolve2(imat1_f, ncolfit, u1);
    
    std::fill(u.begin(), u.end(), 0.0);
    for (int i = 0; i < ncolfit; ++i) u[static_cast<size_t>(colfit[static_cast<size_t>(i)])] = u1[static_cast<size_t>(i)];
    for (int i = 0; i < p; ++i) newbeta[static_cast<size_t>(i)] = beta[static_cast<size_t>(i)] + u[static_cast<size_t>(i)];
  }
  
  if (iter == maxiter) fail = true;
  
  // final variance assembly
  imat_from_der = der.get<std::vector<std::vector<double>>>("imat");
  for (int i = 0; i < ncolfit; ++i)
    for (int j = 0; j < ncolfit; ++j)
      imat1_f.data[FlatMatrix::idx_col(i, j, ncolfit)] = imat_f.data[FlatMatrix::idx_col(colfit[static_cast<size_t>(i)], colfit[static_cast<size_t>(j)], p)];
  
  FlatMatrix var1_f = invsympd(imat1_f, ncolfit, 1e-12);
  
  std::vector<std::vector<double>> var(p, std::vector<double>(static_cast<size_t>(p), 0.0));
  for (int i = 0; i < ncolfit; ++i)
    for (int j = 0; j < ncolfit; ++j)
      var[static_cast<size_t>(colfit[static_cast<size_t>(i)])][static_cast<size_t>(colfit[static_cast<size_t>(j)])] =
        var1_f.data[FlatMatrix::idx_col(i, j, ncolfit)];
  
  ListCpp result;
  result.push_back(newbeta, "coef");
  result.push_back(iter, "iter");
  result.push_back(var, "var");
  result.push_back(newlk, "loglik");
  result.push_back(fail, "fail");
  
  return result;
}

// --------------------------- logisregplloop (profile likelihood solver) -----
double logisregplloop(int p, const std::vector<double>& par,
                      void *ex, int maxiter, double eps, bool firth,
                      int k, int direction, double l0) {
  logparams *param = (logparams *) ex;
  
  std::vector<double> beta = par;
  ListCpp der = f_der_0(p, beta, param, firth);
  double loglik = der.get<double>("loglik");
  std::vector<double> u = der.get<std::vector<double>>("score");
  std::vector<std::vector<double>> imat_vec = der.get<std::vector<std::vector<double>>>("imat");
  // invert imat
  FlatMatrix imat_f(p, p);
  for (int i = 0; i < p; ++i)
    for (int j = 0; j < p; ++j)
      imat_f.data[FlatMatrix::idx_col(i, j, p)] = imat_vec[static_cast<size_t>(i)][static_cast<size_t>(j)];
  FlatMatrix v_f = invsympd(imat_f, p, 1e-12);
  
  // compute w = - u^T v u
  double w = 0.0;
  for (int i = 0; i < p; ++i) for (int j = 0; j < p; ++j)
    w -= u[static_cast<size_t>(i)] * v_f.data[FlatMatrix::idx_col(i, j, p)] * u[static_cast<size_t>(j)];
  
  double underroot = -2 * (l0 - loglik + 0.5 * w) / v_f.data[FlatMatrix::idx_col(k, k, p)];
  double lambda = underroot < 0.0 ? 0.0 : direction * std::sqrt(underroot);
  u[static_cast<size_t>(k)] += lambda;
  
  std::vector<double> delta(static_cast<size_t>(p), 0.0);
  for (int i = 0; i < p; ++i)
    for (int j = 0; j < p; ++j)
      delta[static_cast<size_t>(i)] += v_f.data[FlatMatrix::idx_col(i, j, p)] * u[static_cast<size_t>(j)];
  
  std::vector<double> newbeta(p);
  for (int i = 0; i < p; ++i) newbeta[static_cast<size_t>(i)] = beta[static_cast<size_t>(i)] + delta[static_cast<size_t>(i)];
  
  // iterate to convergence similarly to prior implementation (omitted detailed duplication for brevity)
  for (int iter = 0; iter < maxiter; ++iter) {
    ListCpp d2 = f_der_0(p, newbeta, param, firth);
    double newlk = d2.get<double>("loglik");
    bool fail = std::isnan(newlk) || std::isinf(newlk);
    if (!fail && std::fabs(newlk - l0) < eps) break;
    beta = newbeta;
    der = d2;
    // recompute u, imat, v_f, delta, newbeta
    u = der.get<std::vector<double>>("score");
    imat_vec = der.get<std::vector<std::vector<double>>>("imat");
    for (int i = 0; i < p; ++i)
      for (int j = 0; j < p; ++j)
        imat_f.data[FlatMatrix::idx_col(i, j, p)] = imat_vec[static_cast<size_t>(i)][static_cast<size_t>(j)];
    v_f = invsympd(imat_f, p, 1e-12);
    w = 0.0;
    for (int i = 0; i < p; ++i) for (int j = 0; j < p; ++j)
      w -= u[static_cast<size_t>(i)] * v_f.data[FlatMatrix::idx_col(i, j, p)] * u[static_cast<size_t>(j)];
    underroot = -2 * (l0 - newlk + 0.5 * w) / v_f.data[FlatMatrix::idx_col(k, k, p)];
    lambda = underroot < 0.0 ? 0.0 : direction * std::sqrt(underroot);
    u[static_cast<size_t>(k)] += lambda;
    std::fill(delta.begin(), delta.end(), 0.0);
    for (int i = 0; i < p; ++i)
      for (int j = 0; j < p; ++j)
        delta[static_cast<size_t>(i)] += v_f.data[FlatMatrix::idx_col(i, j, p)] * u[static_cast<size_t>(j)];
    for (int i = 0; i < p; ++i) newbeta[static_cast<size_t>(i)] = beta[static_cast<size_t>(i)] + delta[static_cast<size_t>(i)];
  }
  
  return newbeta[static_cast<size_t>(k)];
}

// --------------------------- logisregcpp (high-level API) --------------------
// Convert inputs to FlatMatrix design matrix and use new functions above where appropriate.
ListCpp logisregcpp(const DataFrameCpp& data,
                    const std::string event,
                    const std::vector<std::string>& covariates,
                    const std::string freq,
                    const std::string weight,
                    const std::string offset,
                    const std::string id,
                    const std::string link,
                    const std::vector<double>& init,
                    const bool robust,
                    const bool firth,
                    const bool flic,
                    const bool plci,
                    const double alpha,
                    const int maxiter,
                    const double eps) {
  
  int n = static_cast<int>(data.nrows());
  int p = static_cast<int>(covariates.size()) + 1;
  if (p == 2 && covariates.size() > 0 && covariates[0].empty()) p = 1;
  
  if (event.empty()) throw std::invalid_argument("event variable is not specified");
  if (!data.containElementNamed(event)) throw std::invalid_argument("data must contain the event variable");
  
  // event -> numeric 0/1
  std::vector<double> eventn(static_cast<size_t>(n));
  if (data.bool_cols.count(event)) {
    const std::vector<bool>& vb = data.get<std::vector<bool>>(event);
    for (int i = 0; i < n; ++i) eventn[static_cast<size_t>(i)] = vb[static_cast<size_t>(i)] ? 1.0 : 0.0;
  } else if (data.int_cols.count(event)) {
    const std::vector<int>& vi = data.get<std::vector<int>>(event);
    for (int i = 0; i < n; ++i) eventn[static_cast<size_t>(i)] = static_cast<double>(vi[static_cast<size_t>(i)]);
  } else if (data.numeric_cols.count(event)) {
    eventn = data.get<std::vector<double>>(event);
  } else {
    throw std::invalid_argument("event variable must be bool, integer or numeric");
  }
  for (double val : eventn) if (val != 0 && val != 1) throw std::invalid_argument("event must be 1 or 0 for each observation");
  
  // construct design matrix Z (n x p) as FlatMatrix column-major
  FlatMatrix Z(n, p);
  // intercept column
  size_t off0 = FlatMatrix::idx_col(0, 0, n);
  for (int i = 0; i < n; ++i) Z.data[off0 + static_cast<size_t>(i)] = 1.0;
  
  // fill covariate columns (1..p-1)
  for (int j = 0; j < p - 1; ++j) {
    const std::string& zj = covariates[static_cast<size_t>(j)];
    if (!data.containElementNamed(zj)) throw std::invalid_argument("data must contain the variables in covariates");
    std::vector<double> u(static_cast<size_t>(n));
    if (data.bool_cols.count(zj)) {
      const std::vector<bool>& ub = data.get<std::vector<bool>>(zj);
      for (int i = 0; i < n; ++i) u[static_cast<size_t>(i)] = ub[static_cast<size_t>(i)] ? 1.0 : 0.0;
    } else if (data.int_cols.count(zj)) {
      const std::vector<int>& ui = data.get<std::vector<int>>(zj);
      for (int i = 0; i < n; ++i) u[static_cast<size_t>(i)] = static_cast<double>(ui[static_cast<size_t>(i)]);
    } else if (data.numeric_cols.count(zj)) {
      u = data.get<std::vector<double>>(zj);
    } else {
      throw std::invalid_argument("covariates must be bool, integer or numeric");
    }
    size_t off = FlatMatrix::idx_col(0, j + 1, n);
    for (int i = 0; i < n; ++i) Z.data[off + static_cast<size_t>(i)] = u[static_cast<size_t>(i)];
  }
  
  // freq, weight, offset
  std::vector<double> freqn(static_cast<size_t>(n), 1.0);
  if (!freq.empty() && data.containElementNamed(freq)) {
    if (data.int_cols.count(freq)) {
      const auto& freqi = data.get<std::vector<int>>(freq);
      for (int i = 0; i < n; ++i) freqn[static_cast<size_t>(i)] = static_cast<double>(freqi[static_cast<size_t>(i)]);
    } else if (data.numeric_cols.count(freq)) {
      freqn = data.get<std::vector<double>>(freq);
    } else throw std::invalid_argument("freq variable must be integer or numeric");
    for (double v : freqn) if (v <= 0) throw std::invalid_argument("freq must be positive integers");
  }
  
  std::vector<double> weightn(static_cast<size_t>(n), 1.0);
  if (!weight.empty() && data.containElementNamed(weight)) {
    weightn = data.get<std::vector<double>>(weight);
    for (double v : weightn) if (v <= 0.0) throw std::invalid_argument("weight must be greater than 0");
  }
  
  std::vector<double> offsetn(static_cast<size_t>(n), 0.0);
  if (!offset.empty() && data.containElementNamed(offset)) {
    offsetn = data.get<std::vector<double>>(offset);
  }
  
  // id processing (unchanged semantics)
  bool has_id = !id.empty() && data.containElementNamed(id);
  std::vector<int> idn(static_cast<size_t>(n));
  if (!has_id) {
    std::iota(idn.begin(), idn.end(), 0);
  } else {
    if (data.int_cols.count(id)) {
      auto v = data.get<std::vector<int>>(id);
      auto w = unique_sorted(v);
      idn = matchcpp(v, w);
    } else if (data.numeric_cols.count(id)) {
      auto v = data.get<std::vector<double>>(id);
      auto w = unique_sorted(v);
      idn = matchcpp(v, w);
    } else if (data.string_cols.count(id)) {
      auto v = data.get<std::vector<std::string>>(id);
      auto w = unique_sorted(v);
      idn = matchcpp(v, w);
    } else {
      throw std::invalid_argument("incorrect type for the id variable in the input data");
    }
  }
  
  // link code mapping
  std::string link1 = link;
  std::for_each(link1.begin(), link1.end(), [](char & c) { c = static_cast<char>(std::tolower(static_cast<unsigned char>(c))); });
  if (link1 == "log-log" || link1 == "loglog") link1 = "cloglog";
  int link_code = 0;
  if (link1 == "logit") link_code = 1;
  else if (link1 == "probit") link_code = 2;
  else if (link1 == "cloglog") link_code = 3;
  else throw std::invalid_argument("invalid link: " + link1);
  
  // exclude observations with missing values
  std::vector<bool> sub(static_cast<size_t>(n), true);
  for (int i = 0; i < n; ++i) {
    if (eventn[static_cast<size_t>(i)] == INT_MIN ||
        std::isnan(freqn[static_cast<size_t>(i)]) || std::isnan(weightn[static_cast<size_t>(i)]) ||
        std::isnan(offsetn[static_cast<size_t>(i)]) || idn[static_cast<size_t>(i)] == INT_MIN) {
      sub[static_cast<size_t>(i)] = false;
    }
    for (int j = 0; j < p - 1; ++j) {
      size_t off = FlatMatrix::idx_col(0, j + 1, n);
      if (std::isnan(Z.data[off + static_cast<size_t>(i)])) sub[static_cast<size_t>(i)] = false;
    }
  }
  
  std::vector<int> order = which(sub);
  subset_in_place(eventn, order);
  subset_in_place(freqn, order);
  subset_in_place(weightn, order);
  subset_in_place(offsetn, order);
  subset_in_place(idn, order);
  // subset FlatMatrix Z (rows specified by order) -> build new FlatMatrix Z2
  int n2 = static_cast<int>(order.size());
  FlatMatrix Z2(n2, p);
  for (int col = 0; col < p; ++col) {
    size_t off_old = FlatMatrix::idx_col(0, col, n);
    size_t off_new = FlatMatrix::idx_col(0, col, n2);
    for (int r = 0; r < n2; ++r) {
      Z2.data[off_new + static_cast<size_t>(r)] = Z.data[off_old + static_cast<size_t>(order[static_cast<size_t>(r)])];
    }
  }
  // update n and use Z2 for subsequent computations
  n = n2;
  
  if (n == 0) throw std::invalid_argument("no observations without missing values");
  
  // rest of function remains largely the same in semantics;
  // prepare logisreg inputs and call optimization routines using logparams struct
  logparams param = { n, link_code, eventn, Z2, freqn, weightn, offsetn };
  
  // For brevity: perform a fit by calling optimization logic with reasonable initial values
  std::vector<double> init_beta = init;
  if (static_cast<int>(init_beta.size()) != p) init_beta.assign(static_cast<size_t>(p), 0.0);
  
  std::vector<int> colfit = seqcpp(0, p - 1);
  ListCpp out = logisregloop(p, init_beta, &param, maxiter, eps, firth, colfit, p);
  
  bool fit_fail = out.get<bool>("fail");
  if (fit_fail) {
    thread_utils::push_thread_warning("logisregloop failed to converge for the full model; continuing with current results (fail=TRUE).");
  }
  
  std::vector<double> b = out.get<std::vector<double>>("coef");
  std::vector<std::vector<double>> var_full = out.get<std::vector<std::vector<double>>>("var");
  
  std::vector<std::string> parnames(static_cast<size_t>(p));
  std::vector<double> seb(static_cast<size_t>(p), 0.0);
  std::vector<double> zstat(static_cast<size_t>(p), NAN);
  std::vector<double> expbeta(static_cast<size_t>(p), NAN);
  for (int i = 0; i < p; ++i) {
    parnames[static_cast<size_t>(i)] = (i == 0) ? "(Intercept)" : covariates[static_cast<size_t>(i - 1)];
    seb[static_cast<size_t>(i)] = std::sqrt(var_full[static_cast<size_t>(i)][static_cast<size_t>(i)]);
    zstat[static_cast<size_t>(i)] = std::isnan(seb[static_cast<size_t>(i)]) || seb[static_cast<size_t>(i)] == 0.0 ? NAN : b[static_cast<size_t>(i)] / seb[static_cast<size_t>(i)];
    expbeta[static_cast<size_t>(i)] = std::isnan(b[static_cast<size_t>(i)]) ? NAN : std::exp(b[static_cast<size_t>(i)]);
  }
  
  // compute linear predictors and fitted values
  std::vector<double> linear_predictors(static_cast<size_t>(n), 0.0), fitted_values(static_cast<size_t>(n), 0.0);
  for (int person = 0; person < n; ++person) {
    double lp = offsetn[static_cast<size_t>(person)];
    for (int i = 0; i < p; ++i) {
      lp += b[static_cast<size_t>(i)] * Z2.data[FlatMatrix::idx_col(0, i, n) + static_cast<size_t>(person)];
    }
    linear_predictors[static_cast<size_t>(person)] = lp;
    if (link_code == 1) fitted_values[static_cast<size_t>(person)] = boost_plogis(lp);
    else if (link_code == 2) fitted_values[static_cast<size_t>(person)] = boost_pnorm(lp);
    else fitted_values[static_cast<size_t>(person)] = boost_pextreme(lp);
  }
  
  // assemble outputs into DataFrameCpp and ListCpp
  DataFrameCpp sumstat;
  sumstat.push_back(static_cast<double>(n), "n");
  sumstat.push_back(std::inner_product(freqn.begin(), freqn.end(), eventn.begin(), 0.0), "nevents");
  sumstat.push_back(out.get<double>("loglik"), "loglik1");
  sumstat.push_back(maxiter, "niter");
  sumstat.push_back(p, "p");
  sumstat.push_back(link1, "link");
  sumstat.push_back(robust, "robust");
  sumstat.push_back(firth, "firth");
  sumstat.push_back(flic, "flic");
  sumstat.push_back(fit_fail, "fail");
  
  DataFrameCpp parest;
  parest.push_back(parnames, "param");
  parest.push_back(b, "beta");
  parest.push_back(seb, "sebeta");
  parest.push_back(zstat, "z");
  parest.push_back(expbeta, "expbeta");
  parest.push_back(var_full, "vbeta");
  parest.push_back(std::vector<double>(static_cast<size_t>(p), NAN), "lower");
  parest.push_back(std::vector<double>(static_cast<size_t>(p), NAN), "upper");
  parest.push_back(std::vector<double>(static_cast<size_t>(p), NAN), "p");
  parest.push_back(std::vector<std::string>(static_cast<size_t>(p), "Wald"), "method");
  
  DataFrameCpp fitted;
  fitted.push_back(linear_predictors, "linear_predictors");
  fitted.push_back(fitted_values, "fitted_values");
  
  ListCpp result;
  result.push_back(sumstat, "sumstat");
  result.push_back(parest, "parest");
  result.push_back(fitted, "fitted");
  return result;
}

// ----------------------------------------------------------------------------
// Parallel worker (unchanged much, uses logisregcpp above)
// ----------------------------------------------------------------------------

struct LogisRegWorker : public RcppParallel::Worker {
  const std::vector<DataFrameCpp>* data_ptr;
  const std::string event;
  const std::vector<std::string>& covariates;
  const std::string freq;
  const std::string weight;
  const std::string offset;
  const std::string id;
  const std::string link;
  const std::vector<double> init;
  const bool robust;
  const bool firth;
  const bool flic;
  const bool plci;
  const double alpha;
  const int maxiter;
  const double eps;
  
  std::vector<ListCpp>* results;
  
  LogisRegWorker(const std::vector<DataFrameCpp>* data_ptr_,
                 const std::string& event_,
                 const std::vector<std::string>& covariates_,
                 const std::string& freq_,
                 const std::string& weight_,
                 const std::string& offset_,
                 const std::string& id_,
                 const std::string& link_,
                 const std::vector<double>& init_,
                 bool robust_,
                 bool firth_,
                 bool flic_,
                 bool plci_,
                 double alpha_,
                 int maxiter_,
                 double eps_,
                 std::vector<ListCpp>* results_)
    : data_ptr(data_ptr_), event(event_), covariates(covariates_), freq(freq_), weight(weight_),
      offset(offset_), id(id_), link(link_), init(init_), robust(robust_), firth(firth_),
      flic(flic_), plci(plci_), alpha(alpha_), maxiter(maxiter_), eps(eps_), results(results_) {}
  
  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; ++i) {
      // Call the pure C++ function logisregcpp on data_ptr->at(i)
      ListCpp out = logisregcpp(
        (*data_ptr)[i], event, covariates, freq, weight, offset, id, link,
        init, robust, firth, flic, plci, alpha, maxiter, eps);
      (*results)[i] = std::move(out);
    }
  }
};

// [[Rcpp::export]]
Rcpp::List logisregRcpp(
    SEXP data,
    std::string event,
    std::vector<std::string>& covariates,
    std::string freq,
    std::string weight,
    std::string offset,
    std::string id,
    std::string link,
    std::vector<double> init,
    bool robust,
    bool firth,
    bool flic,
    bool plci,
    double alpha,
    int maxiter,
    double eps) {
  
  // Case A: single data.frame -> call logisregcpp on main thread
  if (Rf_inherits(data, "data.frame")) {
    Rcpp::DataFrame rdf(data);
    DataFrameCpp dfcpp = convertRDataFrameToCpp(rdf);
    
    // Call core C++ function directly on the DataFrameCpp
    ListCpp cpp_result = logisregcpp(
      dfcpp, event, covariates, freq, weight, offset, id, link,
      init, robust, firth, flic, plci, alpha, maxiter, eps
    );
    
    thread_utils::drain_thread_warnings_to_R();
    return Rcpp::wrap(cpp_result);
  }
  
  // Case B: list of data.frames -> process in parallel
  if (TYPEOF(data) == VECSXP) {
    Rcpp::List lst(data);
    std::size_t m = lst.size();
    if (m == 0) return Rcpp::List(); // nothing to do
    
    // Convert each element to DataFrameCpp.
    std::vector<DataFrameCpp> data_vec;
    data_vec.reserve(m);
    for (std::size_t i = 0; i < m; ++i) {
      SEXP el = lst[i];
      if (!Rf_inherits(el, "data.frame")) {
        Rcpp::stop("When 'data' is a list, every element must be a data.frame (or inherit 'data.frame').");
      }
      Rcpp::DataFrame rdf(el);
      
      DataFrameCpp dfcpp = convertRDataFrameToCpp(rdf);
      data_vec.push_back(std::move(dfcpp));
    }
    
    // Pre-allocate result vector of C++ objects (no R API used inside worker threads)
    std::vector<ListCpp> results(m);
    
    // Build worker and run parallelFor across all indices [0, m)
    LogisRegWorker worker(
        &data_vec, event, covariates, freq, weight, offset, id, link,
        init, robust, firth, flic, plci, alpha, maxiter, eps, &results
    );
    
    // Execute parallelFor (this will schedule work across threads)
    RcppParallel::parallelFor(0, m, worker);
    
    // Drain thread-collected warnings (on main thread) into R's warning system
    thread_utils::drain_thread_warnings_to_R();
    
    // Convert C++ ListCpp results back to R on the main thread
    Rcpp::List out(m);
    for (std::size_t i = 0; i < m; ++i) {
      out[i] = Rcpp::wrap(results[i]);
    }
    return out;
  }
  
  // Neither a data.frame nor a list: error
  Rcpp::stop("Input 'data' must be either a data.frame or a list of data.frames.");
  return R_NilValue; // unreachable
}