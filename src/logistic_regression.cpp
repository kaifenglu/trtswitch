// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>    // RcppParallel::Worker, parallelFor
#include <RcppThread.h>      // RcppThread::Rcerr
#include <Rcpp.h>

#include "logistic_regression.h"
#include "utilities.h"      // boost_pnorm, boost_plogis, etc.
#include "dataframe_list.h"  // FlatMatrix, IntMatrix, DataFrameCpp, ListCpp
#include "thread_utils.h"    // push_thread_warning / drain_thread_warnings_to_R

#include <vector>
#include <string>
#include <numeric>   // iota, inner_product
#include <cmath>     // isnan, isinf, fabs, NaN, exp, log
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

  // compute linear predictor eta efficiently using column-major storage:
  std::vector<double> eta = param->offset; // initialize with offset
  // add contributions of each coefficient times column
  for (int i = 0; i < p; ++i) {
    double beta = par[i];
    if (beta == 0.0) continue;
    int off = i * n;
    for (int r = 0; r < n; ++r) {
      eta[r] += beta * param->z.data[off + r];
    }
  }
  
  double loglik = 0.0;
  std::vector<double> score(p, 0.0);
  // information matrix in column-major p x p
  FlatMatrix imat(p, p);
  
  std::vector<double> pi, d, a, b;
  if (firth) {
    pi.assign(n, 0.0);
    d.assign(n, 0.0);
    a.assign(n, 0.0);
    b.assign(n, 0.0);
  }
  // --- replace the "per-observation temporaries" / person loop + score/imat updates
  // with this column-major friendly two-phase accumulation:
  
  // assume earlier in f_der_0 we computed eta vector and declared:
  // double loglik = 0.0;
  // std::vector<double> score(pp, 0.0);
  // FlatMatrix imat(p, p);  // column-major result
  
  // Pre-allocate per-observation temporaries
  std::vector<double> rvec(n);     // fitted probabilities
  std::vector<double> c1(n, 0.0);  // contribution for score: f*w*(y - r) or variant per link
  std::vector<double> c2(n, 0.0);  // contribution for information: f*w*var-like term
  
  // short-hands to avoid repeated member access
  const double* zptr = param->z.data_ptr();    // column-major: column j starts at zptr + j*n
  const std::vector<double>& freqv = param->freq;
  const std::vector<double>& wtv = param->weight;
  const std::vector<double>& yv   = param->y;
  
  // 1) single pass over observations: compute r, loglik, c1, c2 (+ firth temporaries)
  for (int person = 0; person < n; ++person) {
    double f = freqv[person];
    double w = wtv[person];
    double y = yv[person];
    double et = eta[person];
    
    if (param->link_code == 1) {                // logit
      double r = boost_plogis(et);
      rvec[person] = r;
      loglik += f * w * (y * et + std::log(1.0 - r));
      c1[person] = f * w * (y - r);
      c2[person] = f * w * r * (1.0 - r);
      if (firth) {
        pi[person] = r;
        d[person]  = 1.0;
        a[person]  = r * (1.0 - r);
        b[person]  = 1.0 - 2.0 * r;
      }
    } else if (param->link_code == 2) {         // probit
      double r = boost_pnorm(et);
      double phi = boost_dnorm(et);
      rvec[person] = r;
      // keep original loglik expression used previously for probit case
      loglik += f * w * (y * std::log(r / (1.0 - r)) + std::log(1.0 - r));
      c1[person] = f * w * (y - r) * (phi / (r * (1.0 - r)));
      c2[person] = f * w * (phi * phi / (r * (1.0 - r)));
      if (firth) {
        double dphi = -et;
        pi[person] = r;
        d[person]  = phi / (r * (1.0 - r));
        a[person]  = phi * phi / (r * (1.0 - r));
        b[person]  = (2.0 * r - 1.0) * phi / (r * (1.0 - r)) + 2.0 * dphi;
      }
    } else {                                    // cloglog / extreme
      double r = boost_pextreme(et);
      double phi = boost_dextreme(et);
      rvec[person] = r;
      loglik += f * w * (y * std::log(r / (1.0 - r)) + std::log(1.0 - r));
      c1[person] = f * w * (y - r) * (phi / (r * (1.0 - r)));
      c2[person] = f * w * (phi * phi / (r * (1.0 - r)));
      if (firth) {
        double dphi = 1 - std::exp(et);
        pi[person] = r;
        d[person]  = phi / (r * (1.0 - r));
        a[person]  = phi * phi / (r * (1.0 - r));
        b[person]  = (2.0 * r - 1.0) * phi / (r * (1.0 - r)) + 2.0 * dphi;
      }
    }
  }
  
  // 2) compute score vector by column dot-products: score[j] = sum_i c1[i] * z_{i,j}
  for (int j = 0; j < p; ++j) {
    const double* zcol = zptr + j * n;   // start of column j
    double s = 0.0;
    // inner loop is contiguous load of zcol[i]
    for (int i = 0; i < n; ++i) s += c1[i] * zcol[i];
    score[j] = s;
  }
  
  // 3) compute information matrix (lower triangle) using c2 as per-observation weight
  for (int i = 0; i < p; ++i) {
    const double* zi = zptr + i * n;
    for (int j = 0; j <= i; ++j) {
      const double* zj = zptr + j * n;
      double sum = 0.0;
      for (int k = 0; k < n; ++k) {
        sum += c2[k] * zi[k] * zj[k];
      }
      imat(i, j) = sum;
    }
  }
  
  // 4) fill upper triangle as before (keep symmetric)
  for (int i = 0; i + 1 < p; ++i) {
    for (int j = i + 1; j < p; ++j) {
      imat(i, j) = imat(j, i);
    }
  }
  
  if (firth) {
    // make a copy for Cholesky so we preserve imat for output
    FlatMatrix imat0 = imat;
    cholesky2(imat0, p); // in-place Cholesky on imat0
    
    double sumlog = 0.0;
    for (int i = 0; i < p; ++i) sumlog += std::log(imat0(i, i));
    
    double penloglik = loglik + 0.5 * sumlog;
    
    // precompute per-person scalar c[k] = f * w * a[k]
    std::vector<double> c(n);
    for (int person = 0; person < n; ++person) {
      c[person] = param->freq[person] * param->weight[person] * a[person];
    }
    
    FlatMatrix xwx(p, p); // data initially zeroed by constructor
    const double* Zptr = param->z.data_ptr(); // column-major: col c starts at Zptr + c*n
    
    // compute lower triangle (i >= j) using column-major access
    for (int i = 0; i < p; ++i) {
      const double* zi = Zptr + i * n;        // Z(:, i)
      for (int j = 0; j <= i; ++j) {
        const double* zj = Zptr + j * n;      // Z(:, j)
        double sum = 0.0;
        // inner loop reads zi[k] and zj[k] contiguously
        for (int k = 0; k < n; ++k) {
          sum += c[k] * zi[k] * zj[k];
        }
        xwx(i, j) = sum;
      }
    }
    
    // fill upper triangle of xwx
    for (int i = 0; i + 1 < p; ++i) {
      for (int j = i + 1; j < p; ++j) {
        xwx(i, j) = xwx(j, i);
      }
    }
    
    // compute inverse of xwx as FlatMatrix
    FlatMatrix var = invsympd(xwx, p);
    
    // 1) compute M = Z * var  (n x p)
    FlatMatrix M = mat_mat_mult(param->z, var); // M: n x p
    const double* Mptr = M.data_ptr();
    
    // 2) compute h0[r] = sum_j M(r,j) * Z(r,j)
    std::vector<double> h0(n, 0.0);
    for (int j = 0; j < p; ++j) {
      const double* Mcol = Mptr + j * n;
      const double* Zcol = Zptr + j * n;
      for (int i = 0; i < n; ++i) {
        h0[i] += Mcol[i] * Zcol[i];
      }
    }
    
    // 3) compute u[r] = f*w*resid*d + 0.5*b*(f*w*a*h0)
    std::vector<double> u(n);
    const std::vector<double>& freqv = param->freq;
    const std::vector<double>& wtv  = param->weight;
    const std::vector<double>& yv   = param->y;
    for (int i = 0; i < n; ++i) {
      double f = freqv[i];
      double w = wtv[i];
      double resid = yv[i] - pi[i];
      double h_scaled = h0[i] * f * w * a[i];
      double part1 = f * w * resid * d[i];
      u[i] = part1 + 0.5 * b[i] * h_scaled;
    }
    
    // 4) compute g = Z^T * u
    std::vector<double> g(p, 0.0); 
    for (int j = 0; j < p; ++j) { 
      const double* Zcol = Zptr + j * n; // Z(:, j) 
      double sum = 0.0; 
      for (int i = 0; i < n; ++i) sum += u[i] * Zcol[i]; 
      g[j] = sum;
    }
    
    ListCpp result;
    result.push_back(penloglik, "loglik");
    result.push_back(g, "score");
    result.push_back(imat, "imat");
    result.push_back(loglik, "regloglik");
    result.push_back(score, "regscore");
    return result;
  } else {
    ListCpp result;
    result.push_back(loglik, "loglik");
    result.push_back(score, "score");
    result.push_back(imat, "imat");       // FlatMatrix
    return result;
  }
}



// --------------------------- f_ressco_0 (score residuals) --------------------
// Returns an n x p FlatMatrix (column-major), where entry (r, c) equals 
// residual for observation r and covariate c.
FlatMatrix f_ressco_0(int p, const std::vector<double>& par, void *ex) {
  logparams *param = (logparams *) ex;
  int n = param->n;

  // compute eta similarly to f_der_0
  std::vector<double> eta = param->offset;
  for (int i = 0; i < p; ++i) {
    double beta = par[i];
    if (beta == 0.0) continue;
    int off = i * n;
    for (int r = 0; r < n; ++r) {
      eta[r] += beta * param->z.data[off + r];
    }
  }
  
  FlatMatrix resid(n, p);
  const double* zptr = param->z.data_ptr();
  double* rptr = resid.data_ptr(); 
  
  switch (param->link_code) {
  case 1: { // logit
    std::vector<double> v(n);
    for (int person = 0; person < n; ++person) {
      double r = boost_plogis(eta[person]);
      v[person] = param->y[person] - r;
    }
    
    for (int i = 0; i < p; ++i) {
      int off = i * n;
      const double* zcolptr = zptr + off;
      double* rcolptr = rptr + off;
      for (int person = 0; person < n; ++person) {
        rcolptr[person] = v[person] * zcolptr[person];
      }
    }
    break;
  }
    
  case 2: { // probit
    std::vector<double> vd(n);
    for (int person = 0; person < n; ++person) {
      double r = boost_pnorm(eta[person]);
      double phi = boost_dnorm(eta[person]);
      double d = phi / (r * (1.0 - r));
      double v = param->y[person] - r;
      vd[person] = v * d;
    }
    
    for (int i = 0; i < p; ++i) {
      int off = i * n;
      const double* zcolptr = zptr + off;
      double* rcolptr = rptr + off;
      for (int person = 0; person < n; ++person) {
        rcolptr[person] = vd[person] * zcolptr[person];
      }
    }
    break;
  }
  case 3: {// cloglog / extreme
    std::vector<double> vd(n);
    for (int person = 0; person < n; ++person) {
      double r = boost_pextreme(eta[person]);
      double phi = boost_dextreme(eta[person]);
      double d = phi / (r * (1.0 - r));
      double v = param->y[person] - r;
      vd[person] = v * d;
    }
    
    for (int i = 0; i < p; ++i) {
      int off = i * n;
      const double* zcolptr = zptr + off;
      double* rcolptr = rptr + off;
      for (int person = 0; person < n; ++person) {
        rcolptr[person] = vd[person] * zcolptr[person];
      }
    }
    break;
  }
  }
  return resid;
}

// --------------------------- logisregloop, logisregplloop and logisregcpp ----
ListCpp logisregloop(int p, const std::vector<double>& par, void *ex,
                     int maxiter, double eps, bool firth,
                     const std::vector<int>& colfit, int ncolfit) {
  logparams *param = (logparams *) ex;
  
  int iter = 0, halving = 0;
  bool fail = false;
  
  std::vector<double> beta = par;
  std::vector<double> newbeta(p);
  double loglik = 0.0, newlk = 0.0;
  std::vector<double> u(p);
  std::vector<double> u1(ncolfit);
  FlatMatrix imat(p, p);
  FlatMatrix imat1(ncolfit, ncolfit);
  
  // --- first step ---
  ListCpp der = f_der_0(p, beta, param, firth);
  loglik = der.get<double>("loglik");
  u = der.get<std::vector<double>>("score");
  imat = der.get<FlatMatrix>("imat");
  
  for (int i = 0; i < ncolfit; ++i) u1[i] = u[colfit[i]];
  
  for (int i = 0; i < ncolfit; ++i)
    for (int j = 0; j < ncolfit; ++j)
      imat1(i,j) = imat(colfit[i], colfit[j]);
  
  cholesky2(imat1, ncolfit);
  chsolve2(imat1, ncolfit, u1);
  
  std::fill(u.begin(), u.end(), 0.0);
  for (int i = 0; i < ncolfit; ++i) u[colfit[i]] = u1[i];
  for (int i = 0; i < p; ++i) newbeta[i] = beta[i] + u[i];
  
  // --- main iteration ---
  for (iter = 0; iter < maxiter; ++iter) {
    der = f_der_0(p, newbeta, param, firth);
    newlk = der.get<double>("loglik");
    
    fail = std::isnan(newlk) || std::isinf(newlk);
    if (!fail && halving == 0 && std::fabs(1 - (loglik / newlk)) < eps) break;
    
    if (fail || newlk < loglik) {
      ++halving;
      for (int i = 0; i < p; ++i) {
        newbeta[i] = 0.5 * (beta[i] + newbeta[i]);
      }
      continue;
    }
    
    halving = 0;
    beta = newbeta;
    loglik = newlk;
    u = der.get<std::vector<double>>("score");
    imat = der.get<FlatMatrix>("imat");
    
    for (int i = 0; i < ncolfit; ++i) u1[i] = u[colfit[i]];
    
    for (int i = 0; i < ncolfit; ++i)
      for (int j = 0; j < ncolfit; ++j)
        imat1(i, j) = imat(colfit[i], colfit[j]);
    
    cholesky2(imat1, ncolfit);
    chsolve2(imat1, ncolfit, u1);
    
    std::fill(u.begin(), u.end(), 0.0);
    for (int i = 0; i < ncolfit; ++i) u[colfit[i]] = u1[i];
    for (int i = 0; i < p; ++i) newbeta[i] = beta[i] + u[i];
  }
  
  if (iter == maxiter) fail = true;
  
  // final variance assembly
  imat = der.get<FlatMatrix>("imat");
  
  for (int i = 0; i < ncolfit; ++i)
    for (int j = 0; j < ncolfit; ++j)
      imat1(i, j) = imat(colfit[i], colfit[j]);
  
  FlatMatrix var1 = invsympd(imat1, ncolfit);
  FlatMatrix var(p, p);
  for (int j = 0; j < ncolfit; ++j)
    for (int i = 0; i < ncolfit; ++i)
      var(colfit[i], colfit[j]) = var1(i, j);
  
  ListCpp result;
  result.push_back(newbeta, "coef");
  result.push_back(iter, "iter");
  result.push_back(var, "var");
  result.push_back(newlk, "loglik");
  result.push_back(fail, "fail");
  
  if (firth) {
    double regloglik = der.get<double>("regloglik");
    result.push_back(regloglik, "regloglik");
  }
  
  return result;
}

// --------------------------- logisregplloop (profile likelihood solver) -----
double logisregplloop(int p, const std::vector<double>& par,
                      void *ex, int maxiter, double eps, bool firth,
                      int k, int direction, double l0) {
  logparams *param = (logparams *) ex;
  int iter;
  bool fail = false;
  
  std::vector<double> beta = par;
  std::vector<double> newbeta(p);
  double loglik = 0.0, newlk = 0.0;
  std::vector<double> u(p);
  std::vector<double> delta(p);
  FlatMatrix imat(p, p);
  FlatMatrix v(p, p);
  
  ListCpp der = f_der_0(p, beta, param, firth);
  loglik = der.get<double>("loglik");
  u = der.get<std::vector<double>>("score");
  imat = der.get<FlatMatrix>("imat");
  v = invsympd(imat, p);
  
  // compute w = - u^T v u
  double w = -quadsym(u, v);
  double underroot = -2 * (l0 - loglik + 0.5 * w) / v(k, k);
  double lambda = underroot < 0.0 ? 0.0 : direction * std::sqrt(underroot);
  u[k] += lambda;
  delta = mat_vec_mult(v, u);
  for (int i = 0; i < p; ++i) newbeta[i] = beta[i] + delta[i];
  
  // iterate to convergence 
  for (iter = 0; iter < maxiter; ++iter) {
    der = f_der_0(p, newbeta, param, firth);
    newlk = der.get<double>("loglik");
    fail = std::isnan(newlk) || std::isinf(newlk);
    if (!fail && std::fabs(newlk - l0) < eps && w < eps) break;
    beta = newbeta;
    loglik = newlk;
    u = der.get<std::vector<double>>("score");
    imat = der.get<FlatMatrix>("imat");
    v = invsympd(imat, p);
    w = -quadsym(u, v);
    underroot = -2 * (l0 - newlk + 0.5 * w) / v(k, k);
    lambda = underroot < 0.0 ? 0.0 : direction * std::sqrt(underroot);
    u[k] += lambda;
    delta = mat_vec_mult(v, u);
    for (int i = 0; i < p; ++i) newbeta[i] = beta[i] + delta[i];
  }
  
  if (iter == maxiter) fail = true;
  if (fail) {
    thread_utils::push_thread_warning(
      "logisregplloop did not converge within the maximum iterations.");
    return NaN;
  }
  
  return newbeta[k];
}

// --------------------------- logisregcpp (high-level API) --------------------
// Convert inputs to FlatMatrix design matrix and use new functions where appropriate.
ListCpp logisregcpp(const DataFrameCpp& data,
                    const std::string& event,
                    const std::vector<std::string>& covariates,
                    const std::string& freq,
                    const std::string& weight,
                    const std::string& offset,
                    const std::string& id,
                    const std::string& link,
                    const std::vector<double>& init,
                    const bool robust,
                    const bool firth,
                    const bool flic,
                    const bool plci,
                    const double alpha,
                    const int maxiter,
                    const double eps) {
  
  int n = data.nrows();
  int p = covariates.size() + 1;
  if (p == 2 && covariates[0].empty()) p = 1;
  
  if (event.empty()) throw std::invalid_argument("event variable is not specified");
  if (!data.containElementNamed(event)) 
    throw std::invalid_argument("data must contain the event variable");
  
  // event -> numeric 0/1
  std::vector<double> eventn(n);
  if (data.bool_cols.count(event)) {
    const std::vector<unsigned char>& vb = data.get<unsigned char>(event);
    for (int i = 0; i < n; ++i) eventn[i] = vb[i] ? 1.0 : 0.0;
  } else if (data.int_cols.count(event)) {
    const std::vector<int>& vi = data.get<int>(event);
    for (int i = 0; i < n; ++i) eventn[i] = static_cast<double>(vi[i]);
  } else if (data.numeric_cols.count(event)) {
    eventn = data.get<double>(event);
  } else {
    throw std::invalid_argument("event variable must be bool, integer or numeric");
  }
  for (double val : eventn) if (val != 0 && val != 1) 
    throw std::invalid_argument("event must be 1 or 0 for each observation");
  
  // construct design matrix zn (n x p) as FlatMatrix column-major
  FlatMatrix zn(n, p);
  // intercept column
  for (int i = 0; i < n; ++i) zn.data[i] = 1.0;
  
  // fill covariate columns (1..p-1)
  for (int j = 0; j < p - 1; ++j) {
    const std::string& zj = covariates[j];
    if (!data.containElementNamed(zj)) 
      throw std::invalid_argument("data must contain the variables in covariates");
    std::vector<double> u(n);
    if (data.bool_cols.count(zj)) {
      const std::vector<unsigned char>& ub = data.get<unsigned char>(zj);
      for (int i = 0; i < n; ++i) u[i] = ub[i] ? 1.0 : 0.0;
    } else if (data.int_cols.count(zj)) {
      const std::vector<int>& ui = data.get<int>(zj);
      for (int i = 0; i < n; ++i) u[i] = static_cast<double>(ui[i]);
    } else if (data.numeric_cols.count(zj)) {
      u = data.get<double>(zj);
    } else {
      throw std::invalid_argument("covariates must be bool, integer or numeric");
    }
    int off = (j + 1) * n;
    for (int i = 0; i < n; ++i) zn.data[off + i] = u[i];
  }
  
  // freq, weight, offset
  std::vector<double> freqn(n, 1.0);
  if (!freq.empty() && data.containElementNamed(freq)) {
    if (data.int_cols.count(freq)) {
      const auto& freqi = data.get<int>(freq);
      for (int i = 0; i < n; ++i) freqn[i] = static_cast<double>(freqi[i]);
    } else if (data.numeric_cols.count(freq)) {
      freqn = data.get<double>(freq);
    } else throw std::invalid_argument("freq variable must be integer or numeric");
    for (double v : freqn) if (v <= 0) 
      throw std::invalid_argument("freq must be positive integers");
  }
  
  std::vector<double> weightn(n, 1.0);
  if (!weight.empty() && data.containElementNamed(weight)) {
    weightn = data.get<double>(weight);
    for (double v : weightn) if (v <= 0.0) 
      throw std::invalid_argument("weight must be greater than 0");
  }
  
  std::vector<double> offsetn(n, 0.0);
  if (!offset.empty() && data.containElementNamed(offset)) {
    offsetn = data.get<double>(offset);
  }
  
  // id processing (unchanged semantics)
  bool has_id = !id.empty() && data.containElementNamed(id);
  std::vector<int> idn(n);
  if (!has_id) {
    std::iota(idn.begin(), idn.end(), 0);
  } else {
    if (data.int_cols.count(id)) {
      auto v = data.get<int>(id);
      auto w = unique_sorted(v);
      idn = matchcpp(v, w);
    } else if (data.numeric_cols.count(id)) {
      auto v = data.get<double>(id);
      auto w = unique_sorted(v);
      idn = matchcpp(v, w);
    } else if (data.string_cols.count(id)) {
      auto v = data.get<std::string>(id);
      auto w = unique_sorted(v);
      idn = matchcpp(v, w);
    } else {
      throw std::invalid_argument("incorrect type for the id variable in the input data");
    }
  }
  
  // link code mapping
  std::string link1 = link;
  std::for_each(link1.begin(), link1.end(), [](char & c) { 
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c))); 
  });
  
  if (link1 == "log-log" || link1 == "loglog") link1 = "cloglog";
  
  int link_code = 0;
  if (link1 == "logit") link_code = 1;
  else if (link1 == "probit") link_code = 2;
  else if (link1 == "cloglog") link_code = 3;
  else throw std::invalid_argument("invalid link: " + link1);
  
  // exclude observations with missing values
  std::vector<unsigned char> sub(n, 1);
  for (int i = 0; i < n; ++i) {
    if (eventn[i] == INT_MIN ||
        std::isnan(freqn[i]) || std::isnan(weightn[i]) ||
        std::isnan(offsetn[i]) || idn[i] == INT_MIN) {
      sub[i] = 0;
      continue;
    }
    for (int j = 0; j < p - 1; ++j) {
      int off = (j + 1) * n;
      if (std::isnan(zn.data[off + i])) {
        sub[i] = 0;
        break;
      }
    }
  }
  
  std::vector<int> order = which(sub);
  subset_in_place(eventn, order);
  subset_in_place(freqn, order);
  subset_in_place(weightn, order);
  subset_in_place(offsetn, order);
  subset_in_place(idn, order);
  subset_in_place_flatmatrix(zn, order);
  n = order.size();
  if (n == 0) throw std::invalid_argument("no observations without missing values");
  
  // sumstat data set
  double nobs, nevents;
  double loglik0, loglik1;
  double regloglik0, regloglik1; // regular loglik
  int niter;
  bool fail;
  
  // parest data set
  std::vector<std::string> par(p);
  std::vector<double> b(p), seb(p), rseb(p);
  std::vector<double> z(p), expbeta(p);
  FlatMatrix vb(p, p), rvb(p, p);
  std::vector<double> lb(p), ub(p), prob(p);
  std::vector<std::string> clparm(p);
  
  // linear predictor and fitted values for all observations
  std::vector<double> linear_predictors(n), fitted_values(n);
  
  double zcrit = boost_qnorm(1-alpha/2);
  double xcrit = zcrit * zcrit;
  
  // number of trials and number of events accounting for frequencies
  nobs = std::accumulate(freqn.begin(), freqn.end(), 0.0);
  nevents = std::inner_product(freqn.begin(), freqn.end(), eventn.begin(), 0.0);
  if (nevents == 0) {
    for (int i=0; i<p; ++i) {
      par[i] = (i == 0) ? "(Intercept)" : covariates[i-1];
      b[i] = NaN;
      seb[i] = NaN;
      rseb[i] = NaN;
      z[i] = NaN;
      expbeta[i] = NaN;
      lb[i] = NaN;
      ub[i] = NaN;
      prob[i] = NaN;
      clparm[i] = "Wald";
    }
    
    for (int j=0; j<p; ++j) {
      for (int i=0; i<p; ++i) {
        vb(i,j) = NaN;    
        rvb(i,j) = NaN;
      }
    }
    
    for (int person = 0; person < n; ++person) {
      linear_predictors[person] = offsetn[person];
    }
    switch (link_code) {
    case 1:
      for (int person = 0; person < n; ++person)
        fitted_values[person] = boost_plogis(linear_predictors[person]);
      break;
    case 2:
      for (int person = 0; person < n; ++person)
        fitted_values[person] = boost_pnorm(linear_predictors[person]);
      break;
    case 3:
      for (int person = 0; person < n; ++person)
        fitted_values[person] = boost_pextreme(linear_predictors[person]);
      break;
    }
    
    loglik0 = NaN;
    loglik1 = NaN;
    regloglik0 = NaN;
    regloglik1 = NaN;
    niter = 0;
    fail = true;
  } else {
    // intercept only model
    double num = 0, den = 0;
    for (int i=0; i<n; ++i) {
      num += freqn[i] * weightn[i] * eventn[i];
      den += freqn[i] * weightn[i];
    }
    if (firth) {
      num += 0.5;
      den += 1.0;
    }
    
    std::vector<double> bint0(p);
    bint0[0] = boost_qlogis(num/den);
    
    std::vector<int> colfit0(1);
    logparams param = {n, link_code, eventn, zn, freqn, weightn, offsetn};
    ListCpp outint = logisregloop(p, bint0, &param, maxiter, eps, firth, colfit0, 1);
    
    std::vector<double> bint = outint.get<std::vector<double>>("coef");
    FlatMatrix vbint = outint.get<FlatMatrix>("var");
    
    ListCpp out;
    
    if (p > 1) {
      // parameter estimates and standard errors for the full model
      std::vector<int> colfit = seqcpp(0,p-1);
      if (static_cast<int>(init.size()) == p && 
          std::none_of(init.begin(), init.end(), [](double val){ 
            return std::isnan(val); })) {
        out = logisregloop(p, init, &param, maxiter, eps, firth, colfit, p);
      } else {
        out = logisregloop(p, bint, &param, maxiter, eps, firth, colfit, p);
      }
      
      bool fail = out.get<bool>("fail");
      if (fail) {
        thread_utils::push_thread_warning(
          "logisregloop failed to converge for the full model; continuing with current results.");
      }
      
      b = out.get<std::vector<double>>("coef");
      vb = out.get<FlatMatrix>("var");
      
      // intercept correction
      if (flic) {
        std::vector<double> lp = offsetn;  // linear predictor excluding intercept
        
        // pointers into data for fastest access
        const double* zptr = zn.data_ptr(); // column-major: column j starts at zptr + j*n
        double* lpptr = lp.data();
        
        for (int col = 1; col < p; ++col) {   // skip intercept column 0
          double coef = b[col];
          if (coef == 0.0) continue;                       // skip zero coefficient
          const double* zcol = zptr + col * n;
          // inner loop is contiguous in memory (good locality)
          for (int row = 0; row < n; ++row) {
            lpptr[row] += coef * zcol[row];
          }
        }
        
        logparams param0 = {n, link_code, eventn, zn, freqn, weightn, lp};
        std::vector<double> bint00(1, bint0[0]);
        ListCpp outint0 = logisregloop(1, bint00, &param0, maxiter, eps, 0, colfit0, 1);
        double a = outint0.get<std::vector<double>>("coef")[0];
        double va = outint0.get<FlatMatrix>("var")(0,0);
        
        // update the intercept estimate
        b[0] = a;
        
        // partial derivative of alpha(beta) with respect to beta
        ListCpp derint = f_der_0(p, b, &param, 0);
        FlatMatrix iflic = derint.get<FlatMatrix>("imat");
        std::vector<double> der(p-1);
        for (int i=0; i<p-1; ++i) {
          der[i] = -iflic(i+1, 0)/iflic(0, 0);
        }
        
        // update variance-covariance matrix via the delta method
        if (p <= 1) { // only intercept
          vb(0,0) = va;
        } else {
          int p1 = p - 1;                  // length of der and of the submatrix
          const double* derp = der.data();
          double* vbptr = vb.data_ptr();    // column-major: column c starts at vbptr[c*p]
          
          // accumulator for vb(1..p-1,0): indexed 0..p1-1 corresponds to rows 1..p-1
          std::vector<double> col0(p1, 0.0);
          
          // Compute col0 = M * der, where M = vb[1..p-1, 1..p-1]
          // columns of M correspond to vb columns 1..p-1; for column j (0..p1-1) start at
          // src = vbptr + (j+1)*nrows and element M(row = r+1, col = j+1) is src[r+1].
          for (int j = 0; j < p1; ++j) {
            double dj = derp[j];
            if (dj == 0.0) continue;                   // skip zero coeffs
            const double* src_col = vbptr + (j + 1) * p;
            // accumulate rows 1..p-1 -> indices r=0..p1-1 map to src_col[r+1]
            for (int r = 0; r < p1; ++r) {
              col0[r] += src_col[r + 1] * dj;
            }
          }
          
          // write back vb(1..p-1,0) and vb(0,1..p-1) and compute der^T * col0
          double dot = 0.0;
          // vb( row, col ) -> vbptr[ col*nrows + row ]
          for (int r = 0; r < p1; ++r) {
            double v = col0[r];
            vbptr[0 * p + (r + 1)] = v;          // vb(r+1, 0)
            vbptr[(r + 1) * p + 0] = v;          // vb(0, r+1)
            dot += derp[r] * v;
          }
          
          // final scalar update
          vbptr[0] = va + dot;                       // vb(0,0)  
        }
      }
    } else {
      b = bint;
      vb = vbint;
      out = outint;
      
      if (flic) {
        out = logisregloop(p, bint0, &param, maxiter, eps, 0, colfit0, 1);
        b = out.get<std::vector<double>>("coef");
        vb = outint.get<FlatMatrix>("var");
      }
    }
    
    for (int i=0; i<p; ++i) {
      seb[i] = std::sqrt(vb(i,i));
    }
    
    for (int i=0; i<p; ++i) {
      par[i] = (i == 0) ? "(Intercept)" : covariates[i-1];
    }
    
    // linear predictors and fitted values
    std::vector<double> eta = offsetn;
    for (int i = 0; i < p; ++i) {
      double beta = b[i];
      if (beta == 0.0) continue;
      int off = i * n;
      for (int r = 0; r < n; ++r) {
        eta[r] += beta * zn.data[off + r];
      }
    }
    linear_predictors = eta;
    
    switch (link_code) {
    case 1: // logit
      for (int person = 0; person < n; ++person) {
        fitted_values[person] = boost_plogis(eta[person]);
      }
      break;
    case 2: // probit
      for (int person = 0; person < n; ++person) {
        fitted_values[person] = boost_pnorm(eta[person]);
      }
      break;
    case 3: // cloglog
      for (int person = 0; person < n; ++person) {
        fitted_values[person] = boost_pextreme(eta[person]);
      }
      break;
    }
    
    
    niter = out.get<int>("iter");
    fail = out.get<bool>("fail");
    
    // robust variance estimates
    if (robust) {
      FlatMatrix ressco = f_ressco_0(p, b, &param);
      
      int nr; // number of rows in the score residual matrix
      std::vector<double> freqr;
      if (!has_id) {
        double* rptr = ressco.data_ptr();         // column-major base pointer
        const double* wptr = weightn.data();      // pointer to weights
        
        for (int j = 0; j < p; ++j) {
          double* colptr = rptr + j * n;    // start of column j
          // inner loop reads/writes contiguous memory -> good locality & vectorization
          for (int i = 0; i < n; ++i) {
            colptr[i] *= wptr[i];
          }
        }
        
        nr = n;
        freqr = freqn;
      } else { // need to sum up score residuals by id
        std::vector<int> order(n);
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(), [&](int i, int j) {
          return idn[i] < idn[j];
        });
        
        std::vector<int> id2 = subset(idn, order);
        std::vector<int> idx(1,0);
        for (int i=1; i<n; ++i) {
          if (id2[i] != id2[i-1]) {
            idx.push_back(i);
          }
        }
        
        int nids = static_cast<int>(idx.size());
        idx.push_back(n);
        
        std::vector<double> weight2 = subset(weightn, order);
        std::vector<double> freq2 = subset(freqn, order);
        std::vector<double> freqr0(nids); // cluster frequency
        
        FlatMatrix ressco2(nids, p);
        // raw pointers for fastest access
        const double* src_ptr = ressco.data_ptr();     // column-major: col*nrows + row
        double* dst_ptr = ressco2.data_ptr();          // column-major: col*nids + row
        const double* wptr = weight2.data();
        const int* order_ptr = order.data();           // if order is vector<int>
        const int* idx_ptr = idx.data();               // idx length must be nids+1
        
        // For each column, walk source column contiguously and accumulate grouped sums
        for (int j = 0; j < p; ++j) {
          const double* src_col = src_ptr + j * n;   // start of ressco(:,j)
          double* dst_col = dst_ptr + j * nids;       // start of ressco2(:,j)
          
          for (int i = 0; i < nids; ++i) {
            double sum = 0.0;
            int kstart = idx_ptr[i];
            int kend = idx_ptr[i + 1];
            // accumulate weight2[k] * ressco[ order[k] , j ]
            for (int k = kstart; k < kend; ++k) {
              sum += wptr[k] * src_col[order_ptr[k]];
            }
            dst_col[i] = sum;
          }
        }
        
        for (int i=0; i<nids; ++i) {
          freqr0[i] = freq2[idx_ptr[i]];
        }
        
        // update the score residuals
        ressco = std::move(ressco2);  
        nr = nids;
        freqr = freqr0;
      }
      
      FlatMatrix D = mat_mat_mult(ressco, vb);

      const double* Dptr = D.data_ptr();      // Dcol_j starts at Dptr + j*nr_sz
      double* rvbptr     = rvb.data_ptr();    // rvbcol_k starts at rvbptr + k*p
      const double* fptr = freqr.data();
      
      // compute only upper triangle j <= k
      // rvb(j,k) corresponds to rvbptr[k*p + j] because column-major: column k, row j
      for (int j = 0; j < p; ++j) {
        const double* Dj = Dptr + j * nr;          // pointer to D(:,j)
        for (int k = j; k < p; ++k) {
          const double* Dk = Dptr + k * nr;        // pointer to D(:,k)
          double sum = 0.0;
          // inner loop reads Dj[i] and Dk[i] contiguously
          for (int i = 0; i < nr; ++i) {
            sum += fptr[i] * Dj[i] * Dk[i];
          }
          // write into rvb(j,k) and mirror
          rvbptr[k * p + j] = sum;                // rvb(j,k)
          if (k != j) {
            rvbptr[j * p + k] = sum;              // rvb(k,j) = rvb(j,k)  (mirror)
          }
        }
      }
      
      for (int i=0; i<p; ++i) {
        rseb[i] = std::sqrt(rvb(i, i));
      }
    }
    
    // profile likelihood confidence interval for regression coefficients
    if (plci) {
      double lmax = out.get<double>("loglik");
      double l0 = lmax - 0.5*xcrit;
      
      if (!(firth && flic)) { // PL CI for all parameters
        for (int k=0; k<p; ++k) {
          lb[k] = logisregplloop(p, b, &param, maxiter, eps, firth, k, -1, l0);
          ub[k] = logisregplloop(p, b, &param, maxiter, eps, firth, k, 1, l0);
          
          std::vector<int> colfit1(p-1);
          for (int i = 0, j = 0; i < p; ++i) {
            if (i == k) continue;
            colfit1[j++] = i;
          }
          
          std::vector<double> b0(p);
          ListCpp out0 = logisregloop(p, b0, &param, maxiter, eps, firth, colfit1, p-1);
          double lmax0 = out0.get<double>("loglik");
          prob[k] = boost_pchisq(-2*(lmax0 - lmax), 1, 0);
          clparm[k] = "PL";
        }
      } else { // Wald CI for intercept and PL CI for slopes
        if (!robust) {
          lb[0] = b[0] - zcrit*seb[0];
          ub[0] = b[0] + zcrit*seb[0];
          prob[0] = boost_pchisq(pow(b[0]/seb[0], 2), 1, 0);
        } else {
          lb[0] = b[0] - zcrit*rseb[0];
          ub[0] = b[0] + zcrit*rseb[0];
          prob[0] = boost_pchisq(pow(b[0]/rseb[0], 2), 1, 0);
        }
        clparm[0] = "Wald";
        
        for (int k=1; k<p; ++k) {
          lb[k] = logisregplloop(p, b, &param, maxiter, eps, firth, k, -1, l0);
          ub[k] = logisregplloop(p, b, &param, maxiter, eps, firth, k, 1, l0);
          
          std::vector<int> colfit1(p-1);
          for (int i=0; i<k; ++i) {
            colfit1[i] = i;
          }
          for (int i=k+1; i<p; ++i) {
            colfit1[i-1] = i;
          }
          
          std::vector<double> b0(p);
          ListCpp out0 = logisregloop(p, b0, &param, maxiter, eps, firth, colfit1, p-1);
          double lmax0 = out0.get<double>("loglik");
          prob[k] = boost_pchisq(-2*(lmax0 - lmax), 1, 0);
          clparm[k] = "PL";
        }
      }
    } else { // Wald confidence interval for all parameters
      for (int k=0; k<p; ++k) {
        if (!robust) {
          lb[k] = b[k] - zcrit*seb[k];
          ub[k] = b[k] + zcrit*seb[k];
          prob[k] = boost_pchisq(pow(b[k]/seb[k], 2), 1, 0);
        } else {
          lb[k] = b[k] - zcrit*rseb[k];
          ub[k] = b[k] + zcrit*rseb[k];
          prob[k] = boost_pchisq(pow(b[k]/rseb[k], 2), 1, 0);
        }
        clparm[k] = "Wald";
      }
    }
    
    // log-likelihoods
    if (p > 0 && firth) {
      loglik0 = outint.get<double>("loglik");
      loglik1 = out.get<double>("loglik");
      regloglik0 = outint.get<double>("regloglik");
      regloglik1 = out.get<double>("regloglik");
    } else {
      loglik0 = outint.get<double>("loglik");
      loglik1 = out.get<double>("loglik");
      regloglik0 = loglik0;
      regloglik1 = loglik1;
    }
  }
  
  // compute exp(beta)
  for (int i=0; i<p; ++i) {
    expbeta[i] = std::exp(b[i]);
  }
  
  // compute z statistics
  if (robust) {
    for (int i=0; i<p; ++i) {
      if (rseb[i] == 0) {
        z[i] = NaN;
      } else {
        z[i] = b[i]/rseb[i];
      }
    }
  } else {
    for (int i=0; i<p; ++i) {
      if (seb[i] == 0) {
        z[i] = NaN;
      } else {
        z[i] = b[i]/seb[i];
      }
    }
  }
  
  // prepare the output data sets
  DataFrameCpp sumstat;
  sumstat.push_back(nobs, "n");
  sumstat.push_back(nevents, "nevents");
  sumstat.push_back(loglik0, "loglik0");
  sumstat.push_back(loglik1, "loglik1");
  sumstat.push_back(niter, "niter");
  sumstat.push_back(p, "p");
  sumstat.push_back(link1, "link");
  sumstat.push_back(robust, "robust");
  sumstat.push_back(firth, "firth");
  sumstat.push_back(flic, "flic");
  sumstat.push_back(fail, "fail");
  
  if (p > 0 && firth) {
    sumstat.push_back(regloglik0, "loglik0_unpenalized");
    sumstat.push_back(regloglik1, "loglik1_unpenalized");
  }
  
  DataFrameCpp parest;
  parest.push_back(par, "param");
  parest.push_back(b, "beta");
  parest.push_back(robust ? rseb : seb, "sebeta");
  parest.push_back(z, "z");
  parest.push_back(expbeta, "expbeta");
  parest.push_back(robust ? rvb : vb, "vbeta");
  parest.push_back(lb, "lower");
  parest.push_back(ub, "upper");
  parest.push_back(prob, "p");
  parest.push_back(clparm, "method");
  
  if (robust) {
    parest.push_back(seb, "sebeta_naive");
    parest.push_back(vb, "vbeta_naive");
  }
  
  DataFrameCpp fitted;
  fitted.push_back(linear_predictors, "linear_predictors");
  fitted.push_back(fitted_values, "fitted_values");
  
  ListCpp result;
  result.push_back(sumstat, "sumstat");
  result.push_back(parest, "parest");
  result.push_back(fitted, "fitted");
  
  return result;
}


// [[Rcpp::export]]
Rcpp::List logisregRcpp(const Rcpp::DataFrame& data,
                        const std::string& event,
                        const std::vector<std::string>& covariates,
                        const std::string& freq,
                        const std::string& weight,
                        const std::string& offset,
                        const std::string& id,
                        const std::string& link,
                        const std::vector<double>& init,
                        const bool robust,
                        const bool firth,
                        const bool flic,
                        const bool plci,
                        const double alpha,
                        const int maxiter,
                        const double eps) {
  
  DataFrameCpp dfcpp = convertRDataFrameToCpp(data);
  
  ListCpp cpp_result = logisregcpp(
    dfcpp, event, covariates, freq, weight, offset, id, link,
    init, robust, firth, flic, plci, alpha, maxiter, eps
  );
  
  thread_utils::drain_thread_warnings_to_R();
  return Rcpp::wrap(cpp_result);
}
