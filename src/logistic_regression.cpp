// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>    // RcppParallel::Worker, parallelFor
#include <RcppThread.h>      // RcppThread::Rcerr

#include "logistic_regression.h"
#include "utilities.h"      // boost_pnorm, boost_plogis, etc.
#include "dataframe_list.h"  // DataFrameCpp, ListCpp
#include "thread_utils.h"    // push_thread_warning / drain_thread_warnings_to_R

#include <vector>
#include <string>
#include <numeric>   // iota, inner_product
#include <cmath>     // isnan, isinf, fabs, NAN, exp, log
#include <stdexcept> // std::invalid_argument / std::runtime_error
#include <algorithm> // sort, none_of, any_of

// structure to hold parameters for logistic regression
struct logparams {
  int n;
  int link_code; // 0: logit, 1: probit, 2: cloglog
  std::vector<double> y;
  std::vector<std::vector<double>> z;
  std::vector<double> freq;
  std::vector<double> weight;
  std::vector<double> offset;
};

// all-in-one function for log-likelihood, score, and information matrix for logistic model
ListCpp f_der_0(int p, const std::vector<double>& par, void *ex, bool firth) {
  logparams *param = (logparams *) ex;
  int n = param->n;
  
  std::vector<double> eta(n);
  for (int person = 0; person < n; ++person) {
    eta[person] = param->offset[person];
    for (int i=0; i<p; ++i) {
      eta[person] += par[i]*param->z[person][i];
    }
  }
  
  double loglik = 0.0;
  std::vector<double> score(p, 0.0);
  std::vector<std::vector<double>> imat(p, std::vector<double>(p, 0.0));
  std::vector<double> pi(n), d(n), a(n), b(n);
  switch (param->link_code) {
  case 1: // logit
    for (int person = 0; person < n; ++person) {
      double f = param->freq[person];
      double w = param->weight[person];
      double y = param->y[person];
      double r = boost_plogis(eta[person]);
      double v = y*eta[person] + std::log(1-r);
      loglik += f*w*v;
      
      v = param->y[person] - r;
      
      std::vector<double> z(p);
      for (int i=0; i<p; ++i) {
        z[i] = param->z[person][i]; 
      }
      
      for (int i=0; i<p; ++i) {
        score[i] += f*w*v*z[i];
      }
      
      v = boost_dlogis(eta[person]);
      for (int i=0; i<p; ++i) {
        for (int j=0; j<=i; ++j) {
          imat[i][j] += f*w*v*z[i]*z[j];
        }
      }
      
      if (firth) {
        pi[person] = r;
        d[person] = 1;
        a[person] = r*(1-r);
        b[person] = 1-2*r;
      }
    }
    break;
  case 2: // probit
    for (int person = 0; person < n; ++person) {
      double f = param->freq[person];
      double w = param->weight[person];
      double y = param->y[person];
      double r = boost_pnorm(eta[person]);
      double v = y*std::log(r/(1-r)) + std::log(1-r);
      loglik += f*w*v;
      
      double phi = boost_dnorm(eta[person]);
      double d0 = phi/(r*(1-r));
      v = param->y[person] - r;
      
      std::vector<double> z(p);
      for (int i=0; i<p; ++i) {
        z[i] = param->z[person][i]; 
      }
      
      for (int i=0; i<p; ++i) {
        score[i] += f*w*v*d0*z[i];
      }
      
      v = phi*phi/(r*(1-r));
      for (int i=0; i<p; ++i) {
        for (int j=0; j<=i; ++j) {
          imat[i][j] += f*w*v*z[i]*z[j];
        }
      }
      
      if (firth) {
        double dphi = -eta[person];
        pi[person] = r;
        d[person] = phi/(r*(1-r));
        a[person] = phi*phi/(r*(1-r));
        b[person] = (2*r-1)*phi/(r*(1-r)) + 2*dphi;
      }
    }
    break;
  case 3: // cloglog
    for (int person = 0; person < n; ++person) {
      double f = param->freq[person];
      double w = param->weight[person];
      double y = param->y[person];
      double r = boost_pextreme(eta[person]);
      double v = y*std::log(r/(1-r)) + std::log(1-r);
      loglik += f*w*v;
      
      double phi = boost_dextreme(eta[person]);
      double d0 = phi/(r*(1-r));
      v = param->y[person] - r;
      
      std::vector<double> z(p);
      for (int i=0; i<p; ++i) {
        z[i] = param->z[person][i]; 
      }
      
      for (int i=0; i<p; ++i) {
        score[i] += f*w*v*d0*z[i];
      }
      
      v = phi*phi/(r*(1-r));
      for (int i=0; i<p; ++i) {
        for (int j=0; j<=i; ++j) {
          imat[i][j] += f*w*v*z[i]*z[j];
        }
      }
      
      if (firth) {
        double dphi = 1 - std::exp(eta[person]);
        pi[person] = r;
        d[person] = phi/(r*(1-r));
        a[person] = phi*phi/(r*(1-r));
        b[person] = (2*r-1)*phi/(r*(1-r)) + 2*dphi;
      }
    }
    break;
  }
  
  // fill in the upper triangle of imat
  for (int i=0; i<p-1; ++i) {
    for (int j=i+1; j<p; ++j) {
      imat[i][j] = imat[j][i];
    }
  }
  
  
  if (firth) {
    
    // obtain the determinant of information matrix
    std::vector<std::vector<double>> imat0 = imat;
    cholesky2(imat0, p);
    
    double v = 0;
    for (int i=0; i<p; ++i) {
      v += std::log(imat0[i][i]);
    }
    
    // penalized log-likelihood adjustment
    double penloglik = loglik + 0.5*v;
    
    // compute the bias adjustment to the score function
    std::vector<std::vector<double>> xwx(p, std::vector<double>(p, 0.0));
    for (int person = 0; person < n; ++person) {
      double f = param->freq[person];
      double w = param->weight[person];
      
      std::vector<double> z(p);
      for (int i=0; i<p; ++i) {
        z[i] = param->z[person][i]; 
      }
      
      for (int i=0; i<p; ++i) {
        for (int j=0; j<=i; ++j) {
          xwx[i][j] += f*w*a[person]*z[i]*z[j];
        }
      }
    }
    
    // fill in the upper triangle of xwx
    for (int i=0; i<p-1; ++i) {
      for (int j=i+1; j<p; ++j) {
        xwx[i][j] = xwx[j][i];
      }
    }
    
    std::vector<std::vector<double>> var = invsympd(xwx, p);
    std::vector<double> g(p, 0.0);
    for (int person = 0; person < n; ++person) {
      double f = param->freq[person];
      double w = param->weight[person];
      
      std::vector<double> z(p);
      for (int i=0; i<p; ++i) {
        z[i] = param->z[person][i]; 
      }
      
      double h = 0; // diagonal of H = W^{1/2}*X*inverse(X^T W X)*X^T W^{1/2}
      for (int i=0; i<p; ++i) {
        for (int j=0; j<p; ++j) {
          h += var[i][j]*z[i]*z[j];
        }
      }
      h *= f*w*a[person];
      
      double resid = param->y[person] - pi[person];
      double u = f*w*resid*d[person] + 0.5*b[person]*h;
      for (int i=0; i<p; ++i) {
        g[i] += u*z[i];
      }
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
    result.push_back(imat, "imat");
    return result;
  }
}



// score residual matrix (without firth, weight and frequency)
std::vector<std::vector<double>> f_ressco_0(
    int p, const std::vector<double>& par, void *ex) {
  logparams *param = (logparams *) ex;
  int n = param->n;
  
  std::vector<double> eta(n);
  for (int person = 0; person < n; ++person) {
    eta[person] = param->offset[person];
    for (int i=0; i<p; ++i) {
      eta[person] += par[i]*param->z[person][i];
    }
  }
  
  std::vector<std::vector<double>> resid(n, std::vector<double>(p));
  switch (param->link_code) {
  case 1: // logit
    for (int person = 0; person < n; ++person) {
      double r = boost_plogis(eta[person]);
      double v = param->y[person] - r;
      for (int i=0; i<p; ++i) {
        resid[person][i] = v*param->z[person][i];
      }
    }
    break;
  case 2: // probit
    for (int person = 0; person < n; ++person) {
      double r = boost_pnorm(eta[person]);
      double phi = boost_dnorm(eta[person]);
      double d = phi/(r*(1-r));
      double v = param->y[person] - r;
      for (int i=0; i<p; ++i) {
        resid[person][i] = v*d*param->z[person][i];
      }
    }
    break;
  case 3: // cloglog
    for (int person = 0; person < n; ++person) {
      double r = boost_pextreme(eta[person]);
      double phi = boost_dextreme(eta[person]);
      double d = phi/(r*(1-r));
      double v = param->y[person] - r;
      for (int i=0; i<p; ++i) {
        resid[person][i] = v*d*param->z[person][i];
      }
    }
    break;
  }
  
  return resid;
}


// underlying optimization algorithm for logisreg
//   colfit: vector of indices of parameters to update
//   ncolfit: number of parameters to update
ListCpp logisregloop(int p, const std::vector<double>& par, void *ex,
                     int maxiter, double eps, bool firth,
                     const std::vector<int>& colfit, int ncolfit) {
  logparams *param = (logparams *) ex;
  
  int iter, halving = 0;
  bool fail = false;
  
  std::vector<double> beta(p), newbeta(p);
  double loglik, newlk = 0.0;
  std::vector<double> u(p);
  std::vector<std::vector<double>> imat(p, std::vector<double>(p));
  std::vector<double> u1(ncolfit);
  std::vector<std::vector<double>> imat1(ncolfit, std::vector<double>(ncolfit));
  
  // --- first step ---
  beta = par;
  
  ListCpp der = f_der_0(p, beta, param, firth);
  loglik = der.get<double>("loglik");
  u = der.get<std::vector<double>>("score");
  imat = der.get<std::vector<std::vector<double>>>("imat");
  
  for (int i = 0; i < ncolfit; ++i) u1[i] = u[colfit[i]];
  
  for (int i = 0; i < ncolfit; ++i)
    for (int j = 0; j < ncolfit; ++j)
      imat1[i][j] = imat[colfit[i]][colfit[j]];
  
  cholesky2(imat1, ncolfit);
  chsolve2(imat1, ncolfit, u1);
  
  std::fill(u.begin(), u.end(), 0.0);
  for (int i=0; i<ncolfit; ++i) u[colfit[i]] = u1[i];
  for (int i=0; i<p; ++i) newbeta[i] = beta[i] + u[i];
  
  // --- main iteration ---
  for (iter = 0; iter < maxiter; ++iter) {
    der = f_der_0(p, newbeta, param, firth);
    newlk = der.get<double>("loglik");
    
    fail = std::isnan(newlk) || std::isinf(newlk);
    if (!fail && halving == 0 && std::fabs(1 - (loglik/newlk)) < eps) break;
    
    if (fail || newlk < loglik) {
      ++halving; // adjust step size if likelihood decreases
      for (int i=0; i<p; ++i) {
        newbeta[i] = 0.5*(beta[i] + newbeta[i]);
      }
      continue;
    }
    
    // --- update ---
    halving = 0;
    beta = newbeta;
    loglik = newlk;
    u = der.get<std::vector<double>>("score");
    imat = der.get<std::vector<std::vector<double>>>("imat");
    
    for (int i=0; i<ncolfit; ++i) u1[i] = u[colfit[i]];
    
    for (int i = 0; i < ncolfit; ++i)
      for (int j = 0; j < ncolfit; ++j)
        imat1[i][j] = imat[colfit[i]][colfit[j]];
    
    cholesky2(imat1, ncolfit);
    chsolve2(imat1, ncolfit, u1);
    
    std::fill(u.begin(), u.end(), 0.0);
    for (int i=0; i<ncolfit; ++i) u[colfit[i]] = u1[i];
    for (int i=0; i<p; ++i) newbeta[i] = beta[i] + u[i];
  }
  
  if (iter == maxiter) fail = true;
  
  imat = der.get<std::vector<std::vector<double>>>("imat");
  for (int i = 0; i < ncolfit; ++i)
    for (int j = 0; j < ncolfit; ++j)
      imat1[i][j] = imat[colfit[i]][colfit[j]];
  
  std::vector<std::vector<double>> var1 = invsympd(imat1, ncolfit);
  std::vector<std::vector<double>> var(p, std::vector<double>(p));
  for (int i = 0; i < ncolfit; ++i)
    for (int j = 0; j < ncolfit; ++j)
      var[colfit[i]][colfit[j]] = var1[i][j];
  
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


// confidence limit of profile likelihood method
//   k: index of the parameter to obtain PL CI
//   direction: direction (-1 for lower limit or 1 for upper limit)
//   l0: boundary of profile loglik = max loglik - 0.5*qchisq(1-alpha, 1)
// refer to SAS PROC LOGISTIC documentation for Likelihood Ratio-Based
// Confidence Intervals for Parameters
double logisregplloop(int p, const std::vector<double>& par, 
                      void *ex, int maxiter, double eps, bool firth,
                      int k, int direction, double l0) {
  logparams *param = (logparams *) ex;
  
  int iter;
  bool fail = false;
  
  std::vector<double> beta(p), newbeta(p);
  double loglik, newlk;
  std::vector<double> u(p);
  std::vector<double> delta(p);
  std::vector<std::vector<double>> imat(p, std::vector<double>(p));
  std::vector<std::vector<double>> v(p, std::vector<double>(p));
  
  // --- first step ---
  beta = par;
  
  ListCpp der = f_der_0(p, beta, param, firth);
  loglik = der.get<double>("loglik");
  u = der.get<std::vector<double>>("score");
  imat = der.get<std::vector<std::vector<double>>>("imat");
  
  // Lagrange multiplier method as used in SAS PROC LOGISTIC
  v = invsympd(imat, p);
  
  double w = 0.0;
  for (int i = 0; i < p; ++i)
    for (int j = 0; j < p; ++j)
      w -= u[i] * v[i][j] * u[j];
  
  double underroot = -2*(l0 - loglik + 0.5*w)/v[k][k];
  double lambda = underroot < 0.0 ? 0.0 : direction*std::sqrt(underroot);
  u[k] += lambda;
  
  std::fill(delta.begin(), delta.end(), 0.0);
  for (int i = 0; i < p; ++i)
    for (int j = 0; j < p; ++j)
      delta[i] += v[i][j] * u[j];
  
  // update beta
  for (int i = 0; i < p; ++i) newbeta[i] = beta[i] + delta[i];
  
  // --- main iteration ---
  for (iter = 0; iter < maxiter; ++iter) {
    der = f_der_0(p, newbeta, param, firth);
    newlk = der.get<double>("loglik");
    
    fail = std::isnan(newlk) || std::isinf(newlk);
    if (!fail && std::fabs(newlk - l0) < eps && w < eps) break;
    
    // --- update ---
    beta = newbeta;
    loglik = newlk;
    u = der.get<std::vector<double>>("score");
    imat = der.get<std::vector<std::vector<double>>>("imat");
    
    // Lagrange multiplier method as used in SAS PROC LOGISTIC
    v = invsympd(imat, p);
    
    double w = 0.0;
    for (int i = 0; i < p; ++i)
      for (int j = 0; j < p; ++j)
        w -= u[i] * v[i][j] * u[j];
    
    double underroot = -2*(l0 - loglik + 0.5*w)/v[k][k];
    double lambda = underroot < 0.0 ? 0.0 : direction*std::sqrt(underroot);
    u[k] += lambda;
    
    std::fill(delta.begin(), delta.end(), 0.0);
    for (int i = 0; i < p; ++i)
      for (int j = 0; j < p; ++j)
        delta[i] += v[i][j] * u[j];
    
    // update beta
    for (int i = 0; i < p; ++i)
      newbeta[i] = beta[i] + delta[i];
  }
  
  if (iter == maxiter) fail = true;
  if (fail) {
    thread_utils::push_thread_warning("logisregplloop did not converge within the maximum iterations.");
    return NAN;
  }
  
  return newbeta[k];
}


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
  if (p == 2 && covariates[0] == "") {
    p = 1;
  }
  
  // process the event variable
  if (event.empty()) {
    throw std::invalid_argument("event variable is not specified");
  } else if (!data.containElementNamed(event)) {
    throw std::invalid_argument("data must contain the event variable");
  }
  
  std::vector<double> eventn(n);
  if (data.bool_cols.count(event) > 0) {
    std::vector<bool> eventb = data.get<bool>(event);
    for (int i=0; i<n; ++i) {
      eventn[i] = eventb[i] ? 1.0 : 0.0;
    }
  } else if (data.int_cols.count(event) > 0) {
    std::vector<int> eventi = data.get<int>(event);
    for (int i=0; i<n; ++i) {
      eventn[i] = static_cast<double>(eventi[i]);
    }
  } else if (data.numeric_cols.count(event) > 0) {
    eventn = data.get<double>(event);
  } else {
    throw std::invalid_argument("event variable must be bool, integer or numeric");
  }
  
  for (double val : eventn) {
    if (val != 1 && val != 0) {
      throw std::invalid_argument("event must be 1 or 0 for each observation");
    }
  }
  
  // construct the design matrix
  std::vector<std::vector<double>> zn(n, std::vector<double>(p));
  for (int i=0; i<n; ++i) {
    zn[i][0] = 1.0; // intercept
  }
  for (int j=0; j<p-1; ++j) {
    std::string zj = covariates[j];
    if (!data.containElementNamed(zj)) {
      throw std::invalid_argument("data must contain the variables in covariates");
    }
    
    std::vector<double> u(n);
    if (data.bool_cols.count(zj) > 0) {
      std::vector<bool> ub = data.get<bool>(zj);
      for (int i=0; i<n; ++i) {
        u[i] = ub[i] ? 1.0 : 0.0;
      }
    } else if (data.int_cols.count(zj) > 0) {
      std::vector<int> ui = data.get<int>(zj);
      for (int i=0; i<n; ++i) {
        u[i] = static_cast<double>(ui[i]);
      }
    } else if (data.numeric_cols.count(zj) > 0) {
      u = data.get<double>(zj);
    } else {
      throw std::invalid_argument("covarates must be bool, integer or numeric");
    }
    
    for (int i=0; i<n; ++i) {
      zn[i][j+1] = u[i];
    }
  }
  
  // process freq, weight, offset variables
  std::vector<double> freqn(n, 1.0);
  if (!freq.empty() && data.containElementNamed(freq)) {
    if (data.int_cols.count(freq) > 0) {
      std::vector<int> freqi = data.get<int>(freq);
      for (int i=0; i<n; ++i) {
        freqn[i] = static_cast<double>(freqi[i]);
      }
    } else if (data.numeric_cols.count(freq) > 0) {
      freqn = data.get<double>(freq);
    } else {
      throw std::invalid_argument("freq variable must be integer or numeric");
    }
    
    for (double val : freqn) {
      if (val <= 0) {
        throw std::invalid_argument("freq must be positive integers");
      }
    }
  }
  
  std::vector<double> weightn(n, 1.0);
  if (!weight.empty() && data.containElementNamed(weight)) {
    std::vector<double> weightn = data.get<double>(weight);
    for (double val : weightn) {
      if (val <= 0.0) {
        throw std::invalid_argument("weight must be greater than 0");
      }
    }
  }
  
  std::vector<double> offsetn(n, 0.0);
  if (!offset.empty() && data.containElementNamed(offset)) {
    offsetn = data.get<double>(offset);
  }
  
  // create the numeric id variable
  bool has_id = !id.empty() && data.containElementNamed(id);
  std::vector<int> idn(n);
  if (!has_id) {
    std::iota(idn.begin(), idn.end(), 0);
  } else {
    if (data.int_cols.count(id)) {
      std::vector<int> v = data.get<int>(id);
      std::vector<int> w = unique_sorted(v);
      idn = matchcpp(v, w);
    } else if (data.numeric_cols.count(id)) {
      std::vector<double> v = data.get<double>(id);
      std::vector<double> w = unique_sorted(v);
      idn = matchcpp(v, w);
    } else if (data.string_cols.count(id)) {
      std::vector<std::string> v = data.get<std::string>(id);
      std::vector<std::string> w = unique_sorted(v);
      idn = matchcpp(v, w);
    } else {
      throw std::invalid_argument("incorrect type for the id variable in the input data");
    }
  }
  
  // process link function
  std::string link1 = link;
  std::for_each(link1.begin(), link1.end(), [](char & c) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  });
  
  if (link1 == "log-log" || link1 == "loglog" || link1 == "cloglog") {
    link1 = "cloglog";
  }
  
  int link_code;
  if (link1 == "logit") link_code = 1;
  else if (link1 == "probit") link_code = 2;
  else if (link1 == "cloglog") link_code = 3;
  else throw std::invalid_argument("invalid link: " + link1);
  
  // exclude observations with missing values
  std::vector<bool> sub(n,1);
  for (int i=0; i<n; ++i) {
    if (eventn[i] == INT_MIN ||
        std::isnan(freqn[i]) || std::isnan(weightn[i]) ||
        std::isnan(offsetn[i]) || idn[i] == INT_MIN) {
      sub[i] = 0;
    }
    for (int j=0; j<p-1; ++j) {
      if (std::isnan(zn[i][j+1])) sub[i] = 0;
    }
  }
  
  std::vector<int> order = which(sub);
  subset_in_place(eventn, order);
  subset_in_place(freqn, order);
  subset_in_place(weightn, order);
  subset_in_place(offsetn, order);
  subset_in_place(idn, order);
  subset_in_place(zn, order);
  n = std::accumulate(sub.begin(), sub.end(), 0);
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
  std::vector<std::vector<double>> vb(p, std::vector<double>(p));
  std::vector<std::vector<double>> rvb(p, std::vector<double>(p));
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
      if (i==0) {
        par[i] = "(Intercept)";
      } else {
        par[i] = covariates[i-1];
      }
      
      b[i] = NAN;
      seb[i] = 0;
      rseb[i] = 0;
      z[i] = NAN;
      expbeta[i] = NAN;
      for (int j=0; j<p; ++j) {
        vb[i][j] = 0;
        rvb[i][j] = 0;
      }
      lb[i] = NAN;
      ub[i] = NAN;
      prob[i] = NAN;
      clparm[i] = "Wald";
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
    
    loglik0 = NAN;
    loglik1 = NAN;
    regloglik0 = NAN;
    regloglik1 = NAN;
    niter = 0;
    fail = true;
  } else {
    // intercept only model
    double num = 0, den = 0;
    for (int i=0; i<n; ++i) {
      num += freqn[i]*weightn[i]*eventn[i];
      den += freqn[i]*weightn[i];
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
    std::vector<std::vector<double>> vbint = outint.get<std::vector<std::vector<double>>>("var");
    
    ListCpp out;
    
    if (p > 1) {
      // parameter estimates and standard errors for the full model
      std::vector<int> colfit = seqcpp(0,p-1);
      if (static_cast<int>(init.size()) == p && std::none_of(init.begin(), init.end(), [](double val){ return std::isnan(val); })) {
        out = logisregloop(p, init, &param, maxiter, eps, firth, colfit, p);
      } else {
        out = logisregloop(p, bint, &param, maxiter, eps, firth, colfit, p);
      }
      
      bool fail = out.get<bool>("fail");
      if (fail) {
        thread_utils::push_thread_warning("logisregloop failed to converge for the full model; continuing with current results (fail=TRUE).");
      }
      
      b = out.get<std::vector<double>>("coef");
      vb = out.get<std::vector<std::vector<double>>>("var");
      
      // intercept correction
      if (flic) {
        std::vector<double> lp(n);  // linear predictor excluding intercept
        for (int person = 0; person < n; ++person) {
          lp[person] = offsetn[person];
          for (int i=1; i<p; ++i) {
            lp[person] += b[i]*zn[person][i];
          }
        }
        
        logparams param0 = {n, link_code, eventn, zn, freqn, weightn, lp};
        std::vector<double> bint00(1, bint0[0]);
        ListCpp outint0 = logisregloop(1, bint00, &param0, maxiter, eps, 0, colfit0, 1);
        double a = outint0.get<std::vector<double>>("coef")[0];
        double va = outint0.get<std::vector<std::vector<double>>>("var")[0][0];
        
        // update the intercept estimate
        b[0] = a;
        
        // partial derivative of alpha(beta) with respect to beta
        ListCpp derint = f_der_0(p, b, &param, 0);
        std::vector<std::vector<double>> iflic = derint.get<std::vector<std::vector<double>>>("imat");
        std::vector<double> der(p-1);
        for (int i=0; i<p-1; ++i) {
          der[i] = -iflic[i+1][0]/iflic[0][0];
        }
        
        // update the variance of alpha
        vb[0][0] = va;
        for (int i=0; i<p-1; ++i) {
          for (int j=0; j<p-1; ++j) {
            vb[0][0] += der[i]*vb[i+1][j+1]*der[j];
          }
        }
        
        // update the covariance between alpha and beta
        for (int i=0; i<p-1; ++i) {
          vb[i+1][0] = 0;
          for (int j=0; j<p-1; ++j) {
            vb[i+1][0] += vb[i+1][j+1]*der[j];
          }
          vb[0][i+1] = vb[i+1][0];
        }
      }
    } else {
      b = bint;
      vb = vbint;
      out = outint;
      
      if (flic) {
        out = logisregloop(p, bint0, &param, maxiter, eps, 0, colfit0, 1);
        b = out.get<std::vector<double>>("coef");
        vb = outint.get<std::vector<std::vector<double>>>("var");
      }
    }
    
    for (int j=0; j<p; ++j) {
      seb[j] = std::sqrt(vb[j][j]);
    }
    
    for (int i=0; i<p; ++i) {
      if (i==0) {
        par[i] = "(Intercept)";
      } else {
        par[i] = covariates[i-1];
      }
    }
    
    // linear predictors and fitted values
    std::vector<double> eta(n);
    for (int person = 0; person < n; ++person) {
      eta[person] = offsetn[person];
      for (int i=0; i<p; ++i) {
        eta[person] += b[i]*zn[person][i];
      }
      linear_predictors[person] = eta[person];
    }
    
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
      std::vector<std::vector<double>> ressco = f_ressco_0(p, b, &param);
      
      int nr; // number of rows in the score residual matrix
      std::vector<double> freqr;
      if (!has_id) {
        for (int i=0; i<n; ++i) {
          for (int j=0; j<p; ++j) {
            ressco[i][j] = weightn[i]*ressco[i][j];
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
        
        std::vector<std::vector<double>> ressco2(nids, std::vector<double>(p, 0.0));
        for (int i=0; i<nids; ++i) {
          for (int j=0; j<p; ++j) {
            for (int k=idx[i]; k<idx[i+1]; ++k) {
              ressco2[i][j] += weight2[k]*ressco[order[k]][j];
            }
          }
          freqr0[i] = freq2[idx[i]];
        }
        
        ressco = ressco2;  // update the score residuals
        nr = nids;
        freqr = freqr0;
      }
      
      std::vector<std::vector<double>> D(nr, std::vector<double>(p, 0.0));
      for (int i=0; i<nr; ++i) {
        for (int j=0; j<p; ++j) {
          for (int k=0; k<p; ++k) {
            D[i][j] += ressco[i][k]*vb[k][j];
          }
        }
      }
      
      for (int j=0; j<p; ++j) {
        for (int k=0; k<p; ++k) {
          for (int i=0; i<nr; ++i) {
            rvb[j][k] += freqr[i]*D[i][j]*D[i][k];
          }
        }
      }
      
      for (int i=0; i<p; ++i) {
        rseb[i] = std::sqrt(rvb[i][i]);
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
        z[i] = NAN;
      } else {
        z[i] = b[i]/rseb[i];
      }
    }
  } else {
    for (int i=0; i<p; ++i) {
      if (seb[i] == 0) {
        z[i] = NAN;
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




// Worker that runs logisregcpp on a range of DataFrameCpp inputs.
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
    : data_ptr(data_ptr_), event(event_), covariates(covariates_),
      freq(freq_), weight(weight_), offset(offset_), id(id_), link(link_),
      init(init_), robust(robust_), firth(firth_), flic(flic_), plci(plci_),
      alpha(alpha_), maxiter(maxiter_), eps(eps_), results(results_) {}
  
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