// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>    // RcppParallel::Worker, parallelFor
#include <RcppThread.h>      // RcppThread::Rcerr
#include <Rcpp.h>

#include "survival_analysis.h"
#include "utilities.h"      // boost_pnorm, boost_plogis, etc.
#include "dataframe_list.h"  // FlatMatrix, IntMatrix, DataFrameCpp, ListCpp
#include "thread_utils.h"    // push_thread_warning / drain_thread_warnings_to_R

#include <vector>
#include <string>
#include <numeric>   // iota, inner_product
#include <cmath>     // isnan, isinf, fabs, NaN, exp, log
#include <stdexcept> // exceptions
#include <algorithm> // sort, none_of, any_of
#include <cstring>   // memcpy

struct coxparams {
  int nused;
  std::vector<int> strata;
  std::vector<double> tstart;
  std::vector<double> tstop;
  std::vector<double> event;
  std::vector<double> weight;
  std::vector<double> offset;
  FlatMatrix z;
  std::vector<int> order1;
  int method; // 1: breslow, 2: efron
};

// all-in-one function for log-likelihood, score, and information matrix
// for the Cox model with or without Firth's correction
ListCpp f_der_2(int p, const std::vector<double>& par, void* ex, bool firth) {
  coxparams* param = (coxparams*) ex;
  const int n = param->z.nrow;
  const int nused = param->nused;
  const int method = param->method;
  
  const std::vector<int>& strata = param->strata; 
  const std::vector<double>& tstart = param->tstart; 
  const std::vector<double>& tstop = param->tstop; 
  const std::vector<double>& event = param->event; 
  const std::vector<double>& weight = param->weight; 
  const std::vector<double>& offset = param->offset; 
  const std::vector<int>& order1 = param->order1;
  const FlatMatrix& z = param->z;

  // Precompute eta and exp(eta)
  std::vector<double> eta(nused);
  for (int person = 0; person < nused; ++person) {
    eta[person] = offset[person];
  }
  for (int i = 0; i < p; ++i) {
    double beta = par[i];
    if (beta == 0.0) continue;
    const int off = i * n;
    for (int person = 0; person < nused; ++person) {
      eta[person] += beta * z.data[off + person];
    }
  }
  std::vector<double> exp_eta(nused);
  for (int person = 0; person < nused; ++person) {
    exp_eta[person] = std::exp(eta[person]);
  }
  
  double loglik = 0.0;        // log-likelihood
  std::vector<double> u(p);   // score vector
  FlatMatrix imat(p,p);       // information matrix
  FlatArray dimat(p,p,p);     // tensor for third order derivatives
  std::vector<double> a(p);   // s1(beta,k,t)
  std::vector<double> a2(p);  // sum of w*exp(zbeta)*z for the deaths
  FlatMatrix cmat(p,p);       // s2(beta,k,t)
  FlatMatrix cmat2(p,p);      // sum of w*exp(zbeta)*z*z' for the deaths
  FlatArray dmat(p,p,p);      // q2(beta,k,t)
  FlatArray dmat2(p,p,p);     // sum of w*exp(zbeta)*z*z*z' for the deaths
  double denom = 0.0;         // s0(beta,k,t)
  double denom2 = 0.0;        // sum of weighted risks for deaths
  double deadwt = 0.0;        // sum of weights for the deaths
  double ndead = 0.0;         // number of deaths at this time point
  
  int istrata = strata[0];
  int i1 = 0; // index for removing out-of-risk subjects
  
  // Loop through subjects
  for (int person = 0; person < nused; ) {
    // Reset when entering a new stratum
    if (strata[person] != istrata) {
      istrata = strata[person];
      i1 = person;
      denom = 0.0;
      
      for (int i = 0; i < p; ++i) {
        a[i] = 0.0;
        for (int j = 0; j <= i; ++j) {
          cmat(i,j) = 0.0;
          if (firth) {
            for (int k = 0; k <= j; ++k) {
              dmat(i,j,k) = 0.0;
            }
          }
        }
      }
    }
    
    const double dtime = tstop[person];
    
    // Process all persons tied at this dtime
    for (; person < nused && tstop[person] == dtime && strata[person] == istrata; ++person) {
      
      const double w = weight[person];
      const double r = w * exp_eta[person];
      
      if (event[person] == 0) {
        denom += r;
        
        for (int i = 0; i < p; ++i) {
          const double zi = z(person,i);
          a[i] += r * zi;
          for (int j = 0; j <= i; ++j) {
            const double zj = z(person,j);
            cmat(i,j) += r * zi * zj;
            if (firth) {
              for (int k = 0; k <= j; ++k) {
                dmat(i,j,k) += r * zi * zj * z(person,k);
              }
            }
          }
        }
      } else {
        ++ndead;
        deadwt += w;
        denom2 += r;
        loglik += w * eta[person];
        for (int i = 0; i < p; ++i) {
          const double zi = z(person,i);
          a2[i] += r * zi;
          u[i] += w * zi;
          for (int j = 0; j <= i; ++j) {
            const double zj = z(person,j);
            cmat2(i,j) += r * zi * zj;
            if (firth) {
              for (int k = 0; k <= j; ++k) {
                dmat2(i,j,k) += r * zi * zj * z(person,k);
              }
            }
          }
        }
      }
    }
    
    // Remove subjects leaving risk set
    for (; i1 < nused; ++i1) {
      const int p1 = order1[i1];
      if (tstart[p1] < dtime || strata[p1] != istrata) break;
      
      const double r = weight[p1] * exp_eta[p1];
      denom -= r;
      for (int i = 0; i < p; ++i) {
        const double zi = z(p1,i);
        a[i] -= r * zi;
        for (int j = 0; j <= i; ++j) {
          const double zj = z(p1,j);
          cmat(i,j) -= r * zi * zj;
          if (firth) {
            for (int k = 0; k <= j; ++k) {
              dmat(i,j,k) -= r * zi * zj * z(p1,k);
            }
          }
        }
      }
    }
    
    // Add contributions for deaths at this time
    if (ndead > 0) {
      if (method == 0 || ndead == 1) { // Breslow or single event
        denom += denom2;
        loglik -= deadwt * std::log(denom);
        for (int i = 0; i < p; ++i) {
          a[i] += a2[i];
          const double xbar = a[i] / denom;
          u[i] -= deadwt * xbar;
          for (int j = 0; j <= i; ++j) {
            cmat(i,j) += cmat2(i,j);
            imat(i,j) += deadwt * (cmat(i,j) - xbar * a[j]) / denom;
            if (firth) {
              for (int k = 0; k <= j; ++k) {
                dmat(i,j,k) += dmat2(i,j,k);
                dimat(i,j,k) += deadwt * (dmat(i,j,k) -
                  (cmat(i,j)*a[k] + cmat(i,k)*a[j] + cmat(j,k)*a[i]) / denom +
                  2.0 * a[i] * a[j] * a[k] / (denom * denom)) / denom;
              }
            }
          }
        }
      } else { // Efron method
        const double meanwt = deadwt / ndead;
        const double increment = denom2 / ndead;
        for (int l = 0; l < ndead; ++l) {
          denom += increment;
          loglik -= meanwt * std::log(denom);
          for (int i = 0; i < p; ++i) {
            a[i] += a2[i] / ndead;
            const double xbar = a[i] / denom;
            u[i] -= meanwt * xbar;
            for (int j = 0; j <= i; ++j) {
              cmat(i,j) += cmat2(i,j) / ndead;
              imat(i,j) += meanwt * (cmat(i,j) - xbar * a[j]) / denom;
              if (firth) {
                for (int k = 0; k <= j; ++k) {
                  dmat(i,j,k) += dmat2(i,j,k) / ndead;
                  dimat(i,j,k) += meanwt * (dmat(i,j,k) -
                    (cmat(i,j)*a[k] + cmat(i,k)*a[j] + cmat(j,k)*a[i])/denom +
                    2.0 * a[i] * a[j] * a[k] / (denom * denom)) / denom;
                }
              }
            }
          }
        }
      }
      
      // Reset after processing deaths
      ndead = deadwt = denom2 = 0.0;
      for (int i = 0; i < p; ++i) {
        a2[i] = 0;
        for (int j = 0; j <= i; ++j) {
          cmat2(i,j) = 0;
          if (firth) {
            for (int k = 0; k <= j; ++k) {
              dmat2(i,j,k) = 0;
            }
          }
        }
      }
    }
  }
  
  
  // fill the symmetric elements of the information matrix
  for (int i = 0; i < p - 1; ++i)
    for (int j = i+1; j < p; ++j)
      imat(i,j) = imat(j,i);
  
  
  // fill the symmetric elements of the tensor array
  if (firth) {
    for (int i = 0; i < p-1; ++i)
      for (int j = i+1; j < p; ++j)
        for (int k = 0; k <= i; ++k)
          dimat(i,j,k) = dimat(j,i,k);
    
    for (int j = 0; j < p-1; ++j)
      for (int i = j; i < p; ++i)
        for (int k = j+1; k < p; ++k)
          dimat(i,j,k) = dimat(i,k,j);
    
    for (int i = 0; i < p-1; ++i)
      for (int j = i+1; j < p; ++j)
        for (int k = i+1; k < p; ++k)
          dimat(i,j,k) = dimat(k,j,i);
  }
  
  
  ListCpp result;
  
  // Firth adjustment
  if (p > 0 && firth) {
    // obtain the determinant of information matrix
    FlatMatrix imat0 = imat;
    cholesky2(imat0, p);
    
    double v = 0.0;
    for (int i=0; i<p; ++i) {
      v += std::log(imat0(i,i));
    }
    
    // penalized log-likelihood adjustment
    double penloglik = loglik + 0.5*v;
    
    // compute the bias adjustment to the score function
    FlatMatrix y(p,p);
    std::vector<double> g(p);
    
    for (int k = 0; k < p; ++k) {
      // partial derivative of the information matrix w.r.t. beta[k]
      for (int i = 0; i < p; ++i) {
        for (int j = 0; j < p; ++j) {
          y(i,j) = dimat(i,j,k);
        }
      }
      
      // solve(imat, y)
      for (int h = 0; h < p; ++h) {
        for (int i = 0; i < p; ++i) {
          double temp = y(i,h);
          for (int j = 0; j < i; ++j)
            temp -= y(j,h) * imat0(j,i);
          y(i,h) = temp;
        }
        
        for (int i = p-1; i >= 0; --i) {
          if (imat0(i,i) == 0) y(i,h) = 0;
          else {
            double temp = y(i,h) / imat0(i,i);
            for (int j = i+1; j < p; ++j)
              temp -= y(j,h)*imat0(i,j);
            y(i,h) = temp;
          }
        }
      }
      
      // trace
      for (int i = 0; i < p; ++i) g[k] += y(i,i);
      
      g[k] = u[k] + 0.5 * g[k];
    }
    
    result.push_back(penloglik, "loglik");
    result.push_back(std::move(g), "score");
    result.push_back(std::move(imat), "imat");
    result.push_back(loglik, "regloglik");
    result.push_back(std::move(u), "regscore");
  } else {
    result.push_back(loglik, "loglik");
    if (p > 0) {
      result.push_back(std::move(u), "score");
      result.push_back(std::move(imat), "imat");
    }
  }
  
  return result;
}


// underlying optimization algorithm for phreg
ListCpp phregloop(int p, const std::vector<double>& par, void *ex,
               int maxiter, double eps, bool firth,
               const std::vector<int>& colfit, int ncolfit) {
  coxparams *param = (coxparams *) ex;

  int iter = 0, halving = 0;
  bool fail = false;

  std::vector<double> beta = par;
  std::vector<double> newbeta(p);
  double loglik = 0.0, newlk = 0.0;
  std::vector<double> u(p);
  std::vector<double> u1(ncolfit);
  FlatMatrix imat(p,p);
  FlatMatrix imat1(ncolfit, ncolfit);
  FlatMatrix z1 = param->z;

  // --- first step ---
  ListCpp der = f_der_2(p, beta, param, firth);
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
    der = f_der_2(p, newbeta, param, firth);
    newlk = der.get<double>("loglik");
    
    fail = std::isnan(newlk) || std::isinf(newlk);
    if (!fail && halving == 0 && std::fabs(1 - (loglik / newlk)) < eps) break;

    if (fail || newlk < loglik) {
      ++halving; // adjust step size if likelihood decreases
      for (int i = 0; i < p; ++i) {
        newbeta[i] = 0.5 * (beta[i] + newbeta[i]);
      }
      continue;
    }

    // --- update ---
    halving = 0;
    beta = newbeta;
    loglik = newlk;
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
  }

  if (iter == maxiter) fail = true;

  // --- final variance calculation ---
  imat = der.get<FlatMatrix>("imat");
  for (int i = 0; i < ncolfit; ++i)
    for (int j = 0; j < ncolfit; ++j)
      imat1(i,j) = imat(colfit[i], colfit[j]);

  FlatMatrix var1 = invsympd(imat1, ncolfit);
  FlatMatrix var(p, p);
  for (int i = 0; i < ncolfit; ++i)
    for (int j = 0; j < ncolfit; ++j)
      var(colfit[i], colfit[j]) = var1(i,j);

  ListCpp result;
  result.push_back(std::move(newbeta), "coef");
  result.push_back(iter, "iter");
  result.push_back(std::move(var), "var");
  result.push_back(newlk, "loglik");
  result.push_back(fail, "fail");
  
  if (firth) {
    double regloglik = der.get<double>("regloglik");
    result.push_back(regloglik, "regloglik");
  }
  
  return result;
}


// confidence limit of profile likelihood method
double phregplloop(int p, const std::vector<double>& par, void *ex,
                   int maxiter, double eps, bool firth,
                   int k, int direction, double l0) {
  coxparams *param = (coxparams *) ex;
  int iter;
  bool fail = false;

  std::vector<double> beta = par;
  std::vector<double> newbeta(p);
  double loglik = 0.0, newlk = 0.0;
  std::vector<double> u(p);
  std::vector<double> delta(p);
  FlatMatrix imat(p, p);
  FlatMatrix v(p, p);
  
  ListCpp der = f_der_2(p, beta, param, firth);
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
  
  // --- main iteration ---
  for (iter = 0; iter < maxiter; ++iter) {
    der = f_der_2(p, newbeta, param, firth);
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
  if (fail) thread_utils::push_thread_warning("phregplloop did not converge.");

  return newbeta[k];
}


// baseline hazard estimates
ListCpp f_basehaz(int p, const std::vector<double>& par, void *ex) {
  coxparams *param = (coxparams *) ex;
  const int n = param->z.nrow;
  const int nused = param->nused;
  const int method = param->method;
  
  const std::vector<int>& strata = param->strata; 
  const std::vector<double>& tstart = param->tstart; 
  const std::vector<double>& tstop = param->tstop; 
  const std::vector<double>& event = param->event; 
  const std::vector<double>& weight = param->weight; 
  const std::vector<double>& offset = param->offset; 
  const std::vector<int>& order1 = param->order1;
  const FlatMatrix& z = param->z;
  
  // Precompute eta and exp(eta)
  std::vector<double> eta(nused);
  for (int person = 0; person < nused; ++person) {
    eta[person] = offset[person];
  }
  for (int i = 0; i < p; ++i) {
    double beta = par[i];
    if (beta == 0.0) continue;
    const int off = i * n;
    for (int person = 0; person < nused; ++person) {
      eta[person] += beta * z.data[off + person];
    }
  }
  std::vector<double> exp_eta(nused);
  for (int person = 0; person < nused; ++person) {
    exp_eta[person] = std::exp(eta[person]);
  }
  
  std::vector<double> a(p);   // s1(beta,k,t)
  std::vector<double> a2(p);  // sum of w*exp(zbeta)*z for the deaths
  double deadwt = 0.0;        // sum of weights for the deaths
  double denom = 0.0;         // s0(beta,k,t)
  double denom2 = 0.0;        // sum of weighted risks for the deaths
  double natrisk = 0;         // number at risk at this time point
  double ndead = 0;           // number of deaths at this time point
  double ncens = 0;           // number of censored at this time point

  // locate the first observation within each stratum
  std::vector<int> istratum(1,0);
  for (int i = 1; i < nused; ++i) {
    if (strata[i] != strata[i-1]) {
      istratum.push_back(i);
    }
  }

  int nstrata = static_cast<int>(istratum.size());
  istratum.push_back(nused);

  // add time 0 to each stratum
  int J = nstrata;
  for (int i = 0; i < nstrata; ++i) {
    std::vector<int> idx = seqcpp(istratum[i], istratum[i+1] - 1);
    std::vector<double> utime = subset(tstop, idx);
    utime = unique_sorted(utime);
    J += static_cast<int>(utime.size());
  }

  std::vector<int> stratum(J);
  std::vector<double> time(J), nrisk(J), nevent(J), ncensor(J), haz(J), varhaz(J);
  FlatMatrix gradhaz(J,p);

  int istrata = strata[0];
  int i1 = 0; // index for removing out-of-risk subjects
  int j = J;  // index the unique time in ascending order

  // Loop through subjects
  for (int person = 0; person < nused; ) {
    if (strata[person] != istrata) { // hit a new stratum
      // add time 0 at the start of a new stratum
      j--;
      stratum[j] = istrata;
      time[j] = 0.0;
      nrisk[j] = natrisk;

      istrata = strata[person]; // reset temporary variables
      i1 = person;
      natrisk = 0;
      denom = 0.0;
      std::fill(a.begin(), a.end(), 0.0);
    }

    const double dtime = tstop[person];

    // Process all persons tied at this dtime
    bool first = true;
    for (; person < nused && tstop[person] == dtime && strata[person] == istrata; ++person) {

      if (first) { // first incidence at this time
        j--;
        stratum[j] = strata[person];
        time[j] = dtime;
        first = false;
      }

      const double w = weight[person];
      const double r = w * exp_eta[person];

      ++natrisk;
      if (event[person] == 0) {
        ++ncens;
        denom += r;
        for (int i = 0; i < p; ++i) {
          a[i] += r * z(person,i);
        }
      } else {
        ++ndead;
        deadwt += w;
        denom2 += r;
        for (int i = 0; i < p; ++i) {
          a2[i] += r * z(person,i);
        }
      }
    }

    // Remove subjects leaving risk set
    for (; i1 < nused; ++i1) {
      const int p1 = order1[i1];
      if (tstart[p1] < dtime || strata[p1] != istrata) break;

      const double r = weight[p1] * exp_eta[p1];

      natrisk--;
      denom -= r;
      for (int i = 0; i < p; ++i) {
        a[i] -= r * z(p1,i);
      }
    }

    // Add contributions for deaths at this time
    nrisk[j] = natrisk;
    nevent[j] = ndead;
    ncensor[j] = ncens;
    ncens = 0; // reset for the next time point
    if (ndead > 0) {
      if (method == 0 || ndead == 1) { // Breslow or single event
        denom += denom2;
        const double temp = deadwt / denom;
        haz[j] = temp;
        varhaz[j] = temp / denom;
        for (int i = 0; i < p; ++i) {
          a[i] += a2[i];
          gradhaz(j,i) = temp * a[i] / denom;
        }
      } else { // Efron method
        const double meanwt = deadwt / ndead;
        for (int k = 0; k < ndead; ++k) {
          denom += denom2 / ndead;
          const double temp = meanwt / denom;
          haz[j] += temp;
          varhaz[j] += temp / denom;
          for (int i = 0; i < p; ++i) {
            a[i] += a2[i] / ndead;
            gradhaz(j,i) += temp * a[i] / denom;
          }
        }
      }

      // reset for the next death time
      ndead = deadwt = denom2 = 0.0;
      std::fill(a2.begin(), a2.end(), 0.0);
    }
  }

  // add time 0 for the first stratum
  stratum[0] = istrata;
  time[0] = 0.0;
  nrisk[0] = natrisk;

  ListCpp result;
  result.push_back(stratum, "stratum");
  result.push_back(time, "time");
  result.push_back(nrisk, "nrisk");
  result.push_back(nevent, "nevent");
  result.push_back(ncensor, "ncensor");
  result.push_back(haz, "haz");
  result.push_back(varhaz, "varhaz");
  
  if (p > 0) result.push_back(gradhaz, "gradhaz");

  return result;
}


// martingale residuals
std::vector<double> f_resmart(int p, const std::vector<double>& par, void *ex) {
  coxparams *param = (coxparams *) ex;
  const int n = param->z.nrow;
  const int nused = param->nused;
  const int method = param->method;
  
  const std::vector<int>& strata = param->strata; 
  const std::vector<double>& tstart = param->tstart; 
  const std::vector<double>& tstop = param->tstop; 
  const std::vector<double>& event = param->event; 
  const std::vector<double>& weight = param->weight; 
  const std::vector<double>& offset = param->offset; 
  const std::vector<int>& order1 = param->order1;
  const FlatMatrix& z = param->z;
  
  // Precompute eta and exp(eta)
  std::vector<double> eta(nused);
  for (int person = 0; person < nused; ++person) {
    eta[person] = offset[person];
  }
  for (int i = 0; i < p; ++i) {
    double beta = par[i];
    if (beta == 0.0) continue;
    const int off = i * n;
    for (int person = 0; person < nused; ++person) {
      eta[person] += beta * z.data[off + person];
    }
  }
  std::vector<double> exp_eta(nused);
  for (int person = 0; person < nused; ++person) {
    exp_eta[person] = std::exp(eta[person]);
  }
  
  double denom = 0.0;         // s0(beta,k,t)
  double denom2 = 0.0;        // sum of weighted risks for deaths
  double deadwt = 0.0;        // sum of weights for the deaths
  double ndead = 0.0;         // number of deaths at this time point

  // initialize the residuals to the event indicators
  std::vector<double> resid(n);
  for (int person = 0; person < nused; ++person) {
    resid[person] = event[person];
  }

  int istrata = strata[0];
  int i1 = 0; // index for removing out-of-risk subjects
  int j0 = 0; // first person in the stratum

  // Loop through subjects
  for (int person = 0; person < nused; ) {
    if (strata[person] != istrata) { // hit a new stratum
      istrata = strata[person]; // reset temporary variables
      i1 = person;
      j0 = person;
      denom = 0.0;
    }

    const double dtime = tstop[person];

    // process all persons tied at this dtime
    int j1 = person;   // first person in the stratum with the tied time
    for (; person < nused && tstop[person] == dtime && strata[person] == istrata; ++person) {

      const double w = weight[person];
      const double r = w * exp_eta[person];

      if (event[person] == 0) {
        denom += r;
      } else {
        ++ndead;
        deadwt += w;
        denom2 += r;
      }
    }

    int j2 = person - 1; // last person in the stratum with the tied time

    // Remove subjects leaving risk set
    for (; i1 < nused; ++i1) {
      const int p1 = order1[i1];
      if (tstart[p1] < dtime || strata[p1] != istrata) break;
      denom -= weight[p1] * exp_eta[p1];
    }

    // Add contributions for deaths at this time
    if (ndead > 0) {
      denom += denom2;

      for (int j = j0; j <= j2; ++j) {
        if (tstart[j] < dtime) {
          double hazard;
          if (method == 0 || ndead == 1) {
            hazard = deadwt / denom;
          } else {
            hazard = 0.0;
            const double meanwt = deadwt / ndead;
            if (j < j1 || event[j] == 0) {
              for (int i = 0; i < ndead; ++i) {
                hazard += meanwt /(denom - i/ndead * denom2);
              }
            } else {
              for (int i = 0; i < ndead; ++i) {
                hazard += (1 - i/ndead) * meanwt / (denom - i/ndead * denom2);
              }
            }
          }
          resid[j] -= hazard * exp_eta[j];
        }
      }

      // reset for the next death time
      ndead = deadwt = denom2 = 0.0;
    }
  }

  return resid;
}


// score residual matrix
FlatMatrix f_ressco_2(int p, const std::vector<double>& par, void *ex) {
  coxparams *param = (coxparams *) ex;
  const int n = param->z.nrow;
  const int nused = param->nused;
  const int method = param->method;
  
  const std::vector<int>& strata = param->strata; 
  const std::vector<double>& tstart = param->tstart; 
  const std::vector<double>& tstop = param->tstop; 
  const std::vector<double>& event = param->event; 
  const std::vector<double>& weight = param->weight; 
  const std::vector<double>& offset = param->offset; 
  const std::vector<int>& order1 = param->order1;
  const FlatMatrix& z = param->z;
  
  // Precompute eta and exp(eta)
  std::vector<double> eta(nused);
  for (int person = 0; person < nused; ++person) {
    eta[person] = offset[person];
  }
  for (int i = 0; i < p; ++i) {
    double beta = par[i];
    if (beta == 0.0) continue;
    const int off = i * n;
    for (int person = 0; person < nused; ++person) {
      eta[person] += beta * z.data[off + person];
    }
  }
  std::vector<double> exp_eta(nused);
  for (int person = 0; person < nused; ++person) {
    exp_eta[person] = std::exp(eta[person]);
  }
  
  FlatMatrix resid(n,p);      // residual matrix
  std::vector<double> a(p);   // s1(beta,k,t)
  std::vector<double> a2(p);  // sum of w*exp(zbeta)*z for the deaths
  double denom = 0.0;         // s0(beta,k,t)
  double denom2 = 0.0;        // sum of weighted risks for deaths
  double deadwt = 0.0;        // sum of weights for the deaths
  double ndead = 0.0;         // number of deaths at this time point
  double cumhaz = 0.0;        // cumulative hazard

  std::vector<double> xhaz(p), mh1(p), mh2(p), mh3(p); // temp vectors

  int istrata = strata[0];
  int i1 = 0; // index for removing out-of-risk subjects

  // Loop through subjects
  for (int person = 0; person < nused; ) {
    // Reset when entering a new stratum
    if (strata[person] != istrata) {
      istrata = strata[person];

      // first obs of a new stratum, finish off the prior stratum
      for (; i1 < nused && order1[i1] < person; ++i1) {
        const int p1 = order1[i1];
        for (int i = 0; i < p; ++i) {
          resid(p1,i) -= exp_eta[p1] * (z(p1,i) * cumhaz - xhaz[i]);
        }
      }

      denom = 0.0; // reset temporary variables
      cumhaz = 0.0;
      std::fill(a.begin(), a.end(), 0.0);
      std::fill(xhaz.begin(), xhaz.end(), 0.0);
    }

    const double dtime = tstop[person];

    // process all persons tied at this dtime
    for (; person < nused && tstop[person] == dtime && strata[person] == istrata; ++person) {

      // initialize residuals to score[i] * (x[i] * cumhaz - xhaz), before
      // updating cumhaz and xhaz
      for (int i = 0; i < p; ++i) {
        resid(person,i) = exp_eta[person] * (z(person,i) * cumhaz - xhaz[i]);
      }

      const double w = weight[person];
      const double r = w * exp_eta[person];

      if (event[person] == 0) {
        denom += r;
        for (int i = 0; i < p; ++i) {
          a[i] += r * z(person,i);
        }
      } else {
        ++ndead;
        deadwt += w;
        denom2 += r;
        for (int i = 0; i < p; ++i) {
          a2[i] += r * z(person,i);
        }
      }
    }

    // Remove subjects leaving risk set
    for (; i1 < nused; ++i1) {
      const int p1 = order1[i1];
      if (tstart[p1] < dtime || strata[p1] != istrata) break;

      const double r = weight[p1] * exp_eta[p1];
      denom -= r;
      for (int i = 0; i < p; ++i) {
        // finish the residual by subtracting score[i] * (x[i] * cumhaz - xhaz)
        resid(p1,i) -= exp_eta[p1] * (z(p1,i) * cumhaz - xhaz[i]);
        a[i] -= r * z(p1,i);
      }
    }

    // Add contributions for deaths at this time
    if (ndead > 0) {
      if (method == 0 || ndead == 1) { // Breslow or single event
        denom += denom2;
        const double hazard = deadwt/denom;
        cumhaz += hazard;
        for (int i = 0; i < p; ++i) {
          a[i] += a2[i];
          const double xbar = a[i]/denom;
          xhaz[i] += xbar*hazard;
          for (int j = person-1; j >= person - ndead; j--) {
            resid(j,i) += z(j,i) - xbar;
          }
        }
      } else {  // Efron method
        for (int i = 0; i < p; ++i) {
          mh1[i] = 0.0;
          mh2[i] = 0.0;
          mh3[i] = 0.0;
        }

        const double meanwt = deadwt / ndead;
        const double increment = denom2 / ndead;

        for (int k = 0; k < ndead; ++k) {
          denom += increment;
          const double hazard = meanwt/denom;
          cumhaz += hazard;
          const double downwt = (ndead - k - 1.0)/ndead;
          for (int i = 0; i < p; ++i) {
            a[i] += a2[i] / ndead;
            const double xbar = a[i] / denom;
            xhaz[i] += xbar * hazard;
            mh1[i]  += hazard * downwt;
            mh2[i]  += xbar * hazard * downwt;
            mh3[i]  += xbar / ndead;
          }
        }

        for (int j = person-1; j >= person - ndead; j--) {
          for (int i = 0; i < p; ++i) {
            resid(j,i) += (z(j,i) - mh3[i]) +
              exp_eta[j] * (z(j,i) * mh1[i] - mh2[i]);
          }
        }
      }

      // Reset after processing deaths
      ndead = deadwt = denom2 = 0.0;
      std::fill(a2.begin(), a2.end(), 0.0);
    }
  }

  // finish those remaining in the final stratum
  for (; i1 < nused; ++i1) {
    const int p1 = order1[i1];
    for (int i = 0; i < p; ++i)
      resid(p1,i) -= exp_eta[p1] * (z(p1,i) * cumhaz - xhaz[i]);
  }

  return resid;
}


// main function for phreg
ListCpp phregcpp(const DataFrameCpp& data,
                 const std::vector<std::string>& stratum,
                 const std::string& time,
                 const std::string& time2,
                 const std::string& event,
                 const std::vector<std::string>& covariates,
                 const std::string& weight,
                 const std::string& offset,
                 const std::string& id,
                 const std::string& ties,
                 const std::vector<double>& init,
                 const bool robust,
                 const bool est_basehaz,
                 const bool est_resid,
                 const bool firth,
                 const bool plci,
                 const double alpha,
                 const int maxiter,
                 const double eps) {
  
  int n = data.nrows();
  int p = static_cast<int>(covariates.size());
  if (p == 1 && covariates[0] == "") p = 0;
  
  // --- handle strata (bygroup) ---
  bool has_stratum = false;
  std::vector<int> stratumn(n);
  DataFrameCpp u_stratum;
  int p_stratum = stratum.size();
  if (!(p_stratum == 0 || (p_stratum == 1 && stratum[0] == ""))) {
    ListCpp out = bygroup(data, stratum);
    has_stratum = true;
    stratumn = out.get<std::vector<int>>("index");
    u_stratum = out.get<DataFrameCpp>("lookup");
  }
  
  // --- time / time2 existence and checks ---
  if (!data.containElementNamed(time)) 
    throw std::invalid_argument("data must contain the time variable");
  std::vector<double> timen(n);
  if (data.int_cols.count(time)) {
    const std::vector<int>& vi = data.get<int>(time);
    for (int i = 0; i < n; ++i) timen[i] = static_cast<double>(vi[i]);
  } else if (data.numeric_cols.count(time)) {
    timen = data.get<double>(time);
  } else {
    throw std::invalid_argument("time variable must be integer or numeric");
  }
  for (int i = 0; i < n; ++i) {
    if (!std::isnan(timen[i]) && timen[i] < 0.0)
      throw std::invalid_argument("time must be nonnegative");
  }
  
  bool has_time2 = !time2.empty() && data.containElementNamed(time2);
  std::vector<double> time2n(n);
  if (has_time2) {
    if (data.int_cols.count(time2)) {
      const std::vector<int>& vi = data.get<int>(time2);
      for (int i = 0; i < n; ++i) time2n[i] = static_cast<double>(vi[i]);
    } else if (data.numeric_cols.count(time2)) {
      time2n = data.get<double>(time2);
    } else {
      throw std::invalid_argument("time2 variable must be integer or numeric");
    }    
    for (int i = 0; i < n; ++i) {
      if (!std::isnan(timen[i]) && !std::isnan(time2n[i]) && time2n[i] <= timen[i])
        throw std::invalid_argument("time2 must be greater than time");
    }
  }
  
  // --- event variable ---
  bool has_event = !event.empty() && data.containElementNamed(event);
  if (!has_time2 && !has_event) {
    throw std::invalid_argument("data must contain the event variable for right censored data"); 
  }
  std::vector<double> eventn(n);
  if (has_event) {
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
  }
  
  // --- build design matrix zn (n x p) column-major FlatMatrix ---
  FlatMatrix zn(n, p);
  for (int j = 0; j < p; ++j) {
    const std::string& zj = covariates[j];
    if (!data.containElementNamed(zj)) 
      throw std::invalid_argument("data must contain the variables in covariates");
    if (data.bool_cols.count(zj)) {
      const std::vector<unsigned char>& vb = data.get<unsigned char>(zj);
      int off = j * n;
      for (int i = 0; i < n; ++i) zn.data[off + i] = vb[i] ? 1.0 : 0.0;
    } else if (data.int_cols.count(zj)) {
      const std::vector<int>& vi = data.get<int>(zj);
      int off = j * n;
      for (int i = 0; i < n; ++i) zn.data[off + i] = static_cast<double>(vi[i]);
    } else if (data.numeric_cols.count(zj)) {
      const std::vector<double>& vd = data.get<double>(zj);
      int off = j * n;
      for (int i = 0; i < n; ++i) zn.data[off + i] = vd[i];
    } else {
      throw std::invalid_argument("covariates must be bool, integer or numeric");
    }
  }
  
  // --- weight and offset ---
  std::vector<double> weightn(n, 1.0);
  if (!weight.empty() && data.containElementNamed(weight)) {
    if (data.int_cols.count(weight)) {
      const std::vector<int>& vi = data.get<int>(weight);
      for (int i = 0; i < n; ++i) weightn[i] = static_cast<double>(vi[i]);
    } else if (data.numeric_cols.count(weight)) {
      weightn = data.get<double>(weight);
    } else {
      throw std::invalid_argument("weight variable must be integer or numeric");
    }
    for (double w : weightn) if (std::isnan(w) || w <= 0.0) 
      throw std::invalid_argument("weight must be greater than 0");
  }
  
  std::vector<double> offsetn(n, 0.0);
  if (!offset.empty() && data.containElementNamed(offset)) {
    if (data.int_cols.count(offset)) {
      const std::vector<int>& vi = data.get<int>(offset);
      for (int i = 0; i < n; ++i) offsetn[i] = static_cast<double>(vi[i]);
    } else if (data.numeric_cols.count(offset)) {
      offsetn = data.get<double>(offset);
    } else {
      throw std::invalid_argument("offset variable must be integer or numeric");
    }
  }
  
  // --- id mapping ---
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
    } else throw std::invalid_argument("incorrect type for the id variable in the input data");
  }
  
  if (robust && has_time2 && !has_id) {
    throw std::invalid_argument("id is needed for counting process data with robust variance");
  }
  
  std::string meth = ties;
  std::for_each(meth.begin(), meth.end(), [](char & c) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  });
  
  int method = meth == "efron" ? 1 : 0;
  
  // unify right censored data with counting process data
  std::vector<double> tstartn(n), tstopn(n);
  if (!has_time2) {
    tstopn = timen;
  } else {
    tstartn = timen;
    tstopn = time2n;
  }
  
  
  // exclude observations with missing values
  std::vector<unsigned char> sub(n,1);
  for (int i=0; i<n; ++i) {
    if (stratumn[i] == INT_MIN || idn[i] == INT_MIN ||
        std::isnan(tstartn[i]) || std::isnan(tstopn[i]) ||
        std::isnan(eventn[i]) || std::isnan(weightn[i]) ||
        std::isnan(offsetn[i])) {
      sub[i] = 0;
      continue;
    }
    for (int j = 0; j < p; ++j) {
      if (std::isnan(zn(i,j))) { sub[i] = 0; break; }
    }
  }
  std::vector<int> keep = which(sub);
  if (keep.empty()) throw std::invalid_argument("no observations without missing values");
  
  subset_in_place(stratumn, keep);
  subset_in_place(tstartn, keep);
  subset_in_place(tstopn, keep);
  subset_in_place(eventn, keep);
  subset_in_place(weightn, keep);
  subset_in_place(offsetn, keep);
  subset_in_place(idn, keep);
  subset_in_place_flatmatrix(zn, keep);
  n = keep.size();
  
  // sumstat data set
  double loglik0, loglik1;
  double regloglik0, regloglik1;
  double scoretest;
  int niter;
  bool fail;
  
  // parest data set
  std::vector<std::string> par(p);
  std::vector<double> b(p), seb(p), rseb(p);
  std::vector<double> z(p), expbeta(p);
  FlatMatrix vb(p, p), rvb(p, p);
  std::vector<double> lb(p), ub(p), prob(p);
  std::vector<std::string> clparm(p);
  
  
  // baseline hazards data set
  std::vector<int> stratumn1 = unique_sorted(stratumn);
  int nstrata = static_cast<int>(stratumn1.size());
  int N = n + nstrata; // add a time 0 row for each stratum
  std::vector<int> dstratum;
  std::vector<double> dtime, dnrisk, dnevent, dncensor;
  std::vector<double> dhaz, dvarhaz;
  FlatMatrix dgradhaz(N,p);
  dstratum.reserve(N); 
  dtime.reserve(N); dnrisk.reserve(N); dnevent.reserve(N); dncensor.reserve(N);
  dhaz.reserve(N); dvarhaz.reserve(N);
  
  // martingale residuals
  std::vector<double> resmart(n);
  
  // linear predictors
  std::vector<double> linear_predictors(n);
  
  double zcrit = boost_qnorm(1.0 - alpha / 2.0);
  double xcrit = zcrit * zcrit;
  
  
  double nobs = n;
  double nevents = std::accumulate(eventn.begin(), eventn.end(), 0.0); 
  
  if (nevents == 0) {
    if (p > 0) {
      for (int i=0; i<p; ++i) {
        par[i] = covariates[i];
        b[i] = NaN;
        seb[i] = 0;
        rseb[i] = 0;
        z[i] = NaN;
        expbeta[i] = NaN;
        lb[i] = NaN;
        ub[i] = NaN;
        prob[i] = NaN;
        clparm[i] = "Wald";
      }
      
      for (int j=0; j<p; ++j) {
        for (int i=0; i<p; ++i) {
          vb(i,j) = 0;
          rvb(i,j) = 0;
        }
      }
    }
    
    // baseline hazard
    if (est_basehaz) {
      dstratum[0] = 0;
      dtime[0] = 0;
      dnrisk[0] = n;
      dnevent[0] = 0;
      dncensor[0] = 0;
      dhaz[0] = 0;
      dvarhaz[0] = 0;
      if (p > 0) {
        for (int i=0; i<p; ++i) {
          dgradhaz(0,i) = 0;
        }
      }
    }
    
    // martingale residuals
    if (est_resid) {
      for (int i=0; i<n; ++i) {
        resmart[i] = 0;
      }
    }
    
    // linear predictors
    for (int i=0; i<n; ++i) {
      linear_predictors[i] = offsetn[i];
    }
    
    loglik0 = NaN;
    loglik1 = NaN;
    regloglik0 = NaN;
    regloglik1 = NaN;
    scoretest = NaN;
    niter = 0;
    fail = true;
  } else {
    // sort by stratum
    std::vector<int> order0 = seqcpp(0, n-1);
    std::sort(order0.begin(), order0.end(), [&](int i, int j) {
      return stratumn[i] < stratumn[j];
    });
    
    std::vector<int> stratumnz = subset(stratumn, order0);
    std::vector<double> tstartnz = subset(tstartn, order0);
    std::vector<double> tstopnz = subset(tstopn, order0);
    std::vector<double> eventnz = subset(eventn, order0);
    
    // locate the first observation within each stratum
    std::vector<int> istratum(1,0);
    for (int i=1; i<n; ++i) {
      if (stratumnz[i] != stratumnz[i-1]) {
        istratum.push_back(i);
      }
    }
    
    istratum.push_back(n);
    
    // ignore subjects not at risk for any event time
    std::vector<int> ignorenz(n);
    for (int i = 0; i < nstrata; ++i) {
      std::vector<int> q0 = seqcpp(istratum[i], istratum[i+1]-1);
      std::vector<double> tstart0 = subset(tstartnz, q0);
      std::vector<double> tstop0 = subset(tstopnz, q0);
      std::vector<double> event0 = subset(eventnz, q0);
      int n0 = static_cast<int>(q0.size());
      
      // unique event times
      std::vector<double> etime; etime.reserve(n0);
      for (int j = 0; j < n0; ++j) {
        if (event0[j] == 1) etime.push_back(tstop0[j]);
      }
      etime = unique_sorted(etime);
      
      std::vector<int> index1 = findInterval3(tstart0, etime);
      std::vector<int> index2 = findInterval3(tstop0, etime);
      for (int j = istratum[i]; j < istratum[i+1]; ++j) {
        int j0 = j-istratum[i];
        if (index1[j0] == index2[j0]) { // no event in (tstart, tstop]
          ignorenz[j] = 1;
        } else {
          ignorenz[j] = 0;
        }
      }
    }
    
    std::vector<int> ignoren(n); // back to the original order
    for (int i = 0; i < n; ++i) {
      ignoren[order0[i]] = ignorenz[i];
    }
    
    int nignore = std::accumulate(ignoren.begin(), ignoren.end(), 0);
    int nused = n - nignore; // number of used observations
    
    // sort by stopping time in descending order within each stratum
    std::vector<int> order2 = seqcpp(0, n-1);
    std::sort(order2.begin(), order2.end(), [&](int i, int j) {
      if (ignoren[i] != ignoren[j]) return ignoren[i] < ignoren[j];
      if (stratumn[i] != stratumn[j]) return stratumn[i] < stratumn[j];
      if (tstopn[i] != tstopn[j]) return tstopn[i] > tstopn[j];
      return eventn[i] < eventn[j];
    });
    
    std::vector<int> stratumna = subset(stratumn, order2);
    std::vector<double> tstartna = subset(tstartn, order2);
    std::vector<double> tstopna = subset(tstopn, order2);
    std::vector<double> eventna = subset(eventn, order2);
    std::vector<double> weightna = subset(weightn, order2);
    std::vector<double> offsetna = subset(offsetn, order2);
    std::vector<int> idna = subset(idn, order2);
    std::vector<int> ignorena = subset(ignoren, order2);
    FlatMatrix zna;
    if (p > 0) zna = subset_flatmatrix(zn, order2);
    
    // sort by starting time in descending order within each stratum
    std::vector<int> orderna = seqcpp(0, n-1);
    std::sort(orderna.begin(), orderna.end(), [&](int i, int j) {
      if (ignorena[i] != ignorena[j]) return ignorena[i] < ignorena[j];
      if (stratumna[i] != stratumna[j]) return stratumna[i] < stratumna[j];
      return tstartna[i] > tstartna[j];
    });
    
    coxparams param = {nused, stratumna, tstartna, tstopna, eventna,
                       weightna, offsetna, zna, orderna, method};
    
    std::vector<double> bint(p);
    ListCpp derint = f_der_2(p, bint, &param, firth);
    
    ListCpp out;
    
    if (p > 0) {
      std::vector<int> colfit = seqcpp(0, p - 1);
      if  (static_cast<int>(init.size()) == p && 
           std::none_of(init.begin(), init.end(), [](double val){ 
             return std::isnan(val); })) {
        out = phregloop(p, init, &param, maxiter, eps, firth, colfit, p);
      } else {
        out = phregloop(p, bint, &param, maxiter, eps, firth, colfit, p);
      }
      
      niter = out.get<int>("iter");
      fail = out.get<bool>("fail");
      if (fail) {
        thread_utils::push_thread_warning(
          "phregloop failed to converge for the full model; continuing with current results."); 
      }
      
      b = out.get<std::vector<double>>("coef");
      vb = out.get<FlatMatrix>("var");
      
      for (int j=0; j<p; ++j) {
        seb[j] = std::sqrt(vb(j,j));
      }
      
      for (int i=0; i<p; ++i) {
        par[i] = covariates[i];
      }
      
      // score statistic
      std::vector<double> scorebint;
      if (firth) scorebint = derint.get<std::vector<double>>("regscore");
      else scorebint = derint.get<std::vector<double>>("score");
      FlatMatrix infobint = derint.get<FlatMatrix>("imat");
      FlatMatrix vbint = invsympd(infobint, p);
      scoretest = quadsym(scorebint, vbint);
      
      
      // robust variance estimates
      if (robust) {
        FlatMatrix ressco = f_ressco_2(p, b, &param);
        
        int nr; // number of rows in the score residual matrix
        if (!has_id) {
          for (int j = 0; j < p; ++j) {
            const int off = j * n;
            for (int i = 0; i < n; ++i) {
              ressco.data[off + i] *= weightna[i];
            }
          }
          nr = n;
        } else { // need to sum up score residuals by id
          std::vector<int> order = seqcpp(0, n-1);
          std::sort(order.begin(), order.end(), [&](int i, int j) {
            return idna[i] < idna[j];
          });
          
          std::vector<int> id1 = subset(idna, order);
          std::vector<int> idx(1,0);
          for (int i=1; i<n; ++i) {
            if (id1[i] != id1[i-1]) {
              idx.push_back(i);
            }
          }
          
          int nids = static_cast<int>(idx.size());
          idx.push_back(n);
          
          FlatMatrix ressco1(nids,p);
          for (int j = 0; j < p; ++j) {
            const int off = j * nids;
            for (int i = 0; i < nids; ++i) {
              double sum = 0.0;
              for (int k = idx[i]; k < idx[i+1]; ++k) {
                int row = order[k];
                sum  += weightna[row] * ressco(row,j);
              }
              ressco1.data[off + i] = sum;
            }
          }
          
          ressco = std::move(ressco1);  // update the score residuals
          nr = nids;
        }
        
        FlatMatrix D = mat_mat_mult(ressco, vb); // DFBETA
        
        const double* Dptr = D.data_ptr();
        double* rvbptr = rvb.data_ptr();
        for (int j = 0; j < p; ++j) {
          const double* Dj = Dptr + j * nr; // pointer to D(:,j)
          for (int k = 0; k <= j; ++k) {
            const double* Dk = Dptr + k * nr; // pointer to D(:,k)
            double sum = 0.0;
            for (int i = 0; i < nr; ++i) {
              sum += Dj[i] * Dk[i];
            }
            rvbptr[k * p + j] = sum;
            if (j != k) rvbptr[j * p + k] = sum;
          }
        }
        
        for (int i=0; i<p; ++i) {
          rseb[i] = std::sqrt(rvb(i,i));
        }
      }
      
      // profile likelihood confidence interval for regression coefficients
      if (plci) {
        double lmax = out.get<double>("loglik");
        double l0 = lmax - 0.5 * xcrit;
        
        for (int k=0; k<p; ++k) {
          lb[k] = phregplloop(p, b, &param, maxiter, eps, firth, k, -1, l0);
          ub[k] = phregplloop(p, b, &param, maxiter, eps, firth, k, 1, l0);
          
          std::vector<int> colfit1(p-1);
          for (int i = 0, j = 0; i < p; ++i) {
            if (i == k) continue;
            colfit1[j++] = i;
          }
          
          std::vector<double> b0(p);
          ListCpp out0 = phregloop(p, b0, &param, maxiter, eps, firth, colfit1, p-1);
          double lmax0 = out0.get<double>("loglik");
          prob[k] = boost_pchisq(-2.0 * (lmax0 - lmax), 1, 0);
          clparm[k] = "PL";
        }
      } else {
        for (int k=0; k<p; ++k) {
          if (!robust) {
            lb[k] = b[k] - zcrit*seb[k];
            ub[k] = b[k] + zcrit*seb[k];
            prob[k] = boost_pchisq(sq(b[k] / seb[k]), 1, 0);
          } else {
            lb[k] = b[k] - zcrit*rseb[k];
            ub[k] = b[k] + zcrit*rseb[k];
            prob[k] = boost_pchisq(sq(b[k] / rseb[k]), 1, 0);
          }
          clparm[k] = "Wald";
        }
      }
      
      if (firth) {
        loglik0 = derint.get<double>("loglik");
        loglik1 = out.get<double>("loglik");
        regloglik0 = derint.get<double>("regloglik");
        regloglik1 = out.get<double>("regloglik");
      } else {
        loglik0 = derint.get<double>("loglik");
        loglik1 = out.get<double>("loglik");
        regloglik0 = loglik0;
        regloglik1 = loglik1;
      }
    } else {
      fail = false;
      loglik0 = derint.get<double>("loglik");
      loglik1 = derint.get<double>("loglik");
      regloglik0 = loglik0;
      regloglik1 = loglik1;
      scoretest = 0.0;
      niter = 0;
    }
    
    // estimate baseline hazard
    if (est_basehaz) {
      // prepare the data for estimating baseline hazards at all time points
      
      // sort by stopping time in descending order within each stratum
      std::vector<int> order3 = seqcpp(0, n-1);
      std::sort(order3.begin(), order3.end(), [&](int i, int j) {
        if (stratumn[i] != stratumn[j]) return stratumn[i] > stratumn[j];
        return tstopn[i] > tstopn[j];
      });
      
      std::vector<int> stratumnb = subset(stratumn, order3);
      std::vector<double> tstartnb = subset(tstartn, order3);
      std::vector<double> tstopnb = subset(tstopn, order3);
      std::vector<double> eventnb = subset(eventn, order3);
      std::vector<double> weightnb = subset(weightn, order3);
      std::vector<double> offsetnb = subset(offsetn, order3);
      FlatMatrix znb;
      if (p > 0) znb = subset_flatmatrix(zn, order3);
      
      // sort by starting time in descending order within each stratum
      std::vector<int> ordernb = seqcpp(0, n-1);
      std::sort(ordernb.begin(), ordernb.end(), [&](int i, int j) {
        if (stratumnb[i] != stratumnb[j]) return stratumnb[i] > stratumnb[j];
        return tstartnb[i] > tstartnb[j];
      });
      
      coxparams paramb = {n, stratumnb, tstartnb, tstopnb, eventnb,
                          weightnb, offsetnb, znb, ordernb, method};
      
      ListCpp basehazn = f_basehaz(p, b, &paramb);
      
      dstratum = basehazn.get<std::vector<int>>("stratum");
      dtime = basehazn.get<std::vector<double>>("time");
      dnrisk = basehazn.get<std::vector<double>>("nrisk");
      dnevent = basehazn.get<std::vector<double>>("nevent");
      dncensor = basehazn.get<std::vector<double>>("ncensor");
      dhaz = basehazn.get<std::vector<double>>("haz");
      dvarhaz = basehazn.get<std::vector<double>>("varhaz");
      if (p > 0) dgradhaz = basehazn.get<FlatMatrix>("gradhaz");
    }
    
    // martingale residuals
    if (est_resid) {
      std::vector<double> resid = f_resmart(p, b, &param);
      
      for (int i=0; i<n; ++i) {
        resmart[order2[i]] = resid[i];
      }
    }
    
    // linear predictors
    for (int i=0; i<n; ++i) {
      linear_predictors[order2[i]] = offsetna[i];
    }
    
    if (p > 0) {
      for (int j=0; j<p; ++j) {
        double beta = b[j];
        if (beta == 0.0) continue;
        const int off = j * n;
        for (int i=0; i<n; ++i) {
          linear_predictors[order2[i]] += beta * zna.data[off + i];
        }
      }
    }
    
    // compute exp(beta)
    for (int i = 0; i < p; ++i) {
      expbeta[i] = std::exp(b[i]);
    }
    
    // compute z statistics
    if (robust) {
      for (int i = 0; i < p; ++i) {
        if (rseb[i] == 0) {
          z[i] = NaN;
        } else {
          z[i] = b[i]/rseb[i];
        }
      }
    } else {
      for (int i = 0; i < p; ++i) {
        if (seb[i] == 0) {
          z[i] = NaN;
        } else {
          z[i] = b[i]/seb[i];
        }
      }
    }
  }
  
  if (est_basehaz) {
    std::vector<int> sub = seqcpp(0, dstratum.size() - 1);
    if (p > 0) dgradhaz = subset_flatmatrix(dgradhaz, sub);
  }
  
  // prepare the output data sets
  DataFrameCpp sumstat;
  sumstat.push_back(nobs, "n");
  sumstat.push_back(nevents, "nevents");
  sumstat.push_back(loglik0, "loglik0");
  sumstat.push_back(loglik1, "loglik1");
  sumstat.push_back(scoretest, "scoretest");
  sumstat.push_back(niter, "niter");
  sumstat.push_back(meth, "ties");
  sumstat.push_back(p, "p");
  sumstat.push_back(robust, "robust");
  sumstat.push_back(firth, "firth");
  sumstat.push_back(fail, "fail");
  
  if (p > 0 && firth) {
    sumstat.push_back(regloglik0, "loglik0_unpenalized");
    sumstat.push_back(regloglik1, "loglik1_unpenalized");
  }
  
  ListCpp result;
  result.push_back(sumstat, "sumstat");
  
  if (p > 0) {
    std::vector<double> sebeta = robust ? rseb : seb;
    FlatMatrix vbeta = robust ? rvb : vb;
    DataFrameCpp parest;
    parest.push_back(std::move(par), "param");
    parest.push_back(std::move(b), "beta");
    parest.push_back(std::move(sebeta), "sebeta");
    parest.push_back(std::move(z), "z");
    parest.push_back(std::move(expbeta), "expbeta");
    parest.push_back(std::move(vbeta), "vbeta");
    parest.push_back(std::move(lb), "lower");
    parest.push_back(std::move(ub), "upper");
    parest.push_back(std::move(prob), "p");
    parest.push_back(std::move(clparm), "method");
    
    if (robust) {
      parest.push_back(std::move(seb), "sebeta_naive");
      parest.push_back(std::move(vb), "vbeta_naive");
    }
    
    result.push_back(std::move(parest), "parest");
  }
  
  if (est_basehaz) {
    ListCpp basehaz;
    basehaz.push_back(std::move(dtime), "time");
    basehaz.push_back(std::move(dnrisk), "nrisk");
    basehaz.push_back(std::move(dnevent), "nevent");
    basehaz.push_back(std::move(dncensor), "ncensor");
    basehaz.push_back(std::move(dhaz), "haz");
    basehaz.push_back(std::move(dvarhaz), "varhaz");
    
    if (p > 0) basehaz.push_back(std::move(dgradhaz), "gradhaz");
    
    if (has_stratum) {
      for (int i = 0; i < p_stratum; ++i) {
        std::string s = stratum[i];
        if (u_stratum.int_cols.count(s)) {
          auto v = u_stratum.get<int>(s);
          subset_in_place(v, dstratum);
          result.push_back(std::move(v), s);
        } else if (u_stratum.numeric_cols.count(s)) {
          auto v = u_stratum.get<double>(s);
          subset_in_place(v, dstratum);
          result.push_back(std::move(v), s);
        } else if (u_stratum.string_cols.count(s)) {
          auto v = u_stratum.get<std::string>(s);
          subset_in_place(v, dstratum);
          result.push_back(std::move(v), s);
        } else {
          throw std::invalid_argument("unsupported type for stratum variable " + s);
        }
      }
    }
    
    result.push_back(std::move(basehaz), "basehaz");
  }
  
  if (est_resid) result.push_back(std::move(resmart), "residuals");
  result.push_back(std::move(linear_predictors), "linear_predictors");
  
  return result;
}


// [[Rcpp::export]]
Rcpp::List phregRcpp(const Rcpp::DataFrame& data,
                     const std::vector<std::string>& stratum,
                     const std::string& time,
                     const std::string& time2,
                     const std::string& event,
                     const std::vector<std::string>& covariates,
                     const std::string& weight,
                     const std::string& offset,
                     const std::string& id,
                     const std::string& ties,
                     const std::vector<double>& init,
                     const bool robust,
                     const bool est_basehaz,
                     const bool est_resid,
                     const bool firth,
                     const bool plci,
                     const double alpha,
                     const int maxiter,
                     const double eps) {
  
  DataFrameCpp dfcpp = convertRDataFrameToCpp(data);
  
  ListCpp cpp_result = phregcpp(
    dfcpp, stratum, time, time2, event, covariates, weight, offset, id, 
    ties, init, robust, est_basehaz, est_resid, firth, plci, alpha, 
    maxiter, eps);
  
  thread_utils::drain_thread_warnings_to_R();
  return Rcpp::wrap(cpp_result);
}


// survival function estimation based on Cox model
DataFrameCpp survfit_phregcpp(const int p,
                              const std::vector<double>& beta,
                              const FlatMatrix& vbeta,
                              const DataFrameCpp& basehaz,
                              const DataFrameCpp& newdata,
                              const std::vector<std::string>& covariates,
                              const std::vector<std::string>& stratum,
                              const std::string& offset,
                              const std::string& id,
                              const std::string& tstart,
                              const std::string& tstop,
                              const bool sefit,
                              const std::string& conftype,
                              const double conflev) {
  
  int n0 = basehaz.nrows();
  int n = newdata.nrows();
  int nvar = static_cast<int>(covariates.size());
  
  std::string ct = conftype;
  std::for_each(ct.begin(), ct.end(), [](char & c) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  });

  if (!(ct=="none" || ct=="plain" || ct=="log" || ct=="log-log" || 
      ct=="logit" || ct=="arcsin")) {
    throw std::invalid_argument("conftype must be none, plain, log, log-log, logit, or arcsin");
  }

  if (conflev <= 0.0 || conflev >= 1.0) {
    throw std::invalid_argument("conflev must lie between 0 and 1");
  }

  double zcrit = boost_qnorm((1.0 + conflev) / 2.0);

  std::vector<int> stratumn0(n0);
  DataFrameCpp u_stratum0;
  std::vector<int> nlevels;
  ListCpp& lookups;
  int p_stratum = static_cast<int>(stratum.size());
  if (!(p_stratum == 0 || (p_stratum == 1 && stratum[0] == ""))) {
    ListCpp out = bygroup(basehaz, stratum);
    stratumn0 = out.get<std::vector<int>>("index");
    u_stratum0 = out.get<DataFrameCpp>("lookup");
    nlevels = out.get<std::vector<int>>("nlevels");
    lookups = out.get_list("lookups_per_variable");
  }

  bool nullmodel = (p == 0 || (nvar == 1 && covariates[0] == ""));

  if (!nullmodel && nvar != p) {
    throw std::invalid_argument("incorrect number of covariates for the Cox model");
  }

  FlatMatrix zn(n,p);
  for (int j = 0; j < p; ++j) {
    const std::string& zj = covariates[j];
    if (!newdata.containElementNamed(zj)) 
      throw std::invalid_argument("newdata must contain the variables in covariates");
    if (newdata.bool_cols.count(zj)) {
      const std::vector<unsigned char>& vb = newdata.get<unsigned char>(zj);
      int off = j * n;
      for (int i = 0; i < n; ++i) zn.data[off + i] = vb[i] ? 1.0 : 0.0;
    } else if (newdata.int_cols.count(zj)) {
      const std::vector<int>& vi = newdata.get<int>(zj);
      int off = j * n;
      for (int i = 0; i < n; ++i) zn.data[off + i] = static_cast<double>(vi[i]);
    } else if (newdata.numeric_cols.count(zj)) {
      const std::vector<double>& vd = newdata.get<double>(zj);
      int off = j * n;
      for (int i = 0; i < n; ++i) zn.data[off + i] = vd[i];
    } else {
      throw std::invalid_argument("covariates must be bool, integer or numeric");
    }
  }

  bool has_stratum;
  std::vector<int> stratumn(n);
  if (p_stratum == 1 && (stratum[0] == "" || stratum[0] == "none")) {
    has_stratum = false;
  } else {
    has_stratum = true;
    // match stratum in newdata to stratum in basehaz
    int orep = u_stratum0.nrows();
    for (int i = 0; i < p_stratum; ++i) {
      orep /= nlevels[i];
      std::string s = stratum[i];
      std::vector<int> idx;
      
      if (newdata.bool_cols.count(s) || newdata.int_cols.count(s)) {
        std::vector<int> v;
        std::vector<int> w;
        if (newdata.bool_cols.count(s)) {
          auto vb = newdata.get<unsigned char>(s);
          auto wb = lookups.get<std::vector<unsigned char>>(s);
          v.resize(n);
          w.resize(n);
          for (int j = 0; j < n; ++j) {
            v[j] = vb[j] ? 1 : 0;
            w[j] = wb[j] ? 1 : 0;
          }
        } else {
          v = newdata.get<int>(s);
          w = lookups.get<std::vector<int>>(s);
        }
        idx = matchcpp(v, w);
      } else if (newdata.numeric_cols.count(s)) {
        auto v = newdata.get<double>(s);
        auto w = lookups.get<std::vector<double>>(s);
        idx = matchcpp(v, w);
      } else if (newdata.string_cols.count(s)) {
        auto v = newdata.get<std::string>(s);
        auto w = lookups.get<std::vector<std::string>>(s);
        idx = matchcpp(v, w);
      } else {
        throw std::invalid_argument("Unsupported type for stratum variable: " + s);
      }
      
      for (int person = 0; person < n; ++person) {
        stratumn[person] += idx[person] * orep; 
      }
    }
  }
  
  std::vector<double> offsetn(n, 0.0);
  if (!offset.empty() && newdata.containElementNamed(offset)) {
    if (newdata.int_cols.count(offset)) {
      const std::vector<int>& vi = newdata.get<int>(offset);
      for (int i = 0; i < n; ++i) offsetn[i] = static_cast<double>(vi[i]);
    } else if (newdata.numeric_cols.count(offset)) {
      offsetn = newdata.get<double>(offset);
    } else {
      throw std::invalid_argument("offset variable must be integer or numeric");
    }
  }
  
  std::vector<double> time0 = basehaz.get<double>("time");
  std::vector<double> nrisk0 = basehaz.get<double>("nrisk");
  std::vector<double> nevent0 = basehaz.get<double>("nevent");
  std::vector<double> ncensor0 = basehaz.get<double>("ncensor");
  std::vector<double> haz0 = basehaz.get<double>("haz");
  std::vector<double> vhaz0 = basehaz.get<double>("varhaz");
  FlatMatrix ghaz0(n0,p);
  for (int j = 0; j < p; ++j) {
    std::string col_name = "gradhaz";
    if (p>1) col_name += "." + std::to_string(j+1);
    std::vector<double> u = basehaz.get<double>(col_name);
    for (int i = 0; i < n0; ++i) {
      ghaz0(i,j) = u[i];
    }
  }

  // create the numeric id variable
  bool has_id = !id.empty() && newdata.containElementNamed(id);
  std::vector<int> idn(n);
  std::vector<int> idwi;
  std::vector<double> idwn;
  std::vector<std::string> idwc;
  if (!has_id) {
    idn = seqcpp(0, n-1);
  } else { // input data has the counting process style of input
    if (newdata.int_cols.count(id)) {
      auto v = newdata.get<int>(id);
      idwi = unique_sorted(v);
      idn = matchcpp(v, idwi);
    } else if (newdata.numeric_cols.count(id)) {
      auto v = newdata.get<double>(id);
      idwn = unique_sorted(v);
      idn = matchcpp(v, idwn);
    } else if (newdata.string_cols.count(id)) {
      auto v = newdata.get<std::string>(id);
      idwc = unique_sorted(v);
      idn = matchcpp(v, idwc);
    } else throw std::invalid_argument("incorrect type for the id variable in newdata");
  }
  
  // unify right-censoring data with counting process data
  std::vector<double> tstartn(n), tstopn(n);
  if (!has_id) { // right-censored data
    double maxt0 = *std::max_element(time0.begin(), time0.end()) + 1.0;
    std::fill(tstopn.begin(), tstopn.end(), maxt0);
  } else {
    if (!newdata.containElementNamed(tstart)) 
      throw std::invalid_argument("newdata must contain the tstart variable");
    std::vector<double> tstartn(n);
    if (newdata.int_cols.count(tstart)) {
      const std::vector<int>& vi = newdata.get<int>(tstart);
      for (int i = 0; i < n; ++i) tstartn[i] = static_cast<double>(vi[i]);
    } else if (newdata.numeric_cols.count(tstart)) {
      tstartn = newdata.get<double>(tstart);
    } else {
      throw std::invalid_argument("tstart variable must be integer or numeric");
    }
    for (int i = 0; i < n; ++i) {
      if (!std::isnan(tstartn[i]) && tstartn[i] < 0.0)
        throw std::invalid_argument("tstart must be nonnegative");
    }

    if (!newdata.containElementNamed(tstop)) 
      throw std::invalid_argument("newdata must contain the tstop variable");
    std::vector<double> tstopn(n);
    if (newdata.int_cols.count(tstop)) {
      const std::vector<int>& vi = newdata.get<int>(tstop);
      for (int i = 0; i < n; ++i) tstopn[i] = static_cast<double>(vi[i]);
    } else if (newdata.numeric_cols.count(tstop)) {
      tstopn = newdata.get<double>(tstop);
    } else {
      throw std::invalid_argument("tstop variable must be integer or numeric");
    }
    for (int i = 0; i < n; ++i) {
      if (!std::isnan(tstartn[i]) && !std::isnan(tstopn[i]) && tstopn[i] <= tstartn[i])
        throw std::invalid_argument("tstop must be greater than tstart");
    }
  }

  // order data by id and tstop, assuming consecutive intervals
  std::vector<int> order = seqcpp(0, n-1);
  std::sort(order.begin(), order.end(), [&](int i, int j) {
    if (idn[i] != idn[j]) return idn[i] < idn[j];
    return tstopn[i] < tstopn[j];
  });

  subset_in_place(idn, order);
  subset_in_place(stratumn, order);
  subset_in_place(tstartn, order);
  subset_in_place(tstopn, order);
  subset_in_place(offsetn, order);
  if (p > 0) subset_in_place_flatmatrix(zn, order);

  // exclude observations with missing values
  std::vector<unsigned char> sub(n,1);
  for (int i=0; i<n; ++i) {
    if (stratumn[i] == INT_MIN || idn[i] == INT_MIN ||
        std::isnan(tstartn[i]) || std::isnan(tstopn[i]) ||
        std::isnan(offsetn[i])) {
      sub[i] = 0; continue;
    }
    for (int j=0; j<p; ++j) {
      if (std::isnan(zn(i,j))) { sub[i] = 0; break; }
    }
  }

  std::vector<int> keep = which(sub);
  if (keep.empty()) throw std::invalid_argument("no observations without missing values");
  subset_in_place(stratumn, keep);
  subset_in_place(offsetn, keep);
  subset_in_place(idn, keep);
  subset_in_place(tstartn, keep);
  subset_in_place(tstopn, keep);
  if (p > 0) subset_in_place_flatmatrix(zn, keep);
  n = static_cast<int>(keep.size());
  
  // risk score
  std::vector<double> eta = offsetn;
  for (int j = 0; j < p; ++j) {
    double b = beta[j];
    if (b == 0.0) continue;
    const int off = j * n;
    for (int i = 0; i < n; ++i) {
      eta[i] += b * zn.data[off + i];
    }
  }
  std::vector<double> risk(n);
  for (int i = 0; i < n; ++i) {
    risk[i] = std::exp(eta[i]);
  }

  // count number of observations for each id
  std::vector<int> idx(1,0);
  for (int i=1; i<n; ++i) {
    if (idn[i] != idn[i-1]) {
      idx.push_back(i);
    }
  }

  int nids = static_cast<int>(idx.size());
  idx.push_back(n);

  int N = nids*n0; // upper bound on the number of rows in the output
  std::vector<double> time, nrisk, nevent, ncensor, cumhaz, vcumhaz, secumhaz;
  std::vector<int> strata, ids;
  time.reserve(N); nrisk.reserve(N); nevent.reserve(N); ncensor.reserve(N);
  cumhaz.reserve(N); vcumhaz.reserve(N); secumhaz.reserve(N);
  FlatMatrix z(N,p);
  
  // process by id
  int l = 0;
  for (int h=0; h<nids; ++h) {
    std::vector<int> q1 = seqcpp(idx[h], idx[h+1] - 1);
    int n1 = static_cast<int>(q1.size());

    std::vector<int> id1 = subset(idn, q1);
    std::vector<int> stratum1 = subset(stratumn, q1);
    std::vector<double> tstart1 = subset(tstartn, q1);
    std::vector<double> tstop1 = subset(tstopn, q1);
    std::vector<double> risk1 = subset(risk, q1);
    FlatMatrix z1 = subset_flatmatrix(zn, q1);
    
    std::vector<double> tstop2(n1 + 1);
    std::memcpy(tstop2.data() + 1, tstop1.data(), n1 * sizeof(double));
    
    // match the stratum in basehaz
    std::vector<int> idx1;
    for (int i = 0; i < n0; ++i) {
      if (stratumn0[i] == stratum1[0]) idx1.push_back(i);
    }
    std::vector<double> time01 = subset(time0, idx1);

    // left-open and right-closed intervals containing the event time
    std::vector<int> idx2 = findInterval3(time01, tstop2, 0, 0, 1);
    std::vector<int> sub;
    for (int i = 0; i < static_cast<int>(idx2.size()); ++i) {
      if (idx2[i] >= 1 && idx2[i] <= n1) sub.push_back(i);
    }
    int m1 = sub.size();

    if (m1 != 0) {
      std::vector<int> idx3 = subset(idx1, sub);
      std::vector<double> time1 = subset(time0, idx3);
      std::vector<double> nrisk1 = subset(nrisk0, idx3);
      std::vector<double> nevent1 = subset(nevent0, idx3);
      std::vector<double> ncensor1 = subset(ncensor0, idx3);
      std::vector<double> haz1 = subset(haz0, idx3);

      std::vector<int> idx4 = subset(idx2, sub);
      for (int i = 0; i < m1; ++i) idx4[i] -= 1; // change to 0-1 indexing

      // cumulative hazards
      for (int i = 0; i < m1; ++i) {
        int r = l + i;
        time[r] = time1[i];
        nrisk[r] = nrisk1[i];
        nevent[r] = nevent1[i];
        ncensor[r] = ncensor1[i];

        int k = idx4[i];
        ids[r] = id1[k];
        strata[r] = stratum1[k];
        for (int j = 0; j < p; ++j) {
          z(r,j) = z1(k,j);
        }

        if (i==0) {
          cumhaz[r] = haz1[i] * risk1[k];
        } else {
          cumhaz[r] = cumhaz[r-1] + haz1[i] * risk1[k];
        }
      }

      if (sefit) {
        std::vector<double> vhaz1 = subset(vhaz0, idx3);
        FlatMatrix ghaz1(m1,p);
        for (int j=0; j<p; ++j) {
          for (int i=0; i<m1; ++i) {
            ghaz1(i,j) = ghaz0(idx3[i],j);
          }
        }

        FlatMatrix a(m1,p);
        for (int j=0; j<p; ++j) {
          for (int i=0; i<m1; ++i) {
            int k = idx4[i];
            if (i==0) {
              a(i,j) = (haz1[i] * z1(k,j) - ghaz1(i,j)) * risk1[k];
            } else {
              a(i,j) = a(i-1,j) + (haz1[i] * z1(k,j) - ghaz1(i,j)) * risk1[k];
            }
          }
        }

        // calculate the first component of variance
        for (int i=0; i<m1; ++i) {
          int r = l + i;
          int k = idx4[i];
          if (i==0) {
            vcumhaz[r] = vhaz1[i] * risk1[k] * risk1[k];
          } else {
            vcumhaz[r] = vcumhaz[r-1] + vhaz1[i] * risk1[k] * risk1[k];
          }
        }

        // add the second component of variance
        for (int k=0; k<p; ++k) {
          for (int j=0; j<p; ++j) {
            for (int i=0; i<m1; ++i) {
              int r = l + i;
              vcumhaz[r] += a(i,j) * vbeta(j,k) * a(i,k);
            }
          }
        }
        
        for (int i=0; i<m1; ++i) {
          int r = l + i;
          secumhaz[r] = std::sqrt(vcumhaz[r]);
        }
      }

      l += m1;
    }
  }

  std::vector<double> surv(l);
  for (int i = 0; i < l; ++i) {
    surv[i] = std::exp(-cumhaz[i]);
  }

  DataFrameCpp result;
  result.push_back(std::move(time), "time");
  result.push_back(std::move(nrisk), "nrisk");
  result.push_back(std::move(nevent), "nevent");
  result.push_back(std::move(ncensor), "ncensor");
  result.push_back(std::move(cumhaz), "cumhaz");
  result.push_back(surv, "surv");

  if (sefit) {
    std::vector<double> sesurv(l);
    for (int i = 0; i < l; ++i) {
      sesurv[i] = surv[i] * secumhaz[i];
    }

    std::vector<double> lower(l), upper(l);
    for (int i = 0; i < l; ++i) {
      std::vector<double> ci = fsurvci(surv[i], sesurv[i], ct, zcrit);
      lower[i] = ci[0];
      upper[i] = ci[1];
    }

    result.push_back(std::move(sesurv), "sesurv");
    result.push_back(std::move(lower), "lower");
    result.push_back(std::move(upper), "upper");
    result.push_back(conflev, "conflev");
    result.push_back(ct, "conftype");
  }

  for (int j=0; j<p; ++j) {
    std::string zj = covariates[j];
    std::vector<double> u(l);
    for (int i=0; i<l; ++i) {
      u[i] = z(i,j);
    }
    result.push_back(u, zj);
  }

  if (has_stratum) {
    for (int i=0; i<p_stratum; ++i) {
      std::string s = stratum[i];
      if (u_stratum0.int_cols.count(s)) {
        auto v = u_stratum0.get<int>(s);
        subset_in_place(v, strata);
        result.push_back(std::move(v), s);
      } else if (u_stratum0.numeric_cols.count(s)) {
        auto v = u_stratum0.get<double>(s);
        subset_in_place(v, strata);
        result.push_back(std::move(v), s);
      } else if (u_stratum0.string_cols.count(s)) {
        auto v = u_stratum0.get<std::string>(s);
        subset_in_place(v, strata);
        result.push_back(std::move(v), s);
      } else {
        throw std::invalid_argument("Unsupported type for stratum variable: " + s);
      }
    }
  }

  if (has_id) {
    if (newdata.int_cols.count(id)) {
      auto v = subset(idwi, ids);
      result.push_back(std::move(v), id);
    } else if (newdata.numeric_cols.count(id)) {
      auto v = subset(idwn, ids);
      result.push_back(std::move(v), id);
    } else if (newdata.string_cols.count(id)) {
      auto v = subset(idwc, ids);
      result.push_back(std::move(v), id);
    } else {
      throw std::invalid_argument("incorrect type for the id variable in newdata");
    }
  }

  return result;
}
// 
// // schoenfeld residuals
// List f_ressch(int p, const NumericVector& par, void *ex) {
//   coxparams *param = (coxparams *) ex;
//   
//   const int nused = param->nused;
//   const int method = param->method;
//   
//   // Precompute exp(eta)
//   NumericVector exp_eta(nused);
//   for (int person = 0; person < nused; ++person) {
//     double val = param->offset[person];
//     for (int i = 0; i < p; ++i) {
//       val += par[i] * param->z(person, i);
//     }
//     exp_eta[person] = std::exp(val);
//   }
//   
//   NumericVector xbar(p);      // weighted mean covariate at this time
//   NumericVector a(p);         // s1(beta,k,t)
//   NumericVector a2(p);        // sum of w*exp(zbeta)*z for the deaths
//   double denom = 0.0;         // s0(beta,k,t)
//   double denom2 = 0.0;        // sum of weighted risks for the deaths
//   int ndead = 0;              // number of deaths at this time point
//   
//   int nevent = sum(param->event);   // total number of events
//   NumericMatrix resid(nevent, p);   // residual matrix
//   IntegerVector index(nevent);      // index of residuals
//   
//   int istrata = param->strata[0];
//   int i1 = 0; // index for removing out-of-risk subjects
//   int j = nevent;  // index the events in descending order
//   
//   // Loop through subjects
//   for (int person = 0; person < nused; ) {
//     if (param->strata[person] != istrata) { // hit a new stratum
//       istrata = param->strata[person]; // reset temporary variables
//       i1 = person;
//       denom = 0.0;
//       a.fill(0.0);
//     }
//     
//     const double dtime = param->tstop[person];
//     
//     // process all persons tied at this dtime
//     for (; person < nused && param->tstop[person] == dtime &&
//          param->strata[person] == istrata; ++person) {
//       
//       const double r = param->weight[person] * exp_eta[person];
//       
//       if (param->event[person] == 0) {
//         denom += r;
//         for (int i=0; i<p; ++i) {
//           a[i] += r * param->z(person,i);
//         }
//       } else {
//         j--;
//         resid(j,_) = param->z(person,_);
//         index[j] = person;
//         
//         ++ndead;
//         denom2 += r;
//         for (int i=0; i<p; ++i) {
//           a2[i] += r * param->z(person,i);
//         }
//       }
//     }
//     
//     // remove subjects no longer at risk
//     for (; i1 < nused; ++i1) {
//       const int p1 = param->order1[i1];
//       if (param->tstart[p1] < dtime || param->strata[p1] != istrata) break;
//       const double r = param->weight[p1] * exp_eta[p1];
//       denom -= r;
//       for (int i = 0; i < p; ++i) {
//         a[i] -= r * param->z(p1,i);
//       }
//     }
//     
//     // add to the main terms
//     if (ndead > 0) {
//       if (method == 0 || ndead == 1) {
//         denom += denom2;
//         for (int i = 0; i < p; ++i) {
//           a[i] += a2[i];
//           xbar[i] = a[i] / denom;
//         }
//       } else {
//         xbar.fill(0.0);
//         for (int k = 0; k < ndead; ++k) {
//           denom += denom2 / ndead;
//           for (int i = 0; i < p; ++i) {
//             a[i] += a2[i] / ndead;
//             xbar[i] += a[i] / denom;
//           }
//         }
//         xbar = xbar / ndead;
//       }
//       
//       for (int k = 0; k < ndead; ++k) {
//         for (int i = 0; i < p; ++i) {
//           resid(j+k,i) -= xbar[i];
//         }
//       }
//       
//       // reset for the next death time
//       ndead = 0;
//       denom2 = 0.0;
//       a2.fill(0.0);
//     }
//   }
//   
//   return List::create(
//     Named("resid") = resid,
//     Named("index") = index);
// }
// 
// 
// // [[Rcpp::export]]
// List residuals_phregcpp(const int p,
//                         const NumericVector& beta,
//                         const NumericMatrix& vbeta,
//                         const NumericVector& resmart,
//                         DataFrame data,
//                         const StringVector& stratum = "",
//                         const std::string time = "time",
//                         const std::string time2 = "",
//                         const std::string event = "event",
//                         const StringVector& covariates = "",
//                         const std::string weight = "",
//                         const std::string offset = "",
//                         const std::string id = "",
//                         const std::string ties = "efron",
//                         const std::string type = "schoenfeld",
//                         const bool collapse = false,
//                         const bool weighted = false) {
//   
//   int n = data.nrows();
//   
//   IntegerVector stratumn(n);
//   DataFrame u_stratum;
//   int p_stratum = static_cast<int>(stratum.size());
//   if (p_stratum == 1 && (stratum[0] == "" || stratum[0] == "none")) {
//     stratumn.fill(0);
//   } else {
//     List out = bygroup(data, stratum);
//     stratumn = out["index"];
//     u_stratum = out["lookup"];
//   }
//   
//   bool has_time = hasVariable(data, time);
//   if (!has_time) stop("data must contain the time variable");
//   NumericVector timenz = data[time];
//   NumericVector timen = clone(timenz);
//   if (is_true(any(timen < 0))) {
//     stop("time must be nonnegative for each observation");
//   }
//   
//   bool has_time2 = hasVariable(data, time2);
//   NumericVector time2n(n);
//   if (has_time2) {
//     NumericVector time2nz = data[time2];
//     time2n = clone(time2nz);
//     if (is_true(any(time2n <= timen))) {
//       stop("time2 must be greater than time for each observation");
//     }
//   }
//   
//   bool has_event = hasVariable(data, event);
//   if (!has_event) stop("data must contain the event variable");
//   IntegerVector eventnz = data[event];
//   IntegerVector eventn = clone(eventnz);
//   if (is_true(any((eventn != 1) & (eventn != 0)))) {
//     stop("event must be 1 or 0 for each observation");
//   }
//   
//   NumericMatrix zn(n,p);
//   if (p > 0) {
//     for (int j=0; j<p; ++j) {
//       String zj = covariates[j];
//       if (!hasVariable(data, zj)) {
//         stop("data must contain the variables in covariates");
//       }
//       NumericVector u = data[zj];
//       for (int i=0; i<n; ++i) {
//         zn(i,j) = u[i];
//       }
//     }
//   }
//   
//   bool has_weight = hasVariable(data, weight);
//   NumericVector weightn(n, 1.0);
//   if (has_weight) {
//     NumericVector weightnz = data[weight];
//     weightn = clone(weightnz);
//     if (is_true(any(weightn <= 0))) {
//       stop("weight must be greater than 0");
//     }
//   }
//   
//   bool has_offset = hasVariable(data, offset);
//   NumericVector offsetn(n);
//   if (has_offset) {
//     NumericVector offsetnz = data[offset];
//     offsetn = clone(offsetnz);
//   }
//   
//   // create the numeric id variable
//   bool has_id = hasVariable(data, id);
//   IntegerVector idn(n);
//   if (!has_id) {
//     idn = seq(0, n - 1);
//   } else {
//     SEXP col = data[id];
//     SEXPTYPE col_type = TYPEOF(col);
//     if (col_type == INTSXP) {
//       IntegerVector v = col;
//       IntegerVector w = unique(v);
//       w.sort();
//       idn = match(v, w) - 1;
//     } else if (col_type == REALSXP) {
//       NumericVector v = col;
//       NumericVector w = unique(v);
//       w.sort();
//       idn = match(v, w) - 1;
//     } else if (col_type == STRSXP) {
//       StringVector v = col;
//       StringVector w = unique(v);
//       w.sort();
//       idn = match(v, w) - 1;
//     } else {
//       stop("incorrect type for the id variable in the input data");
//     }
//   }
//   
//   std::string meth = ties;
//   std::for_each(meth.begin(), meth.end(), [](char & c) {
//     c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
//   });
//   
//   int method = meth == "efron" ? 1 : 0;
//   
//   // unify right censored data with counting process data
//   NumericVector tstartn(n), tstopn(n);
//   if (!has_time2) {
//     tstopn = timen;
//   } else {
//     tstartn = timen;
//     tstopn = time2n;
//   }
//   
//   // exclude observations with missing values
//   LogicalVector sub(n,1);
//   for (int i=0; i<n; ++i) {
//     if (stratumn[i] == NA_INTEGER || idn[i] == NA_INTEGER ||
//         std::isnan(tstartn[i]) || std::isnan(tstopn[i]) ||
//         eventn[i] == NA_INTEGER || std::isnan(weightn[i]) ||
//         std::isnan(offsetn[i])) {
//       sub[i] = 0;
//     }
//     for (int j=0; j<p; ++j) {
//       if (std::isnan(zn(i,j))) sub[i] = 0;
//     }
//   }
//   
//   IntegerVector order = which(sub);
//   stratumn = stratumn[order];
//   tstartn = tstartn[order];
//   tstopn = tstopn[order];
//   eventn = eventn[order];
//   weightn = weightn[order];
//   offsetn = offsetn[order];
//   idn = idn[order];
//   if (p > 0) zn = subset_matrix_by_row(zn, order);
//   n = sum(sub);
//   if (n == 0) stop("no observations left after removing missing values");
//   
//   // sort by stratum
//   IntegerVector order0 = seq(0, n-1);
//   std::sort(order0.begin(), order0.end(), [&](int i, int j) {
//     return stratumn[i] < stratumn[j];
//   });
//   
//   IntegerVector stratum1z = stratumn[order0];
//   NumericVector tstart1z = tstartn[order0];
//   NumericVector tstop1z = tstopn[order0];
//   IntegerVector event1z = eventn[order0];
//   
//   // locate the first observation within each stratum
//   IntegerVector istratum(1,0);
//   for (int i=1; i<n; ++i) {
//     if (stratum1z[i] != stratum1z[i-1]) {
//       istratum.push_back(i);
//     }
//   }
//   
//   int nstrata = static_cast<int>(istratum.size());
//   istratum.push_back(n);
//   
//   // ignore subjects not at risk for any event time
//   IntegerVector ignore1z(n);
//   for (int i=0; i<nstrata; ++i) {
//     IntegerVector q0 = Range(istratum[i], istratum[i+1]-1);
//     NumericVector tstart0 = tstart1z[q0];
//     NumericVector tstop0 = tstop1z[q0];
//     IntegerVector event0 = event1z[q0];
//     NumericVector etime = tstop0[event0==1];
//     etime = unique(etime);
//     etime.sort();
//     IntegerVector index1 = findInterval3(tstart0, etime, 0, 0, 0);
//     IntegerVector index2 = findInterval3(tstop0, etime, 0, 0, 0);
//     for (int j=istratum[i]; j<istratum[i+1]; ++j) {
//       int j0 = j - istratum[i];
//       if (index1[j0] == index2[j0]) { // no event in (tstart, tstop]
//         ignore1z[j] = 1;
//       } else {
//         ignore1z[j] = 0;
//       }
//     }
//   }
//   
//   IntegerVector ignore(n);
//   for (int i=0; i<n; ++i) {
//     ignore[order0[i]] = ignore1z[i];
//   }
//   
//   int nused = n - sum(ignore);
//   
//   order = seq(0, n-1);
//   IntegerVector idx(1,0);
//   int nids = n;
//   if (has_id) { // collapse over id
//     std::sort(order.begin(), order.end(), [&](int i, int j) {
//       return idn[i] < idn[j];
//     });
//     
//     IntegerVector id2 = idn[order];
//     for (int i=1; i<n; ++i) {
//       if (id2[i] != id2[i-1]) {
//         idx.push_back(i);
//       }
//     }
//     
//     nids = static_cast<int>(idx.size());
//     idx.push_back(n);
//   }
//   
//   
//   List result;
//   if (type == "martingale") {
//     NumericVector rr = clone(resmart);
//     if (weighted) rr = rr*weightn;
//     
//     if (collapse) { // collapse over id
//       NumericVector rr2(nids);
//       for (int i=0; i<nids; ++i) {
//         for (int j=idx[i]; j<idx[i+1]; ++j) {
//           rr2[i] += rr[order[j]];
//         }
//       }
//       
//       rr = rr2;
//     }
//     
//     result = List::create(Named("resid") = rr);
//   } else if (type == "deviance") {
//     NumericVector rr = clone(resmart);
//     IntegerVector status = clone(eventn);
//     int m = n;
//     
//     if (weighted) rr = rr*weightn;
//     
//     if (collapse) { // collapse over id
//       NumericVector rr2(nids);
//       IntegerVector status2(nids);
//       for (int i=0; i<nids; ++i) {
//         for (int j=idx[i]; j<idx[i+1]; ++j) {
//           int k = order[j];
//           rr2[i] += rr[k];
//           status2[i] += eventn[k];
//         }
//       }
//       
//       rr = rr2;
//       status = status2;
//       m = nids;
//     }
//     
//     for (int i=0; i<m; ++i) {
//       double temp = status[i] == 0 ? 0 : status[i]*std::log(status[i] - rr[i]);
//       rr[i] = ((rr[i]>0) - (rr[i]<0))*std::sqrt(-2*(rr[i] + temp));
//     }
//     
//     result = List::create(Named("resid") = rr);
//   } else if (p == 0) {
//     stop("covariates must be present for score and schoenfeld residuals");
//   } else {
//     NumericMatrix rr(n,p);
//     if (type == "score" || type == "dfbeta" || type == "dfbetas") {
//       // sort by stopping time in descending order within each stratum
//       IntegerVector order2 = seq(0, n-1);
//       std::sort(order2.begin(), order2.end(), [&](int i, int j) {
//         if (ignore[i] != ignore[j]) return ignore[i] < ignore[j];
//         if (stratumn[i] != stratumn[j]) return stratumn[i] < stratumn[j];
//         if (tstopn[i] != tstopn[j]) return tstopn[i] > tstopn[j];
//         return eventn[i] < eventn[j];
//       });
//       
//       IntegerVector stratum1 = stratumn[order2];
//       NumericVector tstart1 = tstartn[order2];
//       NumericVector tstop1 = tstopn[order2];
//       IntegerVector event1 = eventn[order2];
//       NumericVector weight1 = weightn[order2];
//       NumericVector offset1 = offsetn[order2];
//       IntegerVector ignore1 = ignore[order2];
//       NumericMatrix z1 = subset_matrix_by_row(zn, order2);
//       
//       // sort by starting time in descending order within each stratum
//       IntegerVector order1 = seq(0, n-1);
//       std::sort(order1.begin(), order1.end(), [&](int i, int j) {
//         if (ignore1[i] != ignore1[j]) return ignore1[i] < ignore1[j];
//         if (stratum1[i] != stratum1[j]) return stratum1[i] < stratum1[j];
//         return tstart1[i] > tstart1[j];
//       });
//       
//       coxparams param = {nused, stratum1, tstart1, tstop1, event1,
//                          weight1, offset1, z1, order1, method};
//       
//       NumericMatrix ressco = f_ressco_2(p, beta, &param);
//       
//       NumericMatrix score(n,p);
//       for (int i=0; i<n; ++i) {
//         score(order2[i],_) = ressco(i,_); // original order
//       }
//       
//       if (type == "dfbeta" || type == "dfbetas") {
//         for (int i=0; i<n; ++i) {
//           for (int k=0; k<p; ++k) {
//             for (int j=0; j<p; ++j) {
//               rr(i,k) += score(i,j)*vbeta(j,k);
//             }
//             if (type == "dfbetas") {
//               rr(i,k) /= std::sqrt(vbeta(k,k));
//             }
//           }
//         }
//       } else {
//         rr = score;
//       }
//       
//       if (weighted) {
//         for (int i=0; i<n; ++i) {
//           for (int k=0; k<p; ++k) {
//             rr(i,k) = rr(i,k)*weightn[i];
//           }
//         }
//       }
//       
//       if (collapse) { // collapse over id
//         NumericMatrix rr2(nids,p);
//         for (int i=0; i<nids; ++i) {
//           for (int k=0; k<p; ++k) {
//             for (int j=idx[i]; j<idx[i+1]; ++j) {
//               rr2(i,k) += rr(order[j],k);
//             }
//           }
//         }
//         
//         rr = rr2;
//       }
//       
//       result = List::create(Named("resid") = rr);
//     } else if (type == "schoenfeld" || type == "scaledsch") {
//       // sort by stopping time in descending order within each stratum
//       IntegerVector order2 = seq(0, n-1);
//       std::sort(order2.begin(), order2.end(), [&](int i, int j) {
//         if (ignore[i] != ignore[j]) return ignore[i] < ignore[j];
//         if (stratumn[i] != stratumn[j]) return stratumn[i] > stratumn[j];
//         if (tstopn[i] != tstopn[j]) return tstopn[i] > tstopn[j];
//         if (eventn[i] != eventn[j]) return eventn[i] < eventn[j];
//         return idn[i] > idn[j];
//       });
//       
//       IntegerVector stratum1 = stratumn[order2];
//       NumericVector tstart1 = tstartn[order2];
//       NumericVector tstop1 = tstopn[order2];
//       IntegerVector event1 = eventn[order2];
//       NumericVector weight1 = weightn[order2];
//       NumericVector offset1 = offsetn[order2];
//       IntegerVector id1 = idn[order2];
//       IntegerVector ignore1 = ignore[order2];
//       NumericMatrix z1 = subset_matrix_by_row(zn, order2);
//       
//       // sort by starting time in descending order within each stratum
//       IntegerVector order1 = seq(0, n-1);
//       std::sort(order1.begin(), order1.end(), [&](int i, int j) {
//         if (ignore1[i] != ignore1[j]) return ignore1[i] < ignore1[j];
//         if (stratum1[i] != stratum1[j]) return stratum1[i] > stratum1[j];
//         return tstart1[i] > tstart1[j];
//       });
//       
//       coxparams param = {nused, stratum1, tstart1, tstop1, event1,
//                          weight1, offset1, z1, order1, method};
//       
//       List out = f_ressch(p, beta, &param);
//       rr = as<NumericMatrix>(out["resid"]);
//       IntegerVector index = out["index"];
//       IntegerVector stratum2 = stratum1[index];
//       NumericVector time3 = tstop1[index];
//       int ndead = static_cast<int>(index.size());
//       
//       if (weighted) {
//         for (int i=0; i<ndead; ++i) {
//           for (int k=0; k<p; ++k) {
//             rr(i,k) = rr(i,k)*weightn[i];
//           }
//         }
//       }
//       
//       if (type == "scaledsch") {
//         NumericMatrix rr2(ndead,p);
//         for (int i=0; i<ndead; ++i) {
//           for (int k=0; k<p; ++k) {
//             for (int j=0; j<p; ++j) {
//               rr2(i,k) += rr(i,j)*vbeta(j,k);
//             }
//             rr2(i,k) = rr2(i,k)*ndead + beta[k];
//           }
//         }
//         rr = rr2;
//       }
//       
//       result = List::create(
//         Named("resid") = rr,
//         Named("time") = time3);
//       
//       IntegerVector stratum3 = unique(stratum2);
//       if (stratum3.size() > 1) {
//         List strata(p_stratum);
//         for (int i = 0; i < p_stratum; ++i) {
//           std::string s = as<std::string>(stratum[i]);
//           SEXP col = u_stratum[s];
//           SEXPTYPE col_type = TYPEOF(col);
//           if (col_type == INTSXP) {
//             IntegerVector v = col;
//             strata[i] = v[stratum2];
//           } else if (col_type == REALSXP) {
//             NumericVector v = col;
//             strata[i] = v[stratum2];
//           } else if (col_type == STRSXP) {
//             StringVector v = col;
//             strata[i] = v[stratum2];
//           } else {
//             stop("Unsupported type for stratum variable: " + s);
//           }
//         }
//         strata.attr("names") = stratum;
//         result.push_back(as<DataFrame>(strata), "strata");
//       }
//     } else {
//       stop("unknown type of residuals");
//     }
//   }
//   
//   return result;
// }
// 
// 
// // function for individual contributions to score and information matrix
// // for the Cox model
// List f_der_i_2(int p, const NumericVector& par, void* ex) {
//   coxparams* param = (coxparams*) ex;
//   
//   const int nused = param->nused;
//   const int method = param->method;
//   
//   // Precompute eta and exp(eta)
//   NumericVector eta(nused);
//   NumericVector exp_eta(nused);
//   for (int person = 0; person < nused; ++person) {
//     double val = param->offset[person];
//     for (int i = 0; i < p; ++i) {
//       val += par[i] * param->z(person,i);
//     }
//     eta[person] = val;
//     exp_eta[person] = std::exp(val);
//   }
//   
//   NumericMatrix u(nused,p);         // score vector for each individual
//   NumericMatrix imat(nused,p*p);    // information matrix for each individual
//   NumericVector a(p);               // s1(beta,k,t)
//   NumericVector a2(p);              // sum of w*exp(zbeta)*z for the deaths
//   NumericMatrix cmat(p,p);          // s2(beta,k,t)
//   NumericMatrix cmat2(p,p);         // sum of w*exp(zbeta)*z*z' for the deaths
//   double denom = 0.0;               // s0(beta,k,t)
//   double denom2 = 0.0;              // sum of weighted risks for deaths
//   int ndead = 0;                    // number of deaths at this time point
//   
//   
//   int istrata = param->strata[0];
//   int i1 = 0; // index for removing out-of-risk subjects
//   
//   // Loop through subjects
//   for (int person = 0; person < nused; ) {
//     // Reset when entering a new stratum
//     if (param->strata[person] != istrata) {
//       istrata = param->strata[person];
//       i1 = person;
//       denom = 0.0;
//       
//       for (int i = 0; i < p; ++i) {
//         a[i] = 0.0;
//         for (int j = 0; j <= i; ++j) {
//           cmat(i,j) = 0.0;
//         }
//       }
//     }
//     
//     const double dtime = param->tstop[person];
//     
//     // Process all persons tied at this dtime
//     int person1 = person;
//     
//     for (; person < nused && param->tstop[person] == dtime &&
//          param->strata[person] == istrata; ++person) {
//       
//       const double w = param->weight[person];
//       const double r = w * exp_eta[person];
//       
//       if (param->event[person] == 0) {
//         denom += r;
//         
//         for (int i = 0; i < p; ++i) {
//           const double zi = param->z(person,i);
//           a[i] += r * zi;
//           for (int j = 0; j <= i; ++j) {
//             const double zj = param->z(person,j);
//             cmat(i,j) += r * zi * zj;
//           }
//         }
//       } else {
//         ++ndead;
//         denom2 += r;
//         for (int i = 0; i < p; ++i) {
//           const double zi = param->z(person,i);
//           a2[i] += r * zi;
//           u(person,i) = w * zi;
//           for (int j = 0; j <= i; ++j) {
//             const double zj = param->z(person,j);
//             cmat2(i,j) += r * zi * zj;
//           }
//         }
//       }
//     }
//     
//     // Remove subjects leaving risk set
//     for (; i1 < nused; ++i1) {
//       const int p1 = param->order1[i1];
//       if (param->tstart[p1] < dtime || param->strata[p1] != istrata) break;
//       
//       const double r = param->weight[p1] * exp_eta[p1];
//       denom -= r;
//       for (int i = 0; i < p; ++i) {
//         const double zi = param->z(p1,i);
//         a[i] -= r * zi;
//         for (int j = 0; j <= i; ++j) {
//           const double zj = param->z(p1,j);
//           cmat(i,j) -= r * zi * zj;
//         }
//       }
//     }
//     
//     // Add contributions for deaths at this time
//     if (ndead > 0) {
//       if (method == 0 || ndead == 1) { // Breslow or single event
//         denom += denom2;
//         for (int i = 0; i < p; ++i) {
//           a[i] += a2[i];
//           const double xbar = a[i] / denom;
//           for (int human = person1; human < person; ++human) {
//             if (param->event[human] == 1)
//               u(human,i) -= param->weight[human] * xbar;
//           }
//           
//           for (int j = 0; j <= i; ++j) {
//             cmat(i,j) += cmat2(i,j);
//             
//             for (int human = person1; human < person; ++human) {
//               if (param->event[human] == 1)
//                 imat(human, i*p+j) = param->weight[human] *
//                   (cmat(i,j) - xbar * a[j]) / denom;
//             }
//           }
//         }
//       } else { // Efron method
//         const double increment = denom2 / ndead;
//         for (int l = 0; l < ndead; ++l) {
//           denom += increment;
//           for (int i = 0; i < p; ++i) {
//             a[i] += a2[i] / ndead;
//             const double xbar = a[i] / denom;
//             for (int human = person1; human < person; ++human) {
//               if (param->event[human] == 1)
//                 u(human,i) -= param->weight[human]/ndead * xbar;
//             }
//             
//             for (int j = 0; j <= i; ++j) {
//               cmat(i,j) += cmat2(i,j) / ndead;
//               for (int human = person1; human < person; ++human) {
//                 if (param->event[human] == 1)
//                   imat(human,i*p+j) +=  param->weight[human]/ndead *
//                     (cmat(i,j) - xbar * a[j]) / denom;
//               }
//             }
//           }
//         }
//       }
//       
//       // Reset after processing deaths
//       ndead = 0;
//       denom2 = 0.0;
//       for (int i = 0; i < p; ++i) {
//         a2[i] = 0;
//         for (int j = 0; j <= i; ++j) {
//           cmat2(i,j) = 0;
//         }
//       }
//     }
//   }
//   
//   // fill the symmetric elements of the information matrix
//   for (int person = 0; person < nused; ++person)
//     for (int i = 0; i < p - 1; ++i)
//       for (int j = i+1; j < p; ++j)
//         imat(person,i*p+j) = imat(person,j*p+i);
//   
//   List result = List::create(
//     _["score_i"] = u,
//     _["imat_i"] = imat
//   );
//   
//   return result;
// }
// 
// 
// // [[Rcpp::export]]
// List assess_phregcpp(const int p,
//                      const NumericVector& beta,
//                      const NumericMatrix& vbeta,
//                      DataFrame data,
//                      const StringVector& stratum = "",
//                      const std::string time = "time",
//                      const std::string time2 = "",
//                      const std::string event = "event",
//                      const StringVector& covariates = "",
//                      const std::string weight = "",
//                      const std::string offset = "",
//                      const std::string ties = "efron",
//                      const int resample = 1000,
//                      const std::uint32_t seed) {
//   
//   td::mt19937_64 rng(seed);                      // choose generator and seed
//   std::normal_distribution<double> dist(0.0, 1.0);
//   
//   if (seed != NA_INTEGER) set_seed(seed);
//   
//   if (p <= 0) {
//     stop("covariates must be present to test proportional hazards");
//   }
//   
//   int n = data.nrows();
//   
//   IntegerVector stratumn(n);
//   int p_stratum = static_cast<int>(stratum.size());
//   if (p_stratum == 1 && (stratum[0] == "" || stratum[0] == "none")) {
//     stratumn.fill(0);
//   } else {
//     List out = bygroup(data, stratum);
//     stratumn = out["index"];
//   }
//   
//   bool has_time = hasVariable(data, time);
//   if (!has_time) stop("data must contain the time variable");
//   NumericVector timenz = data[time];
//   NumericVector timen = clone(timenz);
//   if (is_true(any(timen < 0))) {
//     stop("time must be nonnegative for each observation");
//   }
//   
//   bool has_time2 = hasVariable(data, time2);
//   NumericVector time2n(n);
//   if (has_time2) {
//     NumericVector time2nz = data[time2];
//     time2n = clone(time2nz);
//     if (is_true(any(time2n <= timen))) {
//       stop("time2 must be greater than time for each observation");
//     }
//   }
//   
//   bool has_event = hasVariable(data, event);
//   if (!has_event) stop("data must contain the event variable");
//   IntegerVector eventnz = data[event];
//   IntegerVector eventn = clone(eventnz);
//   if (is_true(any((eventn != 1) & (eventn != 0)))) {
//     stop("event must be 1 or 0 for each observation");
//   }
//   
//   NumericMatrix zn(n,p);
//   if (p > 0) {
//     for (int j=0; j<p; ++j) {
//       String zj = covariates[j];
//       if (!hasVariable(data, zj)) {
//         stop("data must contain the variables in covariates");
//       }
//       NumericVector u = data[zj];
//       for (int i=0; i<n; ++i) {
//         zn(i,j) = u[i];
//       }
//     }
//   }
//   
//   bool has_weight = hasVariable(data, weight);
//   NumericVector weightn(n, 1.0);
//   if (has_weight) {
//     NumericVector weightnz = data[weight];
//     weightn = clone(weightnz);
//     if (is_true(any(weightn <= 0))) {
//       stop("weight must be greater than 0");
//     }
//   }
//   
//   bool has_offset = hasVariable(data, offset);
//   NumericVector offsetn(n);
//   if (has_offset) {
//     NumericVector offsetnz = data[offset];
//     offsetn = clone(offsetnz);
//   }
//   
//   std::string meth = ties;
//   std::for_each(meth.begin(), meth.end(), [](char & c) {
//     c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
//   });
//   
//   int method = meth == "efron" ? 1 : 0;
//   
//   // unify right censored data with counting process data
//   NumericVector tstartn(n), tstopn(n);
//   if (!has_time2) {
//     tstopn = timen;
//   } else {
//     tstartn = timen;
//     tstopn = time2n;
//   }
//   
//   // exclude observations with missing values
//   LogicalVector sub(n,1);
//   for (int i=0; i<n; ++i) {
//     if (stratumn[i] == NA_INTEGER || std::isnan(tstartn[i]) ||
//         std::isnan(tstopn[i]) || eventn[i] == NA_INTEGER ||
//         std::isnan(weightn[i]) || std::isnan(offsetn[i])) {
//       sub[i] = 0;
//     }
//     for (int j=0; j<p; ++j) {
//       if (std::isnan(zn(i,j))) sub[i] = 0;
//     }
//   }
//   
//   IntegerVector order = which(sub);
//   stratumn = stratumn[order];
//   tstartn = tstartn[order];
//   tstopn = tstopn[order];
//   eventn = eventn[order];
//   weightn = weightn[order];
//   offsetn = offsetn[order];
//   zn = subset_matrix_by_row(zn, order);
//   n = sum(sub);
//   if (n == 0) stop("no observations left after removing missing values");
//   
//   // sort by stratum
//   IntegerVector order0 = seq(0, n-1);
//   std::sort(order0.begin(), order0.end(), [&](int i, int j) {
//     return stratumn[i] < stratumn[j];
//   });
//   
//   IntegerVector stratum1z = stratumn[order0];
//   NumericVector tstart1z = tstartn[order0];
//   NumericVector tstop1z = tstopn[order0];
//   IntegerVector event1z = eventn[order0];
//   
//   // locate the first observation within each stratum
//   IntegerVector istratum(1,0);
//   for (int i=1; i<n; ++i) {
//     if (stratum1z[i] != stratum1z[i-1]) {
//       istratum.push_back(i);
//     }
//   }
//   
//   int nstrata = static_cast<int>(istratum.size());
//   istratum.push_back(n);
//   
//   // ignore subjects not at risk for any event time
//   IntegerVector ignore1z(n);
//   for (int i=0; i<nstrata; ++i) {
//     IntegerVector q0 = Range(istratum[i], istratum[i+1]-1);
//     NumericVector tstart0 = tstart1z[q0];
//     NumericVector tstop0 = tstop1z[q0];
//     IntegerVector event0 = event1z[q0];
//     NumericVector etime = tstop0[event0==1];
//     etime = unique(etime);
//     etime.sort();
//     IntegerVector index1 = findInterval3(tstart0, etime, 0, 0, 0);
//     IntegerVector index2 = findInterval3(tstop0, etime, 0, 0, 0);
//     for (int j=istratum[i]; j<istratum[i+1]; ++j) {
//       int j0 = j - istratum[i];
//       if (index1[j0] == index2[j0]) { // no event in (tstart, tstop]
//         ignore1z[j] = 1;
//       } else {
//         ignore1z[j] = 0;
//       }
//     }
//   }
//   
//   IntegerVector ignoren(n); // final ignore vector
//   for (int i=0; i<n; ++i) {
//     ignoren[order0[i]] = ignore1z[i];
//   }
//   
//   int nused = n - sum(ignoren);
//   
//   // sort by stopping time in descending order within each stratum
//   IntegerVector order1 = seq(0, n-1);
//   std::sort(order1.begin(), order1.end(), [&](int i, int j) {
//     if (ignoren[i] != ignoren[j]) return ignoren[i] < ignoren[j];
//     if (stratumn[i] != stratumn[j]) return stratumn[i] < stratumn[j];
//     if (tstopn[i] != tstopn[j]) return tstopn[i] > tstopn[j];
//     return eventn[i] < eventn[j];
//   });
//   
//   IntegerVector stratum1 = stratumn[order1];
//   NumericVector tstart1 = tstartn[order1];
//   NumericVector tstop1 = tstopn[order1];
//   IntegerVector event1 = eventn[order1];
//   NumericVector weight1 = weightn[order1];
//   NumericVector offset1 = offsetn[order1];
//   IntegerVector ignore1 = ignoren[order1];
//   NumericMatrix z1 = subset_matrix_by_row(zn, order1);
//   
//   // sort by starting time in descending order within each stratum
//   IntegerVector order10 = seq(0, n-1);
//   std::sort(order10.begin(), order10.end(), [&](int i, int j) {
//     if (ignore1[i] != ignore1[j]) return ignore1[i] < ignore1[j];
//     if (stratum1[i] != stratum1[j]) return stratum1[i] < stratum1[j];
//     return tstart1[i] > tstart1[j];
//   });
//   
//   coxparams param = {nused, stratum1, tstart1, tstop1, event1,
//                      weight1, offset1, z1, order10, method};
//   
//   List out = f_der_i_2(p, beta, &param);
//   NumericMatrix score_i = as<NumericMatrix>(out["score_i"]);
//   NumericMatrix imat_i = as<NumericMatrix>(out["imat_i"]);
//   
//   // order data by ascending tstop
//   IntegerVector order2 = seq(0, n-1);
//   std::sort(order2.begin(), order2.end(), [&](int i, int j) {
//     return tstop1[i] < tstop1[j];
//   });
//   
//   NumericVector tstop2 = tstop1[order2];
//   IntegerVector event2 = event1[order2];
//   NumericMatrix score_i2 = subset_matrix_by_row(score_i, order2);
//   NumericMatrix imat_i2 = subset_matrix_by_row(imat_i, order2);
//   
//   // only consider event times
//   IntegerVector q = which(event2 == 1);
//   NumericVector tstop3 = tstop2[q];
//   NumericMatrix score_i3 = subset_matrix_by_row(score_i2, q);
//   NumericMatrix imat_i3 = subset_matrix_by_row(imat_i2, q);
//   int nq = static_cast<int>(q.size());
//   
//   // identify unique event times
//   NumericVector t = unique(tstop3);
//   t.sort();
//   
//   int nt = static_cast<int>(t.size());
//   NumericMatrix score_i4(nt, p);
//   NumericMatrix imat_i4(nt, p*p);
//   
//   tstop3.push_front(-1.0); // to facilitate the loop
//   int k = -1;
//   for (int i=1; i<=nq; ++i) {
//     if (tstop3[i] != tstop3[i-1]) k++;
//     
//     for (int j=0; j<p; ++j) {
//       score_i4(k,j) += score_i3(i-1,j);
//       for (int l=0; l<p; ++l) {
//         imat_i4(k,j*p+l) += imat_i3(i-1,j*p+l);
//       }
//     }
//   }
//   
//   // cumulative sums of score and information processes
//   for (int i=1; i<nt; ++i) {
//     for (int j=0; j<p; ++j) {
//       score_i4(i,j) += score_i4(i-1,j);
//       for (int l=0; l<p; ++l) {
//         imat_i4(i,j*p+l) += imat_i4(i-1,j*p+l);
//       }
//     }
//   }
//   
//   // standardize the score processes
//   for (int i=0; i<nt; ++i) {
//     for (int j=0; j<p; ++j) {
//       score_i4(i,j) *= std::sqrt(vbeta(j,j));
//     }
//   }
//   
//   // resampling to approximate the null distribution
//   List G_list(resample);
//   for (int r=0; r<resample; ++r) {
//     NumericVector G(nq);
//     for (int i=0; i<nq; ++i) {
//       G[i] = dist(rng);
//     }
//     G_list[r] = G;
//   }
//   
//   // sum(u*G)
//   List U_list(resample);
//   for (int r=0; r<resample; ++r) {
//     NumericVector U(p);
//     NumericVector G = G_list[r];
//     for (int i=0; i<nq; ++i) {
//       for (int j=0; j<p; ++j) {
//         U[j] += score_i3(i,j)*G[i];
//       }
//     }
//     U_list[r] = U;
//   }
//   
//   // sum(I(X <= t)*U*G)
//   List score_i_list(resample);
//   for (int r=0; r<resample; ++r) {
//     NumericMatrix score_i4(nt, p);
//     NumericVector G = G_list[r];
//     
//     int k = -1;
//     for (int i=1; i<=nq; ++i) {
//       if (tstop3[i] != tstop3[i-1]) k++;
//       
//       for (int j=0; j<p; ++j) {
//         score_i4(k,j) += score_i3(i-1,j)*G[i-1];
//       }
//     }
//     
//     // cumulative sums of score processes
//     for (int i=1; i<nt; ++i) {
//       for (int j=0; j<p; ++j) {
//         score_i4(i,j) += score_i4(i-1,j);
//       }
//     }
//     
//     score_i_list[r] = score_i4;
//   }
//   
//   // piece together the result
//   List score_i_list_2(resample);
//   for (int r=0; r<resample; ++r) {
//     NumericMatrix score_i4 = score_i_list[r];
//     NumericMatrix score_i4_2(nt, p);
//     NumericVector U = U_list[r];
//     for (int i=0; i<nt; ++i) {
//       for (int j=0; j<p; ++j) {
//         score_i4_2(i,j) = score_i4(i,j);
//         for (int l=0; l<p; ++l) {
//           for (int s=0; s<p; ++s) {
//             score_i4_2(i,j) -= imat_i4(i,j*p+l)*vbeta(l,s)*U[s];
//           }
//         }
//       }
//     }
//     score_i_list_2[r] = score_i4_2;
//   }
//   
//   // standardize the resampled score processes
//   for (int r=0; r<resample; ++r) {
//     NumericMatrix score_i4 = score_i_list_2[r];
//     for (int i=0; i<nt; ++i) {
//       for (int j=0; j<p; ++j) {
//         score_i4(i,j) *= std::sqrt(vbeta(j,j));
//       }
//     }
//     score_i_list_2[r] = score_i4;
//   }
//   
//   // calculate individual p-values for proportional hazards assumption test
//   NumericVector p_values(p+1);
//   NumericVector maxabs(p+1);
//   for (int j=0; j<p; ++j) {
//     NumericVector obs(nt);
//     for (int i=0; i<nt; ++i) {
//       obs[i] = fabs(score_i4(i,j));
//     }
//     double tobs = max(obs);
//     maxabs[j] = tobs;
//     
//     int count = 0;
//     for (int r=0; r<resample; ++r) {
//       NumericMatrix score_i4 = score_i_list_2[r];
//       NumericVector sim(nt);
//       for (int i=0; i<nt; ++i) {
//         sim[i] = fabs(score_i4(i,j));
//       }
//       double tsim = max(sim);
//       if (tsim >= tobs) count++;
//     }
//     p_values[j] = static_cast<double>(count) / static_cast<double>(resample);
//   }
//   
//   // p-value for the global test
//   NumericVector obs(nt);
//   for (int i=0; i<nt; ++i) {
//     for (int j=0; j<p; ++j) {
//       obs[i] += fabs(score_i4(i,j));
//     }
//   }
//   double tobs = max(obs);
//   
//   int count = 0;
//   for (int r=0; r<resample; ++r) {
//     NumericMatrix score_i4 = score_i_list_2[r];
//     NumericVector sim(nt);
//     for (int i=0; i<nt; ++i) {
//       for (int j=0; j<p; ++j) {
//         sim[i] += fabs(score_i4(i,j));
//       }
//     }
//     double tsim = max(sim);
//     if (tsim >= tobs) count++;
//   }
//   p_values[p] = static_cast<double>(count) / static_cast<double>(resample);
//   
//   maxabs[p] = tobs;
//   List result = List::create(
//     Named("time") = t,
//     Named("score_t") = score_i4,
//     Named("score_t_list") = score_i_list_2,
//     Named("max_abs_value") = maxabs,
//     Named("p_value") = p_values
//   );
//   
//   return result;
// }
// 
// 
// // [[Rcpp::export]]
// List zph_phregcpp(const int p,
//                   const NumericVector& beta,
//                   const NumericMatrix& vbeta,
//                   const NumericVector& resmart,
//                   DataFrame data,
//                   const StringVector& stratum = "",
//                   const std::string time = "time",
//                   const std::string time2 = "",
//                   const std::string event = "event",
//                   const StringVector& covariates = "",
//                   const std::string weight = "",
//                   const std::string offset = "",
//                   const std::string ties = "efron",
//                   const std::string transform= "km") {
//   
//   if (p <= 0) {
//     stop("covariates must be present to test proportional hazards");
//   }
//   
//   int n = data.nrows();
//   
//   bool has_stratum;
//   IntegerVector stratumn(n);
//   int p_stratum = static_cast<int>(stratum.size());
//   if (p_stratum == 1 && (stratum[0] == "" || stratum[0] == "none")) {
//     has_stratum = false;
//     stratumn.fill(0);
//   } else {
//     has_stratum = true;
//     List out = bygroup(data, stratum);
//     stratumn = out["index"];
//   }
//   
//   bool has_time = hasVariable(data, time);
//   if (!has_time) stop("data must contain the time variable");
//   NumericVector timenz = data[time];
//   NumericVector timen = clone(timenz);
//   if (is_true(any(timen < 0))) {
//     stop("time must be nonnegative for each observation");
//   }
//   
//   bool has_time2 = hasVariable(data, time2);
//   NumericVector time2n(n);
//   if (has_time2) {
//     NumericVector time2nz = data[time2];
//     time2n = clone(time2nz);
//     if (is_true(any(time2n <= timen))) {
//       stop("time2 must be greater than time for each observation");
//     }
//   }
//   
//   bool has_event = hasVariable(data, event);
//   if (!has_event) stop("data must contain the event variable");
//   IntegerVector eventnz = data[event];
//   IntegerVector eventn = clone(eventnz);
//   if (is_true(any((eventn != 1) & (eventn != 0)))) {
//     stop("event must be 1 or 0 for each observation");
//   }
//   
//   NumericMatrix zn(n,p);
//   if (p > 0) {
//     for (int j=0; j<p; ++j) {
//       String zj = covariates[j];
//       if (!hasVariable(data, zj)) {
//         stop("data must contain the variables in covariates");
//       }
//       NumericVector u = data[zj];
//       for (int i=0; i<n; ++i) {
//         zn(i,j) = u[i];
//       }
//     }
//   }
//   
//   bool has_weight = hasVariable(data, weight);
//   NumericVector weightn(n, 1.0);
//   if (has_weight) {
//     NumericVector weightnz = data[weight];
//     weightn = clone(weightnz);
//     if (is_true(any(weightn <= 0))) {
//       stop("weight must be greater than 0");
//     }
//   }
//   
//   bool has_offset = hasVariable(data, offset);
//   NumericVector offsetn(n);
//   if (has_offset) {
//     NumericVector offsetnz = data[offset];
//     offsetn = clone(offsetnz);
//   }
//   
//   std::string meth = ties;
//   std::for_each(meth.begin(), meth.end(), [](char & c) {
//     c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
//   });
//   
//   int method = meth == "efron" ? 1 : 0;
//   
//   // unify right censored data with counting process data
//   NumericVector tstartn(n), tstopn(n);
//   if (!has_time2) {
//     tstopn = timen;
//   } else {
//     tstartn = timen;
//     tstopn = time2n;
//   }
//   
//   // exclude observations with missing values
//   LogicalVector sub(n,1);
//   for (int i=0; i<n; ++i) {
//     if (stratumn[i] == NA_INTEGER || std::isnan(tstartn[i]) ||
//         std::isnan(tstopn[i]) || eventn[i] == NA_INTEGER ||
//         std::isnan(weightn[i]) || std::isnan(offsetn[i])) {
//       sub[i] = 0;
//     }
//     for (int j=0; j<p; ++j) {
//       if (std::isnan(zn(i,j))) sub[i] = 0;
//     }
//   }
//   
//   IntegerVector order = which(sub);
//   stratumn = stratumn[order];
//   tstartn = tstartn[order];
//   tstopn = tstopn[order];
//   eventn = eventn[order];
//   weightn = weightn[order];
//   offsetn = offsetn[order];
//   zn = subset_matrix_by_row(zn, order);
//   n = sum(sub);
//   if (n == 0) stop("no observations left after removing missing values");
//   
//   // sort by stratum
//   IntegerVector order0 = seq(0, n-1);
//   std::sort(order0.begin(), order0.end(), [&](int i, int j) {
//     return stratumn[i] < stratumn[j];
//   });
//   
//   IntegerVector stratum1z = stratumn[order0];
//   NumericVector tstart1z = tstartn[order0];
//   NumericVector tstop1z = tstopn[order0];
//   IntegerVector event1z = eventn[order0];
//   
//   // locate the first observation within each stratum
//   IntegerVector istratum(1,0);
//   for (int i=1; i<n; ++i) {
//     if (stratum1z[i] != stratum1z[i-1]) {
//       istratum.push_back(i);
//     }
//   }
//   
//   int nstrata = static_cast<int>(istratum.size());
//   istratum.push_back(n);
//   
//   // ignore subjects not at risk for any event time
//   IntegerVector ignore1z(n);
//   for (int i=0; i<nstrata; ++i) {
//     IntegerVector q0 = Range(istratum[i], istratum[i+1]-1);
//     NumericVector tstart0 = tstart1z[q0];
//     NumericVector tstop0 = tstop1z[q0];
//     IntegerVector event0 = event1z[q0];
//     NumericVector etime = tstop0[event0==1];
//     etime = unique(etime);
//     etime.sort();
//     IntegerVector index1 = findInterval3(tstart0, etime, 0, 0, 0);
//     IntegerVector index2 = findInterval3(tstop0, etime, 0, 0, 0);
//     for (int j=istratum[i]; j<istratum[i+1]; ++j) {
//       int j0 = j - istratum[i];
//       if (index1[j0] == index2[j0]) { // no event in (tstart, tstop]
//         ignore1z[j] = 1;
//       } else {
//         ignore1z[j] = 0;
//       }
//     }
//   }
//   
//   IntegerVector ignoren(n); // final ignore vector
//   for (int i=0; i<n; ++i) {
//     ignoren[order0[i]] = ignore1z[i];
//   }
//   
//   int nused = n - sum(ignoren);
//   
//   // sort by stopping time in descending order within each stratum
//   IntegerVector order1 = seq(0, n-1);
//   std::sort(order1.begin(), order1.end(), [&](int i, int j) {
//     if (ignoren[i] != ignoren[j]) return ignoren[i] < ignoren[j];
//     if (stratumn[i] != stratumn[j]) return stratumn[i] < stratumn[j];
//     if (tstopn[i] != tstopn[j]) return tstopn[i] > tstopn[j];
//     return eventn[i] < eventn[j];
//   });
//   
//   IntegerVector stratum1 = stratumn[order1];
//   NumericVector tstart1 = tstartn[order1];
//   NumericVector tstop1 = tstopn[order1];
//   IntegerVector event1 = eventn[order1];
//   NumericVector weight1 = weightn[order1];
//   NumericVector offset1 = offsetn[order1];
//   IntegerVector ignore1 = ignoren[order1];
//   NumericMatrix z1 = subset_matrix_by_row(zn, order1);
//   
//   // sort by starting time in descending order within each stratum
//   IntegerVector order10 = seq(0, n-1);
//   std::sort(order10.begin(), order10.end(), [&](int i, int j) {
//     if (ignore1[i] != ignore1[j]) return ignore1[i] < ignore1[j];
//     if (stratum1[i] != stratum1[j]) return stratum1[i] < stratum1[j];
//     return tstart1[i] > tstart1[j];
//   });
//   
//   coxparams param = {nused, stratum1, tstart1, tstop1, event1,
//                      weight1, offset1, z1, order10, method};
//   
//   List out = f_der_i_2(p, beta, &param);
//   NumericMatrix score_i = as<NumericMatrix>(out["score_i"]);
//   NumericMatrix imat_i = as<NumericMatrix>(out["imat_i"]);
//   
//   // order data by ascending tstop
//   IntegerVector order2 = seq(0, n-1);
//   std::sort(order2.begin(), order2.end(), [&](int i, int j) {
//     return tstop1[i] < tstop1[j];
//   });
//   
//   IntegerVector stratum2 = stratum1[order2];
//   NumericVector tstop2 = tstop1[order2];
//   IntegerVector event2 = event1[order2];
//   NumericMatrix score_i2 = subset_matrix_by_row(score_i, order2);
//   NumericMatrix imat_i2 = subset_matrix_by_row(imat_i, order2);
//   
//   // only consider event times
//   IntegerVector q = which(event2 == 1);
//   IntegerVector stratum3 = stratum2[q];
//   NumericVector tstop3 = tstop2[q];
//   NumericMatrix score_i3 = subset_matrix_by_row(score_i2, q);
//   NumericMatrix imat_i3 = subset_matrix_by_row(imat_i2, q);
//   int nq = static_cast<int>(q.size());
//   
//   // identify unique event times
//   NumericVector t = unique(tstop3);
//   t.sort();
//   
//   int nt = static_cast<int>(t.size());
//   NumericMatrix score_i4(nt, p);
//   NumericMatrix imat_i4(nt, p*p);
//   
//   NumericVector tstop3a = clone(tstop3);
//   tstop3a.push_front(-1.0); // to facilitate the loop
//   int k = -1;
//   for (int i=1; i<=nq; ++i) {
//     if (tstop3a[i] != tstop3a[i-1]) k++;
//     
//     for (int j=0; j<p; ++j) {
//       score_i4(k,j) += score_i3(i-1,j);
//       for (int l=0; l<p; ++l) {
//         imat_i4(k,j*p+l) += imat_i3(i-1,j*p+l);
//       }
//     }
//   }
//   
//   // transformed time points
//   NumericVector g(nt);
//   if (transform == "identity") {
//     for (int i=0; i<nt; ++i) {
//       g[i] = t[i];
//     }
//   } else if (transform == "log") {
//     for (int i=0; i<nt; ++i) {
//       g[i] = std::log(t[i]);
//     }
//   } else if (transform == "rank") {
//     for (int i=0; i<nt; ++i) {
//       g[i] = static_cast<double>(i+1);
//     }
//   } else if (transform == "km") {
//     DataFrame temp = DataFrame::create(
//       Named("tstart") = tstart1,
//       Named("tstop") = tstop1,
//       Named("event") = event1
//     );
//     
//     DataFrame df_km = kmest(temp, "", "", "tstart", "tstop", "event",
//                             "", "none", 0.95, false);
//     
//     NumericVector surv = df_km["surv"];
//     surv.push_front(1.0); // for left-continuous step function
//     for (int i=0; i<nt; ++i) {
//       g[i] = 1.0 - surv[i];
//     }
//   } else {
//     Function f(transform);
//     for (int i=0; i<nt; ++i) {
//       g[i] = as<double>(f(t[i]));
//     }
//   }
//   
//   // score for theta
//   NumericVector u_theta(p);
//   for (int i=0; i<nt; ++i) {
//     for (int j=0; j<p; ++j) {
//       u_theta[j] += score_i4(i,j)*g[i];
//     }
//   }
//   
//   // covariance between u_theta and u_beta
//   NumericMatrix imat_theta_beta(p, p);
//   for (int i=0; i<nt; ++i) {
//     for (int j=0; j<p; ++j) {
//       for (int l=0; l<p; ++l) {
//         imat_theta_beta(j,l) += imat_i4(i,j*p+l)*g[i];
//       }
//     }
//   }
//   
//   // covariance for u_theta
//   NumericMatrix imat_theta(p, p);
//   for (int i=0; i<nt; ++i) {
//     for (int j=0; j<p; ++j) {
//       for (int l=0; l<p; ++l) {
//         imat_theta(j,l) += imat_i4(i,j*p+l)*g[i]*g[i];
//       }
//     }
//   }
//   
//   // conditional variance for u_theta given u_beta
//   NumericMatrix vtheta(p, p);
//   for (int j=0; j<p; ++j) {
//     for (int l=0; l<p; ++l) {
//       vtheta(j,l) = imat_theta(j,l);
//       for (int s=0; s<p; ++s) {
//         for (int r=0; r<p; ++r) {
//           vtheta(j,l) -= imat_theta_beta(j,s)*vbeta(s,r)*imat_theta_beta(l,r);
//         }
//       }
//     }
//   }
//   
//   // individual score test for theta = 0
//   NumericVector score_test(p);
//   for (int j=0; j<p; ++j) {
//     score_test[j] = u_theta[j]*u_theta[j]/vtheta(j,j);
//   }
//   
//   // global score test for theta = 0
//   double score_test_global = 0.0;
//   NumericMatrix vtheta_inv = invsympd(vtheta, p, 1e-12);
//   for (int j=0; j<p; ++j) {
//     for (int l=0; l<p; ++l) {
//       score_test_global += u_theta[j]*vtheta_inv(j,l)*u_theta[l];
//     }
//   }
//   
//   score_test.push_back(score_test_global);
//   NumericVector df(p+1,1.0);
//   df[p] = p; // degrees of freedom for the global test
//   
//   // p-values
//   NumericVector p_values(p+1);
//   for (int j=0; j<=p; ++j) {
//     p_values[j] = 1.0 - boost_pchisq(score_test[j], df[j], true, false);
//   }
//   
//   // assemble the table
//   NumericMatrix table(p+1,3);
//   for (int j=0; j<=p; ++j) {
//     table(j,0) = score_test[j];
//     table(j,1) = df[j];
//     table(j,2) = p_values[j];
//   }
//   
//   StringVector covariates_vec(p+1);
//   for (int j=0; j<p; ++j) covariates_vec[j] = covariates[j];
//   covariates_vec[p] = "GLOBAL";
//   rownames(table) = covariates_vec;
//   colnames(table) = StringVector::create("chisq", "df", "p");
//   
//   // obtain scaled schoenfeld residuals
//   List resid_list = residuals_phregcpp(
//     p, beta, vbeta, resmart, data, stratum, time, time2, event,
//     covariates, weight, offset, "", ties, "scaledsch", false, true);
//   
//   NumericMatrix sresid = as<NumericMatrix>(resid_list["resid"]);
//   colnames(sresid) = covariates;
//   
//   // repeat the g values according to the number of events
//   NumericVector g_rep(nq);
//   int idx = 0;
//   for (int i=0; i<nt; ++i) {
//     int count = 0;
//     for (int j=0; j<nq; ++j) {
//       if (tstop3[j] == t[i]) {
//         count++;
//       }
//     }
//     for (int j=0; j<count; ++j) {
//       g_rep[idx] = g[i];
//       idx++;
//     }
//   }
//   
//   List result = List::create(
//     Named("table") = table,
//     Named("x") = g_rep,
//     Named("time") = tstop3,
//     Named("y") = sresid,
//     Named("var") = vbeta*nq,
//     Named("transform") = transform
//   );
//   
//   if (has_stratum) {
//     result.push_back(stratum3, "strata");
//   }
//   
//   return result;
// }